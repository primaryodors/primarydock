#include <cstring>
#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include <unistd.h>
#include <fstream>
#include <sys/stat.h>
#include <sys/wait.h>
#include "classes/protein.h"
#include "classes/dynamic.h"

using namespace std;

#define dbg_57_contact 0

// TODO: These should be stored in a dedicated config file and made available to the predict/bnav2.pepd script.
const BallesterosWeinstein ebw[8] = { "0.0", "1.40", "2.63", "3.22", "4.60", "5.35", "6.55", "7.32" };
const BallesterosWeinstein cbw[8] = { "0.0", "1.55", "2.40", "3.55", "4.40", "5.58", "6.40", "7.53" };

void save_file(Protein& p, std::string filename, Molecule* ligand = nullptr)
{
    FILE* fp = fopen(filename.c_str(), "wb");
    if (!fp) throw -3;
    p.save_pdb(fp, ligand);
    p.end_pdb(fp);
    fclose(fp);
    cout << "Saved " << filename << endl;
}

float do_template_bend(Protein& p, AminoAcid* aasrc, AminoAcid* aaref, int hxno, SCoord rel_tgt, SCoord x50, AminoAcid* aaopposite_term = nullptr)
{
    int ressrc = aasrc->get_residue_no();
    int resref = aaref->get_residue_no();
    int resoppo = aaopposite_term ? aaopposite_term->get_residue_no() : ressrc;

    std::string region = (std::string)"TMR" + std::to_string(hxno);
    int resterm = (resref < ressrc) ? p.get_region_start(region) : p.get_region_end(region);

    int sr = min(resoppo, resterm);
    int er = max(resoppo, resterm);

    int sr1=0, er1=0;
    if (hxno<7 && !aaopposite_term)
    {
        std::string name = (std::string)"TMR" + to_string(hxno+1);
        sr1 = p.get_region_start(name);
        er1 = p.get_region_end(name);
    }

    Point source = aasrc->get_CA_location();
    Point was = aaref->get_CA_location().subtract(x50);
    Point target = was.add(rel_tgt);

    float theta = find_3d_angle(was, target, source);
    LocatedVector axis = compute_normal(was, target, source);
    axis.origin = source;

    float theta1 = p.region_can_rotate(sr, er, axis, false, 0, sr1, er1);
    if (theta1 < theta) theta = theta1;

    p.rotate_piece(sr, er, axis.origin, axis, theta);
    return theta;
}

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        cout << "Too few args." << endl;
        return -1;
    }
    int i, j, l, m, n, sr, er;
    char c;

    int fam = 0, mem = 0;
    std::string sub;
    std::string protid;
    std::vector<std::string> constraints;

    bool allow_fyg = true;
    bool allow_rock6 = true;
    bool allow_save = true;

    SCoord xl8[10];
    SCoord exr[10];
    SCoord cyt[10];
    Point pt;

    ////////////////////////////////////////////////////////////////////////////////
    // Read program arguments
    ////////////////////////////////////////////////////////////////////////////////

    for (i=1; i<argc; i++)
    {
        if (argv[i][0] == '-')
        {
            if (!strcmp(argv[i], "--nofyg")) allow_fyg = false;
            else if (!strcmp(argv[i], "--norock")) allow_rock6 = false;
            else if (!strcmp(argv[i], "--nosave")) allow_save = false;
            else if (argv[i][1] == '-'
                && (argv[i][2] == 'x' || argv[i][2] == 'e' || argv[i][2] == 'c')
                && argv[i][3] >= '1' && argv[i][3] <= '7'
                && !argv[i][4]
                )
            {
                c = argv[i][2];
                j = atoi(argv[i]+3);
                i++; if (i >= argc)
                {
                    cout << "Too few args." << endl;
                    return -1;
                }
                pt.x = atof(argv[i]);
                i++; if (i >= argc)
                {
                    cout << "Too few args." << endl;
                    return -1;
                }
                pt.y = atof(argv[i]);
                i++; if (i >= argc)
                {
                    cout << "Too few args." << endl;
                    return -1;
                }
                pt.z = atof(argv[i]);

                if (c == 'x') xl8[j] = pt;
                else if (c == 'e') exr[j] = pt;
                else if (c == 'c') cyt[j] = pt;
            }
            else cout << "Warning: unrecognized cli argument " << argv[i] << endl;

            continue;
        }

        // One required argument, in the form of a protein ID, e.g. OR51E2.
        // If the ID does not conform to the format of OR{family}{subfamily}{member}, return an error code.
        if (argv[i][0] == 'O' && argv[i][1] == 'R')
        {
            l = 2;
            std::string tmp = argv[i];
            fam = atoi(tmp.substr(l, 2).c_str());
            if (fam)
            {
                l++;
                if (fam >= 10) l++;
                sub = tmp.substr(l, 1);
                l++;
                if (argv[i][l] >= 'A' && argv[i][l] <= 'Z')
                {
                    sub = tmp.substr(l-1, 2);
                    l++;
                }
                mem = atoi(&argv[i][l]);
                if (!mem)
                {
                    cout << "Invalid OR." << endl;
                    return -2;
                }
            }
        }
        else protid = argv[i];
    }

    if (cyt[6].r) allow_fyg = false;
    if (exr[6].r || cyt[6].r) allow_rock6 = false;


    ////////////////////////////////////////////////////////////////////////////////
    // Load protein
    ////////////////////////////////////////////////////////////////////////////////

    // Locate and load the .upright.pdb in the folder structure, or throw an error if not found.
    std::string path;
    std::string orid;

    if (fam && sub.size() && mem)
    {
        path = (std::string)"pdbs/OR" + std::to_string(fam) + (std::string)"/";
        orid = (std::string)"OR" + std::to_string(fam) + sub + std::to_string(mem);
    }
    else
    {
        orid = protid;
        path = (std::string)"pdbs/" + protid.substr(0, 4) + (std::string)"/";
    }

    std::string in_filename = path + orid + (std::string)".upright.pdb";
    if (!file_exists(in_filename)) return -3;
    std::string out_filename = path + orid + (std::string)".active.pdb";
    std::string cns_filename = path + orid + (std::string)".params";

    Protein p(orid.c_str());
    FILE* fp = fopen(in_filename.c_str(), "rb");
    if (!fp) return -3;
    p.load_pdb(fp);
    fclose(fp);

    cout << "Internal clashes of starting model: " << p.get_internal_clashes() << endl;


    ////////////////////////////////////////////////////////////////////////////////
    // Set up residue vars
    ////////////////////////////////////////////////////////////////////////////////

    AminoAcid *aa1x32 = p.get_residue_bw("1.32");
    AminoAcid *aa1x50 = p.get_residue_bw("1.50");
    AminoAcid *aa1x58 = p.get_residue_bw("1.58");

    AminoAcid *aa2x38 = p.get_residue_bw("2.38");
    AminoAcid *aa2x50 = p.get_residue_bw("2.50");
    AminoAcid *aa2x66 = p.get_residue_bw("2.66");

    AminoAcid *aa3x21 = p.get_residue_bw("3.21");
    AminoAcid *aa3x33 = p.get_residue_bw("3.33");
    AminoAcid *aa3x34 = p.get_residue_bw("3.34");
    AminoAcid *aa3x36 = p.get_residue_bw("3.36");
    AminoAcid *aa3x37 = p.get_residue_bw("3.37");
    AminoAcid *aa3x40 = p.get_residue_bw("3.40");
    AminoAcid *aa3x50 = p.get_residue_bw("3.50");
    AminoAcid *aa3x56 = p.get_residue_bw("3.56");

    AminoAcid *aa4x49 = p.get_residue_bw("4.49");
    AminoAcid *aa4x60 = p.get_residue_bw("4.60");

    AminoAcid *aa45x51 = p.get_residue_bw("45.51");
    AminoAcid *aa45x52 = p.get_residue_bw("45.52");
    AminoAcid *aa45x53 = p.get_residue_bw("45.53");
    AminoAcid *aa45x54 = p.get_residue_bw("45.54");

    AminoAcid *aa5x33 = p.get_residue_bw("5.33");
    AminoAcid *aa5x47 = p.get_residue_bw("5.47");
    AminoAcid *aa5x50 = p.get_residue_bw("5.50");
    AminoAcid *aa5x51 = p.get_residue_bw("5.51");
    AminoAcid *aa5x58 = p.get_residue_bw("5.58");
    AminoAcid *aa5x68 = p.get_residue_bw("5.68");

    AminoAcid *aa56x50 = p.get_residue_bw("56.50");

    AminoAcid *aa6x28 = p.get_residue_bw("6.28");
    AminoAcid *aa6x40 = p.get_residue_bw("6.40");
    AminoAcid *aa6x44 = p.get_residue_bw("6.44");
    AminoAcid *aa6x48 = p.get_residue_bw("6.48");
    AminoAcid *aa6x49 = p.get_residue_bw("6.49");
    AminoAcid *aa6x51 = p.get_residue_bw("6.51");
    AminoAcid *aa6x55 = p.get_residue_bw("6.55");
    AminoAcid *aa6x59 = p.get_residue_bw("6.59");

    AminoAcid *aa7x31 = p.get_residue_bw("7.31");
    AminoAcid *aa7x49 = p.get_residue_bw("7.49");
    AminoAcid *aa7x53 = p.get_residue_bw("7.53");
    
    AminoAcid *aa8x44 = p.get_residue_bw("8.44");

    int n1x32 = aa1x32->get_residue_no();
    int n1x50 = aa1x50->get_residue_no();
    int n1x58 = aa1x58->get_residue_no();

    int n2x38 = aa2x38->get_residue_no();
    int n2x50 = aa2x50->get_residue_no();
    int n2x66 = aa2x66->get_residue_no();

    int n3x21 = aa3x21->get_residue_no();
    int n3x50 = aa3x50->get_residue_no();
    int n3x56 = aa3x56->get_residue_no();

    int n4x60 = aa4x60->get_residue_no();
    
    int n45x51 = aa45x51->get_residue_no();
    int n45x52 = aa45x52->get_residue_no();
    int n45x53 = aa45x53->get_residue_no();
    int n45x54 = aa45x54->get_residue_no();

    int n5x33 = aa5x33->get_residue_no();
    int n5x50 = aa5x50->get_residue_no();
    int n5x58 = aa5x58->get_residue_no();
    int n5x68 = aa5x68->get_residue_no();

    int n56x50 = aa56x50->get_residue_no();

    int n6x28 = aa6x28->get_residue_no();
    int n6x40 = aa6x40->get_residue_no();
    int n6x48 = aa6x48->get_residue_no();
    int n6x49 = aa6x49->get_residue_no();
    int n6x55 = aa6x55->get_residue_no();
    int n6x59 = aa6x59->get_residue_no();

    int n7x31 = aa7x31->get_residue_no();
    int n7x49 = aa7x49->get_residue_no();
    int n7x53 = aa7x53->get_residue_no();

    int n8x44 = aa8x44->get_residue_no();

    char l3x37 = aa3x37->get_letter();
    char l3x40 = aa3x40->get_letter();

    char l45x51 = aa45x51->get_letter();
    char l45x52 = aa45x52->get_letter();
    char l45x53 = aa45x53->get_letter();

    char l5x47 = aa5x47->get_letter();
    char l5x50 = aa5x50->get_letter();
    char l5x58 = aa5x58->get_letter();

    char l6x44 = aa6x44->get_letter();
    char l6x48 = aa6x48->get_letter();
    char l6x49 = aa6x49->get_letter();
    char l6x55 = aa6x55->get_letter();
    char l6x59 = aa6x59->get_letter();
    
    char l7x53 = aa7x53->get_letter();


    float initial_clash_6 = p.get_internal_clashes(n6x28, n6x49);


    ////////////////////////////////////////////////////////////////////////////////
    // n.50 translations.
    ////////////////////////////////////////////////////////////////////////////////
    for (i=1; i<=7; i++)
    {
        if (xl8[i].r)
        {
            std::string name = (std::string)"TMR" + to_string(i);
            sr = p.get_region_start(name);
            er = p.get_region_end(name);

            int sr1=0, er1=0;
            if (i<7)
            {
                name = (std::string)"TMR" + to_string(i+1);
                sr1 = p.get_region_start(name);
                er1 = p.get_region_end(name);
            }

            SCoord motion = xl8[i];
            float r = p.region_can_move(sr, er, motion, false, sr1, er1);
            if (r < motion.r) motion.r = r;

            cout << "Applying helix " << i << " translation..." << endl;
            p.move_piece(sr, er, motion);
        }
    }


    ////////////////////////////////////////////////////////////////////////////////
    // EXR and CYT bends.
    ////////////////////////////////////////////////////////////////////////////////
    for (i=1; i<=7; i++)
    {
        if (exr[i].r)
        {
            AminoAcid* aasrc = p.get_residue_bw(i, 50);
            AminoAcid* aaterm = p.get_residue(ebw[i]);
            if (!aasrc || !aaterm)
            {
                cout << "Error: helix " << i << " EXR terminus or BW50 residue not found in PDB model." << endl;
                return -1;
            }

            cout << "Applying helix " << i << " extracellular bend..." << endl;
            do_template_bend(p, aasrc, aaterm, i, exr[i], xl8[i]);
        }

        if (cyt[i].r)
        {
            AminoAcid* aasrc = p.get_residue_bw(i, 50);
            AminoAcid* aaterm = p.get_residue(cbw[i]);
            if (!aasrc || !aaterm)
            {
                cout << "Error: helix " << i << " CYT terminus or BW50 residue not found in PDB model." << endl;
                return -1;
            }

            cout << "Applying helix " << i << " cytoplasmic bend..." << endl;
            do_template_bend(p, aasrc, aaterm, i, cyt[i], xl8[i]);
        }
    }



    ////////////////////////////////////////////////////////////////////////////////
    // Parameters and eligibility for rock6 motion.
    ////////////////////////////////////////////////////////////////////////////////

    SCoord rock6_dir(0,0,0);
    bool exr2_bend = false;
    if (l6x55 == 'Y' && (l45x51 == 'D' || l45x51 == 'E'))
    {
        p.bridge(aa45x51->get_residue_no(), aa6x55->get_residue_no());
        Atom* reach6x55 = aa6x55->get_reach_atom();
        Atom* reach45x51 = aa45x51->get_nearest_atom(reach6x55->get_location());
        float r = reach6x55->distance_to(reach45x51);
        if (r >= 2) rock6_dir = reach45x51->get_location().subtract(reach6x55->get_location());
        constraints.push_back("STCR 45.51");
        constraints.push_back("STCR 6.55");
    }
    else if (l6x59 == 'R')
    {
        rock6_dir = aa4x49->get_CA_location().subtract(aa6x55->get_CA_location());
        aa6x59->conform_atom_to_location("NE", aa6x59->get_CA_location().add(rock6_dir));
        exr2_bend = true;
    }


    ////////////////////////////////////////////////////////////////////////////////
    // Bend EXR2 so that TMR6 can rock.
    ////////////////////////////////////////////////////////////////////////////////

    if (exr2_bend)
    {
        LocatedVector axis = (SCoord)aa45x52->get_CA_location().subtract(aa45x54->get_CA_location());
        axis.origin = aa45x52->get_CA_location();
        p.rotate_piece(n45x52, n45x54, axis.origin, axis, -fiftyseventh*40);
        aa45x53->conform_atom_to_location(aa45x53->get_reach_atom()->name, aa45x53->get_CA_location().add(Point(0,-10000,0)));
    }


    ////////////////////////////////////////////////////////////////////////////////
    // Rocking motion, if applicable.
    ////////////////////////////////////////////////////////////////////////////////

    if (allow_rock6 && rock6_dir.r)
    {
        DynamicMotion dyn(&p);
        dyn.type = dyn_wind;
        dyn.start_resno = BallesterosWeinstein("6.56");
        dyn.end_resno = BallesterosWeinstein("6.60");
        dyn.bias = -18;
        dyn.apply_absolute(1);

        cout << "Performing rock6..." << endl;
        float theta = do_template_bend(p, aa6x48, aa6x59, 6, rock6_dir, SCoord(0,0,0), aa6x28);
        cout << "TMR6 rocks " << (theta*fiftyseven) << "deg limited by " << *(p.stop1) << "->" << *(p.stop2) << endl;
    }

    
    ////////////////////////////////////////////////////////////////////////////////
    // FYG motion, if applicable.
    ////////////////////////////////////////////////////////////////////////////////

    if (l6x49 == 'G' && allow_fyg)
    {
        cout << "Performing FYG activation..." << endl;

        // TMR6 motion.
        float theta6_initial = 40;
        LocatedVector lv6 = aa6x49->get_psi_vector();
        p.rotate_piece(n6x28, n6x49, lv6.origin, lv6, -fiftyseventh*theta6_initial);

        // TMR5 translation.
        SCoord move5 = aa6x49->get_CA_location().subtract(aa5x50->get_CA_location());
        move5.r = p.region_can_move(n5x33, n5x68, move5, false, n6x28, n6x59);
        p.move_piece(n5x33, n5x68, move5);
        cout << "TMR5 moves " << move5.r << "A toward TMR6, limited by " << *(p.stop1) << "->" << *(p.stop2) << endl;

        // 5.58 ~ 7.53 contact.
        p.bridge(n5x58, n7x53);

        LocatedVector lv5 = compute_normal(aa5x58->get_CA_location(), aa7x53->get_CA_location(), aa5x33->get_CA_location());
        lv5.origin = aa5x33->get_CA_location();

        float theta5 = p.region_can_rotate(n5x33, n5x68, lv5, false, 0, n6x28, n6x59);
        pt = aa7x53->get_reach_atom()->get_location().subtract(aa5x58->get_reach_atom()->get_location());
        pt.scale(fmax(0, pt.magnitude() - contact_r_5x58_7x53));
        if (pt.magnitude())
        {
            Point pt5x58 = aa5x58->get_reach_atom()->get_location();
            float f = find_3d_angle(pt5x58, pt5x58.add(pt), lv5.origin);
            if (f < theta5) theta5 = f;

            p.rotate_piece(n5x33, n5x68, lv5.origin, lv5, theta5);
            cout << "TMR5 rotates " << theta5*fiftyseven << "deg, limited by " << *(p.stop1) << "->" << *(p.stop2) << endl;
        }
        else cout << "5.58~7.53 bridge already met." << endl;

        // Adjust TMR6.
        float theta6 = p.region_can_rotate(n6x28, n6x49, lv6);
        cout << "TMR6 bends " << (theta6_initial - theta6 * fiftyseven) << "deg limited by " << *(p.stop1) << "->" << *(p.stop2) << endl;
        p.rotate_piece(n6x28, n6x49, lv6.origin, lv6, theta6);
    }


    if (!allow_save) return 0;


    ////////////////////////////////////////////////////////////////////////////////
    // 5-7 H-bond.
    ////////////////////////////////////////////////////////////////////////////////

    float r57, rcm, thcr;
    if (l5x58 == 'Y' && l7x53 == 'Y')
    {
        LocatedVector axis = (SCoord)aa7x53->get_atom_location("OH").subtract(aa5x58->get_atom_location("OH"));
        r57 = axis.r;

        // Move TMR5 toward TMR6.
        if (r57 > contact_r_5x58_7x53)
        {
            axis = (SCoord)aa6x48->get_CA_location().subtract(aa5x50->get_CA_location());
            rcm = p.region_can_move(n5x33, n5x68, axis, true, n6x28, n6x59);
            if (rcm < axis.r) axis.r = rcm;
            // axis.r += 0.25;
            p.move_piece(n5x33, n5x68, axis);
            cout << "TMR5 translation I " << axis.r << "A limited by " << *(p.stop1) << "->" << *(p.stop2) << endl;
            p.bridge(n5x58, n7x53);

            axis = (SCoord)aa7x53->get_atom_location("OH").subtract(aa5x58->get_atom_location("OH"));
            r57 = axis.r;
        }

        // Move TMR5 toward TMR7.
        if (r57 > contact_r_5x58_7x53)
        {
            rcm = p.region_can_move(n5x33, n5x68, axis, true, n6x28, n6x59);
            if (rcm < axis.r) axis.r = rcm;
            // axis.r += 0.25;
            p.move_piece(n5x33, n5x68, axis);
            cout << "TMR5 translation II " << axis.r << "A limited by " << *(p.stop1) << "->" << *(p.stop2) << endl;
            p.bridge(n5x58, n7x53);

            axis = (SCoord)aa7x53->get_atom_location("OH").subtract(aa5x58->get_atom_location("OH"));
            r57 = axis.r;
        }

        // Bend TMR7 toward TMR5.
        if (r57 > contact_r_5x58_7x53)
        {
            axis.r = contact_r_5x58_7x53;
            float theta = find_3d_angle(aa7x53->get_atom_location("OH"),
                aa5x58->get_atom_location("OH").add((SCoord)axis),
                aa7x49->get_CA_location());
            axis = (SCoord)compute_normal(aa7x53->get_atom_location("OH"), aa5x58->get_atom_location("OH"), aa7x49->get_CA_location());
            axis.origin = aa7x49->get_CA_location();
            thcr = p.region_can_rotate(n7x49, n7x53, axis, false, 0, n6x28, n6x59);
            thcr += fiftyseventh*10;
            p.rotate_piece(n7x49, n7x53, axis.origin, axis, fmin(theta, thcr));
            cout << "TMR7 bend " << theta*fiftyseven << "deg limited by " << *(p.stop1) << "->" << *(p.stop2) << endl;
            p.bridge(n5x58, n7x53);

            axis = (SCoord)aa7x53->get_atom_location("OH").subtract(aa5x58->get_atom_location("OH"));
            r57 = axis.r;
        }

        // TODO: TMR7 does this weird unwind-wind thing. The following code does not.
        if (false) // r57 > contact_r_5x58_7x53)
        {
            DynamicMotion dyn7u(&p);
            dyn7u.type = dyn_wind;
            dyn7u.start_resno = BallesterosWeinstein("7.50");
            dyn7u.end_resno = BallesterosWeinstein("7.53");
            dyn7u.bias = -50;
            dyn7u.apply_absolute(1);

            DynamicMotion dyn7w(&p);
            dyn7w.type = dyn_wind;
            dyn7w.start_resno = BallesterosWeinstein("7.53");
            dyn7w.end_resno = BallesterosWeinstein("7.57");
            dyn7w.bias = 50;
            dyn7w.apply_absolute(1);
        }

        // Tilt TMR5 toward TMR7.
        if (r57 > contact_r_5x58_7x53)
        {
            axis.r = contact_r_5x58_7x53;
            float theta = find_3d_angle(aa7x53->get_atom_location("OH").subtract((SCoord)axis),
                aa5x58->get_atom_location("OH"),
                aa5x33->get_CA_location());
            axis = (SCoord)compute_normal(aa5x58->get_atom_location("OH"), aa7x53->get_atom_location("OH"), aa5x33->get_CA_location());
            thcr = p.region_can_rotate(n5x33, n5x68, axis, false, 0, n6x28, n6x59);
            p.rotate_piece(n5x33, n5x68, axis.origin, axis, fmin(theta, thcr));
            cout << "TMR5 pivot " << theta*fiftyseven << "deg limited by " << *(p.stop1) << "->" << *(p.stop2) << endl;
            p.bridge(n5x58, n7x53);

            axis = (SCoord)aa7x53->get_atom_location("OH").subtract(aa5x58->get_atom_location("OH"));
            r57 = axis.r;
        }

        // Bend TMR5 toward TMR7.
        if (r57 > contact_r_5x58_7x53)
        {
            axis.r = contact_r_5x58_7x53;
            float theta = find_3d_angle(aa7x53->get_atom_location("OH").subtract((SCoord)axis),
                aa5x58->get_atom_location("OH"),
                aa5x50->get_CA_location());
            axis = (SCoord)compute_normal(aa5x58->get_atom_location("OH"), aa7x53->get_atom_location("OH"), aa5x50->get_CA_location());
            thcr = p.region_can_rotate(n5x50, n5x68, axis, false, 0, n6x28, n6x59);
            p.rotate_piece(n5x50, n5x68, axis.origin, axis, fmin(theta, thcr));
            cout << "TMR5 bend " << theta*fiftyseven << "deg limited by " << *(p.stop1) << "->" << *(p.stop2) << endl;
            p.bridge(n5x58, n7x53);

            axis = (SCoord)aa7x53->get_atom_location("OH").subtract(aa5x58->get_atom_location("OH"));
            r57 = axis.r;
        }

        // Check the result.
        if (r57 > 1.2 * contact_r_5x58_7x53)
        {
            cout << "WARNING: 5.58...7.53 H-bond FAILED (" << r57 << "A)." << endl;
        }
        else
        {
            cout << "5.58...7.53 contact distance: " << r57 << "A." << endl;
        }
    }


    ////////////////////////////////////////////////////////////////////////////////
    // Minimize internal clashes. TODO: exr and cyt domains of all TM helices.
    ////////////////////////////////////////////////////////////////////////////////
    float clash = p.get_internal_clashes(n6x28, n6x49, true);

    if (clash > initial_clash_6)
    {
        cout << "Minimizing TMR6 cytoplasmic clashes..." << endl;

        if (p.stop1 && p.stop2) pt = aa6x40->get_CA_location().add(p.stop1->get_CA_location().subtract(p.stop2->get_CA_location()));
        else pt = aa6x40->get_CA_location().subtract(p.last_int_clash_dir);
        LocatedVector axis = compute_normal(aa6x40->get_CA_location(), pt, aa6x49->get_CA_location());
        axis.origin = aa6x49->get_CA_location();

        float theta = fiftyseventh * 15 / 50;
        for (i=0; i<50; i++)
        {
            p.rotate_piece(n6x28, n6x49, axis.origin, axis, theta);
            float new_clash = p.get_internal_clashes(n6x28, n6x49);
            // cout << *p.stop1 << ", " << *p.stop2 << " " << theta << " " << new_clash << endl;

            if (new_clash > clash)
            {
                p.rotate_piece(n6x28, n6x49, axis.origin, axis, -theta);
                theta *= -0.75;
            }
            else
            {
                clash = new_clash;
            }
            if (new_clash <= initial_clash_6) break;
            if (fabs(theta) < 1e-5) break;

            axis = compute_normal(p.stop2->get_CA_location(), p.stop1->get_CA_location(), aa6x49->get_CA_location());
            axis.origin = aa6x49->get_CA_location();
        }
    }


    ////////////////////////////////////////////////////////////////////////////////
    // Make room near trip switch.
    // Assume at least one of rock6 / FYG bend has occurred.
    ////////////////////////////////////////////////////////////////////////////////

    Molecule inert;
    
    if (l6x48 != 'W')
    {
        aa5x47->movability = MOV_FLEXONLY;
        pt = aa5x47->get_CA_location().subtract(aa3x33->get_CA_location());
        pt = pt.add(aa5x47->get_CA_location());
        aa5x47->conform_atom_to_location(aa5x47->get_reach_atom()->name, pt);
        aa6x51->movability = MOV_FLEXONLY;
        pt = aa6x51->get_CA_location().subtract(aa3x36->get_CA_location());
        pt = pt.add(aa6x51->get_CA_location());
        aa6x51->conform_atom_to_location(aa6x51->get_reach_atom()->name, pt);

        inert.from_smiles("I[Si](I)(I)I");
        pt = aa6x48->get_CA_location().add(aa3x40->get_CA_location());
        pt = pt.add(aa3x36->get_CA_location());
        pt = pt.add(aa5x47->get_CA_location());
        pt = pt.add(aa6x51->get_CA_location());
        pt.multiply(1.0/5);
        inert.recenter(pt);
        inert.movability = MOV_NORECEN;

        AminoAcid* reachres[256];
        int nearby = p.get_residues_can_clash_ligand(reachres, &inert, pt, Point(8,8,8), nullptr);
        Molecule* mols[nearby+4];
        j = 0;
        mols[j++] = &inert;
        for (i=0; i<nearby; i++)
        {
            reachres[i]->movability = MOV_FLEXONLY;
            mols[j++] = reinterpret_cast<Molecule*>(reachres[i]);
        }
        mols[j] = nullptr;

        cout << "Flexing " << nearby << " side chains away from trip switch area." << endl;
        Molecule::conform_molecules(mols, 20);

        if (l3x37 == 'S' || l3x37 == 'N' || l3x37 == 'Q' || l3x37 == 'K' || l3x37 == 'R' || l3x37 == 'D' || l3x37 == 'E')
            aa3x37->conform_atom_to_location(aa3x37->get_reach_atom()->name, pt);
        
        if (l3x40 == 'S' || l3x40 == 'N' || l3x40 == 'Q' || l3x40 == 'K' || l3x40 == 'R' || l3x40 == 'D' || l3x40 == 'E')
            aa3x40->conform_atom_to_location(aa3x40->get_reach_atom()->name, pt);
    }

    if (l5x47 == 'F' || l5x47 == 'L' || l5x47 == 'I' || l5x47 == 'H')
    {
        pt = aa5x47->get_CA_location().subtract(aa3x34->get_CA_location()).add(aa5x47->get_CA_location());
        aa5x47->conform_atom_to_location(aa5x47->get_reach_atom()->name, pt);
    }

    // This side chain shift is observed in all cryo-EM models of active states of TAARs:
    if (l6x44 == 'F' || l6x44 == 'Y')
    {
        aa6x44->movability = MOV_FLEXONLY;
        aa6x44->conform_atom_to_location(aa6x44->get_reach_atom()->name, aa5x51->get_CA_location());
        constraints.push_back("STCR 6.44");
    }

    if (l6x59 == 'R' || l6x59 == 'K')
    {
        aa6x59->conform_atom_to_location(aa6x59->get_reach_atom()->name, aa4x60->get_CA_location());
    }


    ////////////////////////////////////////////////////////////////////////////////
    // Final active protein. No more changes can be made past this point.
    ////////////////////////////////////////////////////////////////////////////////

    cout << "Internal clashes of final model: " << p.get_internal_clashes() << endl;

    if (allow_save)
    {
        save_file(p, out_filename.c_str());

        ////////////////////////////////////////////////////////////////////////////////
        // Save parameters file.
        ////////////////////////////////////////////////////////////////////////////////

        fp = fopen(cns_filename.c_str(), "wb");
        if (!fp) return -3;
        n = constraints.size();
        for (l=0; l<n; l++) fprintf(fp, "%s\n", constraints[l].c_str());
        fclose(fp);
    }
}
