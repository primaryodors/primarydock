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
    if (hxno<7)
    {
        std::string name = (std::string)"TMR" + to_string(hxno+1);
        if (hxno == 6)
        {
            std::string name = (std::string)"TMR" + to_string(hxno);
            sr1 = p.get_region_start(name)-2;
            er1 = aasrc->get_residue_no();
        }
        else
        {
            sr1 = p.get_region_start(name)-2;
            er1 = p.get_end_resno(); // p.get_region_end(name)+2;
        }
    }

    Point source = aasrc->get_CA_location();
    Point was = aaref->get_CA_location().subtract(x50);
    Point target = was.add(rel_tgt);

    float theta = find_3d_angle(was, target, source);
    LocatedVector axis = compute_normal(was, target, source);
    axis.origin = source;

    float theta1 = p.region_can_rotate(sr, er, axis, true, 0, sr1, er1);
    if (theta1 < theta) theta = theta1;

    p.rotate_piece(sr, er, axis.origin, axis, theta);
    return theta;
}

float reduce_iclash_iter(Protein& p, int& fulcrum, bool cterm, float clash, float& theta, int region_start, int region_end)
{
    AminoAcid* aafulcrum = p.get_residue(fulcrum);
    LocatedVector axis = compute_normal(p.stop2->get_CA_location(), p.stop1->get_CA_location(), aafulcrum->get_CA_location());
    axis.origin = aafulcrum->get_CA_location();

    p.rotate_piece(cterm ? fulcrum : region_start, cterm ? region_end : fulcrum, axis.origin, axis, theta);
    float new_clash = p.get_internal_clashes(cterm ? fulcrum : region_start, cterm ? region_end : fulcrum);

    if (new_clash > clash)
    {
        p.rotate_piece(cterm ? fulcrum : region_start, cterm ? region_end : fulcrum, axis.origin, axis, -theta);
        theta *= -0.8;
    }
    else
    {
        clash = new_clash;
    }

    int nstop1 = p.stop1->get_residue_no();
    if (nstop1 < region_start || nstop1 > region_end)
    {
        // return 0;
    }
    else if (p.last_int_clash_dir.r <= clash_limit_per_aa) fulcrum = nstop1;

    // if (new_clash <= initial_clash) return 0;
    if (fabs(theta) < 1e-6) return 0;

    if (fulcrum == nstop1) return 0;

    return new_clash;
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
    char* override_ofname = nullptr;

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
            else if (!strcmp(argv[i], "-o")) override_ofname = argv[++i];
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
    std::string out_filename = override_ofname ? (std::string)override_ofname : (path + orid + (std::string)".active.pdb");
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
    AminoAcid *aa4x64 = p.get_residue_bw("4.64");

    AminoAcid *aa45x51 = p.get_residue_bw("45.51");
    AminoAcid *aa45x52 = p.get_residue_bw("45.52");
    AminoAcid *aa45x53 = p.get_residue_bw("45.53");
    AminoAcid *aa45x54 = p.get_residue_bw("45.54");

    AminoAcid *aa5x33 = p.get_residue_bw("5.33");
    AminoAcid *aa5x35 = p.get_residue_bw("5.35");
    AminoAcid *aa5x39 = p.get_residue_bw("5.39");
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
    AminoAcid *aa6x50 = p.get_residue_bw("6.50");
    AminoAcid *aa6x51 = p.get_residue_bw("6.51");
    AminoAcid *aa6x55 = p.get_residue_bw("6.55");
    AminoAcid *aa6x58 = p.get_residue_bw("6.58");
    AminoAcid *aa6x59 = p.get_residue_bw("6.59");

    AminoAcid *aa7x31 = p.get_residue_bw("7.31");
    AminoAcid *aa7x49 = p.get_residue_bw("7.49");
    AminoAcid *aa7x53 = p.get_residue_bw("7.53");

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
    int n6x50 = aa6x50->get_residue_no();
    int n6x55 = aa6x55->get_residue_no();
    int n6x59 = aa6x59->get_residue_no();

    int n7x31 = aa7x31->get_residue_no();
    int n7x49 = aa7x49->get_residue_no();
    int n7x53 = aa7x53->get_residue_no();

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


    ////////////////////////////////////////////////////////////////////////////////
    // Make room for TMR5 to shift.
    ////////////////////////////////////////////////////////////////////////////////

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
    }


    ////////////////////////////////////////////////////////////////////////////////
    // Compute initial internal clashes.
    ////////////////////////////////////////////////////////////////////////////////

    float initial_clash_[10];

    for (i=1; i<=7; i++)
    {
        std::string region = (std::string)"TMR" + std::to_string(i);
        initial_clash_[i] = p.get_internal_clashes(p.get_region_start(region), p.get_region_end(region));
    }


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
                sr1 = p.get_region_start(name)-2;
                er1 = p.get_end_resno(); // p.get_region_end(name)+2;
            }

            SCoord motion = xl8[i];
            float r = p.region_can_move(sr, er, motion, true, sr1, er1);
            if (r < motion.r) motion.r = r;

            cout << "Applying helix " << i << " translation " << motion.r << "A of " << xl8[i].r << " limited by " << *(p.stop1) << "->" << *(p.stop2) << endl;
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

            cout << "Applying helix " << i << " extracellular bend ";
            do_template_bend(p, aasrc, aaterm, i, exr[i], xl8[i]);
            cout << " limited by " << *(p.stop1) << "->" << *(p.stop2) << endl;
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

            cout << "Applying helix " << i << " cytoplasmic bend ";
            do_template_bend(p, aasrc, aaterm, i, cyt[i], xl8[i]);
            cout << " limited by " << *(p.stop1) << "->" << *(p.stop2) << endl;
        }
    }



    ////////////////////////////////////////////////////////////////////////////////
    // Parameters and eligibility for rock6 motion.
    ////////////////////////////////////////////////////////////////////////////////

    SCoord rock6_dir(0,0,0);
    bool exr2_bend = false;

    AminoAcid* acid45 = nullptr, *ohbridge6 = nullptr;
    if (l6x55 == 'Y' && (l45x51 == 'D' || l45x51 == 'E'))
    {
        acid45 = aa45x51;
        ohbridge6 = aa6x55;
    }
    else if (l6x55 == 'Y' && (l45x52 == 'D' || l45x52 == 'E'))
    {
        acid45 = aa45x52;
        ohbridge6 = aa6x55;
    }
    else if (l6x55 == 'Y' && (l45x53 == 'D' || l45x53 == 'E' || l45x53 == 'N' || l45x53 == 'Q'))
    {
        acid45 = aa45x53;
        ohbridge6 = aa6x55;
    }
    else if (l6x55 == 'D' && (l45x53 == 'N' || l45x53 == 'Q'))
    {
        acid45 = aa45x53;
        ohbridge6 = aa6x55;
    }

    if (acid45 && ohbridge6)
    {
        p.bridge(acid45->get_residue_no(), ohbridge6->get_residue_no());
        Atom* reach6x55 = ohbridge6->get_reach_atom();
        Atom* reach45x51 = acid45->get_nearest_atom(reach6x55->get_location());
        float r = reach6x55->distance_to(reach45x51);
        if (r >= 2) rock6_dir = reach45x51->get_location().subtract(reach6x55->get_location());
        constraints.push_back("STCR "+std::to_string(acid45->get_residue_no()));
        constraints.push_back("STCR "+std::to_string(ohbridge6->get_residue_no()));
    }
    else if (l6x59 == 'R')
    {
        rock6_dir = aa5x39->get_CA_location().subtract(aa6x55->get_CA_location());
        aa6x59->conform_atom_to_location("NE", Point(5000,10000,0));
        aa6x59->movability = MOV_PINNED;
        exr2_bend = true;

        pt = aa4x64->get_CA_location().add(aa5x35->get_CA_location());
        pt.multiply(0.5);
        aa45x52->conform_atom_to_location(aa45x52->get_reach_atom()->name, pt);
        aa45x52->movability = MOV_PINNED;
    }
    else if ((l6x55 == 'D' || l6x55 == 'E') && (l45x53 == 'N' || l45x53 == 'Q'))
    {
        p.bridge(n45x53, n6x55);
        Atom* reach6x55 = aa6x55->get_reach_atom();
        Atom* reach45x53 = aa45x53->get_nearest_atom(reach6x55->get_location());
        float r = reach6x55->distance_to(reach45x53);
        if (r >= 2) rock6_dir = reach45x53->get_location().subtract(reach6x55->get_location());
        constraints.push_back("STCR 45.53");
        constraints.push_back("STCR 6.55");
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
        if (l6x59 == 'R')
        {
            DynamicMotion dyn(&p);
            dyn.type = dyn_wind;
            dyn.start_resno = BallesterosWeinstein("6.56");
            dyn.end_resno = BallesterosWeinstein("6.60");
            dyn.bias = -12;
            dyn.apply_absolute(1);
        }

        cout << "Performing rock6..." << endl;
        float theta = do_template_bend(p, aa6x48, aa6x59, 6, rock6_dir, SCoord(0,0,0), aa6x28);
        cout << "TMR6 rocks " << (theta*fiftyseven) << "deg limited by " << *(p.stop1) << "->" << *(p.stop2) << endl;

        if (acid45 && ohbridge6) p.bridge(acid45->get_residue_no(), ohbridge6->get_residue_no());
    }

    
    ////////////////////////////////////////////////////////////////////////////////
    // FYG motion, if applicable.
    ////////////////////////////////////////////////////////////////////////////////

    AminoAcid* aapivot = aa6x49;
    int npivot = n6x49;
    bool found_fyg = false;

    if (allow_fyg)
    {
        for (i=n6x49; i>n6x40; i--)
        {
            AminoAcid* aatmp = p.get_residue(i);
            char ltmp = aatmp->get_letter();
            if (ltmp == 'G' || ltmp == 'C' || ltmp == 'S' || ltmp == 'T')
            {
                aapivot = aatmp;
                npivot = i;
                found_fyg = true;
                break;
            }
        }
    }

    if (found_fyg && allow_fyg)
    {
        cout << "Performing FYG activation using " << *aapivot << "..." << endl;

        // TMR6 motion.
        float theta6_phi_initial = 10;
        LocatedVector lv6 = aa6x49->get_phi_vector();
        p.rotate_piece(n6x28, n6x49, lv6.origin, lv6, -fiftyseventh*theta6_phi_initial);
        float theta6_psi_initial = 40;
        lv6 = aa6x49->get_psi_vector();
        p.rotate_piece(n6x28, n6x49, lv6.origin, lv6, -fiftyseventh*theta6_psi_initial);

        // TMR5 translation.
        SCoord move5 = aapivot->get_CA_location().subtract(aa5x50->get_CA_location());
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
        float theta6 = p.region_can_rotate(n6x28, npivot, lv6);
        Point pt6x28 = rotate3D(aa6x28->get_CA_location(), aapivot->get_CA_location(), lv6, theta6);
        float r56 = pt6x28.get_3d_distance(aa5x68->get_CA_location());
        cout << "r56 = " << r56 << endl;
        if (r56 > 4)
        {
            theta6 = fmin(theta6, find_3d_angle(aa5x68->get_CA_location(), aa6x28->get_CA_location(), aapivot->get_CA_location()));
            p.stop1 = aa5x68;
            p.stop2 = aa6x28;
        }
        cout << "TMR6 bends " << (theta6_phi_initial + theta6_psi_initial - theta6 * fiftyseven) << "deg limited by " << *(p.stop1) << "->" << *(p.stop2) << endl;
        p.rotate_piece(n6x28, npivot, lv6.origin, lv6, theta6);
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

        sr = n45x54 + 1;
        er = n45x52 + 10;
        p.delete_residues(sr, er);

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

        p.bridge(n5x58, n7x53);

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
    // Minimize internal clashes.
    // TODO: cyt domains of all TM helices, preserving the 5~7 bridge.
    ////////////////////////////////////////////////////////////////////////////////
    float clash = p.get_internal_clashes(n6x28, npivot, true);

    if (clash > initial_clash_[6])
    {
        cout << "Minimizing TMR6 cytoplasmic clashes..." << endl;

        if (p.stop1 && p.stop2) pt = aa6x40->get_CA_location().add(p.stop1->get_CA_location().subtract(p.stop2->get_CA_location()));
        else pt = aa6x40->get_CA_location().subtract(p.last_int_clash_dir);
        LocatedVector axis = compute_normal(aa6x40->get_CA_location(), pt, aapivot->get_CA_location());
        axis.origin = aapivot->get_CA_location();

        float theta = fiftyseventh * 15 / 20;
        int fulcrum = npivot;
        for (i=0; i<200; i++)
        {
            clash = reduce_iclash_iter(p, fulcrum, false, clash, theta, n6x28, fulcrum);
            if (clash <= initial_clash_[6]) break;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Make room near trip switch.
    // Assume at least one of rock6 / FYG bend has occurred.
    ////////////////////////////////////////////////////////////////////////////////

    if (l5x47 == 'F' || l5x47 == 'L' || l5x47 == 'I' || l5x47 == 'H')
    {
        pt = aa5x47->get_CA_location().subtract(aa3x34->get_CA_location()).add(aa5x47->get_CA_location());
        aa5x47->conform_atom_to_location(aa5x47->get_reach_atom()->name, pt);

        // TODO: 5.44 bend so 5.47 doesn't clash with 5.43.
    }

    // This side chain shift is observed in all cryo-EM models of active states of TAARs:
    if (l6x44 == 'F' || l6x44 == 'Y')
    {
        aa6x44->movability = MOV_FLEXONLY;
        aa6x44->conform_atom_to_location(aa6x44->get_reach_atom()->name, aa5x51->get_CA_location());
        constraints.push_back("STCR 6.44");
        aa6x44->movability = MOV_PINNED;
    }

    if (l6x59 == 'R')
    {
        aa6x59->movability = MOV_FLEXONLY;
        pt = aa6x59->get_CA_location().add(aa5x39->get_CA_location());
        pt.multiply(0.5);
        aa6x59->conform_atom_to_location("NE", pt);
    }
    else if (l6x59 == 'K')
    {
        aa6x59->movability = MOV_FLEXONLY;
        aa6x59->conform_atom_to_location(aa6x59->get_reach_atom()->name, aa4x60->get_CA_location());
    }

    aa5x58->movability  = MOV_FLEXONLY;
    aa7x53->movability  = MOV_FLEXONLY;
    if (acid45 && ohbridge6)
    {
        acid45->movability = MOV_FLEXONLY;
        ohbridge6->movability  = MOV_FLEXONLY;
        p.bridge(acid45->get_residue_no(), ohbridge6->get_residue_no());
        acid45->movability = MOV_PINNED;
        ohbridge6->movability = MOV_PINNED;
    }
    if (l5x58 == 'Y' && l7x53 == 'Y')
    {
        p.bridge(n5x58, n7x53);
        aa5x58->movability = MOV_PINNED;
        aa7x53->movability = MOV_PINNED;
    }

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

        if (!override_ofname)
        {
            fp = fopen(cns_filename.c_str(), "wb");
            if (!fp) return -3;
            n = constraints.size();
            for (l=0; l<n; l++) fprintf(fp, "%s\n", constraints[l].c_str());
            fclose(fp);
        }
    }
}
