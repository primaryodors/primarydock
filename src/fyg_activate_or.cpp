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

void do_template_bend(Protein& p, AminoAcid* aasrc, AminoAcid* aaref, int hxno, SCoord rel_tgt, SCoord x50, AminoAcid* aaopposite_term = nullptr)
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

    Protein p(orid.c_str());
    FILE* fp = fopen(in_filename.c_str(), "rb");
    if (!fp) return -3;
    p.load_pdb(fp);
    fclose(fp);


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
    AminoAcid *aa3x50 = p.get_residue_bw("3.50");
    AminoAcid *aa3x56 = p.get_residue_bw("3.56");

    AminoAcid *aa4x49 = p.get_residue_bw("4.49");
    AminoAcid *aa4x60 = p.get_residue_bw("4.60");

    AminoAcid *aa45x51 = p.get_residue_bw("45.51");
    AminoAcid *aa45x52 = p.get_residue_bw("45.52");
    AminoAcid *aa45x53 = p.get_residue_bw("45.53");
    AminoAcid *aa45x54 = p.get_residue_bw("45.54");

    AminoAcid *aa5x33 = p.get_residue_bw("5.33");
    AminoAcid *aa5x50 = p.get_residue_bw("5.50");
    AminoAcid *aa5x58 = p.get_residue_bw("5.58");
    AminoAcid *aa5x68 = p.get_residue_bw("5.68");

    AminoAcid *aa56x50 = p.get_residue_bw("56.50");

    AminoAcid *aa6x28 = p.get_residue_bw("6.28");
    AminoAcid *aa6x40 = p.get_residue_bw("6.40");
    AminoAcid *aa6x48 = p.get_residue_bw("6.48");
    AminoAcid *aa6x49 = p.get_residue_bw("6.49");
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

    char l45x51 = aa45x51->get_letter();
    char l45x52 = aa45x52->get_letter();
    char l45x53 = aa45x53->get_letter();

    char l5x50 = aa5x50->get_letter();
    char l5x58 = aa5x58->get_letter();

    char l6x48 = aa6x48->get_letter();
    char l6x49 = aa6x49->get_letter();
    char l6x55 = aa6x55->get_letter();
    char l6x59 = aa6x59->get_letter();
    
    char l7x53 = aa7x53->get_letter();


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

            do_template_bend(p, aasrc, aaterm, i, cyt[i], xl8[i]);
        }
    }



    ////////////////////////////////////////////////////////////////////////////////
    // Rocking motion, if applicable.
    ////////////////////////////////////////////////////////////////////////////////

    SCoord rock6_dir(0,0,0);
    if (l6x59 == 'R')
    {
        rock6_dir = aa4x49->get_CA_location().subtract(aa6x55->get_CA_location());
        aa6x59->conform_atom_to_location("NE", aa6x59->get_CA_location().add(rock6_dir));
    }
    else if (l6x55 == 'Y' && (l45x51 == 'D' || l45x51 == 'E'))
    {
        p.bridge(aa45x51->get_residue_no(), aa6x55->get_residue_no());
        Atom* reach6x55 = aa6x55->get_reach_atom();
        Atom* reach45x51 = aa45x51->get_nearest_atom(reach6x55->get_location());
        float r = reach6x55->distance_to(reach45x51);
        if (r >= 2) rock6_dir = reach45x51->get_location().subtract(reach6x55->get_location());
    }

    if (allow_rock6 && rock6_dir.r)
    {
        cout << "Performing rock6..." << endl;
        do_template_bend(p, aa6x48, aa6x59, 6, rock6_dir, SCoord(0,0,0), aa6x28);
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
        cout << "TMR6 rotated " << (theta6_initial - theta6 * fiftyseven) << "deg." << endl;
        p.rotate_piece(n6x28, n6x49, lv6.origin, lv6, theta6);
    }

    if (allow_save) save_file(p, out_filename.c_str());
}
