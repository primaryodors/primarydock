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

// TODO: These should be stored in a dedicated config file and made available to the predict/bnav2.pepd script.
const BallesterosWeinstein ebw[8] = { "0.0", "1.40", "2.63", "3.32", "4.58", "5.38", "6.55", "7.42" };
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

Point decode_cartesian_literal(char* c)
{
    if (c[0] == '[') c[0] = ' ';
    char* cy = strchr(c, ',');
    char* cz = cy ? strchr(cy+1, ',') : nullptr;

    Point p(atof(c), cy ? atof(cy+1) : 0, cz ? atof(cz+1) : 0);
    return p;
}


int main(int argc, char** argv)
{

    if (argc < 2)
    {
        cout << "Too few args." << endl;
        return -1;
    }
    int i, j, l, m, n, sr, er, iter;
    char c;

    int fam = 0, mem = 0;
    std::string sub;
    std::string protid;
    std::vector<std::string> constraints;
    std::string insfn;

    bool allow_save = true;
    char* override_ofname = nullptr;

    SCoord vec[10];
    float  rot[10];
    SCoord xl8[10];
    SCoord exr[10];
    SCoord cyt[10];
    float ptrn[10][10];
    Point pt;

    Point hxc[10];
    Point hxe[10];

    FILE* fp;

    for (i=0; i<10; i++)
    {
        vec[i] = xl8[i] = exr[i] = cyt[i] = SCoord(0,0,0);
        rot[i] = 0;
        for (j=0; j<10; j++) ptrn[i][j] = 0;
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Read program arguments
    ////////////////////////////////////////////////////////////////////////////////

    for (i=1; i<argc; i++)
    {
        if (argv[i][0] == '-')
        {
            if (!strcmp(argv[i], "--nosave")) allow_save = false;
            else if (!strcmp(argv[i], "-o")) override_ofname = argv[++i];
            else cout << "Warning: unrecognized cli argument " << argv[i] << endl;

            continue;
        }

        // Required argument, a file of instructions for bending.
        char* shsh = strstr(argv[i], ".shsh");
        if (file_exists(argv[i]) && shsh && !shsh[5])
        {
            insfn = argv[i];
            continue;
        }

        // Required argument, in the form of a protein ID, e.g. OR51E2.
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


    ////////////////////////////////////////////////////////////////////////////////
    // Load instructions
    ////////////////////////////////////////////////////////////////////////////////
    
    fp = fopen(insfn.c_str(), "rb");
    int tmrno = 0;
    char buffer[1024];
    while (!feof(fp))
    {
        char* rfw = fgets(buffer, 1013, fp);
        char** words = chop_spaced_words(buffer);

        if (words[0][0] == 'T' && words[0][1] == 'M' && words[0][2] == 'R'
            && words[0][3] >= '0' && words[0][3] <= '9'
            && words[0][4] == ':' && !words[0][5]
            )
        {
            tmrno = atoi(&words[0][3]);
        }
        else if (!strcmp(words[0], "XFORM"))
        {
            if (!tmrno)
            {
                cout << "Error in input file: XFORM without TMR number." << endl;
                return 1;
            }

            xl8[tmrno] = decode_cartesian_literal(words[1]);
        }
        else if (!strcmp(words[0], "EXR"))
        {
            if (!tmrno)
            {
                cout << "Error in input file: EXR without TMR number." << endl;
                return 1;
            }

            exr[tmrno] = decode_cartesian_literal(words[1]);
        }
        else if (!strcmp(words[0], "CYT"))
        {
            if (!tmrno)
            {
                cout << "Error in input file: CYT without TMR number." << endl;
                return 1;
            }

            cyt[tmrno] = decode_cartesian_literal(words[1]);
        }
        else if (!strcmp(words[0], "AXIS"))
        {
            if (!tmrno)
            {
                cout << "Error in input file: AXIS without TMR number." << endl;
                return 1;
            }

            vec[tmrno] = decode_cartesian_literal(words[1]);
        }
        else if (!strcmp(words[0], "ROT"))
        {
            if (!tmrno)
            {
                cout << "Error in input file: ROT without TMR number." << endl;
                return 1;
            }

            rot[tmrno] = atof(words[1]) * fiftyseventh;
        }
        else if (!strcmp(words[0], "WIND"))
        {
            if (!tmrno)
            {
                cout << "Error in input file: WIND without TMR number." << endl;
                return 1;
            }

            for (i=1; words[i]; i++)
            {
                if (i >= 10)
                {
                    cout << "Error in input file: too many WIND values." << endl;
                    return 1;
                }
                ptrn[tmrno][i] = atof(words[i]);
            }
        }
    }

    fclose(fp);



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
    fp = fopen(in_filename.c_str(), "rb");
    if (!fp) return -3;
    p.load_pdb(fp);
    fclose(fp);

    cout << "Internal clashes of starting model: " << p.get_internal_clashes() << endl;

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
    // Preliminary info for EXR and CYT bends.
    ////////////////////////////////////////////////////////////////////////////////
    for (i=1; i<=7; i++)
    {
        if (exr[i].r)
        {
            AminoAcid* aaterm = p.get_residue(ebw[i]);
            if (!aaterm)
            {
                cout << "Error: helix " << i << " EXR terminus or BW50 residue not found in PDB model." << endl;
                return -1;
            }

            hxe[i] = aaterm->get_CA_location().add(exr[i]);
        }
        if (cyt[i].r)
        {
            AminoAcid* aaterm = p.get_residue(cbw[i]);
            if (!aaterm)
            {
                cout << "Error: helix " << i << " CYT terminus or BW50 residue not found in PDB model." << endl;
                return -1;
            }

            hxc[i] = aaterm->get_CA_location().add(cyt[i]);
        }
    }

    ////////////////////////////////////////////////////////////////////////////////
    // n.50 rotations.
    ////////////////////////////////////////////////////////////////////////////////
    for (i=1; i<=7; i++)
    {
        if (rot[i])
        {
            std::string name = (std::string)"TMR" + to_string(i);
            sr = p.get_region_start(name);
            er = p.get_region_end(name);
            AminoAcid* pr = p.get_residue_bw(i, 50);
            pt = pr->get_CA_location();

            p.rotate_piece(sr, er, pt, vec[i], -rot[i]);
        }
    }

    save_file(p, "tmp/step1.pdb");


    ////////////////////////////////////////////////////////////////////////////////
    // helix perturns.
    ////////////////////////////////////////////////////////////////////////////////
    for (i=1; i<=7; i++)
    {
        std::string name = (std::string)"TMR" + to_string(i);
        sr = p.get_region_start(name);
        er = p.get_region_end(name);
        int b50 = p.get_bw50(i);

        if (i & 1) er = sr;
        sr = b50;

        for (j=1; ptrn[i][j]; j++)
        {
            int ter = sr + sgn(er-sr) * 4 * j;
            int tsr = ter - sgn(er-sr) * 4;

            // Point ps = p.get_residue(tsr)->get_CA_location();
            // Point pt = p.get_residue(ter)->get_CA_location();

            float t = p.helix_tightness(tsr, ter);
            float f = ptrn[i][j] - t;

            p.wind_helix(tsr, ter, f*5.3, er);

            // Rotation rot = align_points_3d(p.get_residue(er)->get_CA_location(), pt, ps);
            // p.rotate_piece(sr, er, rot, sr);
        }
    }

    save_file(p, "tmp/step2.pdb");


    ////////////////////////////////////////////////////////////////////////////////
    // n.50 translations.
    ////////////////////////////////////////////////////////////////////////////////
    /* for (iter = 0; iter < 10; iter++) */ for (i=1; i<=7; i++)
    {
        if (xl8[i].r >= 0.00001)
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
            // float r = p.region_can_move(sr, er, motion, true, sr1, er1);
            // if (r < motion.r) motion.r = r;

            cout << "Applying helix " << i << " translation " << motion.r << "A of " << xl8[i].r << " limited by " << *(p.stop1) << "->" << *(p.stop2) << endl;
            p.move_piece(sr, er, motion);
            xl8[i].r -= motion.r;
        }
    }

    save_file(p, "tmp/step3.pdb");


    ////////////////////////////////////////////////////////////////////////////////
    // EXR and CYT bends.
    ////////////////////////////////////////////////////////////////////////////////
    for (i=1; i<=7; i++)
    {
        bool odd = (i&1);

        AminoAcid* aasrc = p.get_residue_bw(i, 50);
        Point center = aasrc->get_CA_location();
        int b50 = aasrc->get_residue_no();

        std::string name = (std::string)"TMR" + to_string(i);
        sr = p.get_region_start(name);
        er = p.get_region_end(name);

        if (exr[i].r)
        {
            AminoAcid* aaterm = p.get_residue(ebw[i]);
            if (!aasrc || !aaterm)
            {
                cout << "Error: helix " << i << " EXR terminus or BW50 residue not found in PDB model." << endl;
                return -1;
            }

            cout << "Applying helix " << i << " extracellular bend ";

            Rotation rot = align_points_3d(aaterm->get_CA_location(), hxe[i], center);
            p.rotate_piece(odd ? sr : b50, odd ? b50 : er, center, rot.v, rot.a);

            cout << aaterm->get_CA_location() << " target was " << hxe[i];
            cout << endl;
        }

        if (cyt[i].r)
        {
            AminoAcid* aaterm = p.get_residue(cbw[i]);
            if (!aasrc || !aaterm)
            {
                cout << "Error: helix " << i << " CYT terminus or BW50 residue not found in PDB model." << endl;
                return -1;
            }

            cout << "Applying helix " << i << " cytoplasmic bend ";

            Rotation rot = align_points_3d(aaterm->get_CA_location(), hxc[i], center);
            p.rotate_piece(odd ? b50 : sr, odd ? er : b50, center, rot.v, rot.a);
            cout << (odd ? b50 : sr) << "-" << (odd ? er : b50) << " " << center << ", " << (rot.a*fiftyseven) << " deg. ";

            cout << aaterm->get_CA_location() << " target was " << hxc[i];
            cout << endl;
        }
    }

    save_file(p, "tmp/step4.pdb");


    ////////////////////////////////////////////////////////////////////////////////
    // Final reshaped protein. No more changes can be made past this point.
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

    return 0;
}
