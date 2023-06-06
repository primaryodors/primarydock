#include <cstring>
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include "classes/protein.h"
#include "classes/group.h"

using namespace std;

class Contact
{
    public:
    Protein *prot1 = nullptr, *prot2 = nullptr;
    AminoAcid *aa1 = nullptr, *aa2 = nullptr;
    std::string cfgstr1, cfgstr2;

    void interpret_cfgs()
    {
        // TODO:
    }
}

class MovablePiece
{
    public:
    Protein* prot = nullptr;
    AminoAcid *start_residue = nullptr, *end_residue = nullptr;
    AminoAcid *first_pivot = nullptr, *last_pivot = nullptr;
    std::vector<std::string> cfgstrs;

    void interpret_cfgs()
    {
        // TODO:
    }

    function do_motion(SCoord move_amt)
    {
        if (!prot) throw 0xbadc0de;
        if (!start_residue || !end_residue) throw 0xbadc0de;
        if (!first_pivot && !last_pivot) throw 0xbadc0de;

        // Check the direction of motion and correct if necessary.
        Point rel = move_amt;
        if (first_pivot)
        {
            Point pivf = first_pivot->get_CA_location();
            Point old = start_residue->get_CA_location();
            Point nouion = old.add(rel);
            nouion = nouion.multiply_3d_distance(pivf, old.get_3d_distance(pivf) / nouion.get_3d_distance(pivf) );
            rel = nouion.subtract(old);
        }
        if (last_pivot)
        {
            Point pivf = last_pivot->get_CA_location();
            Point old = end_residue->get_CA_location();
            Point nouion = old.add(rel);
            nouion = nouion.multiply_3d_distance(pivf, old.get_3d_distance(pivf) / nouion.get_3d_distance(pivf) );
            rel = nouion.subtract(old);
        }

        // Do the motion.
        prot->move_piece(start_residue->get_residue_no(), end_residue->get_residue_no(), (SCoord)rel);
        Point srca = prot->get_residue(start_residue->get_residue_no()-1)->get_CA_location().add(rel);
        Point eca  = prot->get_residue(end_residue->get_residue_no()+1)->get_CA_location().add(rel);

        // If first pivot, rotate first pivot thru start residue - 1 about first pivot to align with start residue.
        if (first_pivot)
        {
            prot->rotate_piece(first_pivot->get_residue_no(), start_residue->get_residue_no()-1, start_residue->get_residue_no()-1,
                srca, first_pivot->get_residue_no());
        }

        // If no first pivot, rotate head through start residue - 1 about last pivot to align with start residue.
        else
        {
            prot->rotate_piece(1, start_residue->get_residue_no()-1, start_residue->get_residue_no()-1,
                srca, last_pivot->get_residue_no());
        }

        // If last pivot, rotate end residue + 1 thru last pivot about last pivot to align with end residue.
        if (last_pivot)
        {
            prot->rotate_piece(end_residue->get_residue_no()+1, last_pivot->get_residue_no(), end_residue->get_residue_no()+1,
                eca, last_pivot->get_residue_no());
        }

        // If no last pivot, rotate rotate end residue + 1 thru tail about first pivot to align with end residue.
        else
        {
            prot->rotate_piece(end_residue->get_residue_no()+1, 99999, end_residue->get_residue_no()+1,
                eca, first_pivot->get_residue_no());
        }

        // Realign segment to new locations of start residue - 1 and end residue + 1.

    }

    function do_rotation(SCoord axis, float theta)
    {
        // Check the direction of rotation and correct if necessary.

        // TODO:

    }
}

Protein *g_prot1, *g_prot2;         // Global protein pointers; not necessarily G-proteins.
AminoAcid *pivot1, *pivot2, *pivot3, *pivot4, *pivot5, *pivot6, *pivot7;
Molecule** g_contacts_as_mols;
AminoAcid *sbb = nullptr, *sba = nullptr;

std::string prot1fname, prot2fname;
std::vector<Contact> contacts;
std::vector<MovablePiece> segments;

const float montecarlo_theta = fiftyseventh * 1;
const float montecarlo_xform = 0.5;
const float protein_clash_penalty = 10;

void show_usage()
{
    cout << "Usage:" << endl;
    cout << "couple path/to/GPCR.pdb path/to/G-protein.pdb";
    cout << endl;
}

float residue_energy()
{
    int i;
    float e = 0;

    for (i=0; g_contacts_as_mols[i]; i++)
    {
        e += g_contacts_as_mols[i]->get_intermol_binding(g_contacts_as_mols);
    }

    return e;
}

int interpret_cfg_param(char** words)
{
    int i;

    if (!strcmp(words[0], "PROT1"))
    {
        prot1fname = words[1];
        return 2;
    }
    else if (!strcmp(words[0], "PROT2"))
    {
        prot2fname = words[1];
        return 2;
    }
    else if (!strcmp(words[0], "CONTACT"))
    {
        Contact c;
        c.cfgstr1 = words[1];
        c.cfgstr2 = words[2];
        contacts.push_back(c);
        return 3;
    }
    else if (!strcmp(words[0], "SEGMENT"))
    {
        MovablePiece p;
        for (i=1; words[i]; i++)
            p.cfgstrs.push_back(words[i]);

        return i;
    }

    return 0;
}

void read_cfg_lines(FILE* fp)
{
    char buffer[1024];
    char** words;
    int i;
    while (!feof(fp))
    {
        fgets(buffer, 1015, fp);
        words = chop_spaced_words(buffer);
        if (!words) continue;
        for (i=0; words[i]; i++) if (words[i][0] == '#')
        {
            words[i] = nullptr;
            break;
        }

        if (words[0]) interpret_cfg_param(words);

        delete words;
    }
}

void iteration_callback(int iter)
{
    cout << iter << " " << flush;

    int tmrno = (iter % 7) + 1;
    int sr, er;
    Point pivot, axis(0,0,0);

    char buffer[8];
    sprintf(buffer, "TMR%d", tmrno);
    sr = g_prot1->get_region_start(buffer);
    er = g_prot1->get_region_end(buffer);

    switch (tmrno)
    {
        case 1:        pivot = pivot1->get_CA_location();        break;
        case 2:        pivot = pivot2->get_CA_location();        break;
        case 3:        pivot = pivot3->get_CA_location();        break;
        case 4:        pivot = pivot4->get_CA_location();        break;
        case 5:        pivot = pivot5->get_CA_location();        break;
        case 6:        pivot = pivot6->get_CA_location();        break;
        case 7:        pivot = pivot7->get_CA_location();        break;
    
        default:
        break;
    }

    float e = residue_energy(), e1 = 0, theta;
    e -= g_prot1->get_internal_clashes(sr, er, false)*_kJmol_cuA*protein_clash_penalty;

    int i;

    for (i=0; i<3; i++)
    {
        axis = Point( i==0 ? 1000 : 0, i==1 ? 1000 : 0, i==2 ? 1000 : 0 );

        theta = frand(-montecarlo_theta, montecarlo_theta);
        g_prot1->rotate_piece(sr, er, pivot, axis, theta);
        e1 = residue_energy();
        e1 -= g_prot1->get_internal_clashes(sr, er, false)*_kJmol_cuA*protein_clash_penalty;
        if (e1 >= e)
        {
            e = e1;
        }
        else
        {
            g_prot1->rotate_piece(sr, er, pivot, axis, -theta);
        }

        // Do not move salt bridge acid.
        if (sba)
        {
            int resno = sba->get_residue_no();
            if (resno >= sr && resno <= er) continue;
        }

        // Do not move salt bridge base.
        if (sbb)
        {
            int resno = sbb->get_residue_no();
            if (resno >= sr && resno <= er) continue;
        }

        axis.scale(frand(-montecarlo_xform, montecarlo_xform));
        g_prot1->move_piece(sr, er, (SCoord)axis);
        e1 = residue_energy();
        e1 -= g_prot1->get_internal_clashes(sr, er, false)*_kJmol_cuA*protein_clash_penalty;
        if (e1 >= e)
        {
            e = e1;
        }
        else
        {
            axis.scale(-axis.magnitude());
            g_prot1->move_piece(sr, er, (SCoord)axis);
        }
    }

    float c = g_prot1->get_internal_clashes(sr, er, true, 5);

    // for (i=1; i<tmrno; i++) cout << "     ";
    // cout << c << endl;

}

void do_template_homology(const char* template_path)
{
    Protein tpl("Template");
    fp = fopen(template_path, "rb");
    if (!fp)
    {
        cout << "File not found: " << template_path;
        return -1;
    }
    else
    {
        tpl.load_pdb(fp);
        fclose(fp);

        g_prot1->homology_conform(&tpl);
    }
}

int main(int argc, char** argv)
{
    if (argc < 3)
    {
        show_usage();
        return -1;
    }

    Protein p1("Prot1");
    Protein p2("Prot2");
    FILE* fp;

    g_prot1 = &p1;
    g_prot2 = &p2;

    fp = fopen(argv[1], "rb");
    if (!fp)
    {
        cout << "Please ensure " << argv[1] << " exists and is readable." << endl;
        return -1;
    }
    read_cfg_lines(fp);
    fclose(fp);

    // TODO: Command line args.
    

    int i, n;
    p1.load_pdb(prot1fname);
    p2.load_pdb(prot2fname);

    n = contacts.size();
    for (i=0; i<n; i++)
    {
        contacts[i].prot1 = &p1;
        contacts[i].prot2 = &p2;
        contacts[i].interpret_cfgs();
    }

    n = segments.size();
    for (i=0; i<n; i++)
    {
        segments[i].prot = &p1;
        segments[i].interpret_cfgs();
    }


    // TODO: p1 homology.


    // Now rotate p2 to align as many contacts as possible.

    // Next, move p2 away from p1 (along barycenter-to-barycenter axis) until no clashes.

    // Finally, iteratively wiggle the segments around to get optimal contacts and minimal clashes.
    // Include small transformations and rotations of p2 to search the conformational space.


    // Write the output file.
    fp = fopen(output_fname, "wb");
    if (!fp)
    {
        cout << "Failed to open output file." << endl;
        return -1;
    }

    gpcr.set_pdb_chain('A');
    gpcr.save_pdb(fp);
    gnax.set_pdb_chain('B');
    gnax.save_pdb(fp);
    gnax.end_pdb(fp);
    fclose(fp);
    cout << "Wrote output file." << endl;
}
