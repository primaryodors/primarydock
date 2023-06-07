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

#define _dbg_contacts 0
#define _dbg_segments 0
#define xyz_step 0.1

using namespace std;

Protein *g_prot1, *g_prot2;         // Global protein pointers; not necessarily G-proteins.

class Contact
{
    public:
    Protein *prot1 = nullptr, *prot2 = nullptr;
    AminoAcid *aa1 = nullptr, *aa2 = nullptr;
    std::string cfgstr1, cfgstr2;

    void interpret_cfgs()
    {
        int n, i, j;
        char c[256];
        Protein* which;
        for (n=1; n<=2; n++)
        {
            if (n==1) strcpy(c, cfgstr1.c_str());
            if (n==2) strcpy(c, cfgstr2.c_str());

            #if _dbg_contacts
            cout << "Interpreting contact " << c << endl;
            #endif

            i = 0;
            if (c[1] == ':')
            {
                if (c[0] < '1' || c[0] > '2') throw 0xbadcf6;
                else
                {
                    if (n==1) prot1 = (c[0]=='1') ? g_prot1 : g_prot2;
                    if (n==2) prot2 = (c[0]=='1') ? g_prot1 : g_prot2;
                    i += 2;
                }
            }
            else
            {
                prot1 = g_prot1;
                prot2 = g_prot2;
            }

            #if _dbg_contacts
            cout << "Proteins are " << (prot1 ? prot1->get_name() : "null") << ", " << (prot2 ? prot2->get_name() : "null") << endl;
            #endif

            which = (n==1) ? prot1 : prot2;

            #if _dbg_contacts
            cout << "Processing contact for " << which->get_name() << endl;
            #endif

            int mode = 0;
            std::vector<char> allowed_aa;
            char* has_dot = nullptr;
            int resno=0, tol=0;

            for (; c[i]; i++)
            {
                #if _dbg_contacts
                cout << "c[" << i << "] = " << c[i] << endl;
                #endif

                switch (mode)
                {
                    case 0:     // Getting amino acids.
                    if (c[i] >= '0' && c[i] <= '9')
                    {
                        mode = 1;
                        i--;
                    }
                    else
                    {
                        allowed_aa.push_back(c[i]);

                        #if _dbg_contacts
                        cout << "Allowed amino acid: " << c[i] << endl;
                        #endif
                    }
                    break;

                    case 1:     // Getting residue number.
                    for (j=i; (c[j] >= '0' && c[j] <= '9') || c[j] == '.'; j++)
                    {
                        if (c[j] == '.') has_dot = &c[j];
                    }

                    if (c[j] == '~') mode = 2;
                    c[j] = 0;
                    
                    if (has_dot)
                    {
                        #if _dbg_contacts
                        cout << "Getting BW number " << &c[i] << "... ";
                        #endif

                        *has_dot = 0;
                        resno = which->get_bw50(atoi(&c[i]));
                        *has_dot = '.';
                        resno += atoi(&has_dot[1]) - 50;

                        #if _dbg_contacts
                        cout << "residue is " << resno << endl;
                        #endif
                    }
                    else
                    {
                        if (c[i] >= '0' && c[i] <= '9') resno = atoi(&c[i]);

                        #if _dbg_contacts
                        cout << "Residue " << &c[i] << " is " << resno << endl;
                        #endif
                    }

                    if (mode==2) c[j] = '~';
                    i = j-1;
                    break;

                    case 2:     // Getting tolerance.
                    tol = atoi(&c[i]);

                    #if _dbg_contacts
                    cout << "Tolerance is " << tol << endl;
                    #endif
                }
            }

            if (!resno) throw 0xbadcf6;

            bool found = false;
            for (i = 0; i <= tol; i++)
            {
                AminoAcid *aa = which->get_residue(resno+i);
                if (aa && std::find(allowed_aa.begin(), allowed_aa.end(), aa->get_letter()) != allowed_aa.end())
                {
                    if (n==1) aa1 = aa;
                    if (n==2) aa2 = aa;
                    found = true;

                    #if _dbg_contacts
                    cout << "Found allowed " << aa->get_letter() << " at position " << (resno+i) << endl;
                    #endif

                    break;
                }
                else
                {
                    aa = which->get_residue(resno-i);
                    if (aa && std::find(allowed_aa.begin(), allowed_aa.end(), aa->get_letter()) != allowed_aa.end())
                    {
                        if (n==1) aa1 = aa;
                        if (n==2) aa2 = aa;
                        found = true;

                        #if _dbg_contacts
                        cout << "Found allowed " << aa->get_letter() << " at position " << (resno-i) << endl;
                        #endif

                        break;
                    }
                }
            }

            #if _dbg_contacts
            cout << endl;
            #endif

            if (!found)
            {            // If fail condition, blank both AA pointers and return. Contact cannot be made and will be skipped.
                aa1 = aa2 = nullptr;
                return;
            }
        }           // for n.
    }           // interpret cfgs.
};

class MovablePiece
{
    public:
    Protein* prot = nullptr;
    AminoAcid *start_residue = nullptr, *end_residue = nullptr;
    AminoAcid *first_pivot = nullptr, *last_pivot = nullptr;
    std::vector<std::string> cfgstrs;

    void interpret_cfgs()
    {
        prot = g_prot1;

        int n = cfgstrs.size();
        if (n < 3) throw 0xbadcf6;

        int i;
        for (i=0; i<4; i++)
        {
            if (i >= n) continue;

            #if _dbg_segments
            cout << "Interpreting segment string " << cfgstrs[i] << endl;
            #endif

            int resno;
            const char *c = cfgstrs[i].c_str(), *d;

            if (!strcmp(c, "end") || !strcmp(c, "END"))
            {
                resno = prot->get_end_resno();

                #if _dbg_segments
                cout << "Using last residue " << resno << endl;
                #endif
            }
            else if (d = strchr(c, '.'))             // Assignment, not comparison.
            {
                resno = prot->get_bw50(atoi(c)) + atoi(d+1) - 50;

                #if _dbg_segments
                cout << "Using BW residue " << resno << endl;
                #endif
            }
            else
            {
                resno = atoi(c);

                #if _dbg_segments
                cout << "Using literal residue " << resno << endl;
                #endif
            }

            if (!i) start_residue = prot->get_residue(resno);
            else if (i==1) end_residue = prot->get_residue(resno);
            else if (i==2) first_pivot = prot->get_residue(resno);
            else if (i==3) last_pivot = prot->get_residue(resno);

            #if _dbg_segments
            cout << endl;
            #endif
        }
    }

    void do_motion(SCoord move_amt)
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
            nouion = nouion.multiply_3d_distance(&pivf, old.get_3d_distance(pivf) / nouion.get_3d_distance(pivf) );
            rel = nouion.subtract(old);
        }
        if (last_pivot)
        {
            Point pivf = last_pivot->get_CA_location();
            Point old = end_residue->get_CA_location();
            Point nouion = old.add(rel);
            nouion = nouion.multiply_3d_distance(&pivf, old.get_3d_distance(pivf) / nouion.get_3d_distance(pivf) );
            rel = nouion.subtract(old);
        }

        // Do the motion.
        prot->move_piece(start_residue->get_residue_no(), end_residue->get_residue_no(), (SCoord)rel);
        Point srca, eca;

        if (start_residue->get_residue_no() > 1)
        {
            srca = prot->get_residue(start_residue->get_residue_no()-1)->get_CA_location().add(rel);

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
        }

        if (end_residue->get_residue_no() < prot->get_end_resno())
        {
            eca  = prot->get_residue(end_residue->get_residue_no()+1)->get_CA_location().add(rel);

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
        }

        // TODO: Realign segment to new locations of start residue - 1 and end residue + 1.

    }

    void do_rotation(SCoord axis, float theta)
    {
        // Check the direction of rotation and correct if necessary.

        // TODO:

    }
};

Molecule** g_contacts_as_mols;
AminoAcid *sbb = nullptr, *sba = nullptr;
int iters = 50;

std::string prot1fname, prot2fname, tplname;
char output_fname[256];
std::vector<Contact> contacts;
std::vector<MovablePiece> segments;

const float montecarlo_theta = fiftyseventh * 1;
const float montecarlo_xform = 0.5;
const float protein_clash_penalty = 10;

void show_usage()
{
    cout << "Usage:" << endl;
    cout << "couple path/to/file.cplcfg";
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
        segments.push_back(p);

        return i;
    }
    else if (!strcmp(words[0], "TEMPLATE"))
    {
        tplname = words[1];
        return 2;
    }
    else if (!strcmp(words[0], "ITER"))
    {
        iters = atoi(words[1]);
        return 2;
    }
    else if (!strcmp(words[0], "OUT"))
    {
        strcpy(output_fname, words[1]);
        return 2;
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

void do_template_homology(const char* template_path)
{
    Protein tpl("Template");
    FILE* fp = fopen(template_path, "rb");
    if (!fp)
    {
        cout << "File not found: " << template_path;
        throw 0xbadcf6;
    }
    else
    {
        tpl.load_pdb(fp);
        fclose(fp);

        g_prot1->homology_conform(&tpl);
    }
}

float total_contact_binding()
{
    int i, n;
    float f;

    n = contacts.size();
    for (i=0; i<n; i++)
    {
        f += contacts[i].aa1->get_intermol_binding(contacts[i].aa2);
    }
    return f;
}

void optimize_contacts(int iters = 20)
{
    int i, n;

    n = contacts.size();
    for (i=0; i<n; i++)
    {
        Molecule* cfmols[4];
        cfmols[0] = (Molecule*)contacts[i].aa1;
        cfmols[1] = (Molecule*)contacts[i].aa2;
        cfmols[2] = nullptr;

        Molecule::conform_molecules(cfmols, iters);
    }
}

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        show_usage();
        return -1;
    }

    strcpy(output_fname, "output/coupled.pdb");

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
    

    int i, j, l, m, n;
    fp = fopen(prot1fname.c_str(), "rb");
    if (!fp) throw 0xbadcf6;
    cout << "Reading protein 1..." << endl;
    p1.load_pdb(fp);
    fclose(fp);
    fp = fopen(prot2fname.c_str(), "rb");
    if (!fp) throw 0xbadcf6;
    cout << "Reading protein 2..." << endl;
    p2.load_pdb(fp);
    fclose(fp);

    cout << "Reading parameters..." << endl;
    n = contacts.size();
    Star swap;
    for (i=0; i<n; i++)
    {
        contacts[i].prot1 = &p1;
        contacts[i].prot2 = &p2;
        contacts[i].interpret_cfgs();

        if (contacts[i].prot1 == g_prot2 && contacts[i].prot2 == g_prot1)
        {
            swap.pprot = contacts[i].prot2;
            contacts[i].prot2 = contacts[i].prot1;
            contacts[i].prot1 = swap.pprot;

            swap.paa = contacts[i].aa2;
            contacts[i].aa2 = contacts[i].aa1;
            contacts[i].aa1 = swap.paa;
        }
    }

    for (i=n-1; i>=0; i--) if (!contacts[i].aa1 || !contacts[i].aa2)
    {
        contacts.erase(contacts.begin()+i);
    }

    n = segments.size();
    for (i=0; i<n; i++)
    {
        segments[i].prot = &p1;
        segments[i].interpret_cfgs();
    }

    for (i=n-1; i>=0; i--) if (!segments[i].start_residue || !segments[i].end_residue || (!segments[i].first_pivot && !segments[i].last_pivot))
    {
        segments.erase(segments.begin()+i);
    }


    // TODO: p1 homology.
    if (tplname)
    {
        do_template_homology(tplname.c_str());
    }


    // Now rotate p2 to align as many contacts as possible.
    n = contacts.size();
    j = 0;
    l = -1;
    for (i=0; i<n; i++)
    {
        cout << "Found contact residues " << contacts[i].prot1->get_name() << ":" << *contacts[i].aa1
             << " and " << contacts[i].prot2->get_name() << ":" << *contacts[i].aa2
             << endl;
        if (contacts[i].prot1 != contacts[i].prot2)
        {
            j++;
            if (l < 0) l = i;
        }
    }

    if (j < 3)
    {
        cout << "Not enough contacts were found between the proteins. At least 3 contacts are required for coupling." << endl;
        return -1;
    }

    cout << "Performing rough alignment..." << endl;
    Point rel = contacts[l].aa1->get_CA_location().subtract(contacts[l].aa2->get_CA_location());
    int resno1 = contacts[l].aa2->get_residue_no();
    p2.move_piece(1, p2.get_end_resno(), (SCoord)rel);

    l++;
    while (contacts[l].prot1 == contacts[l].prot2)
    {
        l++;
        if (l >= n) throw 0xbadc0de;
    }

    Point ref = contacts[l].aa1->get_CA_location();
    int resno2 = contacts[l].aa2->get_residue_no();
    p2.rotate_piece(1, p2.get_end_resno(), resno2, ref, resno1);

    l++;
    while (contacts[l].prot1 == contacts[l].prot2)
    {
        l++;
        if (l >= n) throw 0xbadc0de;
    }

    ref = contacts[l].aa1->get_CA_location();
    int resno3 = contacts[l].aa2->get_residue_no();
    p2.rotate_piece(1, p2.get_end_resno(), resno3, ref, resno2);

    cout << "Fine-tuning alignment...";
    for (i=0; i<20; i++)
    {
        for (l=0; l<n; l++)
        {
            if (contacts[l].prot1 == contacts[l].prot2) continue;

            Point pt = contacts[l].aa2->get_CA_location(),
                  algn = contacts[l].aa1->get_CA_location(),
                  cen = contacts[l].prot2->get_region_center(1, contacts[l].prot2->get_end_resno());
            
            Rotation rot = align_points_3d(pt, algn, cen);
            rot.a /= 3;

            p2.rotate_piece(1, p2.get_end_resno(), rot, 0);

            rel = contacts[l].aa1->get_CA_location().subtract(contacts[l].aa2->get_CA_location());
            rel.scale(rel.magnitude()/3);

            p2.move_piece(1, 99999, (SCoord)rel);
        }

        cout << "." << flush;
    }
    cout << endl;

    // Next, iteratively wiggle the segments around to get optimal contacts and minimal clashes.
    // Include small transformations and rotations of p2 to search the conformational space.
    Point pcen = p1.get_region_center(1, p1.get_end_resno());
    std::vector<AminoAcid*> cr = p1.get_contact_residues(&p2);
    n = cr.size();
    g_contacts_as_mols = new Molecule*[n+4];

    for (i=0; i<n; i++) g_contacts_as_mols[i] = (Molecule*)cr[i];
    g_contacts_as_mols[i] = nullptr;

    // Test.
    #if 0
    ref = segments[0].prot->get_region_center(segments[0].start_residue->get_residue_no(), segments[0].end_residue->get_residue_no());
    ref.y = pcen.y;
    rel = ref.multiply_3d_distance(&pcen, 2).subtract(ref);
    segments[0].do_motion(rel);

    ref = segments[1].prot->get_region_center(segments[1].start_residue->get_residue_no(), segments[1].end_residue->get_residue_no());
    ref.y = pcen.y;
    rel = ref.multiply_3d_distance(&pcen, 2).subtract(ref);
    segments[1].do_motion(rel);

    #endif

    cout << "Reshaping...";
    n = segments.size();
    for (i=0; i<iters; i++)
    {
        Molecule::conform_molecules(g_contacts_as_mols, 10);
        optimize_contacts(50);

        float e, f = Molecule::total_intermol_binding(g_contacts_as_mols) + total_contact_binding();

        for (j=0; j<n; j++)
        {
            // Point seg = segments[j].prot->get_region_center(segments[j].start_residue->get_residue_no(), segments[j].end_residue->get_residue_no());
            rel = Point( frand(-xyz_step, xyz_step), frand(-xyz_step, xyz_step), frand(-xyz_step, xyz_step) );
            // seg = seg.add(rel);
            segments[j].do_motion(rel);
            optimize_contacts();
            e = Molecule::total_intermol_binding(g_contacts_as_mols) + total_contact_binding();
            if (e < f)
            {
                rel.scale(-rel.magnitude());
                segments[j].do_motion(rel);
                optimize_contacts();
                cout << "-";
            }
            else
            {
                f = e;
                cout << "+";
            }
        }

        rel = Point( frand(-xyz_step, xyz_step), frand(-xyz_step, xyz_step), frand(-xyz_step, xyz_step) );
        p2.move_piece(1, p2.get_end_resno(), (SCoord)rel);
        optimize_contacts();
        e = Molecule::total_intermol_binding(g_contacts_as_mols) + total_contact_binding();
        if (e < f)
        {
            rel.scale(-rel.magnitude());
            p2.move_piece(1, p2.get_end_resno(), (SCoord)rel);
            optimize_contacts();
        }
        else
        {
            f = e;
        }

        cout << "." << flush;
    }
    cout << endl;

    // Write the output file.
    fp = fopen(output_fname, "wb");
    if (!fp)
    {
        cout << "Failed to open output file." << endl;
        return -1;
    }

    optimize_contacts(50);
    p1.set_pdb_chain('A');
    p1.save_pdb(fp);
    p2.renumber_residues(1, p2.get_end_resno(), 1001);
    p2.set_pdb_chain('B');
    p2.save_pdb(fp);
    p2.end_pdb(fp);
    fclose(fp);
    cout << "Wrote output file." << endl;
}
