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

#define xyz_step 0.1
#define xyz_big_step 1
#define contact_importance 100

#define _dbg_contacts 0
#define _dbg_segments 0

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
            char *has_dot = nullptr, *plusminus, *paren, pm;
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
                    else if (c[i] == '(')
                    {
                        plusminus = strchr(&c[i], '+');
                        if (!plusminus) plusminus = strchr(&c[i], '-');
                        if (!plusminus) throw 0xbadf37;
                        paren = strchr(plusminus, ')');
                        if (!paren) throw 0xbadf37;
                        pm = *plusminus;
                        *plusminus = 0;
                        *paren = 0;

                        #if _dbg_contacts
                        cout << "Searching protein for motif " << &c[i+1] << endl;
                        #endif
                        
                        j = which->search_sequence(1, which->get_end_resno(), &c[i+1]);
                        if (!j) throw 0xbadcf6;

                        #if _dbg_contacts
                        cout << "Found motif at position " << j << endl;
                        #endif

                        if (pm == '+') resno = j + atoi(plusminus+1);
                        if (pm == '-') resno = j - atoi(plusminus+1);

                        #if _dbg_contacts
                        cout << "Resno is " << resno << endl << flush;
                        #endif

                        *plusminus = pm;
                        *paren = ')';
                        i = paren - c;
                        mode = 1;
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

    SCoord vector_to_contact_horizon(int protno = 2)
    {
        AminoAcid *mov, *stat;

        if (protno == 2)
        {
            mov = aa2;
            stat = aa1;
        }
        else if (protno == 1)
        {
            mov = aa1;
            stat = aa2;
        }
        else throw 0xbadc0de;

        SCoord v = stat->get_CA_location().subtract(mov->get_CA_location());
        v.r = fmax(0, v.r - mov->get_reach() - stat->get_reach() - 2);

        return v;
    }
};

std::ostream& operator<<(std::ostream& os, const Contact& c)
{
    if (!&c) return os;

    os << c.prot1->get_name() << ":" << *c.aa1 << " ... " << c.prot2->get_name() << ":" << *c.aa2;

    return os;
}

class MovablePiece
{
    public:
    Protein* prot = nullptr;
    AminoAcid *start_residue = nullptr, *end_residue = nullptr;
    AminoAcid *first_pivot = nullptr, *last_pivot = nullptr;
    std::vector<std::string> cfgstrs;
    Pose** undo_poses = nullptr;

    ~MovablePiece()
    {
        if (undo_poses) delete[] undo_poses;
    }

    void interpret_cfgs()
    {
        int i, n;

        prot = g_prot1;
        n = prot->get_end_resno();
        undo_poses = new Pose*[n+4];
        for (i=0; i<=n+4; i++) undo_poses[i] = nullptr;

        n = cfgstrs.size();
        if (n < 3) throw 0xbadcf6;

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

    void save_undo_poses()
    {
        int i, sr, er;

        sr = first_pivot ? first_pivot->get_residue_no() : start_residue->get_residue_no();
        er = last_pivot ? last_pivot->get_residue_no() : end_residue->get_residue_no();

        for (i=sr; i<=er; i++)
        {
            AminoAcid* aa = prot->get_residue(i);
            if (aa)
            {
                if (!undo_poses[i]) undo_poses[i] = new Pose((Molecule*)aa);
                else undo_poses[i]->copy_state((Molecule*)aa);
            }
        }
    }

    void undo()
    {
        int i, sr, er;

        sr = first_pivot ? first_pivot->get_residue_no() : start_residue->get_residue_no();
        er = last_pivot ? last_pivot->get_residue_no() : end_residue->get_residue_no();

        for (i=sr; i<=er; i++)
        {
            AminoAcid* aa = prot->get_residue(i);
            if (aa)
            {
                if (undo_poses[i]) undo_poses[i]->restore_state((Molecule*)aa);
            }
        }
    }

    void do_motion(SCoord move_amt)
    {
        if (!prot) throw 0xbadc0de;
        if (!start_residue || !end_residue) throw 0xbadc0de;
        if (!first_pivot && !last_pivot) throw 0xbadc0de;

        save_undo_poses();

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
        save_undo_poses();

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
    float f = 0;

    n = contacts.size();
    for (i=0; i<n; i++)
    {
        f += contacts[i].aa1->get_intermol_binding(contacts[i].aa2);
    }
    return f;
}

void optimize_contacts(int iters = 30)
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
    if (tplname.length())
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

    Point rel, ref, avg1, avg2;

    cout << "Performing rough alignment..." << endl;

    n = contacts.size();
    for (i=0; i<n; i++)
    {
        avg1 = avg1.add(contacts[i].aa1->get_CA_location());
    }
    avg1.scale(avg1.magnitude() / n);

    rel = avg1.subtract(p1.get_region_center(1, p1.get_end_resno()));
    rel.scale(50);
    p2.move_piece(1, p2.get_end_resno(), rel);

    for (i=0; i<n; i++)
    {
        avg2 = avg2.add(contacts[i].aa2->get_CA_location());
    }
    avg2.scale(avg2.magnitude() / n);

    Rotation rot = align_points_3d(avg2, avg1, p2.get_region_center(1, p2.get_end_resno()));
    p2.rotate_piece(1, p2.get_end_resno(), rot, 0);

    


    cout << "Finding closest and farthest pairs..." << endl;
    float rbest = 0, rworst = 0;
    Point pcen = p1.get_region_center(1, p1.get_end_resno());
    n = contacts.size();
    for (i=0; i<n; i++)
    {
        float r = contacts[i].aa1->get_CA_location().get_3d_distance(pcen);
        if (!i || r < rbest)
        {
            rbest = r;
            l = i;
        }
        else if (r > rworst)
        {
            rworst = r;
            j = i;
        }
    }

    for (m = 0; m < iters; m++)
    {
        cout << "Pullapart..." << endl;
        rel = p2.get_region_center(1, p2.get_end_resno()).subtract(p1.get_region_center(1, p1.get_end_resno()));
        rel.scale(25);
        // _INTERA_R_CUTOFF = 10;
        p2.move_piece(1, p2.get_end_resno(), (SCoord)rel);

        cout << "Lining up closest pair " << contacts[l] << "..." << endl;
        rel = contacts[l].vector_to_contact_horizon(2);
        p2.move_piece(1, p2.get_end_resno(), (SCoord)rel);

        cout << "Lining up farthest pair " << contacts[j] << "..." << endl;
        SCoord axis = contacts[l].aa2->get_CA_location().subtract(contacts[l].aa1->get_CA_location());
        float theta = find_angle_along_vector(contacts[j].aa2->get_CA_location(), contacts[j].aa1->get_CA_location(), contacts[j].aa1->get_CA_location(), axis);
        rot.v = axis;
        rot.a = theta;
        p2.rotate_piece(1, p2.get_end_resno(), rot, contacts[l].aa2->get_residue_no());


        cout << "Performing global average motion..." << endl;
        avg2 = Point(0,0,0);
        for (i=0; i<n; i++)
        {
            avg2 = avg2.add(contacts[i].vector_to_contact_horizon());
        }
        avg2.scale(avg2.magnitude() / n);

        p2.move_piece(1, p2.get_end_resno(), (SCoord)avg2);


        cout << "Performing global average rotation..." << endl;
        SCoord axisx = Point(1,0,0), axisz = Point(0,0,1);
        float xtheta = 0, ztheta = 0;
        pcen = p2.get_region_center(1, p2.get_end_resno());

        n = contacts.size();
        for (i=0; i<n; i++)
        {
            ref = contacts[i].aa2->get_CA_location();
            rel = contacts[i].vector_to_contact_horizon(2);

            xtheta += find_angle_along_vector(ref, ref.add(rel), pcen, axisx);
            ztheta += find_angle_along_vector(ref, ref.add(rel), pcen, axisz);
        }

        xtheta /= n;
        ztheta /= n;

        cout << "Rotating " << (xtheta*fiftyseven) << "deg about X axis..." << endl;
        rot.v = axisx;
        rot.a = xtheta;
        p2.rotate_piece(1, p2.get_end_resno(), rot, 0);

        cout << "Rotating " << (ztheta*fiftyseven) << "deg about Z axis..." << endl;
        rot.v = axisz;
        rot.a = ztheta;
        p2.rotate_piece(1, p2.get_end_resno(), rot, 0);
    }


    n = contacts.size();
    m = 0;
    for (i=0; i<n; i++)
    {
        ref = contacts[i].vector_to_contact_horizon(2);
        cout << "Vector to contact horizon for " << contacts[i] << ": " << ref << endl;
        if (!ref.magnitude()) m++;
    }
    cout << m << " successful contacts out of " << n << " total." << endl;
    
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
