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
#define contact_importance 200

#define _dbg_contacts 0
#define _dbg_segments 0
#define _dbg_iters 0
#define _dbg_flexion 0
#define _dbg_contact_rel 0

using namespace std;

Protein *g_prot1, *g_prot2;         // Global protein pointers; not necessarily G-proteins.

int interpret_resno(Protein* prot, char* str)
{
    int i, j, mode = 0;
    std::vector<char> allowed_aa;
    char *has_dot = nullptr, *plusminus, *paren, pm;
    int resno=0, tol=0;

    for (i=0; str[i]; i++)
    {
        #if _dbg_contacts
        cout << "str[" << i << "] = " << str[i] << endl;
        #endif

        switch (mode)
        {
            case 0:     // Getting amino acids.
            if (str[i] >= '0' && str[i] <= '9')
            {
                mode = 1;
                i--;
            }
            else if (str[i] == '(')
            {
                plusminus = strchr(&str[i], '+');
                if (!plusminus) plusminus = strchr(&str[i], '-');
                if (!plusminus) throw 0xbadf37;
                paren = strchr(plusminus, ')');
                if (!paren) throw 0xbadf37;
                pm = *plusminus;
                *plusminus = 0;
                *paren = 0;

                #if _dbg_contacts
                cout << "Searching protein for motif " << &str[i+1] << endl;
                #endif
                
                j = prot->search_sequence(1, prot->get_end_resno(), &str[i+1]);
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
                i = paren - str;
                mode = 1;
            }
            else
            {
                allowed_aa.push_back(str[i]);

                #if _dbg_contacts
                cout << "Allowed amino acid: " << str[i] << endl;
                #endif
            }
            break;

            case 1:     // Getting residue number.
            for (j=i; (str[j] >= '0' && str[j] <= '9') || str[j] == '.'; j++)
            {
                if (str[j] == '.') has_dot = &str[j];
            }

            if (str[j] == '~') mode = 2;
            str[j] = 0;
            
            if (has_dot)
            {
                #if _dbg_contacts
                cout << "Getting BW number " << &str[i] << "... ";
                #endif

                *has_dot = 0;
                resno = prot->get_bw50(atoi(&str[i]));
                *has_dot = '.';
                resno += atoi(&has_dot[1]) - 50;

                #if _dbg_contacts
                cout << "residue is " << resno << endl;
                #endif
            }
            else
            {
                if (str[i] >= '0' && str[i] <= '9') resno = atoi(&str[i]);

                #if _dbg_contacts
                cout << "Residue " << &str[i] << " is " << resno << endl;
                #endif
            }

            if (mode==2) str[j] = '~';
            i = j-1;
            break;

            case 2:     // Getting tolerance.
            tol = atoi(&str[i]);

            #if _dbg_contacts
            cout << "Tolerance is " << tol << endl;
            #endif
        }
    }

    if (!resno) throw 0xbadcf6;
    if (!allowed_aa.size()) return resno;

    for (i = 0; i <= tol; i++)
    {
        AminoAcid *aa = prot->get_residue(resno+i);
        if (!aa) continue;
        if (allowed_aa[0] == 'X')
        {
            return resno+i;

            #if _dbg_contacts
            cout << "Found allowed " << aa->get_letter() << " at position " << (resno+i) << endl;
            #endif

            break;
        }
        else if (aa && std::find(allowed_aa.begin(), allowed_aa.end(), aa->get_letter()) != allowed_aa.end())
        {
            return resno+i;

            #if _dbg_contacts
            cout << "Found allowed " << aa->get_letter() << " at position " << (resno+i) << endl;
            #endif

            break;
        }
        else
        {
            aa = prot->get_residue(resno-i);
            if (aa && std::find(allowed_aa.begin(), allowed_aa.end(), aa->get_letter()) != allowed_aa.end())
            {
                return resno+i;

                #if _dbg_contacts
                cout << "Found allowed " << aa->get_letter() << " at position " << (resno-i) << endl;
                #endif

                break;
            }
        }
    }

    return 0;
}

class Contact
{
    public:
    Protein *prot1 = nullptr, *prot2 = nullptr;
    AminoAcid *aa1 = nullptr, *aa2 = nullptr;
    std::string cfgstr1, cfgstr2;
    
    float weight()
    {
        if (aa1->get_charge() && sgn(aa1->get_charge()) == -aa2->get_charge()) return 6.0;
        if (fabs(aa1->hydrophilicity()) > 0.25 && fabs(aa2->hydrophilicity()) > 0.25) return 2.0;
        if (aa1->get_aa_definition()->aromatic && aa2->get_aa_definition()->aromatic) return .5;
        return .1;
    }

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

            int resno = interpret_resno(which, &c[i]);
            AminoAcid *aa = which->get_residue(resno);

            if (!resno || !aa)
            {            // If fail condition, blank both AA pointers and return. Contact cannot be made and will be skipped.
                aa1 = aa2 = nullptr;
                return;
            }

            if (n==1) aa1 = aa;
            if (n==2) aa2 = aa;
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
        v.r = fmax(0, v.r - mov->get_reach() - stat->get_reach() - 0);

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
    AminoAcid *prev_start = nullptr, *next_end = nullptr;
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
        for (i=0; i<=n; i++) undo_poses[i] = nullptr;

        n = cfgstrs.size();
        if (n < 3) throw 0xbadcf6;

        for (i=0; i<6; i++)
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
            else if (i==2) prev_start = first_pivot = prot->get_residue(resno);
            else if (i==3) next_end = last_pivot = prot->get_residue(resno);
            else if (i==4) first_pivot = prot->get_residue(resno);
            else if (i==5) last_pivot = prot->get_residue(resno);

            #if _dbg_segments
            cout << endl;
            #endif
        }
    }

    void save_undo_poses()
    {
        int i, sr, er;

        sr = prev_start ? prev_start->get_residue_no() : start_residue->get_residue_no();
        er = next_end ? next_end->get_residue_no() : end_residue->get_residue_no();

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

        sr = prev_start ? prev_start->get_residue_no() : start_residue->get_residue_no();
        er = next_end ? next_end->get_residue_no() : end_residue->get_residue_no();

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
        if (!prev_start && !next_end) throw 0xbadc0de;

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

            // If prev start, rotate prev start thru start residue - 1 about first pivot to align with start residue.
            if (prev_start)
            {
                prot->rotate_piece(prev_start->get_residue_no(), start_residue->get_residue_no()-1, start_residue->get_residue_no()-1,
                    srca, first_pivot->get_residue_no());
            }

            // If no prev start, rotate head through start residue - 1 about last pivot to align with start residue.
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
            if (next_end)
            {
                prot->rotate_piece(end_residue->get_residue_no()+1, next_end->get_residue_no(), end_residue->get_residue_no()+1,
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

class MovableContact
{
    public:
    Protein *prot = nullptr;
    AminoAcid *aa1 = nullptr, *aa2 = nullptr, *segstart = nullptr, *segend = nullptr, *pivot = nullptr;
    std::string cfgstr;

    void interpret_config_str()
    {
        if (!prot) throw 0xbadc0de;

        int i, n;
        char cfg[256];
        strcpy(cfg, cfgstr.c_str());

        char** fields = chop_spaced_words(cfg);

        if (!fields[0] || strcmp(fields[0], "MAKESURE")) throw 0xbadcf6;

        if (!fields[1]) throw 0xbadcf6;
        aa1 = prot->get_residue(interpret_resno(prot, fields[1]));
        if (!aa1) return;
        
        if (!fields[2]) throw 0xbadcf6;
        aa2 = prot->get_residue(interpret_resno(prot, fields[2]));
        if (!aa2) return;
        
        if (!fields[3]) throw 0xbadcf6;
        i = interpret_resno(prot, fields[3]);
        if (!i) throw 0xbadcf6;
        segstart = prot->get_residue(i);
        if (!segstart)
        {
            n = prot->get_end_resno();
            for (; !segstart; i++)
            {
                if (i > n) return;
                segstart = prot->get_residue(i);
            }
        }
        
        if (!fields[4]) throw 0xbadcf6;
        i = interpret_resno(prot, fields[4]);
        if (!i) throw 0xbadcf6;
        segend = prot->get_residue(i);
        if (!segend)
        {
            n = prot->get_end_resno();
            for (; !segend; i--)
            {
                if (i < 1) return;
                segend = prot->get_residue(i);
            }
        }

        if (!fields[5])
        {
            i = segstart->get_residue_no() + segend->get_residue_no() / 2;
            pivot = prot->get_residue(i);
            if (!pivot) throw 0xbadcf6;
        }
        else
        {
            i = interpret_resno(prot, fields[5]);
            if (!i) throw 0xbadcf6;
            pivot = prot->get_residue(i);
            if (!pivot) throw 0xbadcf6;
        }
    }

    void prime_contact()
    {
        if (!aa1 || !aa2) return;

        Molecule* mols[4];
        mols[0] = (Molecule*)aa1;
        mols[1] = (Molecule*)aa2;
        mols[2] = nullptr;
        Molecule::conform_molecules(mols);
    }

    void ensure_contact()
    {
        if (!aa1 || !aa2 || !segstart || !segend || !pivot) return;

        prime_contact();

        Point ought, ref;
        Rotation rot;
        float r = aa1->get_CA_location().get_3d_distance(aa2->get_CA_location()), rorig;
        rorig = r;
        r -= aa1->get_reach();
        r -= aa2->get_reach();
        r -= 2;
        if (r > 0)
        {
            ought = aa1->get_CA_location();
            ref = aa2->get_CA_location();
            ought = ought.multiply_3d_distance(&ref, r/rorig);

            // get angle
            rot = align_points_3d(aa1->get_CA_location(), ought, pivot->get_CA_location());

            // rotate
            prot->rotate_piece(segstart->get_residue_no(), segend->get_residue_no(), rot, pivot->get_residue_no());

            // optimize
            prot->get_internal_clashes(segstart->get_residue_no(), segend->get_residue_no(), true, 50);
        }

        prime_contact();
    }
};

Molecule** g_contacts_as_mols;
AminoAcid *sbb = nullptr, *sba = nullptr;
int iters = 50;

std::string prot1fname, prot2fname, tplname, tplrfnm;
char output_fname[256];
std::vector<Contact> contacts;
std::vector<MovablePiece> segments;
std::vector<MovableContact> makesure;

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

float unmet_contacts(bool same_prot_only = false)
{
    int i, n = contacts.size();
    float f = 0;
    
    for (i=0; i<n; i++)
    {
        if (same_prot_only)
        {
            if (contacts[i].prot1 != contacts[i].prot2) continue;
        }
        
        Atom* a = contacts[i].aa1->get_nearest_atom(contacts[i].aa2->get_CA_location());
        if (!a) continue;
        Atom* b = contacts[i].aa2->get_nearest_atom(a->get_location());
        if (!b) continue;
        a = contacts[i].aa1->get_nearest_atom(b->get_location());
        if (!a) continue;
        
        float r = fmax(0, a->distance_to(b) - 1.5);
        
        if (r) f += r;
    }
    
    return f;
}

float contact_binding(std::vector<AminoAcid*> vca)
{
    int i, j, n;
    float f = 0;

    n = vca.size();
    for (i=0; i<n; i++)
    {
        for (j=i+1; j<n; j++)
        {
            f += vca[i]->get_intermol_binding(vca[j]);
        }
    }

    return f;
}

Point contact_center()
{    
    int i, j=0, n = contacts.size();
    Point result(0,0,0);
    
    for (i=0; i<n; i++)
    {
        if (contacts[i].aa1)
        {
            result = result.add(contacts[i].aa1->get_CA_location());
            j++;
        }
        if (contacts[i].aa2)
        {
            result = result.add(contacts[i].aa2->get_CA_location());
            j++;
        }
    }

    if (j)
    {
        result.x /= j;
        result.y /= j;
        result.z /= j;
    }

    return result;
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
    else if (!strcmp(words[0], "MAKESURE"))
    {
        MovableContact mc;
        for (i=0; words[i]; i++)
            mc.cfgstr += (std::string)words[i] + (std::string)" ";
        makesure.push_back(mc);
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
        if (words[2])
        {
            tplrfnm = words[2];
            return 3;
        }
        else return 2;
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

void do_template_homology(const char* template_path, const char* reference_path)
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

        if (reference_path)
        {
            fp = fopen(reference_path, "rb");
            if (!fp)
            {
                cout << "File not found: " << template_path;
                throw 0xbadcf6;
            }
            else
            {
                Protein ref("reference");
                ref.load_pdb(fp);
                fclose(fp);
                g_prot1->homology_conform(&tpl, &ref);
            }
        }

        else g_prot1->homology_conform(&tpl, nullptr);
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
    

    int h, i, j, k, l, m, n;
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


    std::vector<std::string> matches;
    for (i=n-1; i>=0; i--) 
    {
        if (!contacts[i].aa1 || !contacts[i].aa2)
        {
            matches.push_back((std::string)"No match found for " + contacts[i].cfgstr1 + (std::string)" - " + contacts[i].cfgstr2);
            contacts.erase(contacts.begin()+i);
        }
        else
        {
            matches.push_back((std::string)"Found contact residues " + contacts[i].prot1->get_name() + (std::string)":" + (std::string)contacts[i].aa1->get_name()
                + (std::string)" and " + contacts[i].prot2->get_name() + (std::string)":" + (std::string)contacts[i].aa2->get_name()
                + (std::string)" for " + contacts[i].cfgstr1 + (std::string)" - " + contacts[i].cfgstr2
            );
        }
    }

    n = matches.size();
    for (i=n-1; i>=0; i--) cout << matches[i] << endl;

    n = segments.size();
    for (i=0; i<n; i++)
    {
        segments[i].prot = &p1;
        segments[i].interpret_cfgs();
    }

    for (i=n-1; i>=0; i--) if (!segments[i].start_residue || !segments[i].end_residue || (!segments[i].prev_start && !segments[i].next_end))
    {
        segments.erase(segments.begin()+i);
    }

    n = makesure.size();
    for (i=0; i<n; i++)
    {
        makesure[i].prot = &p1;
        makesure[i].interpret_config_str();
        makesure[i].prime_contact();
    }


    // p1 homology.
    if (tplname.length())
    {
        if (tplrfnm.length())
            do_template_homology(tplname.c_str(), tplrfnm.c_str());
        else
            do_template_homology(tplname.c_str(), nullptr);
    }


    // Now rotate p2 to align as many contacts as possible.
    n = contacts.size();
    j = 0;
    l = -1;
    for (i=0; i<n; i++)
    {
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
    SCoord sc;

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

    #if 1


    cout << "Finding closest and farthest pairs..." << endl;
    float rbest = 99999, rworst = 0;
    Point pcen = p1.get_region_center(1, p1.get_end_resno());
    n = contacts.size();
    for (i=0; i<n; i++)
    {
        if (contacts[i].prot1 == contacts[i].prot2) continue;

        float r = contacts[i].aa1->get_CA_location().get_3d_distance(pcen); // / contacts[i].weight();
        if (r < rbest)
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
    
    // cout << "Closest pair is " << contacts[l] << " and farthest is " << contacts[j] << endl;

    cout << "Iterating...";
    int m2 = iters / 2;
    for (m = 0; m < iters; m++)
    {
        std::vector<AminoAcid*> vca = p1.get_contact_residues(&p2);
        
        // cout << vca.size() << endl;
        
        if (!m)
        {
            #if _dbg_iters
            cout << "Pullapart..." << endl;
            #endif

            // rel = p2.get_region_center(1, p2.get_end_resno()).subtract(p1.get_region_center(1, p1.get_end_resno()));
            rel = avg1.subtract(p1.get_region_center(1, p1.get_end_resno()));
            rel.scale(25);
            // _INTERA_R_CUTOFF = 10;
            p2.move_piece(1, p2.get_end_resno(), (SCoord)rel);
        }
        else if (m > m2)
        {
            Molecule** lm = new Molecule*[vca.size()+4];
            for (i=0; i<vca.size(); i++)
            {
                lm[i] = (Molecule*)vca[i];
            }
            lm[i] = nullptr;
            
            float clash = 0;
            for (i=0; lm[i]; i++)
            {
                clash += lm[i]->get_intermol_clashes(lm);
            }

            delete lm;

            if (clash > 1000)
            {
                // rel = p2.get_region_center(1, p2.get_end_resno()).subtract(p1.get_region_center(1, p1.get_end_resno()));
                rel = avg1.subtract(p1.get_region_center(1, p1.get_end_resno()));
                rel.scale(fmax(2, clash / 500));
                p2.move_piece(1, p2.get_end_resno(), (SCoord)rel);

                #if _dbg_iters
                cout << "Realigning closest pair " << contacts[l] << "..." << endl;
                #endif
                rot = align_points_3d(contacts[l].aa2->get_CA_location(), contacts[l].aa1->get_CA_location(), p2.get_region_center(1, p2.get_end_resno()));
                rot.a /= 2;
                p2.rotate_piece(1, p2.get_end_resno(), rot, 0);
            }
        }

        if (!m)
        {
            #if _dbg_iters
            cout << "Lining up closest pair " << contacts[l] << "..." << endl;
            #endif
            rel = contacts[l].vector_to_contact_horizon(2);
            p2.move_piece(1, p2.get_end_resno(), (SCoord)rel);

            #if _dbg_iters
            cout << "Lining up farthest pair " << contacts[j] << "..." << endl;
            #endif
            SCoord axis = contacts[l].aa2->get_CA_location().subtract(contacts[l].aa1->get_CA_location());
            float theta = find_angle_along_vector(contacts[j].aa2->get_CA_location(), contacts[j].aa1->get_CA_location(), contacts[l].aa1->get_CA_location(), axis);
            rot.v = axis;
            rot.a = theta;
            p2.rotate_piece(1, p2.get_end_resno(), rot, contacts[l].aa2->get_residue_no());
        }
        
        // continue;

        #if _dbg_iters
        cout << "Performing global average motion..." << endl;
        #endif

        avg2 = Point(0,0,0);
        float w = 0;
        for (i=0; i<n; i++)
        {
            if (contacts[i].prot1 == contacts[i].prot2) continue;
            rel = contacts[i].vector_to_contact_horizon();
            rel.scale(rel.magnitude() * contacts[i].weight());
            avg2 = avg2.add(rel);
            w += contacts[i].weight();
        }
        avg2.scale(avg2.magnitude() / n);

        p2.move_piece(1, p2.get_end_resno(), (SCoord)avg2);


        #if _dbg_iters
        cout << "Performing global average rotation..." << endl;
        #endif

        SCoord axisx = Point(1,0,0), axisy = Point(0,1,0), axisz = Point(0,0,1);
        float xtheta = 0, ytheta = 0, ztheta = 0;
        float total_rotation = 0;
        pcen = contacts[l].aa2->get_CA_location(); // p2.get_region_center(1, p2.get_end_resno());

        n = contacts.size();
        for (i=0; i<n; i++)
        {
            ref = contacts[i].aa2->get_CA_location();
            rel = contacts[i].vector_to_contact_horizon(2);

            if (!rel.magnitude()) continue;
            rot = align_points_3d(ref, ref.add(rel), pcen);
            if (rot.a > hexagonal) continue;

            w = fmin(1, pcen.get_3d_distance(ref) / 10);
            rot.a *= w;

            rot.a /= 2;
            total_rotation += rot.a;

            #if _dbg_iters
            cout << "Rotating " << (rot.a*fiftyseven) << " deg to line up " << contacts[i] << endl;
            #endif

            p2.rotate_piece(1, p2.get_end_resno(), rot, contacts[l].aa2->get_residue_no());
        }

        if (total_rotation < 0.01*fiftyseventh) break;

        #if !_dbg_iters
        cout << "." << flush;
        #endif
    }
    cout << endl;


    cout << "Flexing segment backbones..." << endl;
    std::vector<AminoAcid*> vca = p1.get_contact_residues(&p2);
    m = vca.size();
    n = segments.size();
    int iter;
    for (iter=0; iter<iters; iter++)
    {
        rel = Point( frand(-2, 2), frand(-2, 2), frand(-2, 2) );

        float was = contact_binding(vca) - unmet_contacts()*contact_importance;
        p2.move_piece(1, p2.get_end_resno(), (SCoord)rel);
        optimize_contacts(20);
        float now = contact_binding(vca) - unmet_contacts()*contact_importance;
        if (now < was)
        {
            p2.undo();
            optimize_contacts(20);
        }

        ref = contacts[0].aa1->get_CA_location(); // contact_center();
        rel = Point( frand(-2, 2), frand(-2, 2), frand(-2, 2) );
        float theta = frand(-2, 2) * fiftyseventh;

        was = contact_binding(vca) - unmet_contacts()*contact_importance;
        p2.rotate_piece(1, p2.get_end_resno(), ref, rel, theta);
        optimize_contacts(20);
        now = contact_binding(vca) - unmet_contacts()*contact_importance;
        if (now < was)
        {
            p2.undo();
            optimize_contacts(20);
        }

        // Anchor first contact.
        sc = contacts[0].aa1->get_CA_location().subtract(contacts[0].aa2->get_CA_location());
        float r = contacts[0].aa1->get_reach() + contacts[0].aa2->get_reach();
        if (r < sc.r)
        {
            sc.r -= r;
            p2.move_piece(1, p2.get_end_resno(), sc);
        }

        for (i=0; i<n; i++)
        {
            int sr = (segments[i].first_pivot ? segments[i].first_pivot : segments[i].start_residue)->get_residue_no(),
                er = (segments[i].last_pivot ? segments[i].last_pivot : segments[i].end_residue)->get_residue_no();
            l = er - sr;

            std::vector<AminoAcid*> vcc = p1.get_residues_can_clash(sr, er);
            vcc.insert(vcc.end(), vca.begin(), vca.end());
            int m1 = vcc.size();

            Point seg_motion(0,0,0);
            Point straight_dir(0,0,0);
            
            straight_dir = p1.get_residue(sr)->get_CA_location().add(p1.get_residue(er)->get_CA_location());
            straight_dir.x /= 2; straight_dir.y = 0; straight_dir.z /= 2;
            rel = segments[i].start_residue->get_CA_location().add(segments[i].end_residue->get_CA_location());
            rel.x /= 2; rel.y = 0; rel.z /= 2;
            straight_dir = straight_dir.subtract(rel);
            
            for (j=0; j<l; j++)
            {
                AminoAcid* aa = segments[i].prot->get_residue(sr+j);
                if (!aa) continue;

                for (k=0; k<m1; k++)
                {
                    float r = vcc[k]->distance_to((Molecule*)aa), rc, clash;
                    if (r >= 10) continue;
                    rc = vcc[k]->get_reach() + aa->get_reach();
                    if (r > rc) continue;

                    clash = vcc[k]->get_intermol_clashes((Molecule*)aa);
                    if (clash < 2) continue;

                    h = 0;
                    if (clash)
                    {
                        rel = aa->get_CA_location().subtract(vcc[k]->get_CA_location());
                        rel.scale((fmax(1, fmin(rc, 8) - r)) / 3);

                        #if _dbg_flexion
                        cout << *aa << " is clashing with " << *vcc[k] << " by " << clash << " cu.A." << endl;
                        #endif

                        // TODO: Increase rel according to distance from pivot.

                        seg_motion.x = larger(seg_motion.x, rel.x);
                        seg_motion.y = larger(seg_motion.y, rel.y);
                        seg_motion.z = larger(seg_motion.z, rel.z);
                        
                        h++;
                        if (h > 10) break;
                        
                        r = vcc[k]->distance_to((Molecule*)aa), rc, clash;
                        rc = vcc[k]->get_reach() + aa->get_reach();
                        clash = vcc[k]->get_intermol_clashes((Molecule*)aa);
                    }
                }       // for k
            }       // for j

            if (seg_motion.magnitude() >= 0.1)
            {
                straight_dir.scale(seg_motion.magnitude() * 0.666);
                seg_motion = seg_motion.add(straight_dir);

                seg_motion.scale(seg_motion.magnitude() / 33);
                
                #if _dbg_flexion
                cout << "Moving segment " << i << " by " << seg_motion << "..." << endl;
                #endif

                segments[i].do_motion(seg_motion);
            }
        }           // for i

        cout << "." << flush;
    }           // for iter
    cout << endl;

    n = makesure.size();
    if (n > 0) cout << "Ensuring internal contacts..." << endl;
    for (i=0; i<n; i++)
    {
        makesure[i].ensure_contact();
    }

    Point exrloc[8];
    for (i=1; i<=7; i++)
    {
        exrloc[i] = Point(0,0,0);
        n = 0;
        AminoAcid* aa;
        std::string rgn = (std::string)"TMR" + std::to_string(i);
        if (i & 1)                  // "Descending" sequence.
        {
            j = p1.get_region_start(rgn.c_str());
            for (l=0; l<4; l++)
            {
                aa = p1.get_residue(j+l);
                if (aa)
                {
                    exrloc[i] = exrloc[i].add(aa->get_CA_location());
                    n++;
                }
            }
        }
        else                        // "Ascending" sequence.
        {
            j = p1.get_region_end(rgn.c_str());
            for (l=0; l<4; l++)
            {
                aa = p1.get_residue(j-l);
                if (aa)
                {
                    exrloc[i] = exrloc[i].add(aa->get_CA_location());
                    n++;
                }
            }
        }
        
        if (n)
        {
            exrloc[i].x /= n;
            exrloc[i].y /= n;
            exrloc[i].z /= n;
        }
    }
    
    cout << "TMR spacing:" << endl;
    for (i=1; i<=7; i++)
    {
        j = i - 1;
        if (!j) j = 7;
        l = i + 1;
        if (l > 7) l = 1;
        
        float f = (exrloc[i].get_3d_distance(exrloc[j]) + exrloc[i].get_3d_distance(exrloc[l])) / 2;
        
        cout << f << " ";
    }
    cout << endl << endl;

    n = contacts.size();
    m = 0;
    for (i=0; i<n; i++)
    {
        ref = contacts[i].vector_to_contact_horizon(2);

        #if _dbg_contact_rel
        cout << "Vector to contact horizon for " << contacts[i] << ": " << ref << endl;
        #endif

        if (ref.magnitude() < 2) m++;
    }
    cout << m << " successful contacts out of " << n << " total." << endl;

    cout << "Optimizing contacts..." << endl;
    _INTERA_R_CUTOFF = 10;
    optimize_contacts(50);
    p1.get_internal_clashes(1, p1.get_end_resno(), true, 100);

    #endif

    // Add the contacts as binding residues.
    n = contacts.size();
    for (i=0; i<n; i++)
    {
        std::string bsrrem = "REMARK 800 SITE LIGAND_BINDING " + to_string(contacts[i].aa1->get_residue_no()) + (std::string)"\n";
        p1.add_remark(bsrrem);
        bsrrem = "REMARK 800 SITE LIGAND_BINDING " + to_string(1000+contacts[i].aa2->get_residue_no()) + (std::string)"\n";
        p2.add_remark(bsrrem);
    }

    // Write the output file.
    fp = fopen(output_fname, "wb");
    if (!fp)
    {
        cout << "Failed to open output file." << endl;
        return -1;
    }
    p1.set_pdb_chain('A');
    p1.save_pdb(fp);
    p2.renumber_residues(1, p2.get_end_resno(), 1001);
    p2.set_pdb_chain('B');
    p2.save_pdb(fp);
    p2.end_pdb(fp);
    fclose(fp);
    cout << "Wrote output file." << endl;
}







