#include "classes/protein.h"
#include "classes/dynamic.h"

char configfname[256];
char protfname[256];
char protstrand = '\0';

std::vector<std::string> dyn_motion_strings;
float dynamic_minimum = 0;
float dynamic_initial = 1;
int dynamic_every_iter = 13;
DynamicMotion* dyn_motions[64];
int num_dyn_motions = 0;

std::string origbuff = "";
std::vector<std::string> bridges;
std::vector<std::string> atomto;
std::string outpdb;

bool configset=false, protset=false;

Protein* protein;

int interpret_resno(const char* field)
{
    char buffer[strlen(field)+4];
    strcpy(buffer, field);
    char* dot = strchr(buffer, '.');
    if (dot)
    {
        *(dot++) = 0;
        int b = atoi(buffer);
        int w = atoi(dot);
        int _50 = protein->get_bw50(b);
        if (_50 < 1)
        {
            cout << "Error: unknown BW number " << b << "." << w << ", please ensure PDB file has REMARK 800 SITE BW words." << endl;
            throw 0xbad12e5;
        }
        return _50 + w - 50;
    }
    else return atoi(buffer);
}

int interpret_config_line(char** words)
{
    int i;

    if (0) { ; }
    else if (!strcmp(words[0], "ATOMTO"))
    {
        atomto.push_back(origbuff);
    }
    else if (!strcmp(words[0], "BRIDGE"))
    {
        std::string str = words[1];
        str += (std::string)"|" + (std::string)words[2];
        bridges.push_back(str);
        return 2;
    }
    else if (!strcmp(words[0], "OUTPDB"))
    {
        outpdb = words[2] ?: words[1];
        return 2;
    }
    else if (!strcmp(words[0], "PROT"))
    {
        strcpy(protfname, words[1]);
        char* c = strchr(protfname, ':');
        if (c)
        {
            *c = 0;
            protstrand = *(++c);
        }
        else protstrand = 0;
        protset = true;
        return 1;
    }
    else if (!strcmp(words[0], "DYNAMIC"))
    {
        dyn_motion_strings.push_back(origbuff);
    }
    else if (!strcmp(words[0], "DYNMIN"))
    {
        dynamic_minimum = atof(words[1]);
    }
    else if (!strcmp(words[0], "DYNINIT"))
    {
        dynamic_initial = atof(words[1]);
    }
    else if (!strcmp(words[0], "DYNEVERY"))
    {
        dynamic_every_iter = atoi(words[1]);
    }

    return 0;
}

void read_config_file(FILE* pf)
{
    char buffer[65536];
    int i;

    while (!feof(pf))
    {
        char* wgas = fgets(buffer, 1015, pf);
        origbuff = buffer;
        if (buffer[0] >= ' ' && buffer[0] != '#')
        {
            char** words = chop_spaced_words(buffer);
            if (!words) continue;

            interpret_config_line(words);

            delete words;
        }
        buffer[0] = 0;
    }
}

void freeze_bridged_residues()
{
    int i, l;

    if (bridges.size())
    {
        for (i=0; i<bridges.size(); i++)
        {
            int resno1 = interpret_resno(bridges[i].c_str());
            const char* r2 = strchr(bridges[i].c_str(), '|');
            if (!r2) throw 0xbadc0de;
            r2++;
            int resno2 = interpret_resno(r2);
            
            AminoAcid *aa1 = protein->get_residue(resno1), *aa2 = protein->get_residue(resno2);
            if (aa1)
            {
                aa1->movability = MOV_PINNED;
                aa1->been_flexed = true;
                Bond** bb = aa1->get_rotatable_bonds();
                if (bb)
                {
                    for (l=0; bb[l]; l++)
                    {
                        bb[l]->can_rotate = false;
                    }
                }
            }
            if (aa2)
            {
                aa2->movability = MOV_PINNED;
                aa2->been_flexed = true;
                Bond** bb = aa2->get_rotatable_bonds();
                if (bb)
                {
                    for (l=0; bb[l]; l++)
                    {
                        bb[l]->can_rotate = false;
                    }
                }
            }
        }
    }
}

void reconnect_bridges()
{
    int i;
    for (i=0; i<bridges.size(); i++)
    {
        int resno1 = interpret_resno(bridges[i].c_str());
        const char* r2 = strchr(bridges[i].c_str(), '|');
        if (!r2) throw 0xbadc0de;
        r2++;
        int resno2 = interpret_resno(r2);

        #if _dbg_bridges
        cout << "Bridging " << resno1 << " and " << resno2 << "..." << endl;
        #endif

        protein->bridge(resno1, resno2);

        AminoAcid *aa1 = protein->get_residue(resno1), *aa2 = protein->get_residue(resno2);
        if (aa1) aa1->movability = MOV_PINNED;
        if (aa2) aa2->movability = MOV_PINNED; 

        #if _dbg_bridges
        if (!aa1) cout << resno1 << " not found." << endl;
        if (!aa2) cout << resno2 << " not found." << endl;
        if (aa1 && aa2)
        {
            float tb = -aa1->get_intermol_binding(aa2);
            cout << "Bridge energy " << tb << " kJ/mol." << endl;
        }
        #endif
    }

    freeze_bridged_residues();
}

void apply_protein_specific_settings(Protein* p)
{
    int i, n = dyn_motion_strings.size();

    for (i=0; i<n; i++)
    {
        if (dyn_motions[i]) delete dyn_motions[i];
        dyn_motions[i] = new DynamicMotion(p);
        dyn_motions[i]->read_config_line(dyn_motion_strings[i].c_str(), dyn_motions);
        dyn_motions[i]->minimum = dynamic_minimum;
        dyn_motions[i]->apply_absolute(dynamic_minimum);
    }
    num_dyn_motions = i;
    dyn_motions[i] = nullptr;
}


int main(int argc, char** argv)
{
    dyn_motions[0] = nullptr;

    int i, j;
    for (i=0; i<256; i++)
        configfname[i] = protfname[i] = 0;

    strcpy(configfname, "primarydock.config");
    for (i=1; i<argc; i++)
    {
        if (argv[i][0] == '-' && argv[i][1] == '-')
        {
            argv[i] += 2;
            for (j=0; argv[i][j]; j++) if (argv[i][j] >= 'a' && argv[i][j] <= 'z') argv[i][j] &= 0x5f;
            j = interpret_config_line(&argv[i]);
            argv[i] -= 2;
            i += j;
        }
        else
        {
            strcpy(configfname, argv[i]);
            configset = true;
        }
    }

    FILE* pf = fopen(configfname, "r");
    if (!pf)
    {
        cout << "Config file not found: " << configfname << ", exiting." << endl;
        return 0xbadf12e;
    }

    read_config_file(pf);
    fclose(pf);

    for (i=1; i<argc; i++)
    {
        if (argv[i][0] == '-' && argv[i][1] == '-')
        {
            argv[i] += 2;
            for (j=0; argv[i][j]; j++) if (argv[i][j] >= 'a' && argv[i][j] <= 'z') argv[i][j] &= 0x5f;
            j = interpret_config_line(&argv[i]);
            argv[i] -= 2;
            i += j;
        }
    }

    char protid[255];
    char* slash = strrchr(protfname, '/');
    if (!slash) slash = strrchr(protfname, '\\');
    strcpy(protid, slash ? slash+1 : protfname );
    char* dot = strchr(protid, '.');
    if (dot) *dot = 0;

    protein = new Protein(protid);
    pf = fopen(protfname, "r");
    if (!pf)
    {
        cout << "Error trying to read " << protfname << endl;
        return 0xbadf12e;
    }
    protein->load_pdb(pf, 0, protstrand ?: 'A');
    fclose(pf);
    apply_protein_specific_settings(protein);


    for (i=0; i<num_dyn_motions; i++)
    {
        dyn_motions[i]->apply_absolute(1);
    }

    Molecule* repack_residues[1024];
    MovabilityType repacked_movabilities[1024];
    int n = protein->get_end_resno();
    bool already_set_for_repack[n+1];
    for (i=0; i<n; i++) already_set_for_repack[i] = false;

    n = 0;
    for (i=0; i<num_dyn_motions; i++)
    {
        int j;
        AminoAcid *sr, *er;

        sr = protein->get_residue(dyn_motions[i]->start_resno);
        er = protein->get_residue(dyn_motions[i]->end_resno);
        if (!sr || !er) continue;

        for (j = sr->get_residue_no(); j <= er->get_residue_no(); j++)
        {
            if (already_set_for_repack[j]) continue;
            Star s;
            s.paa = protein->get_residue(j);
            if (s.n)
            {
                repacked_movabilities[n] = s.paa->movability;
                if (s.paa->movability != MOV_PINNED) s.paa->movability = MOV_FLEXONLY;
                repack_residues[n++] = s.pmol;
                already_set_for_repack[j] = true;
            }
        }
    }
    repack_residues[n] = 0;

    reconnect_bridges();

    n = atomto.size();
    for (i=0; i<n; i++)
    {
        char buffer[1024];
        strcpy(buffer, atomto[i].c_str());
        char** words = chop_spaced_words(buffer);

        if (!words[1]) throw -1;
        AminoAcid* aa = protein->get_residue_bw(words[1]);
        if (!words[2]) throw -1;
        char* aname = words[2];
        if (!words[3]) throw -1;
        AminoAcid* target = protein->get_residue_bw(words[3]);
        if (words[4]) throw -1;

        if (!aa)
        {
            cout << "Warning: residue " << words[1] << " not found." << endl;
            continue;
        }

        if (!target)
        {
            cout << "Warning: residue " << words[3] << " not found." << endl;
            continue;
        }

        Atom* a = aa->get_atom(aname);
        if (!strcmp("EXTENT", aname)) a = aa->get_reach_atom();
        if (!a)
        {
            cout << "Warning: atom not found " << *aa << ":" << aname << endl;
            continue;
        }

        aa->movability = MOV_FLEXONLY;
        aa->conform_atom_to_location(a->name, target->get_CA_location());
        aa->movability = MOV_FLXDESEL;
    }

    Molecule::conform_molecules(repack_residues);

    FILE* pfout = fopen(outpdb.c_str(), "wb");
    if (!pfout)
    {
        cout << "Failed to open " << outpdb << " for writing." << endl;
        return -1;
    }

    protein->save_pdb(pfout);
    protein->end_pdb(pfout);
    fclose(pfout);

    return 0;
}
