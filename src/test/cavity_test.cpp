
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "../classes/cavity.h"

using namespace std;

Protein* protein;
Molecule* ligand;

int interpret_resno(const char* field)
{
    char buffer[strlen(field)+4];
    strcpy(buffer, field);

    char* offset = buffer;
    while (*offset >= 'A' && *offset <= 'Z') offset++;

    int retval = 0;
    char* dot = strchr(buffer, '.');
    if (dot)
    {
        *(dot++) = 0;
        int b = atoi(offset);
        int w = atoi(dot);
        int _50 = protein->get_bw50(b);
        if (_50 < 1)
        {
            cout << "Error: unknown BW number " << b << "." << w << ", please ensure PDB file has REMARK 800 SITE BW words." << endl;
            throw 0xbad12e5;
        }
        if (_50 < 1) return 0;
        else retval = _50 + w - 50;
    }
    else retval = atoi(offset);

    if (offset == buffer) return retval;

    AminoAcid* aa = protein->get_residue(retval);
    if (!aa) return 0;

    int i;
    for (i=0; buffer[i] >= 'A' && buffer[i] <= 'Z'; i++)
    {
        if (buffer[i] == aa->get_letter()) aa->priority = true;
    }
    return retval;
}

int main(int argc, char** argv)
{
    bool verbose = false;
    Protein p("TheProtein");
    Molecule m("ligand");
    char buffer[65536];
    protein = &p;
    ligand = &m;

    FILE* fp;

    int i;
    for (i=1; i<argc; i++)
    {
        if (!strcmp(argv[i], "-p") || !strcmp(argv[i], "--prot"))
        {
            i++;
            if (file_exists(argv[i]))
            {
                fp = fopen(argv[i], "rb");
                p.load_pdb(fp);
            }
            else
            {
                cout << "File not found: " << argv[i] << endl;
                return -1;
            }
        }
        else if (!strcmp(argv[i], "-l") || !strcmp(argv[i], "--lig"))
        {
            i++;
            if (file_exists(argv[i]))
            {
                fp = fopen(argv[i], "rb");
                if (fp)
                {
                    fread(buffer, 1, 65536, fp);
                    fclose(fp);
                    m.from_sdf(buffer);
                }
                else
                {
                    cout << "Failed to open " << argv[i] << " for reading." << endl;
                    return -1;
                }
            }
            else
            {
                cout << "File not found: " << argv[i] << endl;
                return -1;
            }
        }
        else if (!strcmp(argv[i], "-b") || !strcmp(argv[i], "--bsr"))           // Binding Site Residues.
        {
            if (!p.get_seq_length())
            {
                cout << "Protein must be specified before binding site residues." << endl;
                return -1;
            }

            i++;
            while (i < argc)
            {
                if (argv[i][0] == '-') break;
                int resno = interpret_resno(argv[i]);
                if (!resno) continue;
                AminoAcid* aa = p.get_residue(resno);
                if (aa) aa->priority = true;
                i++;
            }
        }
    }

    if (!p.get_seq_length())
    {
        cout << "No protein specified." << endl;
        return -1;
    }

    ////////////////////////////////////////////////////////////////////////////
    // Phase I: Cavity Search.                                                //
    ////////////////////////////////////////////////////////////////////////////

    Cavity cavities[1029];
    int qfound = Cavity::scan_in_protein(&p, cavities, 1024);

    cout << "Found " << qfound << " cavit" << (qfound == 1 ? "y." : "ies.") << endl;

    fp = fopen("tmp/cavities.js", "w");
    if (fp)
    {
        for (i=0; i<qfound; i++)
        {
            cavities[i].output_ngl_js(fp);
            fprintf(fp, "\n\n\n\n\n");
        }
        fclose(fp);
    }


    ////////////////////////////////////////////////////////////////////////////
    // Phase II: Fleximer Search.                                             //
    ////////////////////////////////////////////////////////////////////////////

    if (m.get_atom_count())
    {
        // TODO: Generate fleximers of pockets and ligand.
        Pose best(&m);
        float bestviol = Avogadro;
        cout << "Trying ligand in pockets..." << flush;
        for (i=0; i<qfound; i++)
        {
            Point cavcen = cavities[i].get_center();
            if (cavcen.y < 5 || cavcen.y > 18) continue;
            float viol = cavities[i].find_best_containment(&m);
            if (verbose) cout << endl << "Cavity centered at " << cavcen << " has " << viol << " containment violations";
            if (viol < bestviol)
            {
                best.copy_state(&m);
                bestviol = viol;
            }
            cout << "." << flush;
        }
        cout << endl;

        if (bestviol > 5e22)
        {
            cout << "No suitable pockets found." << endl;
            return 0;
        }

        best.restore_state(&m);
        cout << "Best pre-iteration pose has " << bestviol << " containment violations." << endl;
        cout << "Ligand centered at " << m.get_barycenter() << endl;
        fp = fopen("tmp/quickdock.pdb", "wb");
        if (!fp) return -1;
        p.save_pdb(fp, &m);
        cout << "Saved output file in tmp/." << endl;
    }

    return 0;
}


