
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "../classes/cavity.h"

using namespace std;

int main(int argc, char** argv)
{
    Protein p("TheProtein");

    FILE* fp;
    if (argc > 1)
    {
        if (file_exists(argv[1]))
        {
            fp = fopen(argv[1], "rb");
            p.load_pdb(fp);
        }
    }

    if (!p.get_seq_length())
    {
        cout << "No protein specified." << endl;
        return -1;
    }

    Cavity cavities[1029];
    int qfound = Cavity::scan_in_protein(&p, cavities, 1024);

    cout << "Found " << qfound << " cavities." << endl;

    int i;
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

    if (argc > 2)
    {
        if (file_exists(argv[2]))
        {
            fp = fopen(argv[2], "rb");
            if (fp)
            {
                Molecule m("ligand");
                char buffer[65536];
                fread(buffer, 1, 65536, fp);
                fclose(fp);
                m.from_sdf(buffer);

                Pose best(&m);
                float bestviol = Avogadro;
                cout << "Trying ligand in pockets..." << flush;
                for (i=0; i<qfound; i++)
                {
                    Point cavcen = cavities[i].get_center();
                    if (cavcen.y < 0 || cavcen.y > 15) continue;
                    float viol = cavities[i].find_best_containment(&m);
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
        }
    }

    return 0;
}


