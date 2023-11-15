
#include "classes/protein.h"

int main(int argc, char** argv)
{
    if (argc < 1) return -1;

    Protein p("testing");
    int i;
    FILE* fp;
    bool show_colors = true;
    for (i=1; i<argc; i++)
    {
        if (argv[i][0] == '-')
        {
            if (!strcmp(argv[i], "-n") || !strcmp(argv[i], "--numbers")) show_colors = false;
            continue;
        }

        fp = fopen(argv[i], "rb");

        if (argc > (i+1)) p.load_pdb(fp, 0, argv[i+1][0]);
        else p.load_pdb(fp);
    }
    if (!fp) return -1;

    p.set_name_from_pdb_name(argv[1]);
    fclose(fp);

    int er = p.get_end_resno();
    int max=0, x, y;
    int grid[36][36];

    for (x=0; x<36; x++) for (y=0; y<36; y++) grid[x][y] = 0;

    for (i=1; i<=er; i++)
    {
        AminoAcid* aa = p.get_residue(i);
        if (!aa) continue;
        float phi = aa->get_phi(), psi = aa->get_psi();

        phi = fmod(phi, (M_PI*2)); if (phi >= M_PI) phi -= M_PI*2;
        psi = fmod(psi, (M_PI*2)); if (psi >= M_PI) psi -= M_PI*2;

        x = (phi+M_PI) * fiftyseven / 10; if (x<0) x=0; else if (x>35) x=35;
        y = 35 - (psi+M_PI) * fiftyseven / 10; if (y<0) y=0; else if (y>35) y=35;

        grid[x][y]++;
        if (grid[x][y] > max) max = grid[x][y];
    }

    float scale = 1.0 / max;
    for (y=0; y<36; y++)
    {
        for (x=0; x<36; x++)
        {
            if (show_colors)
            {
                float value = scale * grid[x][y];
                colorize(value*50);
                cout << "\u2588\u2588";
                colorless();
            }
            else cout << grid[x][y] << " ";
        }
        cout << endl;
    }
}
