
#include "classes/protein.h"

int main(int argc, char** argv)
{
    if (argc < 1) return -1;

    Protein p("testing");
    int i;
    FILE* fp;
    bool show_colors = true;
    bool skip_helices = false;
    bool include_glycine = true, include_not_glycine = true;
    bool include_proline = true, include_not_proline = true;
    bool include_pre_proline = true, include_not_pre_proline = true;
    for (i=1; i<argc; i++)
    {
        if (argv[i][0] == '-')
        {
            if (!strcmp(argv[i], "-n") || !strcmp(argv[i], "--numbers")) show_colors = false;
            if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--skip-helices")) skip_helices = true;
            if (!strcmp(argv[i], "-g") || !strcmp(argv[i], "--omit-glycine")) include_glycine = false;
            if (!strcmp(argv[i], "-G") || !strcmp(argv[i], "--only-glycine")) include_not_glycine = false;
            if (!strcmp(argv[i], "-p") || !strcmp(argv[i], "--omit-proline")) include_proline = false;
            if (!strcmp(argv[i], "-P") || !strcmp(argv[i], "--only-proline")) include_not_proline = false;
            if (!strcmp(argv[i], "-b") || !strcmp(argv[i], "--omit-pre-proline")) include_pre_proline = false;
            if (!strcmp(argv[i], "-B") || !strcmp(argv[i], "--only-pre-proline")) include_not_pre_proline = false;
            if (!include_glycine && !include_not_glycine) return -1;
            continue;
        }

        fp = fopen(argv[i], "rb");

        if (argc > (i+1)) p.load_pdb(fp, 0, argv[i+1][0]);
        else p.load_pdb(fp);
        p.set_name_from_pdb_name(argv[i]);
    }
    if (!fp) return -1;

    fclose(fp);

    int er = p.get_end_resno();
    int max=0, x, y;
    int grid[36][36];

    for (x=0; x<36; x++) for (y=0; y<36; y++) grid[x][y] = 0;

    for (i=1; i<=er; i++)
    {
        AminoAcid* aa = p.get_residue(i);
        if (!aa) continue;
        if (skip_helices && aa->is_alpha_helix()) continue;
        if (aa->get_letter() == 'G')
        {
            if (!include_glycine || !include_not_proline) continue;
        }
        else if (aa->get_letter() == 'P')
        {
            if (!include_proline || !include_not_glycine) continue;
        }
        else if (!include_not_glycine || !include_not_proline) continue;

        if (!include_pre_proline || !include_not_pre_proline)
        {
            AminoAcid* next = p.get_residue(i+1);
            if (next)
            {
                if (next->get_letter() == 'P')
                {
                    if (!include_pre_proline) continue;
                }
                else if (!include_not_pre_proline) continue;
            }
        }


        float phi = aa->get_phi(), psi = aa->get_psi();

        phi = fmod(phi, (M_PI*2)); if (phi >= M_PI) phi -= M_PI*2;
        psi = fmod(psi, (M_PI*2)); if (psi >= M_PI) psi -= M_PI*2;

        x = (phi+M_PI) * fiftyseven / 10; if (x<0) x=0; else if (x>35) x=35;
        y = 35 - (psi+M_PI) * fiftyseven / 10; if (y<0) y=0; else if (y>35) y=35;

        grid[x][y]++;
        if (grid[x][y] > max) max = grid[x][y];
    }

    float scale = 1.0 / max;
    cout << endl << "Ramachandran plot for " << p.get_name() << endl << endl;
    cout << "-180";
    for (i=0; i<32; i++) cout << " ";
    cout << "\u03c6";
    for (i=0; i<32; i++) cout << " ";
    cout << "+180" << endl;

    for (y=0; y<36; y++)
    {
        if (y == 18)
        {
            for (i=0; i<36; i++) cout << "\u2500";
            cout << "\u253c";
            for (i=0; i<36; i++) cout << "\u2500";
            cout << " \u03c8" << endl;
        }
        for (x=0; x<36; x++)
        {
            if (x == 18) cout << "\u2502";
            if (show_colors)
            {
                float value = scale * grid[x][y];
                colorize(value*50);
                cout << "\u2588\u2588";
                colorless();
            }
            else cout << grid[x][y] << " ";
        }
        if (!y) cout << " +180";
        else if (y == 35) cout << " -180";
        cout << endl;
    }
}
