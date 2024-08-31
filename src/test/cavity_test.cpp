
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

    fp = fopen("tmp/cavities.js", "w");
    if (fp)
    {
        int i;
        for (i=0; i<qfound; i++)
        {
            cavities[i].output_ngl_js(fp);
            fprintf(fp, "\n\n\n\n\n");
        }
        fclose(fp);
    }

    return 0;
}


