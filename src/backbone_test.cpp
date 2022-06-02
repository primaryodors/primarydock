
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstring>
#include "classes/protein.h"

using namespace std;

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        cout << "No input file." << endl;
        return -1;
    }

    Protein p(argv[1]);
    FILE* pf = fopen(argv[1], "r");
    if (!pf)
    {
        cout << "Error trying to read " << argv[1] << endl;
        return 0xbadf12e;
    }
    p.load_pdb(pf);
    fclose(pf);

    // Just before the start of TMR1.
    // p.rotate_backbone(21, N_asc, 180.0*M_PI/180);
    p.delete_residues(21, 309);
    p.make_helix(1, 20, M_PI, M_PI);				// Stretched out.
    // p.make_helix(1, 20, ALPHA_PHI, ALPHA_PSI);		// Alpha helix.
    // p.conform_backbone(1, 20, 50);
    p.delete_sidechains(1, 20);

    pf = fopen("bktest.pdb", "wb");

    p.save_pdb(pf);
    p.end_pdb(pf);
    fclose(pf);

    return 0;
}
