#include <cstring>
#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>
#include "protein.h"

using namespace std;

int main(int argc, char** argv)
{
    int i, seqarg=1, namearg=0;

    if (!strcmp("--dat", argv[seqarg]))
    {
        override_aminos_dat = argv[seqarg+1];
        seqarg += 2;
    }

    if (!strcmp("--name", argv[seqarg]))
    {
        namearg = seqarg+1;
        seqarg += 2;
    }

    Protein p(namearg?argv[namearg]:"Test");
    p.add_sequence(argv[seqarg]);

    const char* outfn = "test.pdb";
    FILE* pf = fopen(outfn, "wb");
    p.save_pdb(pf);
    p.end_pdb(pf);
    fclose(pf);
    cout << "Wrote " << outfn << endl;

    Molecule m("Test2");
    pf = fopen(outfn, "rb");
    m.from_pdb(pf);
    fclose(pf);
    // m.minimize_internal_collisions();

    const char* outfn2 = "test2.sdf";
    pf = fopen(outfn2, "wb");
    m.save_sdf(pf);
    fclose(pf);
    cout << "Wrote " << outfn2 << endl;

    Protein p1("Test3");
    pf = fopen(outfn, "rb");
    int rescount = p1.load_pdb(pf);
    fclose(pf);
    cout << "Read " << rescount << " residue(s)." << endl;


    return 0;
}













