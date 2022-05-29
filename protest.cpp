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
    
    if (argc < 2)
    {
    	cout << "Usage:" << endl << "protest [sequence]" << endl;
    	cout << endl;
    	cout << "Optionally, you can override the default aminos.dat file with --dat [filename], for example if using non-standard amino acids.";
    	cout << endl;
    	return -1;
    }

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
    
    // Stretch it out.
    Point pt = p.get_atom_location(1, "N");
    pt.x += 10000;
    int len = p.get_seq_length();
    p.conform_backbone(1, len, p.get_atom(len, "C"), pt, 50);

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













