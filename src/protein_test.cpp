#include <cstring>
#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>
#include "classes/protein.h"

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
    p.set_clashables();

	char outfn[32];
    p.make_helix(1, p.get_seq_length(), M_PI, M_PI);
    strcpy(outfn, "test.pdb");
    FILE* pf = fopen(outfn, "wb");
    p.save_pdb(pf);
    p.end_pdb(pf);
    fclose(pf);
    cout << "Wrote " << outfn << endl;
    
    p.make_helix(1, p.get_seq_length(), ALPHA_PHI, ALPHA_PSI);
    // p.make_helix(1, p.get_seq_length(), ALPHA_PHI, ALPHA_PSI-ALPHA_PHI);
    strcpy(outfn, "test_alpha.pdb");
    pf = fopen(outfn, "wb");
    p.save_pdb(pf);
    p.end_pdb(pf);
    fclose(pf);
    cout << "Wrote " << outfn << endl;
    
    p.make_helix(1, p.get_seq_length(), BETA_PHI, BETA_PSI);
    strcpy(outfn, "test_beta.pdb");
    pf = fopen(outfn, "wb");
    p.save_pdb(pf);
    p.end_pdb(pf);
    fclose(pf);
    cout << "Wrote " << outfn << endl;
    
    p.make_helix(1, p.get_seq_length(), _310_PHI, _310_PSI);
    strcpy(outfn, "test_310.pdb");
    pf = fopen(outfn, "wb");
    p.save_pdb(pf);
    p.end_pdb(pf);
    fclose(pf);
    cout << "Wrote " << outfn << endl;
    
    p.make_helix(1, p.get_seq_length(), PI_PHI, PI_PSI);
    strcpy(outfn, "test_pi.pdb");
    pf = fopen(outfn, "wb");
    p.save_pdb(pf);
    p.end_pdb(pf);
    fclose(pf);
    cout << "Wrote " << outfn << endl;
    
    

    Molecule m("Test2");
    pf = fopen(outfn, "rb");
    m.from_pdb(pf);
    fclose(pf);
    // m.minimize_internal_clashes();

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













