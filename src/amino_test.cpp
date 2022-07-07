#include <cstring>
#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>
#include "classes/aminoacid.h"

using namespace std;

int main(int argc, char** argv)
{
    char letter = 'A';
    int i, j;
    if (argc>1 && strlen(argv[1]) > 1)
    {
        int seqarg=1;

        if (!strcmp("--dat", argv[seqarg]))
        {
            override_aminos_dat = argv[seqarg+1];
            seqarg += 2;
        }

        AminoAcid* firstaa = 0;
        AminoAcid* prevaa = 0;
        for (i=0; argv[seqarg][i]; i++)
        {
            letter = argv[seqarg][i];
            AminoAcid* aa = new AminoAcid(letter, prevaa);
            if (prevaa == aa) throw 0xbad;
            prevaa = aa;
            if (!i) firstaa = aa;

            Bond** b = aa->get_rotatable_bonds();
            AADef* def = aa->get_aa_definition();

            cout << def->name << " bonds can rotate: " << endl;
            if (!b || !b[0]) cout << "None." << endl;
            else
            {
                for (j=0; b[j]; j++)
                {
                    cout << b[j]->atom->name << " - " << b[j]->btom->name << endl;
                }
            }
            cout << endl;
        }

        const char* outfn = "test.pdb";
        FILE* pf = fopen(outfn, "wb");
        AminoAcid* aa = firstaa;
        do
        {
            aa->save_pdb(pf);
        }
        while (aa = aa->get_next());
        fclose(pf);
        cout << "Wrote " << outfn << endl;

		char outfn2[50];
		AADef* aad = firstaa->get_aa_definition();
		if (aad) sprintf(outfn2, "%s.sdf", aad->name);
		else strcpy(outfn2, "aatest.sdf");
		pf = fopen(outfn2, "wb");
		firstaa->save_sdf(pf);
		fclose(pf);
		cout << "Wrote " << outfn2 << endl;
    }
    else
    {
    	if (argc>1) letter = argv[1][0];
        AminoAcid aa(letter);
        
        cout << "Is tyrosine-like (i.e. has an aromatic ring and a separate H-bond acceptor)? "
        	 << (aa.is_tyrosine_like() ? "Y" : "N") << endl;
        
        Bond** bb = aa.get_rotatable_bonds();
        if (bb && bb[0])
        {
        	for (i=0; bb[i]; i++)
        	{
        		Atom** baa = bb[i]->get_moves_with_btom();
        		cout << bb[i]->atom->name << "-" << bb[i]->btom->name << " can rotate, bringing ";
        		
        		if (baa)
        		{
        			for (j=0; baa[j]; j++)
        			{
        				cout << baa[j]->name << " ";
        			}
        			if (!j) cout << "zero atoms.";
        		}
        		else cout << "no atoms.";
        		
        		cout << endl;
        		delete[] baa;
        	}
        }
        else cout << "No rotatable bonds." << endl;
        
        cout << aa.get_name() << " hydrophilicity = " << aa.hydrophilicity() << endl;
        cout << "Similarity to A " << aa.similarity_to('A') << endl;
        cout << "Similarity to D " << aa.similarity_to('D') << endl;
        cout << "Similarity to R " << aa.similarity_to('R') << endl;
        cout << "Similarity to F " << aa.similarity_to('F') << endl;
        cout << "Similarity to C " << aa.similarity_to('C') << endl;
        cout << "Similarity to S " << aa.similarity_to('S') << endl;

        const char* outfn = "test.pdb";
        FILE* pf = fopen(outfn, "wb");
        aa.save_pdb(pf);
        fclose(pf);
        cout << "Wrote " << outfn << endl;

		const char* outfn2 = "test.sdf";
		pf = fopen(outfn2, "wb");
		aa.save_sdf(pf);
		fclose(pf);
		cout << "Wrote " << outfn2 << endl;
    }

    return 0;
}


