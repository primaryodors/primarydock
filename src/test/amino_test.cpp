#include <cstring>
#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>
#include "../classes/aminoacid.h"

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
                    cout << b[j]->atom1->name << " - " << b[j]->atom2->name << endl;
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
        if (argc>1)
        {
            letter = argv[1][0];
            AminoAcid aa(letter);

            cout << "Charge: " << aa.get_charge() << endl;

            cout << "Is tyrosine-like (i.e. has an aromatic ring and a separate H-bond acceptor)? "
                << (aa.is_tyrosine_like() ? "Y" : "N") << endl;

            Bond** bb = aa.get_rotatable_bonds();
            if (bb && bb[0])
            {
                for (i=0; bb[i]; i++)
                {
                    Atom* baa[bb[i]->count_moves_with_atom2()+1];
                    bb[i]->fetch_moves_with_atom2(baa);
                    cout << bb[i]->atom1->name << "-" << bb[i]->atom2->name << " can rotate, bringing ";

                    for (j=0; baa[j]; j++)
                    {
                        cout << baa[j]->name << " ";
                    }
                    if (!j) cout << "zero atoms.";

                    cout << endl;
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
        else
        {
            const char* all_letters = "ARNDCUEQGHIKLMFPSTWYV";

            int i, j, n = strlen(all_letters);
            AminoAcid* all_aa[n+4];

            cout << "   ";
            for (i=0; i<n; i++)
            {
                all_aa[i] = new AminoAcid(all_letters[i]);
                cout << all_letters[i] << " ";
            }
            cout << endl;

            for (i=0; i<n; i++)
            {
                cout << all_letters[i] << "  ";
                for (j=0; j<n; j++)
                {
                    float s = all_aa[i]->similarity_to(all_aa[j]);
                    /*char buffer[16];
                    sprintf(buffer, "%0.2f", 0.02 * roundf(s*50));
                    cout << " " << buffer;*/

                    if (s < 0.1) cout << "  ";
                    else if (s < 0.25) cout << "\u2591\u2591";
                    else if (s < 0.5) cout << "\u2592\u2592";
                    else if (s < 0.8) cout << "\u2593\u2593";
                    else cout << "\u2588\u2588";
                }
                cout << endl;
            }
        }
    }

    return 0;
}


