
#include "classes/protein.h"

int main(int argc, char** argv)
{
    if (argc < 1) return -1;
    Protein p("testing");
    FILE* fp = fopen(argv[1], "rb");
    p.load_pdb(fp);
    fclose(fp);

    int i, j, n = p.get_end_resno();
    for (i=1; i<=n; i++)
    {
        AminoAcid* aa = p.get_residue(i);
        if (!aa) continue;
        AminoAcid** cl = p.get_residues_can_clash(i);
        if (!cl) continue;

        for (j=0; cl[j]; j++)
        {
            if (cl[j] < aa) continue;
            if (fabs(cl[j]->get_residue_no() - i) < 3) continue;

            if (aa->get_residue_no() == 11 && cl[j]->get_residue_no() == 172)
            {
                n++;
                n--;
            }

            float f = aa->get_intermol_binding(cl[j]);
            if (fabs(f) > 2) cout << *aa << "-" << *cl[j] << ": " << -f << endl;
        }
    }
}