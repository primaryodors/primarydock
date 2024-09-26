#include <cstring>
#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>
#include "../classes/protein.h"

using namespace std;

int main(int argc, char** argv)
{
    Protein p("Test Peptide");

    FILE* fp = fopen("pdbs/OR51/OR51E2.active.pdb", "rb");
    p.load_pdb(fp);
    fclose(fp);

    int i, j, k, l, n, seqlen = p.get_end_resno(), fails = 0;
    for (i=1; i<=seqlen; i++)
    {
        AminoAcid* aa = p.get_residue(i);
        if (!aa) continue;
        n = aa->get_atom_count();
        for (j=0; j<n; j++)
        {
            Atom* a = aa->get_atom(j);
            if (!a) continue;
            int ag = a->get_geometry();
            SCoord* abgeo = a->get_geometry_aligned_to_bonds(true);

            for (k=0; k<ag; k++)
            {
                Bond* b = a->get_bond_by_idx(k);
                if (!b || !b->atom2) continue;
                SCoord AB = b->atom2->get_location().subtract(a->get_location());

                float theta;
                for (l=0; l<=ag; l++)
                {
                    float th = find_3d_angle(AB, abgeo[l], Point(0,0,0));
                    if (!l || th < theta) theta = th;
                }

                if (theta > 15.0*fiftyseventh)
                {
                    cout << aa->get_name() << ":" << a->name << "-" << b->atom2->name
                        << " (" << k << "/" << ag << ")"
                        << " geometry is off by " << theta*fiftyseven << "deg." << endl;
                    fails++;
                }
                /*else cout << aa->get_name() << ":" << a->name << "-" << b->atom2->name
                    << " (" << k << "/" << ag << ")"
                    << " " << theta*fiftyseven << "deg." << endl;*/
            }
        }
    }

    if (!fails) cout << "All good." << endl;
    else cout << fails << " fails." << endl;
    return 0;
}
