
#include <cstring>
#include "../classes/molecule.h"

using namespace std;

int main(int argc, char** argv)
{
    Molecule ethane("ethane");
    ethane.from_smiles("CC");

    Molecule butane("butane");
    butane.from_smiles("CCCC");

    int i, j, k, l, m, n;
    Bond* b = ethane.get_atom("C1")->get_bond_between("C2");
    float theta, step = fiftyseventh*2;
    cout << endl << "ETHANE" << endl;
    for (theta=0; theta<M_PI*2; theta += step)
    {
        cout << (theta*fiftyseven) << "deg";

        float energy = 0;

        Bond *abt[16], *bbt[16];
        b->atom1->fetch_bonds(abt);
        b->atom2->fetch_bonds(bbt);
        m = b->atom1->get_geometry();
        n = b->atom2->get_geometry();
        SCoord axis = b->get_axis();

        for (j=0; j<m; j++)
        {
            if (!abt[j]) continue;
            if (!abt[j]->atom1) continue;
            if (!abt[j]->atom2) continue;
            if (abt[j]->atom2 == b->atom2) continue;

            for (l=0; l<n; l++)
            {
                if (!bbt[l]) continue;
                if (!bbt[l]->atom1) continue;
                if (!bbt[l]->atom2) continue;
                if (bbt[l]->atom2 == b->atom1) continue;

                cout << " (" << abt[j]->atom2->name << "~" << bbt[l]->atom2->name << ")";
                float sigma = abt[j]->atom2->get_vdW_radius() + bbt[l]->atom2->get_vdW_radius();
                float r = abt[j]->atom2->distance_to(bbt[l]->atom2);
                float sigma_r = sigma / r;

                energy += Lennard_Jones_epsilon_x4 * pow(sigma_r, 12) / 2;
            }
        }

        cout << " energy " << energy/4.184;

        b->rotate(step);
        cout << endl;
    }

    b = butane.get_atom("C3")->get_bond_between("C2");
    cout << endl << "BUTANE" << endl;
    for (theta=0; theta<M_PI*2; theta += step)
    {
        cout << (theta*fiftyseven) << "deg";

        float energy = 0;

        Bond *abt[16], *bbt[16];
        b->atom1->fetch_bonds(abt);
        b->atom2->fetch_bonds(bbt);
        m = b->atom1->get_geometry();
        n = b->atom2->get_geometry();
        SCoord axis = b->get_axis();

        for (j=0; j<m; j++)
        {
            if (!abt[j]) continue;
            if (!abt[j]->atom1) continue;
            if (!abt[j]->atom2) continue;
            if (abt[j]->atom2 == b->atom2) continue;

            for (l=0; l<n; l++)
            {
                if (!bbt[l]) continue;
                if (!bbt[l]->atom1) continue;
                if (!bbt[l]->atom2) continue;
                if (bbt[l]->atom2 == b->atom1) continue;

                cout << " (" << abt[j]->atom2->name << "~" << bbt[l]->atom2->name << ")";
                float sigma = abt[j]->atom2->get_vdW_radius() + bbt[l]->atom2->get_vdW_radius();
                float r = abt[j]->atom2->distance_to(bbt[l]->atom2);
                float sigma_r = sigma / r;

                energy += Lennard_Jones_epsilon_x4 * pow(sigma_r, 12) / 2;
            }
        }

        cout << " energy " << energy/4.184;

        b->rotate(step);
        cout << endl;
    }
}
