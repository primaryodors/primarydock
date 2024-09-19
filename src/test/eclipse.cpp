
#include <cstring>
#include "../classes/molecule.h"

using namespace std;

int main(int argc, char** argv)
{
    Molecule butane("butane");
    butane.from_smiles("CCCC");
    butane.movability = MOV_ALL;

    int i, j, k, l, m, n;
    Bond* b = butane.get_atom("C3")->get_bond_between("C2");
    Atom *C1 = butane.get_atom("C1"), *C4 = butane.get_atom("C4");
    float theta, step = fiftyseventh*2;

    cout << endl << "BUTANE" << endl;
    for (theta=0; theta<M_PI*2; theta += step)
    {
        cout << (theta*fiftyseven) << "deg";

        float energy = butane.total_eclipses();
        cout << " energy " << energy/4.184;

        b->rotate(step);
        cout << endl;
    }
}
