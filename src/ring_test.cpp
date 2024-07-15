
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstring>
#include "classes/molecule.h"

using namespace std;

int main(int argc, char** argv)
{
    Molecule m("Test");
    char buffer[65536];

    if (argc > 1)
    {
        FILE* pf = fopen(argv[1], "rb");
        if (!pf)
        {
            m.from_smiles(argv[1]);
        }
        else
        {
            size_t rfw = fread(buffer, 1, 65535, pf);
            fclose(pf);
            m.from_sdf(buffer);
        }
    }
    else
    {
        cout << "Usage:" << endl << endl << "test/ring_test {path/to/structure.sdf or SMILES string}" << endl << endl;
        return -1;
    }

    m.hydrogenate();
    m.minimize_internal_clashes();

    int i, n = m.get_num_rings();
    cout << "Molecule has " << n << " ring(s)." << endl;

    for (i=0; i<n; i++)
    {
        cout << "Ring " << i;
        if (m.ring_is_aromatic(i)) cout << " is aromatic and";
        else if (m.ring_is_coplanar(i)) cout << " is coplanar and";
        cout << " consists of ";
        Atom** a = m.get_ring_atoms(i);
        Atom::dump_array(a);
    }
}
