
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "classes/molecule.h"

using namespace std;

Molecule mol("TheMolecule");

int main(int argc, char** argv)
{
    if (argc > 1)
        mol.from_smiles(argv[1]);
    else
        mol.from_smiles("c1ccccc1C=O");
    
    Bond** b = mol.get_rotatable_bonds();
    if (!b || !b[0])
    {
        cout << "No rotatable bonds." << endl;
        return 0;
    }

    int i, n=0;
    for (i=0; b[i]; i++)
    {
        if (b[i]->can_rotate)
        {
            cout << "Bond " << *b[i] << " can rotate." << endl;
            n++;
        }
        else if (b[i]->can_flip)
        {
            cout << "Bond " << *b[i] << " can flip." << endl;
            n++;
        }
    }

    if (!n) cout << "No rotatable bonds." << endl;

    return 0;
}

