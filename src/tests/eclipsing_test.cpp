
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstring>
#include "../classes/molecule.h"

using namespace std;

int main(int argc, char** argv)
{
    Molecule m("Test");
    m.from_smiles("CC");
    float e = m.total_eclipses();
    cout << "Initial eclipses: " << e << endl;

    FILE* fp = fopen("testdata/staggered.sdf", "wb");
    if (!fp) return -1;
    m.save_sdf(fp);
    fclose(fp);

    Atom *C1 = m.get_atom("C1"), *C2 = m.get_atom("C2");
    if (!C1 || !C2) return -1;
    Bond *b = C1->get_bond_between(C2);
    if (!b) return -1;

    b->rotate(hexagonal);
    e = m.total_eclipses();
    cout << "Final eclipses: " << e << endl;

    fp = fopen("testdata/eclipsed.sdf", "wb");
    if (!fp) return -1;
    m.save_sdf(fp);
    fclose(fp);
}
