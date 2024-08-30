
#include <iostream>
#include <math.h>
#include "../classes/atom.h"

using namespace std;

void mol_shape_has_changed(Molecule* mol)
{
    // This stub is necessary whenever we use the Atom class without the Molecule class.
    ;
}

int main (int argc, char** argv)
{
    if (argc < 2)
    {
        cout << "No element symbol given.\n";
        return -1;
    }

    Point pt;
    try
    {
        Atom a(argv[1], &pt);
        cout << "Atom has Z=" << a.get_Z() << " with valence of " << a.get_valence() << " and geometry of " << a.get_geometry()
             << " and belongs to family " << a.get_family()
             << ".\n";
    }
    catch (int ex)
    {
        cout << "Exception! ";
        if (ex == VALENCE_EXCEEDS_GEOMETRY) cout << "Atom's valence should never exceed its geometry. ";
        else cout << "Number: " << ex;
        cout << "\n";
    }

    cout << "Creating a hydroxyl radical...\n";

    Atom O("O");
    Atom OH("H");
    O.bond_to(&OH, 1);

    cout << "Atom O has polarity of " << O.is_polar() << " and its H has polarity " << OH.is_polar() << ".\n";

    cout << "Creating a sulfhydryl radical...\n";

    Atom S("S");
    Atom SH("H");
    SH.bond_to(&S, 1);

    cout << "Atom S has thio of " << S.is_thio() << " and its H has thio " << SH.is_thio() << ".\n";
}
