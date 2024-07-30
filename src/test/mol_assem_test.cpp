
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstring>
#include "../classes/molecule.h"

using namespace std;

#define _use_generational_algorithm 1

int main(int argc, char** argv)
{
    Molecule m("Test");
    cout << "Created empty molecule named " << m.get_name() << ".\n";

    if (argc > 1)
    {
        m.from_smiles(argv[1]);
    }
    else
    {
        Atom* C1 = m.add_atom("C", "C1", 0, 0);
        cout << "Added a carbon atom. Its location is " << C1->get_location().printable() << "." << endl;

        Atom* C2 = m.add_atom("C", "C2", C1, 1);
        cout << "Added another carbon atom. Its location is " << C2->get_location().printable() << "." << endl;

        Atom* O3 = m.add_atom("O", "O3", C2, 1);
        cout << "Added an oxygen atom. Its location is " << O3->get_location().printable() << "." << endl;
    }

    m.hydrogenate();
    srand(time(nullptr));

    #if _use_generational_algorithm
    m.minimize_internal_clashes();
    cout << "# Evolving structure; this may take some time..." << endl << flush;
    float anomaly = m.evolve_structure();
    cout << "# Post-evolution per-atom anomaly: " << anomaly << endl;
    #endif

    cout << "# Internal clashes: " << m.get_internal_clashes() << " cu.A." << endl;

    cout << "Minimizing internal clashes..." << endl;
    m.minimize_internal_clashes();

    /*if (argc > 4)
    {
    	Atom* a = m.get_atom(argv[2]);
    	Atom* b = m.get_atom(argv[3]);
    	float theta = atof(argv[4])/fiftyseven;

    	if (a && b)
    	{
    		Bond* bn = a->get_bond_between(b);
    		if (bn) bn->rotate(theta);
    	}
    }*/

    float int_clsh = m.get_internal_clashes();
    cout << "Adjusted internal clashes: " << int_clsh << " cu.A." << endl;

    int i, n = m.get_atom_count();

    for (i=0; i<n; i++)
    {
        Atom* a = m.get_atom(i);
        if (a && a->get_Z() > 1)
        {
            cout << "Atom " << a->name << " has bond angles:";
            a->print_bond_angles();
            cout << endl;
        }
    }

    char tstoutf[1024];
    if (argc > 2) strcpy(tstoutf, argv[2]);
    else strcpy(tstoutf, "test.sdf");
    FILE* pf = fopen(tstoutf, "wb");
    if (m.save_sdf(pf)) cout << "Saved " << tstoutf << endl;
    else cout << "Failed to save " << tstoutf << endl;
    fclose(pf);

    /*const char* tstpdbf = "test.pdb";
    pf = fopen(tstpdbf, "wb");
    m.save_pdb(pf);
    cout << "Saved " << tstpdbf << endl;
    fclose(pf);*/

    cout << endl;

    return 0;
}

