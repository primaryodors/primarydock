
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <sys/stat.h>
#include "classes/group.h"

using namespace std;

int main(int argc, char** argv)
{
    Molecule m("TheLigand");
    Protein p("TheProtein");
    FILE* fp;
    std::string pdbfn = "testdata/OR12D2.propionic_acid.tester.pdb";

    fp = fopen(pdbfn.c_str(), "rb");
    p.load_pdb(fp);

    fseek(fp, 0, SEEK_SET);
    m.from_pdb(fp, true);

    AminoAcid* aa = p.get_residue(107);
    float before = -m.get_intermol_binding(aa);
    cout << "Initial energy: " << before << endl;

    Atom* a = aa->get_atom("CB");
    Atom* b = aa->get_atom("OG");
    Bond* bt = a->get_bond_between(b);
    float theta=0, during, step = hexagonal/3;

    for (theta=step; theta<M_PI*2; theta += step)
    {
        bt->rotate(step);
        during = -m.get_intermol_binding(aa);
        cout << (theta*fiftyseven) << "deg: " << during << endl;
    }

    aa->movability = MOV_FORCEFLEX;
    m.movability = MOV_PINNED;

    Molecule* mols[4];
    mols[0] = &m;
    mols[1] = (Molecule*)aa;
    mols[2] = nullptr;

    Molecule::conform_molecules(mols);
    float after = -m.get_intermol_binding(aa);
    cout << "Final energy: " << after << endl;
}
