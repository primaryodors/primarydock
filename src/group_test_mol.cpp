
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "classes/group.h"

using namespace std;

int main(int argc, char** argv)
{
    Molecule m("Testing");

    if (argc > 1) m.from_smiles(argv[1]);
    else m.from_smiles("CCO");

    std::vector<std::shared_ptr<AtomGroup>> mgrp = AtomGroup::get_potential_ligand_groups(&m);

    int i, n;
    n = mgrp.size();

    for (i=0; i<n; i++)
    {
        cout << *mgrp[i] << endl;
    }
}
