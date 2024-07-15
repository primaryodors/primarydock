
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstring>
#include "classes/molecule.h"

using namespace std;

void show_usage()
{
    cout << "Usage:" << endl << "ringflip path/to/file.sdf atom1 [ atom2 [ atom3 [ ... ] ] ]" << endl;
}

int main(int argc, char** argv)
{
    int i, l, n;
    std::string filename;
    std::vector<std::string> atom_names;

    if (argc < 3)
    {
        show_usage();
        return -1;
    }

    for (i=1; i<argc; i++)
    {
        if (argv[i][0] == '-')
        {
            // TODO:
        }
        else if (!filename.length()) filename = argv[i];
        else atom_names.push_back(argv[i]);
    }

    if (!filename.length() || !atom_names.size())
    {
        show_usage();
        return -1;
    }

    FILE* fp = fopen(filename.c_str(), "r");
    char buffer[65536];
    size_t rfw = fread(buffer, sizeof(char), 65535, fp);
    fclose(fp);

    Molecule m("TheMolecule");
    m.from_sdf(buffer);

    n = atom_names.size();
    for (i=0; i<n; i++)
    {
        Atom* a = m.get_atom(atom_names[i].c_str());
        if (!a)
        {
            cout << "Atom " << atom_names[i] << " not found in molecule." << endl;
            return -1;
        }

        l = a->num_rings();
        if (l < 1)
        {
            cout << "Atom " << atom_names[i] << " is not part of a ring." << endl;
            return -1;
        }
        else if (l > 1)
        {
            cout << "Atom " << atom_names[i] << " belongs to multiple rings." << endl;
            return -1;
        }

        if (a->is_pi())
        {
            cout << "Atom " << atom_names[i] << " is double bonded and cannot flip." << endl;
            return -1;
        }

        Ring** rr = a->get_rings();
        float theta = rr[0]->flip_atom(a);
        cout << "Flipped " << a->name << " " << (theta*fiftyseven) << "deg." << endl;
    }

    m.evolve_structure(500);
    m.minimize_internal_clashes();

    fp = fopen(filename.c_str(), "w");
    m.save_sdf(fp);
    fclose(fp);
}
