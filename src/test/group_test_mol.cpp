
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <sys/stat.h>
#include "../classes/group.h"

using namespace std;

int main(int argc, char** argv)
{
    Molecule m("Testing");

    if (argc > 1)
    {
        struct stat buffer;   
        if (stat(argv[1], &buffer) == 0)
        {
            FILE* fp = fopen(argv[1], "rb");

            fseek(fp, 0L, SEEK_END);
            int size = ftell(fp);
            rewind(fp);

            char contents[size+4];
            fread(contents, size, 1, fp);
            fclose(fp);
            m.from_sdf(contents);
        }
        else m.from_smiles(argv[1]);
    }
    else m.from_smiles("CCO");

    AtomGroup** mgrp = AtomGroup::get_potential_ligand_groups(&m);

    int i;

    for (i=0; i<mgrp[i]; i++)
    {
        cout << *mgrp[i] << endl;
    }
    delete mgrp;
}
