#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include "classes/scoring.h"

using namespace std;

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        cout << "Usage:" << endl << endl;
        cout << "score_pdb path/to/filename.pdb [{strand}]" << endl;
        return -1;
    }

    Protein p("TheProtein");
    Molecule m("TheLigand");
    FILE* fp = fopen(argv[1], "rb");
    if (!fp)
    {
        cout << "Please check input file exists and is readable." << endl;
        return -2;
    }
    char strand = 'A';
    if (argc > 2) strand = argv[2][0] & 0xdf;

    p.load_pdb(fp, 0, strand);
    fseek(fp, 0, 0);
    m.from_pdb(fp, true);

    DockResult dr(&p, &m, Point(10000,10000,10000));
}

