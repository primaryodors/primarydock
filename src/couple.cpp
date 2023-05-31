#include <cstring>
#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include "classes/protein.h"
#include "classes/group.h"

using namespace std;

void show_usage()
{
    cout << "Usage:" << endl;
    cout << "couple path/to/GPCR.pdb path/to/G-protein.pdb";
    cout << endl;
}

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        show_usage();
        return -1;
    }

    Protein gpcr("GPCR");
    Protein gnax("GNAX");
    FILE* fp;

    fp = fopen(argv[1], "rb");
    if (!fp)
    {
        cout << "Please ensure " << argv[1] << " exists and is readable." << endl;
        return -1;
    }
    gpcr.load_pdb(fp);
    fclose(fp);
    cout << "GPCR: loaded " << gpcr.get_seq_length() << " residues." << endl;

    fp = fopen(argv[2], "rb");
    if (!fp)
    {
        cout << "Please ensure " << argv[2] << " exists and is readable." << endl;
        return -1;
    }
    gnax.load_pdb(fp);
    fclose(fp);
    cout << "G-protein: loaded " << gnax.get_seq_length() << " residues." << endl;

    //
}