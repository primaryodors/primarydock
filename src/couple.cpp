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

    int bw3_50 = gpcr.get_bw50(3),
        bw4_50 = gpcr.get_bw50(4),
        bw7_50 = gpcr.get_bw50(7);

    int c438, c359, c756;

    if (bw3_50 > 0) c359 = bw3_50 + 9;
    else
    {
        int q = gpcr.search_sequence(100, 200, "DRYXAICXPLXY");
        if (q < 1)
        {
            cout << "Cannot find 3.59. Coupling fail." << endl;
            return -1;
        }
        else c359 = q + 10;
    }

    if (bw4_50 > 0) c438 = bw4_50 - 12;
    else c438 = c359 + 6;

    if (bw7_50 > 0) c756 = bw7_50 + 6;
    else
    {
        int q = gpcr.search_sequence(250, 400, "PMLNPLIYSLRNKD");
        if (q < 1)
        {
            cout << "Cannot find 7.56. Coupling fail." << endl;
            return -1;
        }
        else c756 = q + 10;
    }

    cout << "GPCR contact residues: " << c359 << " " << c438 << " " << c756 << endl;
}