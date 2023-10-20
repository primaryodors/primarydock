
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <sys/stat.h>
#include "classes/moiety.h"

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

    Moiety y;

    if (argc > 2) y.pattern = argv[2];
    else y.pattern = "COH";

    Atom* matches[256];
    int result = y.contained_by(&m, matches);

    cout << "Pattern " << y.pattern << " occurs " << result << " times in molecule." << endl;
    int i;
    for (i=0; matches[i]; i++)
    {
        cout << *(matches[i]) << " ";
    }
    cout << endl;
}
