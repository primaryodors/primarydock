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

    cout << " _" << endl << "|/\\" << endl << "|//\\                   PrimaryDock application suite" << endl << "|///\\" << endl << " \\///\\           ___" << endl << "  \\///\\   ____ _/      SSS                          PPPP  DDD   BBBB" << endl << "    \\/| _/ /^\\\\       S   S                         P   P D  D  B   B" << endl << "    //// \\_\\_//       S       CC    O   R RR    EE  P   P D   D B   B" << endl << "   / ----v \\_/\\        SSS   C  C  O O  RR  R  E  E PPPP  D   D BBBB" << endl << "  |||----v                S C     O   O R     EEEEE P     D   D B   B" << endl << "   \\ ----v     ~\\\\    S   S  C  C  O O  R      E    P     D  D  B   B" << endl << "    \\\\\\\\_________/|    SSS    CC    O   R       EEE P     DDD   BBBB" << endl << "     \\__|__|__|__|/" << endl << "     " << endl << "" << endl << endl;

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
    fclose(fp);

    AminoAcid* aa = p.get_residue_bw(3, 50);
    if (!aa) aa = p.get_residue_bw(6, 48);
    if (!aa) aa = p.get_residue(100);
    if (!aa) aa = p.get_residue(p.get_end_resno()/2);

    if (aa && aa->get_hydrogen_count() < 3)
    {
        int i, n = p.get_end_resno();
        for (i=0; i<n; i++)
        {
            aa = p.get_residue(i);
            if (aa) aa->hydrogenate();
        }
    }

    m.identify_acidbase();
    m.hydrogenate();

    fp = fopen("tmp/primarysuck.sdf", "w");
    m.save_sdf(fp);
    fclose(fp);

    DockResult dr(&p, &m, Point(10000,10000,10000));
    dr.include_pdb_data = false;
    dr.display_clash_atom1 = true;
    dr.display_clash_atom2 = true;
    cout << dr;
}

