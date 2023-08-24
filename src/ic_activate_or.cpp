#include <cstring>
#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include <unistd.h>
#include <fstream>
#include <sys/stat.h>
#include <sys/wait.h>
#include "classes/protein.h"
#include "classes/dynamic.h"

using namespace std;

bool file_exists(std::string fname)
{
    struct stat s;
    if (stat(fname.c_str(), &s) == 0) return true;
    else return false;
}

int main(int argc, char** argv)
{
    if (argc < 2) return -1;
    int i, j, l, m, n;

    // One argument, in the form of a protein ID, e.g. OR51E2.
    // If the ID does not conform to the format of OR{family}{subfamily}{member}, return an error code.
    if (argv[1][0] != 'O' || argv[1][1] != 'R') return -2;
    l = 2;
    std::string tmp = argv[1];
    int fam = atoi(tmp.substr(l, 2).c_str());
    if (!fam) return -2;
    l++;
    if (fam >= 10) l++;
    std::string sub = tmp.substr(l, 1);
    l++;
    if (argv[1][l] >= 'A' && argv[1][l] <= 'Z')
    {
        sub = tmp.substr(l-1, 2);
        l++;
    }
    int mem = atoi(&argv[1][l]);
    if (!mem) return -2;


    // Locate and load the .upright.pdb in the folder structure, or throw an error if not found.
    std::string path = (std::string)"pdbs/OR" + std::to_string(fam) + (std::string)"/";
    std::string orid = (std::string)"OR" + std::to_string(fam) + sub + std::to_string(mem);
    std::string in_filename = path + orid + (std::string)".upright.pdb";
    if (!file_exists(in_filename)) return -3;

    Protein p(orid.c_str());
    FILE* fp = fopen(in_filename.c_str(), "rb");
    if (!fp) return -3;
    p.load_pdb(fp);
    fclose(fp);


    // If there is an R at position 6.59, perform a slight helix unwind, enough to make the CA atoms of 6.59 and 6.58 equidistant to 45.53.
    DynamicMotion unwind6(&p);
    AminoAcid* aa6x59 = p.get_residue_bw("6.59");
    AminoAcid* aa6x58 = p.get_residue_bw("6.58");
    AminoAcid* aa45x53 = p.get_residue_bw("45.53");
    if (aa6x59 && aa6x58 && aa45x53 && (aa6x59->get_letter() == 'R'))
    {
        cout << "R6.59" << endl;
        unwind6.start_resno.from_string("6.57");
        unwind6.end_resno.from_string("6.59");
        unwind6.type = dyn_wind;
        unwind6.bias = -13;
        unwind6.apply_absolute(1);

        float delta = 0.1;
        for (i=0; i<10; i++)
        {
            float r58 = aa6x58->distance_to(aa45x53);
            float r59 = aa6x59->distance_to(aa45x53);
            if (r58 > r59) delta = -fabs(delta);
            else delta = fabs(delta);
            unwind6.apply_incremental(delta);
            delta *= 0.8;
            // cout << r58 << ", " << r59 << endl;
        }
    }

    // If there is no Y6.55 or no D/E45.51, measure how far 6.55-6.59 can move toward 45.51-45.54 without clashing. Call it TMR6ex.

    // If there is Y6.55 and D/E45.51 but they are not in contact, measure how far 6.55 must move toward 45.51 to make contact. Call it TMR6ex.

    // Otherwise set TMR6ex to zero.

    // Measure how far 5.43-5.54 can move toward 6.44-6.55 without clashing. Call it TMR5ez.

    // Perform the TMR5ez slide.

    // If TMR6ex is nonzero: 
    
    // Compute the angle to rotate about 6.48 and perform the rotation. Measure how far 56.50 moved; call it TMR6c.

    // Compute the axis and angle to rotate TMR5 about 5.33 to match 56.49 to TMR6c, and perform the rotation.

    // If there is Y5.58 and Y7.53:
    
    // Attempt to bridge 5.58~7.53 and measure the distance necessary to complete the contact. Call it Bridge57.

    // Move the side chain of 6.40 to face 7.53:CA and ensure that 6.40 is not clashing with 7.53.
    
    // If 6.40 is still clashing with 7.53, compute the distance to rotate 6.48 thru 56.50 about 6.48 to eliminate the clash.
    // Then perform the rotation, then compute the rotation of TMR5 about 5.33 to keep up, and perform that rotation.
    // Then re-form the 5.58~7.53 bridge.
    
    // Re-measure Bridge57. If it is nonzero, determine how far 7.53 can move toward 5.58 without clashing with 3.43. Call it TMR7cz.
    
    // Compute and execute a bend of TMR7 at 7.48 to move 7.53 the minimum distance of TMR7cz and Bridge57 toward 5.58.
    
    // Re-measure Bridge57. If it is nonzero, compute and execute a pivot of TMR5 from 5.33 to move 5.58 the rest of the way to make contact with 7.53.
    // Then compute and execute a y-axis rotation of TMR6 to bring 6.28 as far along horizontally as 5.68 moved.
    // Translate the CYT3 region (BW numbers 56.x) to stay with 5.68 and 6.28 as smoothly as possible.

    // If there is R6.59, adjust its side chain to keep pointing inward while avoiding clashes with other nearby side chains.

    // Now save the output file. It will be the same name as the input file except it will end in .icactive.pdb instead of .upright.pdb.

    // Finally, write all BRIDGE, STCR, and FLXR parameters to a constraints file ending in .params instead of .upright.pdb.
}

