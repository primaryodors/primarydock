#include <cstring>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "classes/protein.h"

// One thousand is suitable for ORs. If desired, this limit can be increased for larger proteins.
#define max_sequence_length 1000

using namespace std;

vector<string> names;
vector<string> sequences;
vector<int[max_sequence_length]> rel_to_first;
int num_seq = 0;

int main(int argc, char** argv)
{
    // Get input filename from command line.
    int l=1;
    string inpfname = "data/sequences.txt";
    if (argv[l]) inpfname = argv[l];

    // Read sequences from input file.
    FILE* pf = fopen(inpfname.c_str(), "rb");
    if (!pf)
    {
        cout << "Failed to open input file " << inpfname << " check file exists and you have permissions." << endl;
        return -1;
    }
    while (!feof(pf))
    {
        char buffer[max_sequence_length+256];
        fgets(buffer, max_sequence_length+256, pf);

        if (buffer[0])
        {
            char** fields = chop_spaced_fields(buffer);
            if (fields[1])
            {
                names.push_back(fields[0]);
                sequences.push_back(fields[1]);
                num_seq++;
            }
        }
    }
    fclose(pf);

    // For each sequence, identify the most similar preexisting sequence.
    // Go letter by letter finding matches of whatever length substrings.

    // Fill rel_to_first with sequence alignments. For example, if 100-110 of sequence 0 matches 120-130 of sequence 1, then
    // rel_to_first[120] thru [130] will contain the values 100 thru 110.
    // If an existing rel_to_first is nonzero, only overwrite it if the new aa is a better match.
    // If the sequence has a deletion, there will be a gap in its numbers (e.g. 231, 232, 234, 235).
    // If it has an insertion, the inserted residue will be zero (e.g. 231, 232, 0, 233, 234).

    // When detect an insertion, renumber all later positions across all alignments.
    // Also try to match them with a different previous sequence.
    // For instance if sequence 3 has 173, 174, 0, 175, then increment all sequences' >=175 and assign 175 to the zero in seq 3.

    // Finally, output the results and return a success.
    return 0;
}