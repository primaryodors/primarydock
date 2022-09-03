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

vector<string> sequences;
vector<int[max_sequence_length]> rel_to_first;

int main(int argc, char** argv)
{
    // Get input filename from command line.

    // Read sequences from input file.

    // For the first five sequences, grab random segments from first sequence and search for them in the other sequences.
    // For remaining sequences, identify recurring motifs from the first five and search for those.

    // Fill rel_to_first with sequence alignments. For example, if 100-110 of sequence 0 matches 120-130 of sequence 1, then
    // rel_to_first[120] thru [130] will contain the values 100 thru 110.
    // If the sequence has a deletion, there will be a gap in its numbers (e.g. 231, 232, 234, 235).
    // If it has an insertion, the inserted residue will be zero (e.g. 231, 232, 0, 233, 234).

    // Fill in any remaining unassigned alignments as best as possible.

    // When all sequences are aligned, determine the maximum length required to accommodate all insertions.

    // Renumber all sequences to fill the new longer alignment length. Note this will probably cause sequence 0 to have gaps
    // in its numbering for the first time.

    // Finally, output the results and return a success.
    ;
}