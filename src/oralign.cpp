#include <cstring>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "classes/protein.h"

// One thousand is suitable for ORs. If desired, this limit can be increased for larger proteins.
#define max_sequence_length 1000
#define seek_len 7

using namespace std;

vector<string> names;
vector<string> sequences;
vector<vector<int>> rel_align;
vector<vector<int>> align_sim;
int num_seq = 0;

void increment_since(int lowest_num_to_incr)
{
    int i, j;

    for (i=0; i<num_seq; i++)
    {
        int ilen = sequences[i].length();
        for (j=0; j<=ilen; j++)
        {
            if (rel_align[i][j] >= lowest_num_to_incr) rel_align[i][j]++;
        }
    }
}

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
                rel_align.push_back(std::vector<int>(sequences[num_seq].length() + 4));
                align_sim.push_back(std::vector<int>(sequences[num_seq].length() + 4));
                num_seq++;
            }
        }
    }
    fclose(pf);

    int i, j, k, m, n;

    for (i=0; i<num_seq; i++)
    {
        int ilen = sequences[i].length();
        for (j=0; j<=ilen; j++)
        {
            rel_align[i][j] = -1;
        }
    }

    for (i=0; i<num_seq; i++)
    {
        // For each sequence, identify the most similar preexisting sequence.
        n = 1e6;
        k = -1;
        for (j=0; j<i; j++)
        {
            m = levenshtein_distance(sequences[i], sequences[j]);
            if (m < n)
            {
                n = m;
                k = j;
            }
        }

        int ilen = sequences[i].length();
        int klen = (k>=0) ? sequences[k].length() : 0;

        if (k >= 0)
        {
            cout << names[i] << "'s closest Levenshtein match is " << names[k] << "." << endl;
        }

        // Go letter by letter finding matches of whatever length substrings.
        for (j=0; j<ilen-seek_len; j++)
        {
            if (k < 0)
            {
                rel_align[i][j] = j;
                align_sim[i][j] = 1000;
            }
            else
            {
                SearchResult sr = find_in_sequence(sequences[i].substr(j, seek_len), sequences[k]);

                // Fill rel_align with sequence alignments. For example, if 100-110 of sequence 0 matches 120-130 of sequence 1, then
                // rel_align[120] thru [130] will contain the values 100 thru 110.
                // If an existing rel_align is nonzero, only overwrite it if the new aa is a better match.
                // If the sequence has a deletion, there will be a gap in its numbers (e.g. 231, 232, 234, 235).
                // If it has an insertion, the inserted residue will be zero (e.g. 231, 232, 0, 233, 234).
                for (l = 0; l < seek_len; l++)
                {
                    int ioff = l + j;
                    int koff = l + sr.position;
                    if (ioff >= ilen || koff >= klen) break;

                    if (rel_align[i][ioff] >= 0)
                    {
                        char oldm = sequences[k].at(in_array(rel_align[i][ioff], rel_align[k]));
                        char newm = sequences[k].at(koff);
                        char cmpr = sequences[i].at(ioff);

                        // TODO: Asses which aa is the closer fit, not simply equals.
                        int incumbent = (cmpr == oldm) ? 1 : 0;
                        int challenge = (cmpr == newm) ? 1 : 0;

                        if (challenge > incumbent)
                        {
                            rel_align[i][ioff] = rel_align[k][koff];
                            align_sim[i][ioff] = sr.similarity;
                        }
                    }
                    else
                    {
                        rel_align[i][ioff] = rel_align[k][koff];
                        align_sim[i][ioff] = sr.similarity;
                    }
                }

                // When detect an insertion, renumber all later positions across all alignments.
                // Also try to match them with a different previous sequence.
                // For instance if sequence 3 has 173, 174, 0, 175, then increment all sequences' >=175 and assign 175 to the zero in seq 3.
                if (j > 0)
                {
                    if (rel_align[i][j] < rel_align[i][j-1])
                    {
                        // TODO: Determine which is the better fit and realign.
                        ;
                    }
                    else if (rel_align[i][j] == rel_align[i][j-1])
                    {
                        // Insertion.
                        if (align_sim[i][j] > align_sim[i][j-1])
                            increment_since(rel_align[i][j]);
                        else
                            increment_since(rel_align[i][j-1]);
                    }
                }
            }       // end if k >= 0
        }       // end for j
    }       // end for i

    int max = 0;
    for (i=0; i<num_seq; i++)
    {
        int ilen = sequences[i].length();
        for (j=0; j<=ilen; j++)
        {
            if (rel_align[i][j] > max) max = rel_align[i][j];
        }
    }

    // Finally, output the results and return a success.
    for (n = 0; n < max; n++)
    {
        for (i=0; i<num_seq; i++)
        {
            cout << names[i] << "\t";
            int ilen = sequences[i].length();
            bool yet = false;
            for (j=0; j<ilen; j++)
            {
                if (rel_align[i][j] == n)
                {
                    cout << sequences[i].at(j);
                    yet = (j >= ilen-1);
                }
                else if (yet) cout << "-";
                else cout << " ";
            }
            cout << endl;
        }
    }

    return 0;
}