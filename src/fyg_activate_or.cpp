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

#define dbg_57_contact 0



void save_file(Protein& p, std::string filename, Molecule* ligand = nullptr)
{
    FILE* fp = fopen(filename.c_str(), "wb");
    if (!fp) throw -3;
    p.save_pdb(fp, ligand);
    p.end_pdb(fp);
    fclose(fp);
    cout << "Saved " << filename << endl;
}

int main(int argc, char** argv)
{
    if (argc < 2) return -1;
    int i, j, l, m, n;

    std::vector<std::string> constraints;

    ////////////////////////////////////////////////////////////////////////////////
    // Read program arguments
    ////////////////////////////////////////////////////////////////////////////////

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


    ////////////////////////////////////////////////////////////////////////////////
    // Load protein
    ////////////////////////////////////////////////////////////////////////////////

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


    ////////////////////////////////////////////////////////////////////////////////
    // Set up residue vars
    ////////////////////////////////////////////////////////////////////////////////

    AminoAcid *aa1x32 = p.get_residue_bw("1.32");
    AminoAcid *aa1x50 = p.get_residue_bw("1.50");
    AminoAcid *aa1x58 = p.get_residue_bw("1.58");

    AminoAcid *aa2x38 = p.get_residue_bw("2.38");
    AminoAcid *aa2x50 = p.get_residue_bw("2.50");
    AminoAcid *aa2x66 = p.get_residue_bw("2.66");

    AminoAcid *aa3x21 = p.get_residue_bw("3.21");
    AminoAcid *aa3x50 = p.get_residue_bw("3.50");
    AminoAcid *aa3x56 = p.get_residue_bw("3.56");

    AminoAcid *aa4x60 = p.get_residue_bw("4.60");

    AminoAcid *aa45x51 = p.get_residue_bw("45.51");
    AminoAcid *aa45x52 = p.get_residue_bw("45.52");
    AminoAcid *aa45x53 = p.get_residue_bw("45.53");
    AminoAcid *aa45x54 = p.get_residue_bw("45.54");

    AminoAcid *aa5x33 = p.get_residue_bw("5.33");
    AminoAcid *aa5x50 = p.get_residue_bw("5.50");
    AminoAcid *aa5x58 = p.get_residue_bw("5.58");
    AminoAcid *aa5x68 = p.get_residue_bw("5.68");

    AminoAcid *aa56x50 = p.get_residue_bw("56.50");

    AminoAcid *aa6x28 = p.get_residue_bw("6.28");
    AminoAcid *aa6x40 = p.get_residue_bw("6.40");
    AminoAcid *aa6x48 = p.get_residue_bw("6.48");
    AminoAcid *aa6x49 = p.get_residue_bw("6.49");
    AminoAcid *aa6x55 = p.get_residue_bw("6.55");
    AminoAcid *aa6x59 = p.get_residue_bw("6.59");

    AminoAcid *aa7x31 = p.get_residue_bw("7.31");
    AminoAcid *aa7x49 = p.get_residue_bw("7.49");
    AminoAcid *aa7x53 = p.get_residue_bw("7.53");
    
    AminoAcid *aa8x44 = p.get_residue_bw("8.44");

    int n1x32 = aa1x32->get_residue_no();
    int n1x50 = aa1x50->get_residue_no();
    int n1x58 = aa1x58->get_residue_no();

    int n2x38 = aa2x38->get_residue_no();
    int n2x50 = aa2x50->get_residue_no();
    int n2x66 = aa2x66->get_residue_no();

    int n3x21 = aa3x21->get_residue_no();
    int n3x50 = aa3x50->get_residue_no();
    int n3x56 = aa3x56->get_residue_no();

    int n4x60 = aa4x60->get_residue_no();
    
    int n45x51 = aa45x51->get_residue_no();
    int n45x52 = aa45x52->get_residue_no();
    int n45x53 = aa45x53->get_residue_no();
    int n45x54 = aa45x54->get_residue_no();

    int n5x33 = aa5x33->get_residue_no();
    int n5x50 = aa5x50->get_residue_no();
    int n5x58 = aa5x58->get_residue_no();
    int n5x68 = aa5x68->get_residue_no();

    int n56x50 = aa56x50->get_residue_no();

    int n6x28 = aa6x28->get_residue_no();
    int n6x40 = aa6x40->get_residue_no();
    int n6x48 = aa6x48->get_residue_no();
    int n6x49 = aa6x49->get_residue_no();
    int n6x55 = aa6x55->get_residue_no();
    int n6x59 = aa6x59->get_residue_no();

    int n7x31 = aa7x31->get_residue_no();
    int n7x49 = aa7x49->get_residue_no();
    int n7x53 = aa7x53->get_residue_no();

    int n8x44 = aa8x44->get_residue_no();

    char l45x51 = aa45x51->get_letter();
    char l45x52 = aa45x52->get_letter();
    char l45x53 = aa45x53->get_letter();

    char l5x50 = aa5x50->get_letter();
    char l5x58 = aa5x58->get_letter();

    char l6x48 = aa6x48->get_letter();
    char l6x55 = aa6x55->get_letter();
    char l6x59 = aa6x59->get_letter();
    
    char l7x53 = aa7x53->get_letter();


    // TODO:

}
