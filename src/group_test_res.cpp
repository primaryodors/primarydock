
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <sys/stat.h>
#include "classes/group.h"

using namespace std;

int main(int argc, char** argv)
{
    Molecule m("Testing");
    Protein p("TheProtein");
    FILE* fp;

    p.add_sequence("AAAAAQGSAYAAVLAAAAAT");
    p.make_helix(1, p.get_seq_length(), M_PI, M_PI);
    p.make_helix(1, p.get_end_resno(), p.get_end_resno(), ALPHA_PHI, ALPHA_PSI);
    p.bridge(8, 10);

    if (argc > 1)
    {
        if (file_exists(argv[1]))
        {
            m.from_sdf(argv[1]);
        }
        else
        {
            m.from_smiles(argv[1]);
        }
    }
    else m.from_smiles("C12=C(C=CN2)C=CC=C1");

    int i, n;

    std::vector<std::shared_ptr<AtomGroup>> ligand_groups = AtomGroup::get_potential_ligand_groups(&m);

    AminoAcid* reach_residues[SPHREACH_MAX];
    Point target = p.get_residue(10)->get_CA_location();
    Point center = p.get_region_center(1, p.get_end_resno());
    SCoord d = target.subtract(center);
    d.r = 3.5;
    target = target.add(d);

    Point size(9,9,9);
    int sphres = p.get_residues_can_clash_ligand(reach_residues, &m, target, size, nullptr);
    cout << "# Reach residues:" << endl;
    cout << "# ";
    for (i=0; i<sphres; i++)
    {
        cout << *reach_residues[i] << " ";
    }
    cout << endl << "# " << endl;

    std::vector<std::shared_ptr<ResidueGroup>> sidechain_groups = ResidueGroup::get_potential_side_chain_groups(reach_residues, target);
    n = sidechain_groups.size();
    cout << "Sidechain groups:" << endl;
    for (i=0; i<n; i++) cout << *sidechain_groups[i] << endl;
    cout << endl;

    std::vector<std::shared_ptr<GroupPair>> pairs = GroupPair::pair_groups(ligand_groups, sidechain_groups, target, 0);
    m.recenter(target);
    GroupPair::align_groups(&m, pairs);
    n = pairs.size();
    cout << "Group pairs:" << endl;
    for (i=0; i<n; i++) cout << *(pairs[i]->ag) << " ~ " << *(pairs[i]->scg) << endl;
    cout << endl;

    fp = fopen("tmp/group_test_res.pdb", "wb");
    p.save_pdb(fp, &m);
    p.end_pdb(fp);
    fclose(fp);
}
