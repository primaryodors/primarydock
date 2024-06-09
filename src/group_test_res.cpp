
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

    if (argc <= 2)
    {
        cout << "Usage:" << endl
            << "test/group_test_res path/to/protein.pdb path/to/molecule.sdf {pocket residues}" << endl
            << "test/group_test_res path/to/protein.pdb {SMILES string} {pocket residues}" << endl;
        return -1;
    }

    /*
    p.add_sequence("AAAAAQGSAYAAVLAAAAAT");
    p.make_helix(1, p.get_seq_length(), M_PI, M_PI);
    p.make_helix(1, p.get_end_resno(), p.get_end_resno(), ALPHA_PHI, ALPHA_PSI);
    p.bridge(8, 10);
    */

    fp = fopen(argv[1], "rb");
    if (!fp) throw 0xffff;
    p.load_pdb(fp);
    fclose(fp);

    if (file_exists(argv[2]))
    {
        FILE* fp = fopen(argv[2], "rb");

        fseek(fp, 0L, SEEK_END);
        int size = ftell(fp);
        rewind(fp);

        char contents[size+4];
        fread(contents, size, 1, fp);
        fclose(fp);
        m.from_sdf(contents);
    }
    else
    {
        m.from_smiles(argv[2]);
    }

    int i, n, l;

    std::vector<std::shared_ptr<AtomGroup>> ligand_groups = AtomGroup::get_potential_ligand_groups(&m, true);

    AminoAcid* reach_residues[SPHREACH_MAX];
    Point target, center;
    std::vector<Point> pocket;

    for (l=3; l<argc; l++)
    {
        AminoAcid* a = p.get_residue_bw(argv[l]);
        if (!a) continue;
        pocket.push_back(a->get_CA_location());
    }

    if (n = pocket.size())
    {
        target = Point(0,0,0);
        for (l=0; l<n; l++)
        {
            target = target.add(pocket[l]);
        }
        target.multiply(1.0/n);
    }
    else
    {
        target = p.get_residue(10)->get_CA_location();
        center = p.get_region_center(1, p.get_end_resno());
        SCoord d = target.subtract(center);
        d.r = 3.5;
        target = target.add(d);
    }

    Point size(7.5,7.5,7.5);
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
