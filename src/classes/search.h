
#include "protein.h"
#include "group.h"

#ifndef _SEARCH
#define _SEARCH

class Search
{
    public:
    static void prepare_constrained_search(Protein* protein, Molecule* ligand, Point l_pocket_cen);
    static void do_constrained_search(Protein* protein, Molecule* ligand);
    static void do_best_binding(Protein* protein, Molecule* ligand, Point l_pocket_cen, AminoAcid** reaches_spheroid);
    static void do_tumble_spheres(Protein* protein, Molecule* ligand, Point l_pocket_cen);
};

extern Point size, loneliest;
extern std::vector<int> exclusion;

extern AtomGroup ligand_groups[3];
extern ResidueGroup sc_groups[3];

#endif
