
#include "point.h"
#include <vector>

#ifndef _NEIGHBOR
#define _NEIGHBOR

#define block_size 8.0
#define block_max_atoms 256

class Atom;
class Molecule;

struct Block
{
    int cx, cy, cz;
    Atom* atoms[block_max_atoms];
    int atom_count = 0;
    Molecule* mols[block_max_atoms];
    int mol_count = 0;
};

class Neighborhood
{
    public:
    void add_atom(Atom* neighborly_atom);
    void update_atom(Atom* neighborly_atom, Point old_location);
    int fetch_atoms_near(Atom** results, const int max_results, const Point location, const int depth = 1);
    int fetch_molecules_near(Molecule** results, const int max_results, const Point location, const int depth = 1);

    protected:
    // TODO: Want some kind of index on cx, cy, cz values.
    std::vector<Block> blocks;
    Block* get_block_from_location(Point location);
};

extern Neighborhood the_neighborhood;

#endif
