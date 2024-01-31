
#include "point.h"
#include <vector>

#ifndef _NEIGHBOR
#define _NEIGHBOR

#define block_size 8.0
#define block_max_atoms 1024

class Atom;
class Molecule;

struct Block
{
    int serial;
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
    void remove_atom(Atom* vacating_atom);
    int fetch_atoms_near(Atom** results, const int max_results, const Point location, const int depth = 1);
    int fetch_molecules_near(Molecule** results, const int max_results, const Point location, const int depth = 1);
    static double total_system_energy();
    void set_active_protein(Protein* p);
    void set_active_ligand(Molecule* m);
    void clear_active_neighbors();

    protected:
    // TODO: Want some kind of index on cx, cy, cz values.
    std::vector<Block> blocks;
    Block* get_block_from_location(Point location);
    std::vector<Atom*> actives;
};

extern Neighborhood the_neighborhood;


#endif
