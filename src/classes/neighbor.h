
#include "point.h"
#include <vector>

#ifndef _NEIGHBOR
#define _NEIGHBOR

#define block_max_atoms 8192

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
    Block* next_door[26];
    int nxdrct = 0;
};

class Neighborhood
{
    public:
    void add_atom(Atom* neighborly_atom);
    void update_atom(Atom* neighborly_atom, Point old_location);
    void remove_atom(Atom* vacating_atom);
    int fetch_atoms_near(Atom** results, const int max_results, const Point location, const int depth = 1);
    int fetch_molecules_near(Molecule** results, const int max_results, const Point location, const int depth = 1);
    double total_system_energy();
    void set_initial_energy();
    double total_energy_delta();
    double total_atom_energy(Atom* atom);
    double total_molecule_energy(Molecule* mol);
    void set_active_protein(Protein* p);
    void set_active_ligand(Molecule* m);
    void clear_active_neighbors();

    double get_worst_clash() { return worst_neighbor_clash; }
    void output_worst_clash(std::ostream& os);

    protected:
    // TODO: Want some kind of index on cx, cy, cz values.
    Block* blocks[65536];
    int blockqty = 0;
    Block* get_block_from_location(Point location);
    std::vector<Atom*> actives;
    double initial_energy = 0;
    Atom *worst_clash_1, *worst_clash_2;
    double worst_neighbor_clash = 0;
};

extern Neighborhood the_neighborhood;


#endif
