
#ifndef _GRPCLS
#define _GRPCLS

#include <memory>
#include <algorithm>
#include "protein.h"
#include "moiety.h"

class AtomGroup
{
    public:
    Atom* atoms[256];
    int atct = 0;
    Point get_center();
    float get_pi();
    float get_polarity();
    float get_ionic();
    float get_mcoord();
    float get_sum();
    float get_avg_elecn();
    float max_potential_binding(intera_type type);
    int contains_element(const char* esym);
    bool contains_atom(Atom* a);
    Atom* get_mcoord_atom();
    void remove_atom(Atom* a);
    float distance_to(Point pt);
    float bounds();
    bool is_bonded_to(Atom* a);
    Molecule* get_ligand() { return ligand; }
    void update_atom_pointers(Molecule* new_ligand);
    int intersecting(AtomGroup* compare_with);
    void merge(AtomGroup* merge_with);
    float average_similarity(AtomGroup* compare_with);
    float hydrophilicity();
    bool has_hbond_acceptors();
    bool has_hbond_donors();

    int heavy_atom_count();
    static AtomGroup** get_potential_ligand_groups(Molecule* mol, bool separate_metal_coord = false);
    static AtomGroup** make_hbond_subgroups(AtomGroup* from_group);
    void remove_duplicates();

    protected:
    Molecule* ligand = nullptr;
};

class ResidueGroup
{
    public:
    AminoAcid* aminos[16];
    int naminos = 0;
    bool metallic = false;
    Atom* metal = nullptr;

    Point get_center();
    Atom* get_nearest_atom(Point pt);
    float distance_to(Point pt);
    float group_reach();
    void conform_to(Molecule* mol);
    float hydrophilicity();
    float pi_stackability();

    static ResidueGroup** get_potential_side_chain_groups(AminoAcid** aalist, Point pocketcen);

    static AminoAcid* disqualified_residues[1024];
    static int ndisreq;
};

class GroupPair
{
    public:
    AtomGroup* ag;
    ResidueGroup* scg;

    float get_potential();
    float get_weighted_potential();
    bool is_priority() { return priority; }

    static GroupPair** pair_groups(AtomGroup** agroups, ResidueGroup** scgroups, Point pocketcen, float rel_stochasticity = 1);
    static void align_groups(Molecule* ligand, GroupPair** group_pairs);
    static void align_groups_noconform(Molecule* ligand, GroupPair** group_pairs);
    static void align_groups(Molecule* ligand, GroupPair** group_pairs, bool do_conforms, float amount=1);    // Assumes the ligand is already centered in the pocket.

    void disqualify();

    protected:
    float potential = 0;
    Point pocketcen;
    bool priority = false;
};

std::ostream& operator<<(std::ostream& os, const AtomGroup& ag);
std::ostream& operator<<(std::ostream& os, const ResidueGroup& scg);

extern MCoord mtlcoords[16];
extern int nmetals;
extern GroupPair* global_pairs[256];
extern int nglobal_pairs;
extern Moiety predef_grp[65536];
extern int npredef_grp;

#endif
