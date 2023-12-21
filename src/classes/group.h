
#ifndef _GRPCLS
#define _GRPCLS

#include <memory>
#include <algorithm>
#include "protein.h"
#include "moiety.h"

struct ResiduePlaceholder
{
    int node = 0;
    int resno = 0;
    std::string bw;

    void set(const char* str);
    void resolve_resno(Protein* prot);
};

struct MCoord
{
    int Z = 29;
    int charge = 2;
    Atom* mtl = nullptr;
    std::vector<ResiduePlaceholder> coordres;
};

class AtomGroup
{
    public:
    std::vector<Atom*> atoms;
    Point get_center();
    float get_pi();
    float get_polarity();
    float get_ionic();
    float get_mcoord();
    float get_sum();
    float get_avg_elecn();
    int contains_element(const char* esym);
    bool contains_atom(Atom* a);
    void remove_atom(Atom* a);
    float distance_to(Point pt);
    float bounds();
    bool is_bonded_to(Atom* a);
    Molecule* get_ligand() { return ligand; }
    int intersecting(AtomGroup* compare_with);
    void merge(AtomGroup* merge_with);
    float average_similarity(AtomGroup* compare_with);
    float hydrophilicity();

    int heavy_atom_count();
    static std::vector<std::shared_ptr<AtomGroup>> get_potential_ligand_groups(Molecule* mol, bool separate_metal_coord = false);
    void remove_duplicates();

    protected:
    Molecule* ligand = nullptr;
};

class ResidueGroup
{
    public:
    std::vector<AminoAcid*> aminos;
    bool metallic = false;
    Atom* metal = nullptr;

    Point get_center();
    Atom* get_nearest_atom(Point pt);
    float distance_to(Point pt);
    float group_reach();
    void conform_to(Molecule* mol);
    float hydrophilicity();
    float pi_stackability();

    static std::vector<std::shared_ptr<ResidueGroup>> get_potential_side_chain_groups(AminoAcid** aalist, Point pocketcen);
};

class GroupPair
{
    public:
    std::shared_ptr<AtomGroup> ag;
    std::shared_ptr<ResidueGroup> scg;

    float get_potential();
    float get_weighted_potential();
    bool is_priority() { return priority; }

    static std::vector<std::shared_ptr<GroupPair>> pair_groups(std::vector<std::shared_ptr<AtomGroup>> agroups, std::vector<std::shared_ptr<ResidueGroup>> scgroups, Point pocketcen, float rel_stochasticity = 1);
    static void align_groups(Molecule* ligand, std::vector<std::shared_ptr<GroupPair>> group_pairs);
    static void align_groups_noconform(Molecule* ligand, std::vector<std::shared_ptr<GroupPair>> group_pairs);
    static void align_groups(Molecule* ligand, std::vector<std::shared_ptr<GroupPair>> group_pairs, bool do_conforms, float amount=1);    // Assumes the ligand is already centered in the pocket.

    protected:
    float potential = 0;
    Point pocketcen;
    bool priority = false;
};

std::ostream& operator<<(std::ostream& os, const AtomGroup& ag);
std::ostream& operator<<(std::ostream& os, const ResidueGroup& scg);

extern std::vector<MCoord> mtlcoords;
extern std::vector<std::shared_ptr<GroupPair>> global_pairs;
extern std::vector<Moiety> predef_grp;

#endif
