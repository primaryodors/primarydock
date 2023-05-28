
#ifndef _GRPCLS
#define _GRPCLS

#include <memory>
#include <algorithm>
#include "protein.h"

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
    float distance_to(Point pt);
    float bounds();
    float compatibility(AminoAcid* aa);
    bool is_bonded_to(Atom* a);
    Molecule* get_ligand() { return ligand; }
    int intersecting(AtomGroup* compare_with);
    void merge(AtomGroup* merge_with);
    float average_similarity(AtomGroup* compare_with);

    static std::vector<std::shared_ptr<AtomGroup>> get_potential_ligand_groups(Molecule* mol);

    protected:
    Molecule* ligand;
};

class ResidueGroup
{
    public:
    std::vector<AminoAcid*> aminos;
    bool metallic = false;
    Atom* metal = nullptr;

    Point get_center();
    float distance_to(Point pt);
    float compatibility(AtomGroup* ag);
    float group_reach();
    void conform_to(Molecule* mol);

    static std::vector<std::shared_ptr<ResidueGroup>> get_potential_side_chain_groups(AminoAcid** aalist, Point pocketcen);
};

class GroupPair
{
    public:
    std::shared_ptr<AtomGroup> ag;
    std::shared_ptr<ResidueGroup> scg;

    float get_potential();
    bool is_priority() { return priority; }

    static std::vector<std::shared_ptr<GroupPair>> pair_groups(std::vector<std::shared_ptr<AtomGroup>> agroups, std::vector<std::shared_ptr<ResidueGroup>> scgroups, Point pocketcen);
    static void align_groups(Molecule* ligand, std::vector<std::shared_ptr<GroupPair>> group_pairs);    // Assumes the ligand is already centered in the pocket.

    protected:
    float potential = 0;
    Point pocketcen;
    bool priority = false;
};

std::ostream& operator<<(std::ostream& os, const AtomGroup& ag);
std::ostream& operator<<(std::ostream& os, const ResidueGroup& scg);

extern std::vector<MCoord> mtlcoords;
extern std::vector<std::shared_ptr<GroupPair>> global_pairs;

#endif
