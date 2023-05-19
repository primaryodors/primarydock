
#ifndef _GLOMCLS
#define _GLOMCLS

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

class AtomGlom
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

    static std::vector<std::shared_ptr<AtomGlom>> get_potential_ligand_gloms(Molecule* mol);
};

class ResidueGlom
{
    public:
    std::vector<AminoAcid*> aminos;
    bool metallic = false;
    Atom* metal = nullptr;

    Point get_center();
    float distance_to(Point pt);
    float compatibility(AtomGlom* ag);
    float glom_reach();

    static std::vector<std::shared_ptr<ResidueGlom>> get_potential_side_chain_gloms(AminoAcid** aalist, Point pocketcen);
};

class GlomPair
{
    public:
    std::shared_ptr<AtomGlom> ag;
    std::shared_ptr<ResidueGlom> scg;

    float get_potential();

    static std::vector<std::shared_ptr<GlomPair>> pair_gloms(std::vector<std::shared_ptr<AtomGlom>> agloms, std::vector<std::shared_ptr<ResidueGlom>> scgloms, Point pocketcen);
    static void align_gloms(Molecule* ligand, std::vector<std::shared_ptr<GlomPair>> glom_pairs);    // Assumes the ligand is already centered in the pocket.

    protected:
    float potential = 0;
    Point pocketcen;
};

std::ostream& operator<<(std::ostream& os, const AtomGlom& ag);
std::ostream& operator<<(std::ostream& os, const ResidueGlom& scg);

extern std::vector<int> extra_wt;
extern std::vector<MCoord> mtlcoords;
extern std::vector<std::shared_ptr<GlomPair>> gp;

#endif
