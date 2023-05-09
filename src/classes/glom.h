
#ifndef _GLOMCLS
#define _GLOMCLS

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

    static std::vector<AtomGlom> get_potential_ligand_gloms(Molecule* mol);
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

    static std::vector<ResidueGlom> get_potential_side_chain_gloms(AminoAcid** aalist);
};

class GlomPair
{
    public:
    AtomGlom* ag;
    ResidueGlom* scg;

    float get_potential();

    static std::vector<GlomPair> pair_gloms(std::vector<AtomGlom> agloms, std::vector<ResidueGlom> scgloms, Point pocketcen);

    protected:
    float potential = 0;
    Point pocketcen;
};

std::ostream& operator<<(std::ostream& os, const AtomGlom& ag);
std::ostream& operator<<(std::ostream& os, const ResidueGlom& scg);

extern std::vector<int> extra_wt;
extern std::vector<MCoord> mtlcoords;

#endif
