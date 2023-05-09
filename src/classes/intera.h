
#include "atom.h"
#include <vector>

#ifndef _INTERATOMIC
#define _INTERATOMIC

class InteratomicForce
{
public:
    InteratomicForce();
    intera_type get_type()
    {
        return type;
    }
    float get_arity()
    {
        return arity;
    }
    float get_distance()
    {
        return distance;
    }
    float get_kJmol()
    {
        return kJ_mol;
    }
    float get_dp()
    {
        return dirprop;
    }

    std::string get_config_string() const;

    static bool atom_is_capable_of(Atom* a, intera_type t);
    static InteratomicForce** get_applicable(Atom* a, Atom* b);
    static float potential_binding(Atom* a, Atom* b);
    static float total_binding(Atom* a, Atom* b);
    static float distance_anomaly(Atom* a, Atom* b);
    static float covalent_bond_radius(Atom* a, Atom* b, float cardinality);
    static float coordinate_bond_radius(Atom* a, Atom* b, intera_type btype);

protected:
    int Za=0;
    int bZa=0;
    int Zb=0;
    int bZb=0;
    intera_type type=vdW;
    float arity=0;
    int aritybZa=0;
    int aritybZb=0;
    float distance=0;
    float kJ_mol=0;
    float dirprop=0;		// Directional Propensity.

    void read_dat_line(char* line);
    static void read_all_forces();

    friend std::ostream& operator<<(std::ostream& os, const InteratomicForce& iff);
};

SCoord* get_geometry_for_pi_stack(SCoord* in_geometry);

static InteratomicForce intertmp;
static InteratomicForce* all_forces[_MAX_NUM_FORCES];
static InteratomicForce** forces_by_Z[36][36];				// Good for H thru Br. Includes BCNOFPSKNaClBrMgCaFeCuZn.
static bool read_forces_dat = false;
static bool reading_forces = false;
std::ostream& operator<<(std::ostream& os, const intera_type& it);
std::ostream& operator<<(std::ostream& os, const InteratomicForce& f);
extern float total_binding_by_type[_INTER_TYPES_LIMIT];

#if active_persistence
extern int active_persistence_resno[active_persistence_limit];
#endif

std::ostream& operator<<(std::ostream& os, const InteratomicForce& iff);

#if _peratom_audit
extern std::vector<std::string> interaudit;
extern bool interauditing;
#endif

#endif






