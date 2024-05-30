
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <memory>
#include "point.h"

#ifndef _ATOM
#define _ATOM

enum intera_type
{
    covalent,
    ionic,
    hbond,
    pi,
    polarpi,
    mcoord,
    vdW
};

enum bond_rotation_fail_reason
{
    bf_none,
    bf_empty_atom,
    bf_terminal_atom,
    bf_bond_not_found,
    bf_limited_rotation,
    bf_disallowed_rotation,
    bf_sidechain_hierarchy,
    bf_not_connected,
    bf_unknown
};

std::ostream& operator<<(std::ostream& os, const bond_rotation_fail_reason& bf);

class Atom;
class Bond
{
public:
    Atom* atom = 0;
    Atom* btom = 0;
    float cardinality=0;			// aromatic bonds = 1.5.
    bool can_rotate=false;
    bool can_flip=false;
    float flip_angle=0;				// signed.
    float angular_momentum=0;
    float total_rotations=0;
    intera_type type = covalent;
    float optimal_radius = 1;
    bond_rotation_fail_reason last_fail = bf_none;

    #if _debug_active_bond_rot
    bool echo_on_rotate = false;
    #endif

    Bond();
    Bond(Atom* a, Atom* b, int card);
    ~Bond();

    bool rotate(float angle_radians, bool allow_backbone = false, bool skip_inverse_check = false);
    Point ring_rotate(float angle_radians, Atom* stop_at);
    void clear_moves_with_cache()
    {
        moves_with_btom = 0;
    }
    void fetch_moves_with_btom(Atom** result);
    int count_moves_with_btom();
    int count_heavy_moves_with_atom();
    int count_heavy_moves_with_btom();
    Bond* get_reversed();
    void swing(SCoord newdir);		// Rotate btom, and all its moves_with atoms, about atom so that the bond points to newdir.

protected:
    void fill_moves_with_cache();
    void enforce_moves_with_uniqueness();
    Atom** moves_with_btom = 0;
    Bond* reversed = nullptr;
};

enum RING_TYPE
{
    AROMATIC,
    ANTIAROMATIC,
    COPLANAR,
    OTHER,
    UNKNOWN
};

class Ring
{
public:
    Ring() { ; }
    Ring(Atom** from_atoms);
    Ring(Atom** from_atoms, RING_TYPE type);

    int get_atom_count()
    {
        return atcount;
    }
    Atom* get_atom(int index);
    Atom** get_atoms() const;
    int get_overlap_count(Ring* ringb);
    RING_TYPE get_type();
    Point get_center();
    SCoord get_normal();
    LocatedVector get_center_and_normal();
    bool is_coplanar();
    bool is_conjugated();
    bool Huckel();						// Compiler doesn't allow ü in an identifier - boo hiss!
    Atom* traverse_ring(Atom* from, Atom* away_from = nullptr);     // If away_from is null, traverse in either direction.
    float flip_atom(Atom* which_atom);
    void aromatize();

protected:
    Atom* atoms[256];
    int atcount = 0;
    RING_TYPE type = UNKNOWN;

    void fill_with_atoms(Atom** from_atoms);
    void determine_type();
    void make_coplanar();
};

class Atom
{
    friend class Bond;
    friend class Ring;

public:
    // Constructors and destructors.
    Atom(const char* elem_sym);
    Atom(const char* elem_sym, const Point* location);
    Atom(const char* elem_sym, const Point* location, const float charge);
    Atom(FILE* instream);
    ~Atom();

    // Basic getters.
    const char* get_elem_sym();
    int get_Z()
    {
        return Z;
    }
    int get_family()
    {
        return family;
    }
    int get_valence()
    {
        return valence;
    }
    int get_geometry()
    {
        return geometry;
    }
    Point get_location();
    float get_vdW_radius()
    {
        return vdW_rad;
    }
    float get_atomic_weight()
    {
        return at_wt;
    }
    float get_electronegativity()
    {
        return elecn;
    }
    bool is_pKa_near_bio_pH() { return is_imidazole_like; }
    float get_acidbase();
    float get_charge();
    float get_orig_charge() { return origchg; }
    float get_max_conj_charge() { return max_localized_charge; }
    float is_polar();						// -1 if atom is H-bond acceptor; +1 if donor.
    bool is_metal();
    int is_thio();							// -1 if atom is S; +1 if atom is H of a sulfhydryl.
    bool is_pi();
    bool is_amide();
    bool is_aldehyde();

    // Setters.
    void set_aa_properties();
    void set_acidbase(float ab)
    {
        acidbase = ab;
    }
    void clear_all_moves_cache();
    void increment_charge(float lcharge)
    {
        charge += lcharge;
        origchg = charge;
    }

    // Bond functions.
    void fetch_bonds(Bond** result);
    int get_bonded_atoms_count();
    int get_bonded_heavy_atoms_count();
    int get_count_pi_bonds();
    float get_sum_pi_bonds();

    bool bond_to(Atom* btom, float cardinality);
    void unbond(Atom* btom);
    void unbond_all();
    void consolidate_bonds();

    float is_bonded_to(Atom* btom);			// If yes, return the cardinality.
    Atom* is_bonded_to(const char* element);
    Atom* is_bonded_to(const char* element, const int cardinality);
    Atom* is_bonded_to(const int family);
    Atom* is_bonded_to(const int family, const int cardinality);
    Atom* is_bonded_to_pi(const int family, const bool other_atoms_pi);

    int num_bonded_to(const char* element);
    int num_bonded_to_in_ring(const char* element, Ring* member_of);

    bool shares_bonded_with(Atom* btom);
    bool check_Greek_continuity();
    Atom* get_heaviest_bonded_atom_that_isnt(Atom* excluded);

    Bond* get_bond_between(Atom* btom);
    Bond* get_bond_between(const char* bname);
    Bond* get_bond_by_idx(int bidx);
    Bond* get_bond_closest_to(Point target);
    int get_idx_bond_between(Atom* btom);

    float hydrophilicity_rule();

    bool is_conjugated_to(Atom* a, Atom* break_if_reach = nullptr, Atom* caller = nullptr);
    float is_conjugated_to_charge(Atom* break_if_reach = nullptr, Atom* caller = nullptr);
    std::vector<Atom*> get_conjugated_atoms(Atom* break_if_reach = nullptr, Atom* caller = nullptr);

    // Ring membership.
    int num_rings();
    int num_conj_rings();
    Ring** get_rings();
    bool is_in_ring(Ring* ring);
    Ring* closest_arom_ring_to(Point target);
    Ring* in_same_ring_as(Atom* b);
    void aromatize()
    {
        geometry=3;
        // if (valence>3) valence--;
        if (bonded_to)
        {
            int i;
            for (i=0; i<geometry; i++)
            {
                if (bonded_to[i].btom && in_same_ring_as(bonded_to[i].btom))
                {
                    if (bonded_to[i].cardinality > 1
                            ||
                            (	bonded_to[i].cardinality == 1
                                && bonded_to[i].btom->get_Z() > 1
                                && bonded_to[i].btom->get_bonded_atoms_count() < 4
                            )
                    )
                    {
                        bonded_to[i].cardinality = 1.5;
                    }
                }
            }
        }
        geov=0;
    }

    // Serialization
    void save_pdb_line(FILE* pf, unsigned int atomno);
    void stream_pdb_line(ostream& os, unsigned int atomno);

    // Spatial functions.
    bool move(Point* pt);
    bool move(Point pt)
    {
        return move(&pt);
    }
    bool move_rel(SCoord* v);
    int move_assembly(Point* pt, Atom* excluding);			// Return number of atoms moved. Note excluding must be a bonded atom.
    SCoord* get_basic_geometry();
    SCoord* get_geometry_aligned_to_bonds(bool prevent_infinite_loop = false);
    float get_geometric_bond_angle();
    float get_bond_angle_anomaly(SCoord v, Atom* ignore = nullptr);	// Assume v is centered on current atom.
    float distance_to(Atom* btom)
    {
        if (!btom) return -1;
        else return location.get_3d_distance(&btom->location);
    };
    float similarity_to(Atom* btom);
    SCoord get_next_free_geometry(float lcard);
    int get_idx_next_free_geometry();
    void rotate_geometry(Rotation rot);			// Necessary for bond rotation.
    void clear_geometry_cache()
    {
        geov=0;
    }
    void swing_all(int startat=0);
    void swap_chirality()
    {
        swapped_chirality = !swapped_chirality;
        chirality_unspecified = false;
    }

    void print_bond_angles();                   // For unit tests.

    // Static fuctions.
    static int Z_from_esym(const char* elem_sym);
    static const char* esym_from_Z(const int lZ)
    {
        if (!lZ || lZ >= _ATOM_Z_LIMIT) return 0;
        else if (lZ == any_element) return asterisk;
        else return elem_syms[lZ];
    }
    static void dump_array(Atom** aarr);
    static float electronegativity_from_Z(int atom_Z) { return electronegativities[atom_Z]; }

    // Public member vars.
    float pK = nanf("n/a");         // To be managed and used by external classes.
    int pdbidx=0;                   // "
    int residue=0;					// "
    char aaletter;					// "
    char aa3let[4];					// "
    char* region;					// "
    bool is_backbone=false;			// "
    char* name;						// "
    bool used = false;      		// Required for certain algorithms such as Molecule::identify_rings().
    int mirror_geo=-1;				// If >= 0, mirror the geometry of the btom of bonded_to[mirror_geo].
    bool flip_mirror=false;			// If true, do trans rather than cis bond conformation.
    bool dnh=false;					// Do Not Hydrogenate. Used for bracketed atoms in SMILES conversion.
    bool EZ_flip = false;
    float last_bind_energy = 0;
    float strongest_bind_energy = 0;
    Atom* strongest_bind_atom = nullptr;
    float shielding_angle = 0;
    char pdbchain = ' ';
    bool doing_ring_closure = false;
    Conjugation* conjugation = nullptr;

    #if debug_break_on_move
    bool break_on_move = false;		// debugging feature.
    #endif

protected:
    int Z=0;
    Point location;
    int valence=0;
    int geometry=0;						// number of vertices, so 4 = tetrahedral; 6 = octahedral; etc.
    bool geometry_dirty = true;
    int origgeo=0;
    SCoord* geov=0;
    float at_wt = 0;
    float vdW_rad = 0;
    float elecn = 0;
    float Eion = 0;
    float Eaffin = 0;
    float charge = 0;					// can be partial.
    float origchg = 0;
    float max_localized_charge = 0;     // for conjugated charged systems.
    float acidbase = 0;					// charge potential; negative = acid / positive = basic.
    float polarity = 0;					// maximum potential relative to -OH...H-.
    bool polar_calcd = false;
    int thiol = 0;
    Bond* bonded_to = 0;
    bool reciprocity = false;
    int family=0;
    // InteratomicForce** Zforces;			// Non-covalent bond types where the atom's Z = either Za or Zb.
    Rotation geo_rot_1, geo_rot_2;
    bool swapped_chirality = false;
    bool chirality_unspecified = true;
    Ring** member_of = nullptr;
    int recursion_counter = 0;
    bool is_imidazole_like = false;     // Rings having a pKa near the biological pH of 7.4, that aromatic pnictogens can protonate.

    static void read_elements();
    void figure_out_valence();

    static char* elem_syms[_ATOM_Z_LIMIT];
    static float vdW_radii[_ATOM_Z_LIMIT];
    static float electronegativities[_ATOM_Z_LIMIT];
    static float ioniz_energies[_ATOM_Z_LIMIT];
    static float elec_affinities[_ATOM_Z_LIMIT];
    static float atomic_weights[_ATOM_Z_LIMIT];
    static int valences[_ATOM_Z_LIMIT];
    static int geometries[_ATOM_Z_LIMIT];
};

bool atoms_are_conjugated(Atom** atoms);

static bool read_elem_syms = false;

std::ostream& operator<<(std::ostream& os, const Atom& a);
std::ostream& operator<<(std::ostream& os, const Bond& b);
std::ostream& operator<<(std::ostream& os, const Ring& r);




#endif

