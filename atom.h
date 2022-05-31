
#include <iostream>
#include <fstream>
#include <iomanip>
#include "point.h"

#ifndef _ATOM
#define _ATOM

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

    Bond();
    Bond(Atom* a, Atom* b, int card);
    ~Bond();

    bool rotate(float angle_radians, bool allow_backbone = false);
    void clear_moves_with_cache()
    {
        moves_with_btom = 0;
    }
    Atom** get_moves_with_btom();
    int count_moves_with_btom();
    void swing(Vector newdir);		// Rotate btom, and all its moves_with atoms, about atom so that the bond points to newdir.

protected:
    void fill_moves_with_cache();
    Atom** moves_with_btom = 0;
};

class Atom
{
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
    float get_charge();
    float get_acidbase();
    float is_polar();						// -1 if atom is H-bond acceptor; +1 if donor.
    int is_thio();							// -1 if atom is S; +1 if atom is H of a sulfhydryl.
    bool is_metal();
    bool is_pi();

    // Setters.
    void set_acidbase(float ab)
    {
        acidbase = ab;
    }
    void set_aa_properties();
    void clear_all_moves_cache();
    void increment_charge(float lcharge)
    {
        charge += lcharge;
    }

    // Bond functions.
    Bond** get_bonds();
    int get_bonded_atoms_count();
    bool bond_to(Atom* btom, float cardinality);
    void unbond(Atom* btom);
    float is_bonded_to(Atom* btom);			// If yes, return the cardinality.
    Atom* is_bonded_to(const char* element);
    Atom* is_bonded_to(const char* element, const int cardinality);
    Atom* is_bonded_to(const int family);
    bool shares_bonded_with(Atom* btom);
    Bond* get_bond_between(Atom* btom);
    Bond* get_bond_between(const char* bname);
    Bond* get_bond_by_idx(int bidx);
    int get_idx_bond_between(Atom* btom);
    void aromatize()
    {
        geometry=3;
        if (valence>3) valence--;
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
    bool move_rel(Vector* v);
    int move_assembly(Point* pt, Atom* excluding);			// Return number of atoms moved. Note excluding must be a bonded atom.
    Vector* get_basic_geometry();
    Vector* get_geometry_aligned_to_bonds();
    float distance_to(Atom* btom)
    {
        if (!btom) return -1;
        else return location.get_3d_distance(&btom->location);
    };
    Vector get_next_free_geometry(float lcard);
    int get_idx_next_free_geometry();
    void rotate_geometry(Rotation rot);			// Necessary for bond rotation.
    void clear_geometry_cache()
    {
        geov=0;
    }
    void swing_all(int startat=0);

    // Static fuctions.
    static int Z_from_esym(const char* elem_sym);
    static char* esym_from_Z(const int lZ)
    {
        if (!lZ || lZ >= _ATOM_Z_LIMIT) return 0;
        else return elem_syms[lZ];
    }
    static void dump_array(Atom** aarr);

    // Public member vars.
    int residue=0;				// To be managed and used by the AminoAcid class.
    char aaletter;				// "
    char aa3let[4];				// "
    char* region;				// "
    bool is_backbone=false;		// "
    char* name;					// "
    bool used;					// Required for certain algorithms such as Molecule::identify_rings().
    int mirror_geo=-1;			// If >= 0, mirror the geometry of the btom of bonded_to[mirror_geo].
    bool flip_mirror=false;		// If true, do trans rather than cis bond conformation.
    bool dnh=false;				// Do Not Hydrogenate. Used for bracketed atoms in SMILES conversion.
    Point* arom_center=0;
    bool swap_chirality = false;
    bool EZ_flip = false;
    float last_bind_energy=0;

protected:
    void figure_out_valence();
    int Z=0;
    Point location;
    int valence=0;
    int geometry=0;						// number of vertices, so 4 = tetrahedral; 6 = octahedral; etc.
    int origgeo=0;
    Vector* geov=0;
    float at_wt = 0;
    float vdW_rad = 0;
    float elecn = 0;
    float Eion = 0;
    float Eaffin = 0;
    float charge = 0;					// can be partial.
    float acidbase = 0;					// charge potential; negative = acid / positive = basic.
    float polarity = 0;					// maximum potential relative to -OH...H-.
    int thiol = 0;
    Bond* bonded_to = 0;
    bool reciprocity = false;
    int family=0;
    // InteratomicForce** Zforces;			// Non-covalent bond types where the atom's Z = either Za or Zb.
    Rotation geo_rot_1, geo_rot_2;

    static void read_elements();

    static char* elem_syms[_ATOM_Z_LIMIT];
    static float vdW_radii[_ATOM_Z_LIMIT];
    static float electronegativities[_ATOM_Z_LIMIT];
    static float ioniz_energies[_ATOM_Z_LIMIT];
    static float elec_affinities[_ATOM_Z_LIMIT];
    static float atomic_weights[_ATOM_Z_LIMIT];
    static int valences[_ATOM_Z_LIMIT];
    static int geometries[_ATOM_Z_LIMIT];
};

static bool read_elem_syms = false;





#endif

