
#include "intera.h"

#ifndef _MOLECULE
#define _MOLECULE

#include <vector>

class SMILES_Parenthetical
{
	public:
    Atom* startsfrom=0;
    char* smilesstr=0;
};

enum MovabilityType
{
    MOV_ALL			= 1000,
    MOV_NORECEN		=  200,
    MOV_FLEXONLY	=   50,
    MOV_NONE		=    0
};

enum MoleculeType
{
	MOLTYP_UNKNOWN,
	MOLTYP_LIGAND,
	MOLTYP_WATER,
	MOLTYP_AMINOACID
};

class Molecule
{
	public:
    Molecule(const char* name);
    Molecule(const char* name, Atom** collection);
    virtual ~Molecule();

    // Load and save functions.
    int from_sdf(const char* sdf_dat);		// returns number of atoms loaded.
    bool save_sdf(FILE* outf);
    bool save_sdf(FILE* outf, Molecule** included_ligands);
    void save_pdb(FILE* outf, int atomno_offset=0);
    int from_pdb(FILE* inf);				// returns number of atoms loaded.
    void identify_acidbase();				// called within every load.
    bool from_smiles(char const * smilesstr);
    float close_loop(Atom** path, float closing_bond_cardinality);

    // Getters.
    const char* get_name() const	{	return name;	}
    int get_atom_count() const	{	return atcount;	}
    int get_bond_count(bool unidirectional) const;
    Atom* get_nearest_atom(Point loc) const;
    Atom* get_nearest_atom(Point loc, intera_type capable_of) const;
    Point get_bounding_box() const;				// Return the +x+y+z vertex of a bounding box, including vdW radii, if center={0,0,0}.

    // Spatial functions.
    Point get_barycenter() const;
    virtual void move(SCoord move_amt);
    virtual void move(Point move_amt);
    virtual void recenter(Point new_location);
    virtual void rotate(SCoord* SCoord, float theta);
    void rotate(LocatedVector SCoord, float theta);
    bool shielded(Atom* a, Atom* b) const;
    Point get_atom_location(char const * aname);
    float correct_structure(int iters = 500);

    // Atom functions.
    Atom* add_atom(const char* elemsym, const char* aname, Atom* bond_to, const float bcard);
    Atom* add_atom(const char* elemsym, const char* aname, const Point* location, Atom* bond_to, const float bcard);
    char** get_atom_names() const;
    Atom* get_atom(const char* aname) const;
    Atom* get_atom(const int a_idx) const	{	return atoms[a_idx];	}
    void delete_atom(Atom* a);
    void hydrogenate(bool steric_only = false);
    void clear_atom_binding_energies();

    // Bond functions.
    Bond** get_rotatable_bonds();
    Bond** get_all_bonds(bool unidirectional);
    void clear_all_bond_caches();					// Call this any time you add or remove an atom.
    bool rotate_bond(const Bond* rot8b, const float angle);

    // Ring functions.
    int identify_rings();
    int get_num_rings()	{	return ringcount;	}
    bool ring_is_coplanar(int ringid);
    bool ring_is_aromatic(int ringid);
    Point get_ring_center(int ringid);
    SCoord get_ring_normal(int ringid);
    Atom** get_ring_atoms(int ringid);
    int get_ring_num_atoms(int ringid);
    bool Huckel(int ringid);						// Compiler doesn't allow ü in an identifier - boo hiss!

    // Interaction functions.
    float get_internal_clashes();
    void minimize_internal_clashes();
    float get_intermol_clashes(Molecule* ligand);
    float get_intermol_clashes(Molecule** ligands);
    float get_intermol_binding(Molecule* ligand);
    float get_intermol_binding(Molecule** ligands);

    static void multimol_conform(Molecule** interactors, int iters = 50, void (*iter_callback)(int) = NULL);

    // Returns the sum of all possible atom-molecule interactions if all distances and anisotropies were somehow optimal.
    float get_atom_mol_bind_potential(Atom* a);


    void intermol_conform(Molecule* ligand, int iters = 50);
    void intermol_conform(Molecule* ligand, int iters, Molecule** avoid_clashing_with);
    void intermol_conform(Molecule* ligand, int iters, AminoAcid** avoid_clashing_with);
    void intermol_conform(Molecule** ligands, int iters = 50);
    void intermol_conform_norecen(AminoAcid** ligands, int iters = 50);
    void intermol_conform(Molecule** ligands, int iters, Molecule** avoid_clashing_with);
    void intermol_conform_norecen(Molecule* ligand, int iters, Molecule** avoid_clashing_with);
    void intermol_conform_norecen(Molecule* ligand, int iters, AminoAcid** avoid_clashing_with);
    void intermol_conform_norecen(Molecule** ligands, int iters, Molecule** avoid_clashing_with);
    void intermol_conform_norecen(Molecule** ligands, int iters, AminoAcid** avoid_clashing_with);
    void intermol_conform_flexonly(Molecule* ligand, int iters, Molecule** avoid_clashing_with);
    void intermol_conform_flexonly(Molecule** ligands, int iters, Molecule** avoid_clashing_with);
    void reset_conformer_momenta();
    Atom** get_most_bindable(int max_num = 3);						// Return the atoms with the greatest potential intermol binding.

    bool echo_iters = false;
    MovabilityType movability = MOV_ALL;
    float lastbind = 0;

	protected:
    Molecule();
    
    Atom** atoms = 0;
    int atcount = 0;
    char* name = 0;
    char* smiles = 0;
    int ringcount = 0;
    Atom*** ring_atoms = 0;
    bool* ring_aromatic = 0;
    Bond** rotatable_bonds = 0;
    bool immobile = false;
    float minclash = 0;					// Baseline computed internal clashes due to unavoidably close atoms.
    bool doing_bkbend = false;

    // For intermol conformer optimization:
    float lmx=0,lmy=0,lmz=0;			// Linear momentum xyz.
    float amx=0,amy=0,amz=0;			// Angular momentum xyz.

    bool from_smiles(char const * smilesstr, Atom* ipreva);
    int smlen = 0;
    SMILES_Parenthetical* paren;
    int spnum = 0;
	MoleculeType mol_typ = MOLTYP_UNKNOWN;

    int aidx(Atom* a);
    void reallocate();
    float fsb_lsb_anomaly(Atom* first, Atom* last, float lcard, float bond_length);
    void make_coplanar_ring(Atom** ring_members, int ringid);
    void recenter_ring(int ringid, Point new_ring_cen);
    void rotate_ring(int ringid, Rotation rot);
    bool in_same_ring(Atom* a, Atom* b);
    float get_atom_error(int atom_idx, LocatedVector* best_lv);

    void intermol_conform_norecen(Molecule* ligand, int iters, Molecule** avoid_clashing_with, float lastbind);
    void intermol_conform_norecen(Molecule** ligands, int iters, Molecule** avoid_clashing_with, float lastbind);
    void intermol_conform_flexonly(Molecule* ligand, int iters, Molecule** avoid_clashing_with, float lastbind);
    void intermol_conform_flexonly(Molecule** ligands, int iters, Molecule** avoid_clashing_with, float lastbind);
};

extern float potential_distance;


#endif

