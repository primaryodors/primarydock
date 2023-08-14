
#include "intera.h"

#ifndef _MOLECULE
#define _MOLECULE

#include <vector>
#include <memory>

struct SMILES_Parenthetical
{
    Atom* startsfrom=0;
    char* smilesstr=0;
};

struct HistidineFlip
{
    Atom* H;
    Atom* C;
    Atom* N1;
    Atom* N2;
};

enum MovabilityType
{
    MOV_CAN_RECEN   =   0x8000,           // Molecule can move through space.
    MOV_CAN_AXIAL   =    0x800,           // Whole molecule can rotate in space.
    MOV_MC_AXIAL    =    0x400,           // Monte Carlo whole molecule rotations allowed.
    MOV_MUST_FLEX   =     0x80,           // Molecule's rotatable bonds are guaranteed free to rotate.
    MOV_MC_FLEX     =     0x40,           // Monte Carlo flexion is allowed.
    MOV_CAN_FLEX    =     0x20,           // Molecule's rotatable bonds may rotate if selected for flexion.
    MOV_ALL			=   0xfff0,
    MOV_NORECEN		=   0x0ff0,
    MOV_NOAXIAL     =   0xf0f0,
    MOV_FORCEFLEX   =     0xf0,
    MOV_FLEXONLY	=     0x70,
    MOV_PINNED      =     0x04,
    MOV_FLXDESEL    =     0x02,
    MOV_FORBIDDEN   =     0x0f,
    MOV_NONE		=     0x00,
    MOV_BKGRND      = 0x020000,
    MOV_CAN_CLASH   = 0x100000
};

enum MoleculeType
{
    MOLTYP_UNKNOWN,
    MOLTYP_LIGAND,
    MOLTYP_WATER,
    MOLTYP_AMINOACID
};

class GroupPair;
class AtomGroup;

class Pose
{
public:
    Pose();
    Pose(Molecule* from_mol);
    ~Pose();
    void copy_state(Molecule* from_mol);
    void restore_state(Molecule* to_mol);
    void reset();

protected:
    int sz = 0;
    Point* saved_atom_locs = nullptr;
    Molecule* saved_from = nullptr;
};


class Molecule
{
    friend class Pose;

public:
    Molecule(const char* name);
    Molecule(const char* name, Atom** collection);
    virtual ~Molecule();

    // Load and save functions.
    int from_sdf(const char* sdf_dat);		// returns number of atoms loaded.
    bool save_sdf(FILE* outf);
    bool save_sdf(FILE* outf, Molecule** included_ligands);
    void save_pdb(FILE* outf, int atomno_offset=0, bool endpdb = true);
    int from_pdb(FILE* inf);				// returns number of atoms loaded.
    void identify_acidbase();				// called within every load.
    bool from_smiles(char const * smilesstr);
    void clear_cache();

    // Getters.
    const char* get_name() const
    {
        return name;
    }
    int get_atom_count() const
    {
        return atcount;
    }
    int get_bond_count(bool unidirectional) const;
    Atom* get_nearest_atom(Point loc) const;
    Atom* get_nearest_atom(Point loc, intera_type capable_of) const;
    Point get_bounding_box() const;				// Return the +x+y+z vertex of a bounding box, including vdW radii, if center={0,0,0}.
    float get_charge() const;
    int is_residue();
    bool is_thiol();

    // Spatial functions.
    Point get_barycenter(bool bond_weighted = false) const;
    virtual void move(SCoord move_amt, bool override_residue = false);
    virtual void move(Point move_amt, bool override_residue = false);
    virtual void recenter(Point new_location);
    void rotate(SCoord* SCoord, float theta, bool bond_weighted = false);
    void rotate(LocatedVector SCoord, float theta);
    bool shielded(Atom* a, Atom* b) const;
    float correct_structure(int iters = 500);
    float close_loop(Atom** path, float closing_bond_cardinality);
    void crumple(float theta);					// Randomly rotate all rotatable bonds by +/- the specified angle.
    float distance_to(Molecule* other_mol);
    std::vector<Atom*> longest_dimension();

    // Atom functions.
    Atom* add_atom(const char* elemsym, const char* aname, Atom* bond_to, const float bcard);
    Atom* add_atom(const char* elemsym, const char* aname, const Point* location, Atom* bond_to, const float bcard);
    void add_existing_atom(Atom* to_add);
    char** get_atom_names() const;
    Atom* get_atom(const char* aname) const;
    Atom* get_atom(const int a_idx) const
    {
        return atoms[a_idx];
    }
    int count_atoms_by_element(const char* esym);
    Point get_atom_location(const char* aname);
    int atom_idx_from_ptr(Atom* a);
    void delete_atom(Atom* a);
    int get_hydrogen_count();
    virtual void hydrogenate(bool steric_only = false);
    void clear_atom_binding_energies();
    int has_hbond_donors();
    int has_hbond_acceptors();                    // N+ is not an h-bond acceptor.

    // Bond functions.
    Bond** get_rotatable_bonds(bool include_heavy = true);
    Bond** get_all_bonds(bool unidirectional);
    void clear_all_bond_caches();					// Call this any time you add or remove an atom.
    bool rotate_bond(const Bond* rot8b, const float angle);
    void do_histidine_flip(HistidineFlip* hf);

    // Ring functions.
    int identify_rings();
    int get_num_rings();
    int add_ring(Atom** atoms);
    bool ring_is_coplanar(int ringid);
    bool ring_is_aromatic(int ringid);
    Point get_ring_center(int ringid);
    SCoord get_ring_normal(int ringid);
    Atom** get_ring_atoms(int ringid);
    int get_ring_num_atoms(int ringid);

    // Interaction functions.
    float get_internal_clashes();
    void minimize_internal_clashes();
    float get_intermol_clashes(Molecule* ligand);
    float get_intermol_clashes(Molecule** ligands);
    float get_intermol_binding(Molecule* ligand, bool subtract_clashes = true);
    float get_intermol_binding(Molecule** ligands, bool subtract_clashes = true);
    float get_intermol_potential(Molecule* ligand, bool disregard_distance = false);
    float get_intermol_potential(Molecule** ligands, bool disregard_distance = false);
    float hydrophilicity();
    float get_intermol_polar_sat(Molecule* ligand);
    float get_intermol_contact_area(Molecule* ligand, bool hydrophobic_only = false);

    float get_vdW_repulsion(Molecule* ligand);

    float bindability_by_type(intera_type type, bool include_backbone = false);

    static float total_intermol_binding(Molecule** ligands);
    static void conform_molecules(Molecule** molecules, int iterations = 50, void (*callback)(int, Molecule**) = nullptr, void (*group_realign)(Molecule*, std::vector<std::shared_ptr<GroupPair>>) = nullptr);
    static void conform_molecules(Molecule** molecules, Molecule** background, int iterations = 50, void (*callback)(int, Molecule**) = nullptr, void (*group_realign)(Molecule*, std::vector<std::shared_ptr<GroupPair>>) = nullptr);
    static void conform_molecules(Molecule** molecules, Molecule** background, Molecule** clashables, int iterations = 50, void (*callback)(int, Molecule**) = nullptr, void (*group_realign)(Molecule*, std::vector<std::shared_ptr<GroupPair>>) = nullptr);
    void conform_atom_to_location(int atom_idx, Point target, int iterations = 50);
    void conform_atom_to_location(char* atom_name, Point target, int iterations = 50);
    SCoord motion_to_optimal_contact(Molecule* ligand);

    // Returns the sum of all possible atom-molecule interactions if all distances and anisotropies were somehow optimal.
    float get_atom_mol_bind_potential(Atom* a);

    float get_springy_bond_satisfaction();

    void reset_conformer_momenta();
    Atom** get_most_bindable(int max_num = 3);						// Return the atoms with the greatest potential intermol binding.
    Atom** get_most_bindable(int max_num, Atom* for_atom);

    void allocate_mandatory_connections(int mcmax);
    void add_mandatory_connection(Molecule* addmol);
    void remove_mandatory_connection(Molecule* rmvmol);
    void zero_mandatory_connection_cache();
    void delete_mandatory_connections();

    // Debug stuff.
    #if debug_break_on_move
    void set_atoms_break_on_move(bool break_on_move)
    {
        if (atoms)
        {
            int i;
            for (i=0; atoms[i]; i++) atoms[i]->break_on_move = break_on_move;
        }
    }
    #endif

    bool echo_iters = false;
    MovabilityType movability = MOV_ALL;
    float lastbind = 0;
    float lastbind_history[10];
    float lastshielded = 0;
    HistidineFlip** hisflips = nullptr;
    Bond* springy_bonds = nullptr;
    int springy_bondct = 0;
    bool been_flexed = false;
    bool priority = false;
    std::vector<std::shared_ptr<GroupPair>> agroups;

protected:
    Molecule();

    Atom** atoms = 0;
    int atcount = 0;
    char* name = 0;
    char* smiles = 0;
    Ring** rings = nullptr;
    Bond** rotatable_bonds = nullptr;
    Bond** rotatable_bonds_nh = nullptr;
    bool immobile = false;
    bool doing_bkbend = false;
    float base_internal_clashes = 0;					// Baseline computed internal clashes due to unavoidably close atoms.
    std::string sdfgen_aboutline = "";
    Molecule** mandatory_connection = nullptr;
    float* last_mc_binding = nullptr;

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
    float intermol_bind_for_multimol_dock(Molecule* othermol, bool allow_clash);

    void intermol_conform_norecen(Molecule* ligand, int iters, Molecule** avoid_clashing_with, float lastbind);
    void intermol_conform_norecen(Molecule** ligands, int iters, Molecule** avoid_clashing_with, float lastbind);
    void intermol_conform_flexonly(Molecule* ligand, int iters, Molecule** avoid_clashing_with, float lastbind);
    void intermol_conform_flexonly(Molecule** ligands, int iters, Molecule** avoid_clashing_with, float lastbind);
    static float cfmol_multibind(Molecule* mol, Molecule** nearby_mols);
};

extern float conformer_momenta_multiplier;
extern float conformer_tumble_multiplier;
extern bool allow_ligand_360_tumble;
extern bool allow_ligand_360_flex;
extern bool wet_environment;
extern float _momentum_rad_ceiling;

#endif

