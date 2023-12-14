
#include "aminoacid.h"

#ifndef _PROTEIN
#define _PROTEIN

#include <string>

void ext_mtl_coord_cnf_cb(int iter);

enum region_source
{
    rgn_none,
    rgn_pdb,
    rgn_manual
};

struct Region
{
    int start=0;
    int end=0;
    std::string name="";
};

class BallesterosWeinstein
{
    public:
    int helix_no = 0;
    int member_no = 0;

    BallesterosWeinstein() { ; }
    BallesterosWeinstein(int b, int w) { helix_no = b; member_no = w; }
    BallesterosWeinstein(const char* fromstr) { this->from_string(fromstr); }

    void from_string(const char* inpstr);
};

struct SoftBias
{
    std::string region_name;
    float radial_transform = 0;                 // Motion away from or towards the pocket center.
    float angular_transform = 0;                // Motion towards and away from neighboring helices.
    float vertical_transform = 0;               // Motion in the extracelllar or cytoplasmic direction.
    float helical_rotation = 0;                 // Rotation about the helical axis.
    float radial_rotation = 0;                  // Rotation about the imaginary line to the pocket center.
    float transverse_rotation = 0;              // Rotation about the imaginary line perpendicular to the pocket center.
};

class Protein
{
public:
    // Constructors.
    Protein();
    Protein(const char* name);
    ~Protein();

    // Build functions.
    bool add_residue(const int resno, const char aaletter);
    bool add_sequence(const char* sequence);
    bool add_residue(const char* pdbdata);
    void set_clashables(int resno = -1, bool recursed = false);
    void delete_residue(int resno);
    void delete_sidechain(int resno);
    void delete_residues(int startres, int endres);
    void delete_sidechains(int startres, int endres);
    MetalCoord* coordinate_metal(Atom* metal, int residues, int* resnos, std::vector<string> res_anames);
    void set_region(std::string name, int start, int end);
    void set_bw50(int helixno, int resno);
    void renumber_residues(int startres, int endres, int new_startres);
    bool disulfide_bond(int resno1, int resno2);

    // Serialization.
    void set_name_from_pdb_name(const char* pdb_name);
    int load_pdb(FILE* infile, int resno_offset = 0, char chain = 'A');				// Returns number of residues loaded.
    void save_pdb(FILE* outfile, Molecule* ligand = nullptr);
    void end_pdb(FILE* outfile);
    void revert_to_pdb();

    // Getters.
    int get_seq_length();
    int get_start_resno();
    int get_end_resno();
    std::string get_sequence();
    Molecule* metals_as_molecule();
    int get_metals_count();
    AminoAcid* get_residue(int resno);
    AminoAcid* get_residue(BallesterosWeinstein bw);
    AminoAcid* get_residue_bw(int helixno, int bwno);
    AminoAcid* get_residue_bw(const char* bwno);
    BallesterosWeinstein get_bw_from_resno(int resno);
    Region get_region(std::string name);
    const Region* get_regions() { return regions; }
    int get_region_end(std::string name);
    int get_region_start(std::string name);
    bool aa_ptr_in_range( AminoAcid* aaptr );
    Atom* get_atom(int resno, const char* aname);
    std::string get_name()
    {
        return std::string(name);
    }
    Point get_atom_location(int resno, const char* aname);
    void add_remark(const char* remark);
    void add_remark(std::string new_remark);
    std::vector<std::string> get_remarks(std::string search_for = "");
    int get_bw50(int helixno);
    int search_sequence(const int start_resno, const int end_resno, const char* search_for, const int threshold = -1, int* similarity = nullptr);

    char get_pdb_chain() const { return pdbchain; }
    char set_pdb_chain(char chain);

    // Metrics functions.
    float get_internal_clashes(int start_resno = 0, int end_resno = 0, bool repack = false, int repack_iters = 10);
    float get_rel_int_clashes();
    float get_internal_binding();
    float get_intermol_clashes(Molecule* ligand);
    float get_intermol_binding(Molecule* ligand);
    AminoAcid** get_residues_can_clash(int resno);
    std::vector<AminoAcid*> get_residues_can_clash(int start_resno, int end_resno);
    int get_residues_can_clash_ligand
    (	AminoAcid** reaches_spheroid,
        Molecule* ligand,
        const Point nodecen,
        const Point size,
        const int* addl_resno
    );

    std::vector<AminoAcid*> get_residues_near(Point pt, float max_distance, bool facing=true);
    std::vector<AminoAcid*> get_contact_residues(Protein* other_prot, float contact_distance = 2.5);
    Molecule** all_residues_as_molecules();
    Molecule** all_residues_as_molecules_except(Molecule** mm);
    Point get_region_center(int startres, int endres);
    SCoord get_region_axis(int startres, int endres);
    float get_helix_orientation(int startres, int endres);
    Point find_loneliest_point(Point search_center, Point spheroid_size);
    Point estimate_pocket_size(std::vector<AminoAcid*> ba);
    float binding_to_nearby_residues(int resno);
    void minimize_residue_clashes(int resno);
    float region_can_move(int startres, int endres, SCoord direction, bool repack = false, int ignore_startres = 0, int ignore_endres = 0);
    float region_can_rotate(int startres, int endres, LocatedVector axis, bool repack = false, float extra_clash_allowance = 0, int ignore_startres = 0, int ignore_endres = 0);     // Searches positive theta.
    void region_optimal_positioning(int startres, int endres, SCoord* output_transformation, Rotation* output_rotation, Protein** other_strands = nullptr);
    void set_conditional_basicities();

    // Motion functions
    void upright();
    void move_piece(int start_res, int end_res, Point new_center);
    void move_piece(int start_res, int end_res, SCoord move_amt);
    LocRotation rotate_piece(int start_res, int end_res, int align_res, Point align_target, int pivot_res = 0);		// If no pivot res, rotate about the center.
    LocRotation rotate_piece(int start_res, int end_res, Rotation rot, int pivot_res);
    LocRotation rotate_piece(int start_res, int end_res, Point origin, SCoord axis, float theta);

    void rotate_backbone(int residue_no, bb_rot_dir direction, float angle);
    void conform_backbone(int startres, int endres, Atom* a, Point target, int iters = 50);
    void rotate_backbone_partial(int startres, int endres, bb_rot_dir direction, float angle);
    void conform_backbone(int startres, int endres, int iters = 50, bool backbone_atoms_only = false);
    void conform_backbone(int startres, int endres, Atom* a1, Point target1, Atom* a2, Point target2, Atom* a3, Point target3, int iters = 50);
    void conform_backbone(int startres, int endres, Atom* a1, Point target1, Atom* a2, Point target2, int iters = 50, bool backbone_atoms_only = false);
    void conform_backbone(int startres, int endres,
                          Atom* a1, Point target1,
                          Atom* a2, Point target2,
                          Atom* a3, Point target3,
                          int iters, bool backbone_atoms_only
                         );

    void backconnect(int startres, int endres);
    void find_residue_initial_bindings();
    void undo();

    void make_helix(int startres, int endres, float phi, float psi);
    void make_helix(int startres, int endres, int stopat, float phi, float psi);
    float orient_helix
    (	int startres, int endres,						// Boundaries of helix.
        int stopat,										// Last residue to move with helix.
        float angle,									// 0 = horizontal; positive = ascending (+Y) with increasing resno.
        int iterations
    );

    void homology_conform(Protein* target_structure, Protein* reference_structure);
    void bridge(int resno1, int resno2);
    void soft_iteration(std::vector<Region> l_soft_rgns, Molecule* ligand = nullptr);
    int replace_side_chains_from_other_protein(Protein* other);

    int* mcoord_resnos = NULL;

    SCoord last_uprighted_xform;
    LocRotation last_uprighted_A, last_uprighted_B;
    SCoord last_int_clash_dir;

    AminoAcid *stop1, *stop2;
    Atom *stop1a, *stop2a;
    int last_saved_atom_number = 0;

protected:
    Atom** ca = nullptr;
    int arrlimit = 0;
    std::string name;
    char* sequence = nullptr;
    AminoAcid** residues = nullptr;
    AminoAcid*** res_can_clash = nullptr;
    float* res_reach = nullptr;
    Atom** metals = nullptr;
    int metcount = 0;
    Star aaptrmin, aaptrmax;
    float initial_int_clashes = 0;
    Region regions[PROT_MAX_RGN];
    region_source regions_from = rgn_none;
    char** remarks = nullptr;
    int remarksz = 0;
    MetalCoord** m_mcoord = nullptr;
    int Ballesteros_Weinstein[79];
    std::vector<AABridge> aabridges;
    std::vector<Bond*> connections;
    std::vector<Pose> origpdb_residues;
    char pdbchain = ' ';
    Pose** undo_poses = nullptr;
    bool mass_undoable = false;

    int* get_residues_in_reach(int resno);
    float get_coord_anomaly(Atom* metal, AminoAcid* coord_res);
    friend void ext_mtl_coord_cnf_cb(int iter);
    void mtl_coord_cnf_cb(int iter);
    void allocate_undo_poses();
    void save_undo_state();
};

extern float *g_rgnxform_r, *g_rgnxform_theta, *g_rgnxform_y, *g_rgnrot_alpha, *g_rgnrot_w, *g_rgnrot_u;


#endif

