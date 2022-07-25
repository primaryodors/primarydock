
#include "aminoacid.h"

#ifndef _PROTEIN
#define _PROTEIN

#include <string>

void ext_mtl_coord_cnf_cb(int iter);

struct Region
{
    int start=0;
    int end=0;
    std::string name="";
};

class Protein
{
	public:
    // Constructors.
    Protein(const char* name);

    // Build functions.
    bool add_residue(const int resno, const char aaletter);
    bool add_sequence(const char* sequence);
    bool add_residue(const char* pdbdata);
    void set_clashables();
    void delete_residue(int resno);
    void delete_sidechain(int resno);
    void delete_residues(int startres, int endres);
    void delete_sidechains(int startres, int endres);
    MetalCoord* coordinate_metal(Atom* metal, int residues, int* resnos, std::vector<string> res_anames);
    void set_region(std::string name, int start, int end);

    // Serialization.
    int load_pdb(FILE* infile);				// Returns number of residues loaded.
    void save_pdb(FILE* outfile);
    void end_pdb(FILE* outfile);

    // Getters.
    int get_seq_length();
    int get_start_resno();
    std::string get_sequence();
    Molecule* metals_as_molecule();
    AminoAcid* get_residue(int resno);
    Region get_region(std::string name);
    int get_region_end(std::string name);
    int get_region_start(std::string name);
    bool aa_ptr_in_range( AminoAcid* aaptr );
    Atom* get_atom(int resno, const char* aname);
    std::string get_name() { return std::string(name); }
    Point get_atom_location(int resno, const char* aname);
    std::vector<std::string> get_remarks(std::string search_for = "");

    // Metrics functions.
    float get_internal_clashes();
    float get_intermol_clashes(const Molecule* ligand);
    float get_intermol_binding(const Molecule* ligand);
    AminoAcid** get_residues_can_clash(int resno);
    int get_residues_can_clash_ligand
    (	AminoAcid** reaches_spheroid,
        Molecule* ligand,
        const Point nodecen,
        const Point size,
        const int* mcoord_resno
    );
    std::vector<AminoAcid*> get_residues_near(Point pt, float max_distance, bool facing=true);
    Point get_region_center(int startres, int endres);
    float get_helix_orientation(int startres, int endres);
    Point find_loneliest_point(Point search_center, Point spheroid_size);

    // Motion functions
    void move_piece(int start_res, int end_res, Point new_center);		// After calling this, you should reconnect the broken ends with conform_backbone().
    void rotate_piece(int start_res, int end_res, int align_res, Point align_target, int pivot_res = 0);		// If no pivot res, rotate about the center.
    
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
    
    void make_helix(int startres, int endres, float phi, float psi);
    void make_helix(int startres, int endres, int stopat, float phi, float psi);
    float orient_helix
    (	int startres, int endres,						// Boundaries of helix.
        int stopat,										// Last residue to move with helix.
        float angle,									// 0 = horizontal; positive = ascending (+Y) with increasing resno.
        int iterations
    );

    int* mcoord_resnos = NULL;

	protected:
    Atom** ca = nullptr;
    std::string name;
    char* sequence = nullptr;
    AminoAcid** residues = nullptr;
    AminoAcid*** res_can_clash = nullptr;
    float* res_reach = nullptr;
    Atom** metals = nullptr;
    int metcount = 0;
    Star aaptrmin, aaptrmax;
    Region regions[PROT_MAX_RGN];
    std::vector<string> remarks;
    MetalCoord** m_mcoord = nullptr;

    int* get_residues_in_reach(int resno);
    float get_coord_anomaly(Atom* metal, AminoAcid* coord_res);
    friend void ext_mtl_coord_cnf_cb(int iter);
    void mtl_coord_cnf_cb(int iter);
};


#endif

