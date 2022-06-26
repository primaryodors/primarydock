
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
    bool add_residue(const char* pdbdata);
    bool add_sequence(const char* sequence);
    float coordinate_metal(const char* metal_elemsym, const int* coord_residues);		// Returns coordination anomaly.
    void set_clashables();
    void delete_residue(int resno);
    void delete_sidechain(int resno);
    void delete_residues(int startres, int endres);
    void delete_sidechains(int startres, int endres);
    MetalCoord* coordinate_metal(Atom* metal, int residues, int* resnos, char** res_anames);
    void set_region(std::string name, int start, int end);

    // Serialization.
    int load_pdb(FILE* infile);				// Returns number of residues loaded.
    void save_pdb(FILE* outfile);
    void end_pdb(FILE* outfile);

    // Getters.
    int get_seq_length();
    std::string get_sequence();
    int get_start_resno();
    AminoAcid* get_residue(int resno);
    Molecule* metals_as_molecule();
    bool aa_ptr_in_range(AminoAcid* aaptr);
    Region get_region(std::string name);
    int get_region_start(std::string name);
    int get_region_end(std::string name);
    Atom* get_atom(int resno, const char* aname);
    Point get_atom_location(int resno, const char* aname);
    std::vector<std::string> get_remarks(std::string search_for = "");
    std::string get_name() { return std::string(name); }

    // Metrics functions.
    float get_internal_clashes();
    float get_intermol_clashes(const Molecule* ligand);
    float get_intermol_binding(const Molecule* ligand);
    AminoAcid** get_residues_can_clash(int resno);
    int get_residues_can_clash_ligand
    (	AminoAcid** reaches_spheroid,
        const Molecule* ligand,
        const Point nodecen,
        const Point size,
        const int* mcoord_resno
    );
    float get_helix_orientation(int startres, int endres);
    Point get_region_center(int startres, int endres);

    // Motion functions
    void move_piece(int start_res, int end_res, Point new_center);		// After calling this, you should reconnect the broken ends with conform_backbone().
    void rotate_piece(int start_res, int end_res, int align_res, Point align_target, int pivot_res = 0);		// If no pivot res, rotate about the center.
    
    void rotate_backbone(int residue_no, bb_rot_dir direction, float angle);
    void rotate_backbone_partial(int startres, int endres, bb_rot_dir direction, float angle);
    void conform_backbone(int startres, int endres, int iters = 50, bool backbone_atoms_only = false);
    void conform_backbone(int startres, int endres, Atom* a, Point target, int iters = 50);
    void conform_backbone(int startres, int endres, Atom* a1, Point target1, Atom* a2, Point target2, int iters = 50);
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
    std::string name;
    char* sequence=0;
    AminoAcid** residues=0;
    AminoAcid*** res_can_clash = 0;
    Atom** ca=0;
    float* res_reach=0;
    Atom** metals=0;
    int metcount=0;
    Star aaptrmin, aaptrmax;
    MetalCoord** m_mcoord=0;
    Region regions[PROT_MAX_RGN];
    std::vector<string> remarks;

    int* get_residues_in_reach(int resno);
    float get_coord_anomaly(Atom* metal, AminoAcid* coord_res);
    void conform_backbone(int startres, int endres,
                          Atom* a1, Point target1,
                          Atom* a2, Point target2,
                          int iters, bool backbone_atoms_only
                         );
    void mtl_coord_cnf_cb(int iter);
    friend void ext_mtl_coord_cnf_cb(int iter);
};


#endif

