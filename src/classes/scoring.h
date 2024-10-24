
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include "protein.h"

#ifndef _SCORING
#define _SCORING

class DockResult
{
    public:
    DockResult();
    DockResult(Protein* prot, Molecule* lig, Point size, int* addl_resno = nullptr, int pose = 1, Molecule** waters = nullptr);

    int pose = 0;
    int auth = 0;
    float kJmol = 0;
    float ikJmol = 0;
    char** metric = nullptr;
    float* mkJmol = 0;
    float* imkJmol = 0;
    #if compute_vdw_repulsion
    float* mvdWrepl = 0;
    float* imvdWrepl = 0;
    #endif
    float ligand_self = 0;
    float worst_energy = 0;
    float worst_nrg_aa = 0;
    Atom* worst_clash_1 = nullptr;
    Atom* worst_clash_2 = nullptr;
    #if compute_clashdirs
    float* residue_clash = nullptr;
    SCoord* res_clash_dir = nullptr;
    #endif
    #if compute_missed_connections
    float* missed_connections = 0;
    #endif
    const char** mb_atom1_name = nullptr;
    const char** mb_atom2_name = nullptr;
    const char** mc_atom1_name = nullptr;
    const char** mc_atom2_name = nullptr;
    std::string pdbdat;
    std::string isomer;
    std::string miscdata;
    float bytype[_INTER_TYPES_LIMIT];
    float ibytype[_INTER_TYPES_LIMIT];
    float proximity = 0;                    // How far the ligand center is from the node's center.
    float polsat = 0;
    float protclash = 0;
    float A100 = 0;
    bool do_output_colors = false;
    bool include_pdb_data = true;
    bool display_binding_atoms = false;
    bool display_clash_atoms = false;
    float energy_mult = 1;
    Pose ligpos;

    bool out_per_res_e = true;
    bool out_per_btyp_e = true;
    float out_itemized_e_cutoff = 0.01;
    bool out_lig_int_e = true;
    bool out_lig_pol_sat = false;
    bool out_prox = false;
    bool out_pro_clash = false;
    bool out_mc = false;
    bool out_vdw_repuls = false;

    bool disqualified = false;
};

extern float init_total_binding_by_type[_INTER_TYPES_LIMIT];
extern float fin_total_binding_by_type[_INTER_TYPES_LIMIT];

extern float* initial_binding;
#if compute_vdw_repulsion
extern float* initial_vdWrepl;
#endif

std::ostream& operator<<(std::ostream& output, const DockResult& dr);

#endif
