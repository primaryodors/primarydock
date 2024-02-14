
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
    DockResult(Protein* prot, Molecule* lig, Point size, int* addl_resno = nullptr, int pose = 1, bool differential_dock = false);

    int pose = 0;
    float kJmol = 0;
    float ikJmol = 0;
    char** metric = nullptr;
    float* mkJmol = 0;
    float* imkJmol = 0;
    float* mvdWrepl = 0;
    float* imvdWrepl = 0;
    float ligand_self = 0;
    float worst_energy = 0;
    float worst_nrg_aa = 0;
    Atom* worst_clash_1 = nullptr;
    Atom* worst_clash_2 = nullptr;
    float* residue_clash = nullptr;
    SCoord* res_clash_dir = nullptr;
    float* missed_connections = 0;
    const char** m_atom1_name = nullptr;
    const char** m_atom2_name = nullptr;
    std::string pdbdat;
    std::string isomer;
    std::string softrock;
    std::string miscdata;
    float bytype[_INTER_TYPES_LIMIT];
    float ibytype[_INTER_TYPES_LIMIT];
    float proximity = 0;                    // How far the ligand center is from the node's center.
    #if use_trip_switch
    float tripswitch = 0;                   // Effect of the ligand on the receptor's trip switch.
    #endif
    float polsat = 0;
    float protclash = 0;
    bool do_output_colors = false;
    bool include_pdb_data = true;
    bool display_clash_atom1 = false;
    bool display_clash_atom2 = false;
    float energy_mult = 1;

    bool out_per_res_e = true;
    bool out_per_btyp_e = true;
    float out_itemized_e_cutoff = 0.01;
    bool out_lig_int_e = true;
    bool out_lig_pol_sat = false;
    bool out_prox = false;
    bool out_pro_clash = false;
    bool out_mc = false;
    bool out_vdw_repuls = false;
};

extern float init_total_binding_by_type[_INTER_TYPES_LIMIT];
extern float fin_total_binding_by_type[_INTER_TYPES_LIMIT];

extern float* initial_binding;
extern float* initial_vdWrepl;

std::ostream& operator<<(std::ostream& output, const DockResult& dr);

#endif
