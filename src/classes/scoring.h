
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
    std::string pdbdat;
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
};

extern float init_total_binding_by_type[_INTER_TYPES_LIMIT];
extern float fin_total_binding_by_type[_INTER_TYPES_LIMIT];

extern float* initial_binding;
extern float* initial_vdWrepl;


#endif
