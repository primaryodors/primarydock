
#ifndef _CONSTANTS
#define _CONSTANTS

#define fiftyseven (180.0/M_PI)
#define fiftyseventh (M_PI/180)
#define tetrahedral arccos(−1.0/3)
#define triangular (M_PI/1.5)
#define square (M_PI/2);
#define hexagonal (M_PI/3)

#define _kcal_per_kJ 0.239006
#define _kJmol_cuA 0.5
#define _INTERA_R_CUTOFF 20
#define _INTER_TYPES_LIMIT 10
#define BOND_DEF_NOT_FOUND 0xbadb09d

#define ALKMETAL 1
#define ALKEARTH 2
#define TRIEL 3
#define TETREL 4
#define PNICTOGEN 5
#define CHALCOGEN 6
#define HALOGEN 7
#define NOBLEGAS 8
#define HEAVYMETAL 30
#define LANTHANIDE 60
#define ACTINIDE 90

#define VALENCE_EXCEEDS_GEOMETRY 22

// Allow for future expansion of the periodic table.
#define _ATOM_Z_LIMIT 200

#define _def_atc 100
#define _ALLOW_FLEX_RINGS 0
#define _shield_angle 130.0 * M_PI / 180
#define _fullrot_stepdeg 30
#define _fullrot_steprad M_PI/180*_fullrot_stepdeg
#define _fullrot_every 10
#define _def_lin_momentum 0.5
#define _def_ang_momentum _fullrot_steprad/2
#define _def_bnd_momentum _fullrot_steprad/2

#define _voxel_resolution 0.1
#define USE_VOXEL_ARRAY false

#define ATOM_NOT_OF_AMINO_ACID 0x907aa

#define SPHREACH_MAX 1024

#define PROT_MAX_RGN 40

// Torsion angles for various helices and for beta strands.
#define ALPHA_PHI fiftyseventh*-61
#define ALPHA_PSI fiftyseventh*-43
#define BETA_PHI fiftyseventh*-140
#define BETA_PSI fiftyseventh*130
#define _310_PHI fiftyseventh*-49
#define _310_PSI fiftyseventh*-26
#define PI_PHI fiftyseventh*-55
#define PI_PSI fiftyseventh*-70


#endif











