
#ifndef _CONSTANTS
#define _CONSTANTS

#define fiftyseven (180.0/M_PI)
#define fiftyseventh (M_PI/180)
#define tetrahedral acos(-1.0/3)
#define triangular (M_PI/1.5)
#define square (M_PI/2)
#define hexagonal (M_PI/3)

#define _kcal_per_kJ 0.239006
#define _kJmol_cuA 0.5
#define _INTERA_R_CUTOFF 11
#define _INTER_TYPES_LIMIT 10
#define BOND_DEF_NOT_FOUND 0xbadb09d

#if defined(__linux__) || defined(__sun) || defined(__FreeBSD__) || defined(__NetBSD__) || defined(__OpenBSD__) || defined(__APPLE__)
	#define CMD_CHECK_INSTALLED_3P_SMILES_PARSER "which obabel"
	#define CMD_CALL_3P_SMILES_PARSER "obabel -:'%s' --gen3D -osdf 2> /dev/null"
#elif _WIN32
	#define CMD_CHECK_INSTALLED_3P_SMILES_PARSER "dir \"C:\\Program Files\\OpenBabel*\""
	#define CMD_CALL_3P_SMILES_PARSER "obabel.exe -:'%s' --gen3D -osdf"
#else
	#error It appears your operating system is not supported yet - would you be willing to add it and submit a pull request?
#endif


#define _SANOM_BOND_ANGLE_WEIGHT 25.0
#define _SANOM_CLASHES_WEIGHT 5.0
#define _SANOM_BOND_RAD_WEIGHT 30.0

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
#define _def_lin_momentum 0.1
#define _def_ang_momentum _fullrot_steprad/2
#define _def_bnd_momentum _fullrot_steprad/2

#define _voxel_resolution 0.1
#define USE_VOXEL_ARRAY false

#define ATOM_NOT_OF_AMINO_ACID 0x907aa
#define NOT_ATOM_RECORD 0xb1207e19

#define _MAX_NUM_FORCES 65536

#define SPHREACH_MAX 1024

#define PROT_MAX_RGN 40

#define AMINOACID_HYDROPHILICITY_THRESHOLD 0.2

// Torsion angles for various helices and for beta strands.
#define ALPHA_PHI fiftyseventh*-57.8
#define ALPHA_PSI fiftyseventh*-47.0
#define BETA_PHI fiftyseventh*-140
#define BETA_PSI fiftyseventh*130
#define _310_PHI fiftyseventh*-49
#define _310_PSI fiftyseventh*-26
#define PI_PHI fiftyseventh*-55
#define PI_PSI fiftyseventh*-70
#define POLYPRO2_PHI fiftyseventh*-75
#define POLYPRO2_PSI fiftyseventh*145
#define POLYPRO2_OMEGA fiftyseventh*124
#define POLYPRO1_PHI fiftyseventh*-75
#define POLYPRO1_PSI fiftyseventh*160
#define POLYPRO1_OMEGA fiftyseventh*113


#endif











