
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
#define vdw_clash_allowance 1.0
#define oxytocin 0.007
#define _DEFAULT_INTERA_R_CUTOFF 8
#define _INTER_TYPES_LIMIT 10
#define BOND_DEF_NOT_FOUND 0xbadb09d

#define any_element -5141
#define Avogadro 6.25e+23

// Give the atoms a sort of lookahead to know what kind of potential binding they could have if only they would rotate properly.
#define intermol_ESP 0.05

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
#define _shield_angle (130.0 * fiftyseventh)
#define _shield_angle_pi (100.0 * fiftyseventh)
#define _can_clash_angle (180.0 * fiftyseventh)
#define _fullrot_stepdeg 30
#define _fullrot_steprad (fiftyseventh*_fullrot_stepdeg)
#define _fullrot_every 10
#define _def_lin_momentum 0.1
#define _def_ang_momentum (_fullrot_steprad/2)
#define _def_bnd_momentum (_fullrot_steprad/2)

#define pi_mult_dkytw 264
#define pi_CH_dkytw 0.0766
#define pi_HT_dkytw 0.002657

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

// Warning - increasing these constants significantly above their stable branch values
// will cause docking fails in the Big Three tests.
#define polar_repulsion 35.0
#define charge_repulsion 120.0

#define lmpush 0.0005

#define _enhanced_pi_stacking 1
#define _preflex_alignment_res 1
#define bb_stochastic 0.15
#define enforce_no_bb_pullaway 1
#define bb_pullaway_allowance 0.13

// If using an activation matrix, active_persistence "rewards" the ligand for keeping
// bindings to the same residues post-activation as pre-activation. The noflex option
// prevents rotating the ligand's bonds in the node immediately after activation.
#define active_persistence 1
#define active_persistence_follow 1
#define active_persistence_limit 16
#define active_persistence_noflex 0
#define active_persistence_ratio 5
#define active_persistence_threshold 5

#define redo_tumble_spheres_on_activation 0
#define redo_tumble_spheres_every_node 1

// Output the activation matrix or the transmembrane regions' active rotations so that
// the viewer can update its cartoon backbone.
// Only one of these switches is actually necessary.
// If both sets of data are given, the viewer will ignore the matrix and use only the
// rotation data.
// Also, currently the code interlaces the two datasets which is not ideal.
#define write_activation_matrix 0
#define write_active_rotation 1

#define _use_gloms 1

#define soft_ligand_importance 20
#define soft_bias_overlap 0.1

#define soft_rock_clash_allowance 0
#define soft_rock_clash_penalty 0.1

// Amount to reduce momenta for path nodes beyond zero. Since the point of path based
// docking is to keep as closely as possible the same ligand pose and move it through
// the protein, we want to minimize the ligand's conformational changes from node to
// node. Linear momenta are not affected by this number. Angular momenta are multiplied
// by this number, and bond rotation momenta by the square of this number.
#define internode_momentum_mult 0.25

#define internode_momentum_only_on_activation 1

// Switches for conformational space search.
#define allow_axial_tumble 1
#define allow_bond_rots 1
#define allow_linear_motion 1
#define monte_carlo_axial 0
#define monte_carlo_flex 1

// Drift pulls the ligand towards the loneliest point if it encounters clashes.
// Turning it off can cause the ligand to be ejected from the protein.
#define allow_drift 1
#define initial_drift 0.333
#define drift_decay_rate 0.08

// Allows full 360 degree whole molecule rotations to search for lower energy configurations.
#define allow_mol_fullrot_iter 1
#define sidechain_fullrot_lig_bmult 3

// Extra weight given to sidechain-ligand binding strengths during conformer search.
#define dock_ligand_bias 0.5

// Turns off the 360 degree rotations for all but the zeroth node of a path.
#define nodes_no_ligand_360_tumble 0
#define nodes_no_ligand_360_flex 0
#define prevent_ligand_360_on_activate 1

// Iteration callback function feature.
#define allow_iter_cb 1

// Allows the ligand to "see through" shielding and anisotropy and seek the greatest
// potential binding with residues, irrespective of whether the side chain is in an
// optimum orientation in space.
#define allow_ligand_esp 1
#define shielding_avoidance_factor 2.5

// Uses the ligand's most strongly bound atom as the center for full molecule rotations,
// instead of the barycenter, to prevent "letting go" of the strongest binding.
#define allow_tethered_rotations 1

// Reject any conformation where any single residue's total clashes exceed this threshold.
#define individual_clash_limit 5

// Overwrite the user supplied pocket center with the loneliest point as determined by
// distances to the nearest residue atoms to the supplied pocket center.
#define pocketcen_is_loneliest 1

// Switches whether the best-binding algorithm is active by default, instead of tumble spheres.
#define default_bestbind 1
#define preemptively_minimize_intermol_clashes 0
#define bestbind_springiness 15

// In the best-binding algorithm, how strongly to associate hydrophilic features of the ligand
// with hydrophilicity of the side chains.
#define hydrophilicity_boost 5
#define best_binding_stochastic 0.3

// For differential docking, whether to multimol_conform() all the protein's residues into an
// optimized initial conformation before adding the ligand.
#define preconform_protein 0
#define default_pre_ligand_multimol_radius 15
#define default_pre_ligand_flex_radius 10
#define pre_ligand_iteration_ratio 1

// Force the ligand to the new pocket center (or loneliest point) each node.
#define recenter_ligand_each_node 0

// Generate an output file called tmp/active.pdb containing the active matrix modified protein.
#define save_active_protein 0

// Whether to count van der Waals interactions towards the potential energies of candidate
// starting poses in tumble spheres.
#define tumble_spheres_include_vdW 0

// Enable the EXCL feature of the config file.
#define use_exclusions 1

// Output PDB data to the command line even if an output file was specified.
#define echo_pdb_data 0

// Auto hydroxy makes geraniol fail in OR1A1. So does pre-rotate side chains.
// Auto hydroxy was supposed to always point polar hydrogens towards nearby H-bond acceptors.
// Prerot side chains was supposed to minimize clashes between the ligand and side chains after
// choosing an initial candidate conformer from tumble spheres.
#define allow_auto_hydroxy 0
#define prerot_sidechains_from_ligand 0

// How much variance to allow in the strongest atom-to-atom binding energy when multimol conforming.
// Decreasing this value requires any positional or rotational change to adhere more tightly to the
// strongest interatomic interaction.
#define strongest_loss_tolerance 0.25
#define _slt1 (1.0 - strongest_loss_tolerance)

// Whether to add the indicated partial charge to a neutral pnictogen, within pKa limits,
// if a negatively charged atom is nearby.
#define _ALLOW_PROTONATE_PNICTOGENS 1
#define pnictogen_partial_protonation 0.2
#define pn_protonation_pKa_min 5

#define prealign_iters 50
#define prealign_momenta_mult 0

// Whether to move water molecules around that are clashing or are not forming intermolecular bonds.
#define _teleport_dissatisfied_waters 1
// Threshold is positive for binding, negative for clashes.
#define _water_satisfaction_threshold 5
#define _water_teleport_tries 25

#define polar_sat_influence_for_dock 30
#define polar_sat_influence_for_bb 30
#define polar_sat_influence_for_scoring 0

// How strong an intermolecular bond is required to prevent a histidine hydrogen from flipping
// to a less favorable state.
#define _hisflip_binding_threshold 25

// Debugging stuff.
#define _bb_maxglom 3
#define _dummy_atoms_for_debug 0

#define _DBG_LONELINESS 0
#define _DBG_STEPBYSTEP 0
#define _DBG_TOOLARGE_DIFFNUMS 0
#define _DBG_TUMBLE_SPHERES 0
#define _DBG_MAX_CLASHES 0
#define output_tumble_debug_docs 0
#define debug_break_on_move 0
#define debug_stop_after_tumble_sphere 0
#define _DORESPHRES 0
#define _DBG_RESBMULT 0
#define _debug_active_bond_rot 0
#define _DBG_SPACEDOUT 0
#define _DBG_H2O_TELEPORT 0
#define _DBG_HISFLIP 0
#define _DBG_MOLBB 0
#define _dbg_bb_rots 0
#define _dbg_bb_pullaway 0
#define _dbg_soft 0
#define _dbg_glomsel 0
#define _dbg_polsat 0
#define _dbg_softrock 0
#define _dbg_rock_pic 0

#endif











