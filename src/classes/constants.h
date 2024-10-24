
#ifndef _CONSTANTS
#define _CONSTANTS

#define fiftyseven (180.0/M_PI)
#define fiftyseventh (M_PI/180)
#define tetrahedral acos(-1.0/3)
#define triangular (M_PI/1.5)
#define square (M_PI/2)
#define hexagonal (M_PI/3)
#define circle M_PI*2

#define _kcal_per_kJ 0.239006
// #define coplanar_threshold 0.5
#define coplanar_threshold 2.5
#define default_cavity_stuffing 0.03
#define _DEFAULT_INTERA_R_CUTOFF 6
#define _INTER_TYPES_LIMIT 10
#define BOND_DEF_NOT_FOUND 0xbadb09d

#define any_element -5141
#define Avogadro 6.02214076e+23

#define helix_hbond_cutoff 2.8
#define hx_tight_cutoff 3.2
#define peptide_bond_length 1.32
#define unconnected_residue_mindist 4.82
#define clash_limit_per_aa 7.0
#define clash_limit_per_atom 3.5
#define iteration_additional_clash_coefficient 10.0
#define contact_r_5x58_7x53 3.93

#define include_eclipses 1
#define include_residue_eclipses 0

#define memsanity 0x10000000

#define pH 6.0
#define auto_pK_protonation 0
#define hydrophilicity_cutoff 0.25
#define group_simil_threshold 0.83

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


#define _SANOM_BOND_RADIUS_WEIGHT 50.0
#define _SANOM_BOND_ANGLE_WEIGHT 10.0
#define _SANOM_CLASHES_WEIGHT 10.0
#define _SANOM_BOND_RAD_WEIGHT 30.0

#define _evolution_atom_displacement 1.0
#define _default_mutation_rate 0.05
#define _default_population_size 80
#define _evolution_aromatic_rigidity 50
#define _evolution_default_generations 5000

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
#define _ALLOW_FLEX_RINGS 1
#define _shield_angle (130.0 * fiftyseventh)
#define _shield_angle_pi (100.0 * fiftyseventh)
#define _can_clash_angle (180.0 * fiftyseventh)
#define _allow_fullrot 0
#define _fullrot_stepdeg 20
#define _fullrot_steprad (fiftyseventh*_fullrot_stepdeg)
#define _fullrot_every 7
#define _def_lin_momentum 0.1
#define _def_ang_momentum (fiftyseventh*5)
#define _def_bnd_momentum (fiftyseventh*15)
#define speed_limit 0.666
#define warn_orphan_atoms 0

#define pi_mult_dkytw 264
#define pi_CH_dkytw 0.0766
#define pi_HT_dkytw 0.002657

#define _voxel_resolution 0.1
#define USE_VOXEL_ARRAY false

#define ATOM_NOT_OF_AMINO_ACID 0x907aa
#define NOT_ATOM_RECORD 0xb1207e19

#define _MAX_NUM_FORCES 65536

#define SPHREACH_MAX 256

#define PROT_MAX_RGN 40

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

// Warning - increasing these constants significantly above the original 35.0, 60.0 values
// will cause docking fails in the unit tests.
#define polar_repulsion 10.0
#define charge_repulsion 1.0

#define atom_object_forces 1
#define charge_attraction 60.0

#define global_clash_allowance 0.4
#define double_hydrogen_clash_allowance_multiplier 1.5
#define ignore_double_hydrogen_clashes 0
#define ignore_nonpolar_hydrogen_clashes 0
#define Lennard_Jones_epsilon 1.0
#define Lennard_Jones_epsilon_x4 Lennard_Jones_epsilon*4
#define lmpush 1.0
#define lmpull 1.0
#define lmsteps 3
#define recapture_ejected_ligand 0
#define summed_missed_connections 1

#define amide_zwitterionic_amount 0.1

#define priority_weight_group 10
#define ts_priority_coefficient 10

#define cs_max_tries 10
#define cs_bondweight_exponent 2.0
#define cs_360_step M_PI / 6.7
#define cs_ligand_rotation 0.333
#define cs_keep_ratio 0.8
#define iter_lostreturns_threshold 0.05
#define max_iters_without_ligand_change 5
#define _enhanced_pi_stacking 0
#define _preflex_alignment_res 1
#define bb_stochastic 0.25
#define bb_group_distance_cutoff 4.3
#define bb_enable_residue_disqualifications 0
#define bb_disqualification_energy 1000
#define enforce_no_bb_pullaway 0
#define bb_pullaway_allowance 1.0
#define bb_realign_amount 0.1
#define bb_realign_mol 0
#define bb_realign_flexions 1
#define bb_realign_threshold_distance 2.8
#define bb_realign_threshold_angle (fiftyseventh * 22.5)
#define enable_bb_scooch 1
#define bb_scooch_threshold_distance 1.5
#define flexion_selection 1
#define no_zero_flexions 1
#define flexion_probability_multiplier 1.0
#define ignore_invalid_partial 1
#define use_best_binding_iteration 0
#define compute_missed_connections 1

// Mandatory coordination bindings that meet this kJ/mol threshold may freely fluctuate above it.
#define mandatory_coordination_threshold 5

// Obsolete stuff.
#define compute_vdw_repulsion 0
#define compute_clashdirs 0

#define redo_tumble_spheres_on_activation 0
#define redo_tumble_spheres_every_node 0

// Output the activation matrix or the transmembrane regions' active rotations so that
// the viewer can update its cartoon backbone.
// Only one of these switches is actually necessary.
// If both sets of data are given, the viewer will ignore the matrix and use only the
// rotation data.
// Also, currently the code interlaces the two datasets which is not ideal.
#define write_activation_matrix 0
#define write_active_rotation 1

// Amount to reduce momenta for path nodes beyond zero. Since the point of path based
// docking is to keep as closely as possible the same ligand pose and move it through
// the protein, we want to minimize the ligand's conformational changes from node to
// node. Linear momenta are not affected by this number. Angular momenta are multiplied
// by this number, and bond rotation momenta by the square of this number.
#define internode_momentum_mult 0.25

#define internode_momentum_only_on_activation 1

// Switches for conformational space search.
#define allow_linear_motion 1
#define allow_axial_tumble 1
#define allow_bond_rots 1
#define monte_carlo_axial 0
#define monte_carlo_flex 1
#define multiflex 0
#define multiflex_iter 25
#define _multiflex_stepdiv 8
#define _multiflex_stepdeg (360.0/_multiflex_stepdiv)
#define _multiflex_steprad (fiftyseventh*_multiflex_stepdeg)
#define repack_on_hxr false

// Drift pulls the ligand towards the loneliest point if it encounters clashes.
// Turning it off can cause the ligand to be ejected from the protein.
#define allow_drift 1
#define initial_drift 0.3
#define drift_decay_rate 0.01
#define drift_energy_threshold 1000

// Allows full 360 degree whole molecule rotations to search for lower energy configurations.
#define allow_mol_fullrot_iter 0
#define sidechain_fullrot_lig_bmult 3

// Extra weight given to sidechain-ligand binding strengths during conformer search.
#define dock_ligand_bias 1.0

// Turns off the 360 degree rotations.
#define nodes_no_ligand_360_tumble 1
#define nodes_no_ligand_360_flex 1
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

// Overwrites the user supplied pocket center with the loneliest point as determined by
// distances to the nearest residue atoms to the supplied pocket center.
#define pocketcen_is_loneliest 0
#define pocketcen_from_reach_atoms 0

#define preemptively_minimize_intermol_clashes 0
#define bestbind_springiness 15

// In the best-binding algorithm, how strongly to associate hydrophilic features of the ligand
// with hydrophilicity of the side chains.
#define hydrophilicity_boost 5
#define best_binding_stochastic 0.3
#define limit_interactions_by_hydrophilicity 0

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
#define strongest_loss_tolerance 0.05
#define _slt1 (1.0 - strongest_loss_tolerance)

// Whether to add the indicated partial charge to a neutral pnictogen, within pKa limits,
// if a negatively charged atom is nearby.
#define _ALLOW_PROTONATE_PNICTOGENS 0
#define _allow_conditional_basicity 1
#define _allow_conditional_basicity_with_acid_ligand 0
#define cond_bas_hbond_energy_threshold 10.0
#define cond_bas_hbond_distance_threshold 2.9

// Whether to move water molecules around that are clashing or are not forming intermolecular bonds.
#define _teleport_dissatisfied_waters 1
// Threshold is positive for binding, negative for clashes.
#define _water_satisfaction_threshold 5
#define _water_teleport_tries 200

#define polar_sat_influence_for_dock 30
#define polar_sat_influence_for_bb 30
#define polar_sat_influence_for_scoring 0

// How strong an intermolecular bond is required to prevent a histidine hydrogen from flipping
// to a less favorable state.
#define _hisflip_binding_threshold 25

#define hydrogenate_add_missing_heavy_atoms 1

#define homology_phi_psi_rotations 0
#define homology_long_axis_rotations 0
#define homology_region_optimization 0

//////////////////////////////////////////////////////////////////////////////////////////////

// Debugging stuff.

// Should normally be 3:
#define _bb_max_grp 3

// For auditing binding energies between individual atoms:
#define _peratom_audit 0
#define _peratom_audit_nans 0

// A short-term feature that can be used in case the bond reciprocity problem recurs. See issue #423.
#define bond_reciprocity_fix 0

// Should normally be false or zero:
#define _dbg_259 0
#define _dbg_51e2_ionic 0
#define _dbg_A100 0
#define _dbg_allow_excessive_aa_clashes 0
#define _dbg_anemia 0
#define _dbg_asunder_atoms 0
#define _dbg_bb_pullaway 0
#define _dbg_bb_realign 0
#define _dbg_bb_rots 0
#define _dbg_bridges 0
#define _dbg_can_rotate 0
#define _dbg_cond_basic 0
#define _dbg_cond_basic_acd_lig 0
#define _dbg_conj_chg 0
#define _dbg_conjugation 0
#define _dbg_cvty_pose_filter 0
#define _dbg_eclipses 0
#define _dbg_fitness_plummet 0
#define _dbg_flexion_selection 0
#define _dbg_groupsalign 0
#define _dbg_groupsel 0
#define _DBG_H2O_TELEPORT 0
#define _DBG_HISFLIP 0
#define _dbg_homology 0
#define _dbg_Huckel 0
#define _dbg_hxrax 0
#define _dbg_hydrogenate 0
#define _dbg_identify_rings 0
#define _dbg_imidazole_check 0
#define _dbg_improvements_only_rule 0
#define _dbg_interatomic_forces 0
#define _dbg_internal_energy 0
#define _dbg_linear_motion 0
#define _DBG_LONELINESS 0
#define _dbg_mand_conn 0
#define _DBG_MAX_CLASHES 0
#define _dbg_moieties 0
#define _DBG_MOLBB 0
#define _dbg_mol_flexion 0
#define _dbg_mol_frames 0
#define _dbg_molstruct_evolutions 0
#define _dbg_molstruct_evolution_bond_lengths 0
#define _dbg_multiflex 0
#define _dbg_null_flexions 0
#define _dbg_path_search 0
#define _dbg_point_avg 0
#define _dbg_polsat 0
#define _dbg_repack 0
#define _DBG_RESBMULT 0
#define _dbg_residue_poses 0
#define _dbg_rock_pic 0
#define _dbg_soft 0
#define _dbg_stays_enforce 0
#define _DBG_SPACEDOUT 0
#define _DBG_STEPBYSTEP 0
#define _DBG_TOOLARGE_DIFFNUMS 0
#define _DBG_TUMBLE_SPHERES 0
#define _dbg_unreciprocated_bonds 0
#define _dbg_worst_energy 0
#define _debug_active_bond_rot 0
#define debug_stop_after_tumble_sphere 0
#define _DORESPHRES 0
#define _dummy_atoms_for_debug 0
#define output_tumble_debug_docs 0
#define _show_cond_bas_hbond_energy 0
#define _show_final_group_pairs 0

#endif











