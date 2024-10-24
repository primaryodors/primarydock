<?php

// method_directmdl.php
//
// Performs a dock of an odorant inside inactive and active conformer PDB files.
// The active PDB may be an experimental structure from a crystallography or cryo-EM
// measurement, or it may be a homology model based on experimental models.
//
// Example call syntax:
// php -f predict/method_directmdl.php prot=OR1A1 lig=d-limonene
//

chdir(__DIR__);
// require_once("algorithm_delta.php");
require_once("algorithm_affinity.php");
chdir(__DIR__);
require_once("methods_common.php");
chdir(__DIR__);
require_once("../data/protutils.php");
chdir(__DIR__);

// Configurable variables
$flex = 1;                      // Flexion (0 or 1).
$pose = 15;
$iter = 30;
$elim = 1e3;                    // Energy limit for poses. (Not the tailor/spy from the space station.)
$num_std_devs = 2.0;            // How many standard deviations to move the helices for active clash compensation.

prepare_outputs();

$metrics_to_process["BEST"] = "Pose1";
$metrics_to_process["A100"] = "A100";
$metrics_to_process["occlusion"] = "occlusion";

chdir(__DIR__);
chdir("..");

$pdbfname_inactive = $pdbfname;
$pdbfname_active = str_replace(".upright.pdb", ".active.pdb", $pdbfname);
$paramfname = str_replace(".upright.pdb", ".params", $pdbfname);

if (!file_exists($pdbfname_active)
    /* || filemtime($pdbfname_active) < filemtime("hm/allgpcr.ali")
    || filemtime($pdbfname_active) < filemtime("hm/build_alignment_file.php")
    || filemtime($pdbfname_active) < filemtime("hm/dohm.php")
    || filemtime($pdbfname_active) < filemtime("hm/experimental.ali")
    */
    )
{
    if (filemtime("hm/experimental.ali") > filemtime("hm/allgpcr.ali")) passthru("php -f hm/build_alignment_file.php");
    passthru("php -f hm/dohm.php $protid");
}

if (!file_exists($pdbfname_active)) die("No active model.\n");

$flex_constraints = "";
if (file_exists($paramfname)) $flex_constraints = file_get_contents($paramfname);

$fam = family_from_protid($protid);
$pdbfname = $pdbfname_inactive;
$outfname = "output/$fam/$protid/$protid.$ligname.inactive.dock";
$cenres = $cenres_inactive;

prepare_receptor($pdbfname, "$flxr $iflxr");

// Convert ligand as well.
prepare_ligand($ligname);

chdir(__DIR__);
chdir("..");
if (!file_exists("output/$fam")) mkdir("output/$fam");
if (!file_exists("output/$fam/$protid")) mkdir("output/$fam/$protid");

if (!@$_REQUEST["acvonly"]) process_dock("i");


$pdbfname = $pdbfname_active;
$outfname = "output/$fam/$protid/$protid.$ligname.active.dock";
$cenres = $cenres_active;

prepare_receptor($pdbfname, "$flxr $aflxr");
$poses = process_dock("a");
