<?php

// method_fygactive.php
//
// Performs a dock of an odorant inside inactive and FYG-motif active conformer PDB files.
//
// Example call syntax:
// php -f predict/method_fygactive.php prot=OR1A1 lig=d-limonene
//

// Configurable variables
$flex = 1;                      // Flexion (0 or 1) for active dock.
$flxi = 0;                      // Flexion for inactive dock.
$pose = 25;
$iter = 20;

chdir(__DIR__);
require_once("methods_common.php");
chdir(__DIR__);
require_once("../data/protutils.php");
chdir(__DIR__);

prepare_outputs();

$metrics_to_process =
[
    "BENERG" => "BindingEnergy",
    // "BENERG.rgn" => "BindingEnergy.rgn",
    "BEST" => "Pose1"
];

function make_prediction($data)
{
    global $protid, $ligname;

    if (isset($data["a_Pose1"]))
    {
        $ae = floatval(@$data['a_BindingEnergy']);
        $ie = floatval(@$data["i_BindingEnergy"]);
        $a1 = floatval( $data['a_Pose1']);
        $i1 = floatval(@$data['i_Pose1']);
        if ($a1 < 0 && $a1 < $i1)
        {
            $data['Predicted'] = 'Agonist';
            $data['DockScore'] = (min($i1, 0) - $a1) / 2;
        }
        else if ($ae < 0 && $ae < $ie)
        {
            $data['Predicted'] = 'Agonist';
            $data['DockScore'] = (min($ie, 0) - $ae) / 2;
        }
        else if ($a1 < 0 && $a1 > $i1)
        {
            $data['Predicted'] = 'Inverse Agonist';
            $data['DockScore'] = (min($i1, 0) - $a1) / 2;
        }
        else
        {
            $data['Predicted'] = 'Non-agonist';
            $data['DockScore'] = 0.0;
        }

        echo "\nProtein: $protid\nLigand: $ligname";
        echo "\nResult: " . print_r($data, true) . "\n";
    }
    else if (isset($data["a_POSES"]) && !$data["a_POSES"])
    {
        $data['Predicted'] = 'Non-agonist';
        $data['DockScore'] = 0.0;

        echo "\nProtein: $protid\nLigand: $ligname";
        echo "\nResult: " . print_r($data, true) . "\n";
    }

    return $data;
}


chdir(__DIR__);
chdir("..");

$pdbfname_active = str_replace(".upright.pdb", ".active.pdb", $pdbfname);
$paramfname = str_replace(".upright.pdb", ".params", $pdbfname);

if (!file_exists($pdbfname_active) || filemtime($pdbfname_active) < filemtime("bin/fyg_activate_or"))
{
    $cryoem = json_decode(file_get_contents("data/cryoem_motions.json"), true);

    $args = "$protid";
    $template = [];

    if (substr($protid, 0, 4) == "TAAR")
    {
        $template = $cryoem["mTAAR9"];
        // TODO: Blend mTAAR9 with TAAR1 depending on sequence similarity.
    }
    else if (substr($protid, 0, 2) == "OR")
    {
        $fam = intval(substr($protid, 2, 2));
        if ($fam >= 50)
        {
            $template = $cryoem["OR51E2"];
        }
        else
        {
            foreach ($cryoem["OR51E2"] as $hxno => $metrics)
            {
                foreach ($metrics as $metric => $dimensions)
                {
                    foreach (array_keys($dimensions) as $dimension)
                    {
                        // TODO: Make the proportions dependent on sequence similarity.
                        $template[$hxno][$metric][$dimension]
                            = 0.66 * $cryoem["OR51E2"][$hxno][$metric][$dimension]
                            + 0.25 * $cryoem["mTAAR9"][$hxno][$metric][$dimension]
                            + 0.09 * $cryoem["mTAAR9"][$hxno][$metric][$dimension];
                    }
                }
            }
        }
    }

    // Get typology. Do a no-save run of fyg_activate_or to determine whether protein has FYG or rock6 capabilities.
    $output = [];
    exec("bin/fyg_activate_or --nosave $protid", $output);
    $has_rock6 = $has_fyg = false;
    foreach ($output as $line)
    {
        if (false!==strpos($line, "Performing rock6")) $has_rock6 = true;
        if (false!==strpos($line, "Performing FYG activation")) $has_fyg = true;
    }

    foreach ($template as $hxno => $metrics)
    {
        foreach ($metrics as $metric => $dimensions)
        {
            if ($hxno == 6)
            {
                if ($metric == "cyt" && ($has_fyg || $has_rock6)) continue;
                else if ($metric == "exr" && $has_rock6) continue;
            }

            $cmdarg = "--" . substr($metric, 0, 1) . $hxno;
            $args .= " $cmdarg {$dimensions['x']} {$dimensions['y']} {$dimensions['z']}";
        }
    }

    if (false) // $protid == "OR51E2")
    {
        $pepd = <<<heredoc
LOAD pdbs/OR51/OR51E2.8f76.pdb
CENTER
UPRIGHT
HYDRO
SAVE $pdbfname_active


heredoc;
        $fp = fopen("tmp/8f76.pepd", "wb");
        fwrite($fp, $pepd);
        fclose($fp);
        passthru("bin/pepteditor tmp/8f76.pepd");
    }
    else
    {
        $cmd = "bin/fyg_activate_or $args";
        echo "$cmd\n";
        passthru($cmd);
    }
}

$flex_constraints = "";
if (file_exists($paramfname)) $flex_constraints = file_get_contents($paramfname);

$fam = family_from_protid($protid);
$outfname = "output/$fam/$protid/$protid.$ligname.inactive.dock";

$configf = <<<heredoc

PROT $pdbfname
LIG sdf/$ligname.sdf

$cenres_inactive
SIZE $size
# H2O 5
$mcoord
$atomto
$stcr
$flxr

EXCL 1 56		# Head, TMR1, and CYT1.

SEARCH $search
POSE $pose
ELIM 10000
$flex_constraints
ITERS $iter
PROGRESS

FLEX $flxi
WET

OUT $outfname
OUTPDB 1 output/$fam/$protid/%p.%l.inactive.model%o.pdb


heredoc;

chdir(__DIR__);
chdir("..");
if (!file_exists("output/$fam")) mkdir("output/$fam");
if (!file_exists("output/$fam/$protid")) mkdir("output/$fam/$protid");

if (!@$_REQUEST["acvonly"]) process_dock("i");


$pdbfname = $pdbfname_active;

$outfname = "output/$fam/$protid/$protid.$ligname.active.dock";

$configf = <<<heredoc

PROT $pdbfname
LIG sdf/$ligname.sdf

$cenres_active
SIZE $size
# H2O 5
$mcoord
$atomto
$stcr
$flxr

EXCL 1 56		# Head, TMR1, and CYT1.

# SOFT 4.53 4.64
# SOFT 5.50 5.33
# SOFT 6.48 6.59

SEARCH $search
POSE $pose
ELIM 10000
$flex_constraints
ITERS $iter
PROGRESS

FLEX $flex
WET

OUT $outfname
OUTPDB 1 output/$fam/$protid/%p.%l.active.model%o.pdb


heredoc;

$poses = process_dock("a");
