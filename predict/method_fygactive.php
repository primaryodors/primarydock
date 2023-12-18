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
$flxi = 1;                      // Flexion for inactive dock.
$pose = 10;
$iter = 50;
$elim = 1e4;                    // Energy limit for poses. (Not the tailor/spy from the space station.)
$num_std_devs = 1.0;            // How many standard deviations to move the helices for active clash compensation.

chdir(__DIR__);
require_once("methods_common.php");
chdir(__DIR__);
require_once("../data/protutils.php");
chdir(__DIR__);

prepare_outputs();

$metrics_to_process =
[
    "BENERG" => "BindingEnergy",
    "BENERG.rgn" => "BindingEnergy.rgn",
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

function build_template()
{
    global $template, $protid, $cryoem, $has_rock6, $has_fyg, $args;

    if (substr($protid, 0, 4) == "TAAR")
    {
        $template = $cryoem["mTAAR9"];
        // TODO: Blend mTAAR9 with TAAR1 depending on sequence similarity.
    }
    else if (substr($protid, 0, 2) == "OR")
    {
        $fam = intval(substr($protid, 2, 2));
        if ($fam == 51)
        {
            $template = $cryoem["OR51E2"];
        }
        else if ($fam == 52)
        {
            foreach ($cryoem["OR52"] as $hxno => $metrics)
            {
                foreach ($metrics as $metric => $dimensions)
                {
                    $template[$hxno][$metric] = $cryoem["OR52"][$hxno][$metric];

                    $template[$hxno][$metric]["sigma"]
                        + 0.4 * $cryoem["TAAR1"][$hxno][$metric]["sigma"]
                        + 0.6 * $cryoem["mTAAR9"][$hxno][$metric]["sigma"];
                }
            }
        }
        else
        {
            foreach ($cryoem["OR51E2"] as $hxno => $metrics)
            {
                foreach ($metrics as $metric => $dimensions)
                {
                    foreach (array_keys($dimensions) as $dimension)
                    {
                        if ($dimension == "sigma")
                        {
                            $template[$hxno][$metric][$dimension]
                                = 0.4 * $cryoem["TAAR1"][$hxno][$metric][$dimension]
                                + 0.6 * $cryoem["mTAAR9"][$hxno][$metric][$dimension];
                        }
                        else
                        {
                            // TODO: Make the proportions dependent on sequence similarity.
                            $template[$hxno][$metric][$dimension]
                                = 0.40 * $cryoem["OR51E2"][$hxno][$metric][$dimension]
                                + 0.30 * $cryoem["OR52"][$hxno][$metric][$dimension]
                                + 0.10 * $cryoem["TAAR1"][$hxno][$metric][$dimension]
                                + 0.20 * $cryoem["mTAAR9"][$hxno][$metric][$dimension];
                        }
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
}


chdir(__DIR__);
chdir("..");

$pdbfname_active = str_replace(".upright.pdb", ".active.pdb", $pdbfname);
$paramfname = str_replace(".upright.pdb", ".params", $pdbfname);
$template = [];
$args = "$protid";

$cryoem = json_decode(file_get_contents("data/cryoem_motions.json"), true);

if (!file_exists($pdbfname_active) || filemtime($pdbfname_active) < filemtime("bin/fyg_activate_or"))
{
    build_template();

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

    $cmd = "bin/fyg_activate_or $args";
    echo "$cmd\n";
    passthru($cmd);
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
ELIM $elim
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

SEARCH $search
POSE $pose
ELIM $elim
$flex_constraints
ITERS $iter
PROGRESS

FLEX $flex
WET

OUT $outfname
OUTPDB 1 output/$fam/$protid/%p.%l.active.model%o.pdb


heredoc;

$poses = process_dock("a");

if ((!$poses || $best_energy >= 0) && count($clashcomp))
{
    if (!count($template)) build_template();

    // Ensure template contains standard deviations. If not, average the ones from the templates that do.
    foreach ($template as $hxno => $metrics)
    {
        foreach ($metrics as $metric => $xyz)
        {
            if (!isset($xyz["sigma"]))
            {
                $averages = [];
                $counts = [];

                foreach ($cryoem as $lrecep => $lhelixes)
                {
                    foreach ($lhelixes as $lhelix => $lmetrics)
                    {
                        foreach ($lmetrics as $lmetric => $lxyz)
                        {
                            if (isset($lxyz["sigma"]))
                            {
                                if (!isset($averages[$lhelix][$lmetric])) $averages[$lhelix][$lmetric] = floatval($lxyz["sigma"]);
                                else $averages[$lhelix][$lmetric] += floatval($lxyz["sigma"]);

                                if (!isset($counts[$lhelix][$lmetric])) $counts[$lhelix][$lmetric] = 1;
                                else $counts[$lhelix][$lmetric]++;
                            }
                        }
                    }
                }

                foreach ($averages as $lhelix => $lmetrics)
                {
                    foreach ($lmetrics as $lmetric => $total)
                    {
                        if (!isset($template[$lhelix][$lmetric]["sigma"]))
                        {
                            if ($counts[$lhelix][$lmetric]) $total /= $counts[$lhelix][$lmetric];
                            $template[$lhelix][$lmetric]["sigma"] = $total;
                        }
                    }
                }

                break;
            }
        }
    }

    foreach ($clashcomp as $tmrno => $cc)
    {
        foreach ($cc as $segment => $xyz)
        {
            // Check each clash compensation is within the SD limit. If not, reduce the magnitude to the deviation limit.
            $sigma = floatval(@$template[$tmrno][$segment]["sigma"]);
            $rlimit = $sigma*$num_std_devs;
            if (!$sigma) continue;
            $tmparr = []; foreach (['x','y','z'] as $idx => $var) $tmparr[$var] = floatval($xyz[$idx]); extract($tmparr);
            $r = sqrt($x*$x+$y*$y+$z*$z);

            if ($r > $rlimit)
            {
                $divisor = $rlimit / $r;
                $x *= $divisor; $y *= $divisor; $z *= $divisor;
            }

            // Apply the corrected compensation.
            foreach (['x','y','z'] as $var)
            {
                $$var = round($$var, 3);
                $template[$tmrno][$segment][$var] = floatval($template[$tmrno][$segment][$var]) + $$var;
                echo "Compensating TMR$tmrno.$segment.$var to {$template[$tmrno][$segment][$var]}...\n";
            }
        }
    }

    // Create a temporary custom output PDB and retry the active dock.
    $args = "$protid";
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

    $tmpoutpdb = "tmp/$protid.".getmypid().".pdb";
    $args .= " -o $tmpoutpdb";

    $cmd = "bin/fyg_activate_or $args";
    echo "$cmd\n";
    passthru($cmd);


    $configf = str_replace($pdbfname, $tmpoutpdb, $configf);
    $poses = process_dock("a");

    // Delete the tmp PDB.
    unlink($tmpoutpdb);
}
