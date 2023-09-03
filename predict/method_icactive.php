<?php

// method_icactive.php
//
// Performs a dock of an odorant inside inactive and internal-contact active conformer PDB files.
//
// Example call syntax:
// php -f predict/method_icactive.php prot=OR1A1 lig=d-limonene
//

// Configurable variables
$dock_metals = false;

chdir(__DIR__);
require("methods_common.php");
chdir(__DIR__);
$binding_pockets = json_decode(file_get_contents("../data/binding_pocket.json"), true);

prepare_outputs();

$size = "5.0 6.0 5.0";
$search = "BB";
$atomto = [];
$stcr = "";
$flxr = "";
$mcoord = "";

$mbp = false;                       // Matched binding pocket.

if (isset($binding_pockets[$protid])) $mbp = $binding_pockets[$protid];
else foreach ($binding_pockets as $pocketid => $pocket)
{
    if (substr($pocketid, -1) == '*' && substr($pocketid, 0, -1) == substr($protid, 0, strlen($pocketid)-1))
    {
        $mbp = $pocket;
        echo "Matched $pocketid via wildcard.\n";
        break;
    }
    else if (preg_match("/^$pocketid\$/", $protid))
    {
        $mbp = $pocket;
        echo "Matched $pocketid via regex.\n";
        break;
    }
}

if ($mbp)
{
    if (isset($mbp["size"])) $size = $mbp["size"];
    if (isset($mbp["search"])) $search = $mbp["search"];
    if (isset($mbp["mcoord"])) $mcoord = "MCOORD {$mbp["mcoord"]}";
    if (isset($mbp["stcr"])) $stcr = "STCR {$mbp["stcr"]}";
    if (isset($mbp["flxr"])) $flxr = "FLXR {$mbp["flxr"]}";

    if (isset($mbp["atomto"]))
    {
        foreach ($mbp["atomto"] as $a2)
        {
            $atomto[] = "ATOMTO $a2";
        }
    }
}

if ($mbp && isset($mbp["pocket"]))
{
    $cenres_active = $cenres_inactive = "CEN RES {$mbp["pocket"]}";
}
else if ($mbp && isset($mbp["active_pocket"]) && isset($mbp["inactive_pocket"]))
{
    $cenres_active = isset($mbp["active_pocket"]);
    $cenres_inactive = isset($mbp["inactive_pocket"]);
}
else
{
    if (substr($fam, 0, 2) == "OR")
    {
        $cenres_active = $cenres_inactive = "CEN RES 3.37 5.47 6.55 7.41";
    }
    else if (substr($fam, 0, 4) == "TAAR")
    {
        die("There is not yet an internal contacts activation app for TAARs.\n");
        $cenres_active = $cenres_inactive = "CEN RES 3.32 3.37 5.43 6.48 7.43";
    }
    else die("Unsupported receptor family.\n");
}

$atomto = implode("\n", $atomto);

$metrics_to_process =
[
    "BENERG" => "BindingEnergy",
    "BENERG.rgn" => "BindingEnergy.rgn",
    "BEST" => "Pose1"
];

function make_prediction($data)
{
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

        echo "\nResult: " . print_r($data, true) . "\n";
    }

    return $data;
}


chdir(__DIR__);
chdir("..");

$pdbfname_active = str_replace(".upright.pdb", ".active.pdb", $pdbfname);
$paramfname = str_replace(".upright.pdb", ".params", $pdbfname);

if (!file_exists($pdbfname_active) || filemtime($pdbfname_active) < filemtime("bin/ic_activate_or"))
{
    passthru("bin/ic_activate_or $protid");
}

$flex_constraints = "";
if (file_exists($paramfname)) $flex_constraints = file_get_contents($paramfname);


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
POSE 10
ELIM 5000
$flex_constraints
ITERS 50
PROGRESS

FLEX 1
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
POSE 10
ELIM 5000
$flex_constraints
ITERS 50
PROGRESS

FLEX 1
WET

OUT $outfname
OUTPDB 1 output/$fam/$protid/%p.%l.active.model%o.pdb


heredoc;

process_dock("a");
