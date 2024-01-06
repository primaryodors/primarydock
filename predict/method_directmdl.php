<?php

// method_fygactive.php
//
// Performs a dock of an odorant inside inactive and FYG-motif active conformer PDB files.
//
// Example call syntax:
// php -f predict/method_fygactive.php prot=OR1A1 lig=d-limonene
//

chdir(__DIR__);
require_once("methods_common.php");
chdir(__DIR__);
require_once("../data/protutils.php");
chdir(__DIR__);

prepare_outputs();

function make_prediction($data)
{
    global $protid, $ligname;

    if (isset($data["a_BENERG"]))
    {
        $ascore = min(0, floatval( $data['a_BENERG']));
        $iscore = min(0, floatval(@$data['i_BENERG']));

        if ($ascore < 0 && $ascore < $iscore)
        {
            $data['Predicted'] = 'Agonist';
            $data['DockScore'] = (min($iscore, 0) - $ascore);
        }
        else if ($iscore < 0 && $iscore < $ascore)
        {
            $data['Predicted'] = 'Inverse Agonist';
            $data['DockScore'] = (min($iscore, 0) - $ascore);
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

$pdbfname_inactive = str_replace(".upright.pdb", ".apo.pdb", $pdbfname);
$pdbfname_active = str_replace(".upright.pdb", ".bound.pdb", $pdbfname);
$paramfname = str_replace(".upright.pdb", ".params", $pdbfname);

if (!file_exists($pdbfname_inactive) && file_exists($pdbfname_active))
    $pdbfname_inactive = $pdbfname;

if (!file_exists($pdbfname_active))
{
    exec("php -f predict/cryoem_motions.php");
    exec("bin/pepteditor data/OR51.pepd");
    exec("bin/pepteditor data/OR52.pepd");
}

$flex_constraints = "";
if (file_exists($paramfname)) $flex_constraints = file_get_contents($paramfname);

$fam = family_from_protid($protid);
$pdbfname = $pdbfname_inactive;
$outfname = "output/$fam/$protid/$protid.$ligname.inactive.model1.pdb";

// Filter out everything except ATOM records.
// exec("cat $pdbfname_inactive | grep ATOM > tmp/prot.pdb");
$lines = explode("\n", file_get_contents($pdbfname_inactive));
$rf = split_pdb_to_rigid_and_flex($protid, $lines, explode(" ", "$flxr $iflxr"));
$fp = fopen("tmp/prot.pdb", "w");
if (!$fp) die("Failed to write to tmp/prot.pdb.\n");
fwrite($fp, implode("\n",$rf[0]));
fclose($fp);
$fp = fopen("tmp/flex.pdb", "w");
if (!$fp) die("Failed to write to tmp/flex.pdb.\n");
fwrite($fp, implode("\n",$rf[1]));
fclose($fp);

// Convert to PDBQT format.
exec("obabel -i pdb tmp/prot.pdb -xr -o pdbqt -O tmp/prot.pdbqt");
exec("obabel -i pdb tmp/flex.pdb -xs -o pdbqt -O tmp/flex.pdbqt");

// Convert ligand as well.
exec("obabel -i sdf \"sdf/$ligname.sdf\" -o pdbqt -O tmp/lig.pdbqt");

chdir(__DIR__);
chdir("..");
if (!file_exists("output/$fam")) mkdir("output/$fam");
if (!file_exists("output/$fam/$protid")) mkdir("output/$fam/$protid");

if (!@$_REQUEST["acvonly"]) process_dock("i");


$pdbfname = $pdbfname_active;
$outfname = "output/$fam/$protid/$protid.$ligname.active.model1.pdb";

// Filter out everything except ATOM records.
// exec("cat $pdbfname_active | grep ATOM > tmp/prot.pdb");
$lines = explode("\n", file_get_contents($pdbfname_active));
$rf = split_pdb_to_rigid_and_flex($protid, $lines, explode(" ", "$flxr $aflxr"));
$fp = fopen("tmp/prot.pdb", "w");
if (!$fp) die("Failed to write to tmp/prot.pdb.\n");
fwrite($fp, implode("\n",$rf[0]));
fclose($fp);
$fp = fopen("tmp/flex.pdb", "w");
if (!$fp) die("Failed to write to tmp/flex.pdb.\n");
fwrite($fp, implode("\n",$rf[1]));
fclose($fp);

// Convert to PDBQT format.
exec("obabel -i pdb tmp/prot.pdb -xr -o pdbqt -O tmp/prot.pdbqt");
exec("obabel -i pdb tmp/flex.pdb -xs -o pdbqt -O tmp/flex.pdbqt");

$poses = process_dock("a");
