<?php
header("Access-Control-Allow-Origin: *");

header("Cache-Control: no-cache, no-store, must-revalidate");
header("Pragma: no-cache");
header("Expires: 0");

chdir(__DIR__);
require_once("../data/protutils.php");

$protid = @$_REQUEST['p'];
if (!$protid || !isset($prots[$protid]) )
{
    header("Location: receptors.php");
    exit;
}

$prot = $prots[$prot];

$fam = family_from_protid($protid);

chdir(__DIR__);
chdir("..");

$mod = "upright";
if (@$_REQUEST['mod'] == 'm') $mod = "metal";
else if (@$_REQUEST['mod'] == 'a') $mod = "active";
$pdbfn = "pdbs/$fam/$protid.$mod.pdb";
if (!file_exists($pdbfn)) $pdbfn = "pdbs/$fam/$protid.upright.pdb";

if (!file_exists($pdbfn))
{
    header("Location: receptors.php");
    exit;
}

echo file_get_contents($pdbfn);

