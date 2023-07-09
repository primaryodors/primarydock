<?php
header("Access-Control-Allow-Origin: *");

header("Cache-Control: no-cache, no-store, must-revalidate");
header("Pragma: no-cache");
header("Expires: 0");

chdir(__DIR__);
require_once("../predict/protutils.php");

$protid = @$_REQUEST['p'];
$protid = str_replace("-", "_", $protid);
if (!$protid || !isset($prots[explode("_",$protid)[0]]) )
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
if (file_exists("pdbs/coupled/$fam/$protid.pdb")) $pdbfn = "pdbs/coupled/$fam/$protid.pdb";
else $pdbfn = "pdbs/$fam/$protid.$mod.pdb";
if (!file_exists($pdbfn)) $pdbfn = "pdbs/$fam/$protid.upright.pdb";

if (!file_exists($pdbfn))
{
    header("Location: receptors.php");
    exit;
}

echo file_get_contents($pdbfn);

