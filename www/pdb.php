<?php
chdir(__DIR__);
require_once("../predict/protutils.php");

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

$pdbfn = "pdbs/$fam/$protid.upright.pdb";

if (!file_exists($pdbfn))
{
    header("Location: receptors.php");
    exit;
}

echo file_get_contents($pdbfn);

