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

$result = file_get_contents($pdbfn);

$cavfn = str_replace(".pdb", ".cav", $pdbfn);
if (file_exists($cavfn))
{
    $lines = explode("\n", $result);
    foreach ($lines as $i => $ln)
    {
        if (substr($ln, 0, 10) != "REMARK 800") continue;
        $next = @$lines[$i+1];
        if (substr($next, 0, 10) != "REMARK 800")
        {
            $ln .= "\nREMARK 821";
            $cavs = explode("\n", file_get_contents($cavfn));
            foreach ($cavs as $cav)
            {
                if (!trim($cav)) continue;
                $ln .= "\nREMARK 821 $cav";
            }
            $ln .= "\nREMARK 821";
            $lines[$i] = $ln;
            break;
        }
    }

    $result = implode("\n", $lines);
}

echo $result;
