<?php
header("Access-Control-Allow-Origin: *");

header("Cache-Control: no-cache, no-store, must-revalidate");
header("Pragma: no-cache");
header("Expires: 0");

chdir(__DIR__);
require_once("../data/protutils.php");

$protid = @$_REQUEST['recep'];

if (!isset($prots[$protid])) die("$protid not found.");

$rcp = $prots[$protid];
echo "id={$rcp['id']}|";
echo "id={$rcp['uniprot_id']}|";
echo "id={$rcp['sequence']}|";
echo "\n";

foreach ($rcp['region'] as $rgname => $se)
{
    echo "REGION|";
    echo "$rgname|";
    echo "{$se['start']}|";
    echo "{$se['end']}|";
    echo "\n";
}

if (@$rcp["bw"])
{
    foreach ($rcp["bw"] as $tmr => $bw50)
    {
        echo "BW|";
        echo "$tmr|";
        echo "$bw50|";
        echo "\n";
    }
}

// TODO: Metal coordination.

$bsr = array_flip(binding_site($protid));
$fam = family_from_protid($protid);

$px = $py = $pz = 0.0;
$divisor = 0;
chdir(__DIR__);
$pdbname = "../pdbs/$fam/$protid.upright.pdb";
if (file_exists($pdbname))
{
    $fp = fopen($pdbname, "r");
    while (!feof($fp))
    {
        $ln = fgets($fp);

        //           1111111111222222222233333333334444444444555555555566666666667777777777
        // 01234567890123456789012345678901234567890123456789012345678901234567890123456789
        // ATOM   1734  CA  ASN   109      -1.539  -0.510  -6.991  1.00001.00           C
        if (trim(substr($ln, 0, 5)) == "ATOM")
        {
            $resno = intval(substr($ln, 22, 5));
            if (isset($bsr[$resno]))
            {
                $px += floatval(substr($ln, 31, 8));
                $py += floatval(substr($ln, 39, 8));
                $pz += floatval(substr($ln, 47, 8));
                $divisor++;
            }
        }
    }
    fclose($fp);
}

if ($divisor)
{
    $px /= $divisor;
    $py /= $divisor;
    $pz /= $divisor;
}

echo "POCKET|type=A|";
echo "center_x=$px|";
echo "center_y=$py|";
echo "center_z=$pz|";
echo "\n";

