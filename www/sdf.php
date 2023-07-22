<?php
header("Access-Control-Allow-Origin: *");

header("Cache-Control: no-cache, no-store, must-revalidate");
header("Pragma: no-cache");
header("Expires: 0");

chdir(__DIR__);
require_once("../data/odorutils.php");

if (@$_REQUEST['m'] == "rand") $odor = $odors[array_keys($odors)[rand(0,count($odors)-1)]];
else
{
	$odor = find_odorant(@$_REQUEST['m']);
	if (!$odor)
	{
		header("Location: odorants.php");
		exit;
	}
}

chdir(__DIR__);
chdir("..");
$fullname = $odor['full_name'];
$sdfname = "sdf/$fullname.sdf";
if (!file_exists($sdfname))
{
    $smilesu = urlencode($odor['smiles']);
    $url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/SDF";

    $ch = curl_init( $url );
    curl_setopt( $ch, CURLOPT_POST, 1);
    curl_setopt( $ch, CURLOPT_POSTFIELDS, $arg = "record_type=3d&smiles=$smilesu");
    curl_setopt( $ch, CURLOPT_FOLLOWLOCATION, 1);
    curl_setopt( $ch, CURLOPT_HEADER, 0);
    curl_setopt( $ch, CURLOPT_RETURNTRANSFER, 1);

    $sdfdat = curl_exec( $ch );
    $sdfln = explode("\n", $sdfdat);
    if (false!==strpos($sdfln[0], "Status: ")) die("$arg\n\n$sdfdat");

    $fp = fopen($sdfname, "wb");
    if ($fp)
    {
        fwrite($fp, $sdfdat);
        fclose($fp);
    }
    else error_log("Warning: Unable to write to $sdfname");
}
else
{
    $sdfdat = file_get_contents($sdfname);
    $sdfln = explode("\n", $sdfdat);
    if (false!==strpos($sdfln[0], "Status: ")) unlink($sdfname);
}

$sdfdat = explode("\n", $sdfdat);
$sdfdat[0] = $fullname;
$sdfdat = implode("\n", $sdfdat);

echo $sdfdat;
