<?php

chdir(__DIR__);
chdir('..');
require_once("predict/protutils.php");

define("_ALPHAFOLD", 22);
define("_ZHANGLAB",  17);
define("_DEFSRC", _ALPHAFOLD);

$source = _DEFSRC;

foreach (@$argv as $a)
{
	$a = explode('=',$a,2);
	$_REQUEST[$a[0]] = (count($a)>1) ? $a[1] : true;
}

if (@$_REQUEST['src'])
{
    $c = strtoupper(substr($_REQUEST['src'], 0, 1));
    switch ($c)
    {
        case 'A':   $source = _ALPHAFOLD;   break;
        case 'Z':   $source = _ZHANGLAB;    break;
        default: break;
    }
}

foreach ($prots as $protid => $p)
{
    if (!isset($p['region'])) continue;
    $rgns = "";
    foreach ($p['region'] as $rgname => $se)
    {
        $rgns .= "REGION $rgname {$se['start']} {$se['end']}\n";
    }

    $rems = "";
    if (isset($p['bw']))
    {
        foreach ($p['bw'] as $b => $w)
        {
            $rems .= "REMARK 800 SITE BW $b $w\n";
        }
    }

    $bsr = binding_site($protid);
    foreach ($bsr as $resno)
    {
        $rems .= "REMARK 800 SITE LIGAND_BINDING $resno\n";
    }

    $uid = $p['uniprot_id'];

    switch ($source)
    {
        case _ALPHAFOLD:
        $url = "https://alphafold.ebi.ac.uk/files/AF-$uid-F1-model_v4.pdb";
        break;

        case _ZHANGLAB:
        $zhang = 1;
        if ($protid == 'OR5K1') $zhang = 2;
        $url = "https://zhanggroup.org/GPCR-EXP/pdb/hgmod/$uid/{$uid}_results/model$zhang.pdb";

        default:
        die("Unknwon source $source.\n");
    }

    $fam = family_from_protid($protid);
    if (!file_exists("pdbs")) mkdir("pdbs");
    if (!file_exists("pdbs/$fam")) mkdir("pdbs/$fam");
    if (!file_exists("pdbs/$fam/import")) mkdir("pdbs/$fam/import");
    $infname = "pdbs/$fam/import/".substr($url, strrpos($url, "/")+2);
    $outfname = "pdbs/$fam/$protid.upright.pdb";

    $antitax_server = true;
    if (!file_exists($infname))
    {
        echo ($cmd = "wget -O \"$infname\" \"$url\"");
        echo "\n";
        passthru($cmd);
    }
    else $antitax_server = false;

    if (filesize($infname) < 100000) die("File too short: $infname\n");

    $editf = <<<heredoc
LOAD $infname
REMARK   1
REMARK   1 This file has been modified for PrimaryDock compatibility.
REMARK   1
REMARK 600
$rgns
UPRIGHT
HYDRO
DISULF
REMARK 800
$rems
REMARK 800
SAVE $outfname

heredoc;

    $pepdname = "tmp/get.pepd";
    $fp = fopen($pepdname, "wb") or die("FAILED to open $pepdname for writing.\n");
    fwrite($fp, $editf);
    fclose($fp);

    passthru("bin/pepteditor $pepdname");

    if ($antitax_server)
    {
        sleep(1);
        usleep(1048576 + rand(0, 524288));
    }
    exit;
}