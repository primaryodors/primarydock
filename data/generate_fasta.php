<?php

chdir(__DIR__);
chdir("..");
require_once("data/protutils.php");

$seqa = [];
if (@$_REQUEST['aligned'])
{
    foreach (explode("\n", file_get_contents("data/sequences_aligned.txt")) as $ln)
    {
        $rcp = explode(" ", $ln)[0];
        if (isset($prots[$rcp])) $seqa[$rcp] = str_replace(" ", "-", substr($ln, 8));
    }
}

function get_sequence($protid)
{
    global $prots, $seqa;

    if (@$_REQUEST['aligned']) return $seqa[$protid];
    else if (@$_REQUEST['byrgn'])
    {
        $ret = "";
        $p = $prots[$protid];
        $s = $p['sequence'];
        $j = 0;
        foreach ($p['region'] as $rgn => $regionse)
        {
            $ret .= substr($s, $j, $regionse['start']-$j-1);
            $ret .= "----------";
            $ret .= substr($s, $regionse['start']-1, $regionse['end'] - $regionse['start']+1);
            $ret .= "----------";
            $j = $regionse['end'];
        }
        if ($j) $ret .= substr($s, $j);
        return $ret;
    }
    else if (@$_REQUEST['bybw'])
    {
        $ret = "";
        $p = $prots[$protid];
        $s = $p['sequence'];
        $j = 0;
        foreach ($p['region'] as $rgn => $regionse)
        {
            // Make each rgend-rgend take up 200 places, with the x.50 residue at position 150.
            $seg = substr($s, $j, $regionse['end']-$j);

            $rgno = intval(preg_replace("/[^0-9]/", "", $rgn));
            if ($rgno >= 10) continue;
            try
            {
                $n50 = resno_from_bw($protid, "$rgno.50");
            }
            catch (Exception $ex)
            {
                continue;
            }

            $seg = str_repeat("-", 150 - ($n50 - $j)).$seg;
            $seg .= str_repeat("-", 200-strlen($seg));

            $ret .= $seg;
            $j = $regionse['end'];
        }
        if ($j) $ret .= substr($s, $j);
        return $ret;
    }
    else return $prots[$protid]['sequence'];
}

function make_fasta($protid, $sequence)
{
    $ret = ">$protid";
    while ($sequence)
    {
        $ret .= "\n".substr($sequence, 0, 80);
        $sequence = substr($sequence, 80);
    }
    $ret .= "*\n\n";
    return $ret;
}

foreach (@$argv as $a)
{
	$a = explode('=',$a,2);
	$_REQUEST[$a[0]] = (count($a)>1) ? $a[1] : true;
}

$rcp = @$_REQUEST['rcp'];
$gp = @$_REQUEST['gprot'];

if (!$rcp && !$gp)
{
    die("No protein selected. Please specify a GPCR and/or a G protein, e.g. generate_fasta.php rcp=OR51E2 gprot=hGNAS2\n");
}

if ($rcp)
{
    if ($rcp == "all")
    {
        foreach ($prots as $rcp => $p) echo make_fasta($rcp, get_sequence($rcp));
    }
    else if ($rcp == "allgpcr")
    {
        foreach ($prots as $rcp => $p) if (substr($rcp, 0, 4) != "MS4A") echo make_fasta($rcp, get_sequence($rcp));
    }
    else
    {
        if (!isset($prots[$rcp])) die("Unable to locate $rcp in data/receptor.json; please check spelling.\n");
        echo make_fasta($rcp, get_sequence($rcp));
    }
}

if ($gp)
{
    if (!isset($gprots[$gp])) die("Unable to locate $gp in data/gprot.json; please check spelling.\n");
    echo ">$gp\n{$gprots[$gp]['sequence']}\n\n";
}
