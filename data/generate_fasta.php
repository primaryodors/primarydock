<?php

chdir(__DIR__);
chdir("..");
require_once("predict/protutils.php");

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
    if (!isset($prots[$rcp])) die("Unable to locate $rcp in data/receptor.json; please check spelling.\n");
    echo ">$rcp\n{$prots[$rcp]['sequence']}\n\n";
}

if ($gp)
{
    if (!isset($gprots[$gp])) die("Unable to locate $gp in data/gprot.json; please check spelling.\n");
    echo ">$gp\n{$gprots[$gp]['sequence']}\n\n";
}