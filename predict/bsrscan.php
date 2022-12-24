<?php

require_once("protutils.php");

foreach (@$argv as $a)
{
	$a = explode('=',$a,2);
	$_REQUEST[$a[0]] = (count($a)>1) ? $a[1] : true;
}

$search = @$_REQUEST['search'] ?: false;
$result = [];
foreach ($prots as $protid => $p)
{
    if (!isset($p['bw'])) continue;


    $bsrs = binding_site($protid);
    $string = "";
    $count = 0;

    foreach ($bsrs as $resno)
    {
        $c = substr($p['sequence'], $resno-1, 1);
        $string .= $c;
        if ($c == $search) $count++;
    }
    
    $result[] = ($search ? "$count " : "").str_pad("$protid: ", 13).$string;
}

if ($search) natsort($result);
echo implode("\n", $result)."\n";