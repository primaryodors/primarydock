<?php

// protutils.php
//
// Loads receptor data into memory and provides useful protein functions. 
//

global $prots;

function resno_from_bw($protid, $bw)
{
	global $prots;
	if (!isset($prots[$protid])) die("Protein not found: $protid.\n");
	
	$pettia = explode(".", $bw);
	$tmrno = intval($pettia[0]);
	$offset = intval($pettia[1]);
	
	$res50 = intval(@$prots[$protid]["bw"]["$tmrno.50"]) or die("Unknown Ballesteros-Weinstein number: $bw.\n");
	
	return $res50 + $offset - 50;
}

function family_from_protid($protid)
{
	if (substr($protid, 0, 2) == "OR") return "OR".intval(substr($protid, 2, 2));
	else return substr($protid, 0, 4);
}

function json_encode_pretty($array)
{
	return preg_replace("/([ \t]*)([^\\s]*) {\n/", "\$1\$2\n\$1{\n", json_encode($array, JSON_PRETTY_PRINT));
}

chdir(__DIR__);
chdir("..");
$prots = json_decode(file_get_contents("data/receptor.json"), true);


