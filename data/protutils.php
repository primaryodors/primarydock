<?php

// protutils.php
//
// Loads receptor data into memory and provides useful protein functions. 
//

global $prots, $aminos;

function bw_insdel($prot, $tmrno, $offset)
{
	$insdel = 0;
	if (isset($prot['deletion']))
		foreach ($prot['deletion'] as $del)
		{
			$pettia = explode(".", $del);
			$dtmr = intval($pettia[0]);
			$doff = intval($pettia[1]);

			if ($dtmr == $tmrno)
			{
				if ($doff <  50 && $doff >= $offset) $insdel++;
				if ($doff >= 50 && $doff <  $offset) $insdel--;
			}
		}
	if (isset($prot['insertion']))
		foreach ($prot['insertion'] as $ins)
		{
			$pettia = explode(".", $ins);
			$itmr = intval($pettia[0]);
			$ioff = intval($pettia[1]);

			if ($itmr == $tmrno)
			{
				if ($ioff <  50 && $ioff >= $offset) $insdel--;
				if ($ioff >= 50 && $ioff <  $offset) $insdel++;
			}
		}
	
	return $insdel;
}

function resno_from_bw($protid, $bw)
{
	global $prots;
	if (!isset($prots[$protid])) die("Protein not found: $protid.\n");
	
	$pettia = explode(".", $bw);
	$tmrno = intval($pettia[0]);
	$offset = intval($pettia[1]);

	$insdel = bw_insdel($prots[$protid], $tmrno, $offset);
	
	$res50 = intval(@$prots[$protid]["bw"]["$tmrno.50"]);
	if (!$res50) throw new Exception("Unknown Ballesteros-Weinstein number $bw");
	
	return $res50 + $offset - 50 + $insdel;
}

function bw_from_resno($protid, $resno)
{
	global $prots;
	if (!isset($prots[$protid])) die("Protein not found: $protid.\n");

	$prot = $prots[$protid];
	
	foreach ($prot['region'] as $rgn => $se)
	{
		if (substr($rgn, 0, 3) == 'TMR' || substr($rgn, 0, 3) == 'HXR')
		{
			$tmrno = intval(substr($rgn, -1));
			if ($resno >= $se['start'] && $resno <= $se['end'])
			{
				$res50 = intval(@$prot["bw"]["$tmrno.50"]) or die("Unknown Ballesteros-Weinstein number: $bw.\n");
				$offset = $resno - $res50 + 50;

				$insdel = bw_insdel($prot, $tmrno, $offset);
				$offset -= $insdel;

				return "$tmrno.$offset";
			}
			else if ($tmrno > 1
					 &&
					 $resno < $se['start']
					 &&
					 $resno > (@$prot['region']['TMR'.($tmrno-1)]['end'] ?: @$prot['region']['HXR'.($tmrno-1)]['end']) 
					)
			{
				$tmr_1 = $tmrno-1;
				if (isset($prot['bw']["{$tmr_1}{$tmrno}.50"]))
				{
					$tmrno = intval("{$tmr_1}{$tmrno}");
					$res50 = intval(@$prot["bw"]["$tmrno.50"]) or die("Unknown Ballesteros-Weinstein number: $bw.\n");
					$offset = $resno - $res50 + 50;

					$insdel = bw_insdel($prot, $tmrno, $offset);
					$offset -= $insdel;

					return "$tmrno.$offset";
				}
			}
		}
	}

	return "-";
}

function family_from_protid($protid)
{
	if (substr($protid, 0, 2) == "OR") return "OR".intval(substr($protid, 2, 2));
	else return substr($protid, 0, 4);
}

function filename_protid($protid)
{
	$fam = family_from_protid($protid);
	return "pdbs/$fam/$protid.upright.pdb";
}

function binding_site($protid)
{
	$retval = [];

	$retval['bsr2a'] = resno_from_bw($protid, "2.53");
	$retval['bsr3a'] = resno_from_bw($protid, "3.29");
	$retval['bsr3b'] = resno_from_bw($protid, "3.32");
	$retval['bsr3c'] = resno_from_bw($protid, "3.33");
	$retval['bsr3d'] = resno_from_bw($protid, "3.36");
	$retval['bsr3e'] = resno_from_bw($protid, "3.37");
	$retval['bsr3f'] = resno_from_bw($protid, "3.40");
	$retval['bsr3g'] = resno_from_bw($protid, "3.41");
	$retval['bsr4a'] = resno_from_bw($protid, "4.53");
	$retval['bsr4b'] = resno_from_bw($protid, "4.57");
	$retval['bsr4c'] = resno_from_bw($protid, "4.60");
	$retval['bsr5a'] = resno_from_bw($protid, "5.39");
	$retval['bsr5b'] = resno_from_bw($protid, "5.43");
	$retval['bsr5c'] = resno_from_bw($protid, "5.46");
	$retval['bsr5d'] = resno_from_bw($protid, "5.47");
	$retval['bsr6a'] = resno_from_bw($protid, "6.48");
	$retval['bsr6b'] = resno_from_bw($protid, "6.51");
	$retval['bsr7a'] = resno_from_bw($protid, "7.38");
	$retval['bsr7b'] = resno_from_bw($protid, "7.39");
	$retval['bsr7c'] = resno_from_bw($protid, "7.42");

	return $retval;
}

function json_encode_pretty($array)
{
	return preg_replace("/([ \t]*)([^\\s]*) ([{\\[])\n/", "\$1\$2\n\$1\$3\n", json_encode($array, JSON_PRETTY_PRINT));
}

function split_pdb_to_rigid_and_flex($protid, $pdblines, $flxr_array)
{
	$rigid = [];
	$flex  = [];

	$flxr_res = [];
	foreach ($flxr_array as $bw)
	{
		if (false===strpos($bw, "."))
		{
			if (intval($bw)) $flxr_res[] = intval($bw);
			continue;
		}
		$resno = resno_from_bw($protid, $bw);
		if ($resno) $flxr_res[] = $resno;
	}

	foreach ($pdblines as $ln)
	{
		if (substr($ln, 0, 5) != "ATOM ") continue;
		$y = floatval(substr($ln, 38, 8));
		$resno = intval(substr($ln, 22, 4));
		if (in_array($resno, $flxr_res)) $flex[] = $ln;
		else $rigid[] = $ln;
	}

	return [$rigid,$flex];
}

$cwd = getcwd();
chdir(__DIR__);
chdir("..");
$prots = json_decode(file_get_contents("data/receptor.json"), true);
$gprots = json_decode(file_get_contents("data/gprot.json"), true);
$treenodes = json_decode(file_get_contents("data/tree_nodes.json"), true);
chdir($cwd);

$aminos = 
[
	'A' => 'Ala',
	'R' => 'Arg',
	'N' => 'Asn',
	'D' => 'Asp',
	'C' => 'Cys',
	'E' => 'Glu',
	'Q' => 'Gln',
	'G' => 'Gly',
	'H' => 'His',
	'I' => 'Ile',
	'L' => 'Leu',
	'K' => 'Lys',
	'M' => 'Met',
	'F' => 'Phe',
	'P' => 'Pro',
	'S' => 'Ser',
	'T' => 'Thr',
	'W' => 'Trp',
	'Y' => 'Tyr',
	'V' => 'Val',
	'O' => 'Pyl',
	'U' => 'Sec',
];
