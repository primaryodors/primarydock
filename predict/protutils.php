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
	
	$res50 = intval(@$prots[$protid]["bw"]["$tmrno.50"]) or die("Unknown Ballesteros-Weinstein number: $bw.\n");
	
	return $res50 + $offset - 50 + $insdel;
}

function bw_from_resno($protid, $resno)
{
	global $prots;
	if (!isset($prots[$protid])) die("Protein not found: $protid.\n");

	$prot = $prots[$protid];
	
	foreach ($prot['region'] as $rgn => $se)
	{
		if (substr($rgn, 0, 3) == 'TMR')
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
			else if ($tmrno > 1 && ($tmrno == 7 || ($resno < @$prot['region']['TMR'.($tmrno+1)]['start']) ) )
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

	return 0;
}

function family_from_protid($protid)
{
	if (substr($protid, 0, 2) == "OR") return "OR".intval(substr($protid, 2, 2));
	else return substr($protid, 0, 4);
}

function json_encode_pretty($array)
{
	return preg_replace("/([ \t]*)([^\\s]*) ([{\\[])\n/", "\$1\$2\n\$1\$3\n", json_encode($array, JSON_PRETTY_PRINT));
}

$cwd = getcwd();
chdir(__DIR__);
chdir("..");
$prots = json_decode(file_get_contents("data/receptor.json"), true);
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