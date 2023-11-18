<?php

// protutils.php
//
// Loads receptor data into memory and provides useful protein functions. 
//

global $prots, $aminos;

function find_poid($id)
{
	global $prots;
	if (isset($prots[$id])) return $id;
	foreach ($prots as $k => $p) if ($p['id'] == $id) return $k;
	return false;
}

function find_prot($id)
{
	global $prots;
	if (isset($prots[$id])) return $prots[$id];
	foreach ($prots as $p) if ($p['id'] == $id) return $p;
	return false;
}

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
	$prot = find_prot($protid);
	if (!$prot) die("Protein not found: $protid.\n");
	$protid = $prot['id'];
	
	$pettia = explode(".", $bw);
	$tmrno = intval($pettia[0]);
	$offset = intval($pettia[1]);

	$insdel = bw_insdel($prot, $tmrno, $offset);
	
	$res50 = intval(@$prot["bw"]["$tmrno.50"]);
	if (!$res50) throw new Exception("Unknown Ballesteros-Weinstein number $bw");
	
	return $res50 + $offset - 50 + $insdel;
}

function bw_from_resno($protid, $resno)
{
	global $prots;
	$prot = find_prot($protid);
	if (!$prot) die("Protein not found: $protid.\n");
	$protid = $prot['id'];
	
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
	else if (substr(strtoupper($protid), 0, 3) == "ODR") return "Odr".intval(substr($protid, 3, 2));
	else return substr($protid, 0, 4);
}

function rcpid_cmp($a, $b)
{
	if (substr($a, 0, 4) == "MS4A" && substr($b, 0, 4) != "MS4A") return 1;
	if (substr($a, 0, 4) != "MS4A" && substr($b, 0, 4) == "MS4A") return -1;

	if (substr($a, 0, 2) > substr($b, 0, 2)) return 1;
	if (substr($a, 0, 2) < substr($b, 0, 2)) return -1;

	$fam1 = intval(preg_replace("/[^0-9]/", "", family_from_protid($a)));
	$fam2 = intval(preg_replace("/[^0-9]/", "", family_from_protid($b)));

	if ($fam1 > 20 && $fam1 < 30) $fam1 = floatval($fam1) / 10;
	if ($fam2 > 20 && $fam2 < 30) $fam2 = floatval($fam2) / 10;

	if ($fam1 >= 16 && $fam1 <= 18) $fam1 = floatval($fam1-10) / 10 + 6;
	if ($fam2 >= 16 && $fam2 <= 18) $fam2 = floatval($fam2-10) / 10 + 6;

	if ($fam1 < $fam2) return -1;
	else if ($fam1 > $fam2) return 1;

	$fam1 = family_from_protid($a);
	$fam2 = family_from_protid($b);

	$sub1 = preg_replace("/[0-9]/", "", substr($a, strlen($fam1)));
	$sub2 = preg_replace("/[0-9]/", "", substr($b, strlen($fam2)));

	if (strlen($sub1) < strlen($sub2)) return -1;
	else if (strlen($sub1) > strlen($sub2)) return 1;

	if ($sub1 < $sub2) return -1;
	else if ($sub1 > $sub2) return 1;

	$a = intval(substr($a, strlen($fam1) + strlen($sub1)));
	$b = intval(substr($b, strlen($fam2) + strlen($sub2)));

	if ($a < $b) return -1;
	else if ($a > $b) return 1;
	else return 0;
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

$cwd = getcwd();
chdir(__DIR__);
chdir("..");
$prots = json_decode(file_get_contents("data/receptor.json"), true);
$gprots = json_decode(file_get_contents("data/gprot.json"), true);
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
