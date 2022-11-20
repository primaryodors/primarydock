<?php

global $odors, $refs;

$cwd = getcwd();
chdir(__DIR__);
chdir("..");
$odors = json_decode(file_get_contents("data/odorant.json"), true);
$refs = json_decode(file_get_contents("data/reference.json"), true);
chdir($cwd);

$types =
[
	0 => "na",
	1 => "vwa",
	2 => "wa",
	3 => "ma",
	4 => "sa",
	5 => "vsa",
   -1 => "ia",
   -999 => "?",
];

$sepyt = array_flip($types);

function best_empirical_pair($protein, $aroma, $as_object = false)
{
	global $odors, $sepyt;
	
	$btyp = $as_object ? false : $sepyt["?"];
	
	foreach ($odors as $oid => $o)
	{
		if ($oid == $aroma
			||
			$o['full_name'] == $aroma
			||
			str_replace(" ", "_", $o['full_name']) == $aroma
		)
		{
			if (@$o['activity']) foreach ($o['activity'] as $ref => $acv)
			{
				if (isset($acv[$protein]))
				{
					// echo "Reference $ref says $aroma is a {$acv[$protein]['type']} for $protein.\n";

					// TODO: This is horribly inadequate for $as_object=true.
					if ($sepyt[$acv[$protein]['type']] > $btyp || $btyp == $sepyt["?"])
					{
						$btyp = $as_object ? $acv[$protein] : $sepyt[trim($acv[$protein]['type'])];
						if ($as_object) $btyp['ref'] = $ref;
						// echo "Ligand type set to $btyp.\n";
					}
				}
			}
			
			return $btyp;
		}
	}
}

function empirical_response($protein, $odorobj)
{
	global $sepyt;
	
	$retval = false;
	foreach ($odorobj['activity'] as $ref => $acv)
	{
		if (isset($acv[$protein]))
		{
			if (!$retval) $retval = $acv[$protein];
			else
			{
				if (isset($acv[$protein]['adjusted_curve_top']) && $acv[$protein]['adjusted_curve_top'] > @$retval['adjusted_curve_top'])
					$retval['adjusted_curve_top'] = $acv[$protein]['adjusted_curve_top'];
				if (isset($acv[$protein]['ec50']) && $acv[$protein]['ec50'] < @$retval['ec50'])
					$retval['ec50'] = $acv[$protein]['ec50'];
				if (isset($acv[$protein]['type']) && isset($retval['type']) && $sepyt[$acv[$protein]['type']] > $sepyt[$retval['type']])
					$retval['type'] = $acv[$protein]['type'];
			}
		}
	}
	
	return $retval;
}

function is_agonist($response)
{
	if (isset($response['adjusted_curve_top']))
	{
		if ($response['adjusted_curve_top'] > 0) return 1;
		else if ($response['adjusted_curve_top'] < 0) return -1;
		else return 0;
	}
	if (isset($response['type']) && $response['type'] != '?')
	{
		if ($response['type'] == 'ia') return -1;
		if ($response['type'] == 'na') return 0;
		return 1;
	}
	if (isset($response['ec50']))
	{
		if ($response['ec50'] < 0) return 1;
		else return 0;
	}
	
	return 0;
}

function has_antagonists($protein)
{
	global $odors;
	
	foreach ($odors as $oid => $o)
	{
		if (isset($o['activity'])) foreach ($o['activity'] as $ref => $acv)
		{
			if (isset($acv[$protein]))
			{
				if (@$acv[$protein]['antagonist']) return true;
			}
		}
	}

	return false;
}

function all_empirical_pairs_for_receptor($protein, $return_1dim = false)
{
	global $odors;

	$array = [];
	$sortable = [];
	
	foreach ($odors as $oid => $o)
	{
		if (isset($o['activity'])) foreach ($o['activity'] as $ref => $acv)
		{
			if (isset($acv[$protein]))
			{
				if (!isset($array[$oid]))
				{
					if ($acv[$protein]['type'] == 'na') $acv[$protein]['adjusted_curve_top'] = 0;
					$array[$oid] = $acv[$protein];
					if (isset($array[$oid]['adjusted_curve_top'])) $array[$oid]['top_ref'] = $ref;
					if (isset($array[$oid]['ec50'])) $array[$oid]['ec50_ref'] = $ref;
				}
				else
				{
					if (@$acv[$protein]['adjusted_curve_top'] > $array[$oid]['adjusted_curve_top'])
					{
						$array[$oid]['adjusted_curve_top'] = $acv[$protein]['adjusted_curve_top'];
						$array[$oid]['top_ref'] = $ref;
					}
					if (@$acv[$protein]['ec50'] < $array[$oid]['ec50'])
					{
						$array[$oid]['ec50'] = $acv[$protein]['ec50'];
						$array[$oid]['ec50_ref'] = $ref;
					}
				}

				$value = 0.0;
				$samples = 0;
				if (@$array[$oid]['adjusted_curve_top'])
				{
					$value += $array[$oid]['adjusted_curve_top'];
					$samples++;
				}
				if (@$array[$oid]['ec50'] && $value >= 0)
				{
					$value -= $array[$oid]['ec50']*1.666;
					$samples++;					
				}

				if ($samples) $value /= $samples;
				$sortable[$oid] = $value;
			}
		}
	}

	arsort($sortable);
	if ($return_1dim) return $sortable;

	$retval = [];

	foreach ($sortable as $oid => $value)
	{
		$retval[$oid] = $array[$oid];
	}

	return $retval;
}

function ensure_sdf_exists($ligname)
{
	global $odors;

	foreach ($odors as $o)
	{
		$full_name = str_replace(" ", "_", $o['full_name']);
		if ($ligname && $ligname != $full_name) continue;
		
		$sdfname = "sdf/$ligname.sdf";
		if (!file_exists($sdfname) || filesize($sdfname) < 10)
		{
			$obresult = [];
			exec("which obabel", $obresult);
			if (trim(@$obresult[0]))
			{
				exec("obabel -:'{$o['smiles']}' --gen3D -osdf -O$sdfname");
			}
			else
			{
				$f = fopen($sdfname, "wb");
				if (!$f) die("Unable to create $sdfname, please ensure write access.\n");
				$sdfdat = file_get_contents("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{$o['smiles']}/SDF?record_type=3d");
				fwrite($f, $sdfdat);
				fclose($f);
			}
		}
		return;
	}

	die("Odorant not found $ligname.\n");
}

function find_odorant($aroma)
{
	global $odors;
	if (!$aroma) return false;

	if (isset($odors[$aroma]))
	{
		$retval = $odors[$aroma];
		$retval['oid'] = $aroma;
		return $retval;
	}

	$aroma1 = preg_replace( "/[^a-z0-9]/", "", strtolower($aroma) );
	foreach ($odors as $oid => $o)
	{
		if ( $o['smiles'] == $aroma || preg_replace( "/[^a-z0-9]/", "", $o['full_name'] ) == $aroma1 )
		{
			$retval = $o;
			$retval['oid'] = $oid;
			return $retval;
		}
	}
	return false;
}

