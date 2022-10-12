<?php

global $odors;

$cwd = getcwd();
chdir(__DIR__);
chdir("..");
$odors = json_decode(file_get_contents("data/odorant.json"), true);
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
			foreach ($o['activity'] as $ref => $acv)
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

function all_empirical_pairs_for_receptor($protein)
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
				if ($array[$oid]['adjusted_curve_top'])
				{
					$value += $array[$oid]['adjusted_curve_top'];
					$samples++;
				}
				if ($array[$oid]['ec50'] && $value >= 0)
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

		if (!file_exists("sdf/$ligname.sdf"))
		{
			$obresult = [];
			exec("which obabel", $obresult);
			if (trim(@$obresult[0]))
			{
				exec("obabel -:'{$o['smiles']}' --gen3D -osdf -Osdf/$ligname.sdf");
			}
			else
			{
				$f = fopen("sdf/$ligname.sdf", "wb");
				if (!$f) die("Unable to create sdf/$ligname.sdf, please ensure write access.\n");
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

	$aroma1 = preg_replace( "/^[a-z0-9]/", "", strtolower($aroma) );
	foreach ($odors as $oid => $o)
	{
		if ( $o['smiles'] == $aroma || preg_replace( "/^[a-z0-9]/", "", $o['full_name'] ) == $aroma1 )
		{
			$retval = $o;
			$retval['oid'] = $oid;
			return $retval;
		}
	}
	return false;
}

