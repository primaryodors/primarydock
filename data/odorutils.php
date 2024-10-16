<?php

global $odors, $refs;

$cwd = getcwd();
chdir(__DIR__);
chdir("..");
$odors = json_decode(file_get_contents("data/odorant.json"), true);
$refs = json_decode(file_get_contents("data/reference.json"), true);
chdir($cwd);

foreach ($odors as $oid => $odor)
{
	$odors[$oid]["oid"] = $oid;
	if (!empty($odor["activity"])) foreach ($odor["activity"] as $url => $acv)
	{
		if (@$refs[$url]["hidden"]) unset($odors[$oid]["activity"][$url]);
	}

	if (!empty($odor["aroma"])) foreach ($odor["aroma"] as $url => $aroma)
	{
		if (count($aroma) == 1 && false!==strpos($aroma[0], " ")) $odors[$oid]['aroma'][$url] = explode(" ", $aroma[0]);
	}
}

function trim_prefixes($what)
{
	return
		preg_replace("/[()]/", "",
		preg_replace(
		"/^(([0-9-]*[(][ERSZ0-9, +-]+[)]?[-,]?)|(([0-9A-Z,]+|[a-z])-)|((alpha|beta|gamma|delta|epsilon|zeta|eta|theta|iota|kappa|lambda|mu|nu|xi|omicron|pi|rho|sigma|tau|upsilon|phi|chi|psi|omega|o|ortho|m|meta|p|para|cis|trans)-))+/"
		, "", $what));
}

function odorcmp($a, $b)
{
	$A = strtolower(trim_prefixes($a["full_name"]));
	$B = strtolower(trim_prefixes($b["full_name"]));

	if ($A == $B)
	{
		$A = $a["full_name"];
		$B = $b["full_name"];
	}

	if ($A == $B) return 0;
	else return ($A > $B) ? 1 : -1;
}

usort($odors, "odorcmp");

foreach ($odors as $k => $odor)
{
	$oid = $odor["oid"];
	$odors[$oid] = $odor;
	unset($odors[$k]);
}

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

$typenames =
[
	0 => "non-agonist",
	1 => "v. weak agonist",
	2 => "weak agonist",
	3 => "moderate agonist",
	4 => "strong agonist",
	5 => "v. strong agonist",
   -1 => "inverse agonist",
   -999 => "unknown",
];

$sepyt = array_flip($types);

function hellenicize($what)
{
	$what = preg_replace("/alpha([,-])/",   "&#x3B1;$1", $what);
    $what = preg_replace("/beta([,-])/",    "&#x3B2;$1", $what);
    $what = preg_replace("/gamma([,-])/",   "&#x3B3;$1", $what);
    $what = preg_replace("/delta([,-])/",   "&#x3B4;$1", $what);
    $what = preg_replace("/epsilon([,-])/", "&#x3B5;$1", $what);
    $what = preg_replace("/omega([,-])/",   "&#x3C9;$1", $what);
	return $what;
}

function best_empirical_pair($protein, $aroma, $as_object = false)
{
	global $odors, $sepyt;

	$btyp = $sepyt["?"];
	$obj = [];

	$o = find_odorant($aroma);

	if ($o)
	{
		if (@$o['activity']) foreach ($o['activity'] as $ref => $acv)
		{
			if (isset($acv[$protein]))
			{
				// echo "<!-- Reference $ref says {$o['full_name']} is ec50 {$acv[$protein]['ec50']} for $protein. -->\n";

				if (isset($acv[$protein]['adjusted_curve_top']) &&
					(!isset($obj['adjusted_curve_top']) || $obj['adjusted_curve_top'] < floatval($acv[$protein]['adjusted_curve_top'])))
				{
					$obj['adjusted_curve_top'] = floatval($acv[$protein]['adjusted_curve_top']);
					$obj['reft'] = $ref;
					$btyp = @$sepyt[trim($acv[$protein]['type'])];
				}
				if (isset($acv[$protein]['ec50']) &&
					(!isset($obj['ec50']) || ($obj['ec50'] > floatval($acv[$protein]['ec50']))))
				{
					$obj['ec50'] = floatval($acv[$protein]['ec50']);
					$obj['refe'] = $ref;
					$btyp = @$sepyt[trim($acv[$protein]['type'])];
				}
			}
		}
	}

	$obj['ref'] = @$obj['refe'] ?: @$obj['reft'];
	unset($obj['reft']);
	unset($obj['refe']);
	if (@$obj['ec50'] < 0 && isset($obj['adjusted_curve_top']) && @$obj['adjusted_curve_top'] <= 0) unset($obj['adjusted_curve_top']);

	return $as_object ? $obj : $btyp;
}

function empirical_response($protein, $odorobj)
{
	global $sepyt;
	
	$retval = false;
	if (isset($odorobj['activity'])) foreach ($odorobj['activity'] as $ref => $acv)
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

function make_clickable_notes($notes)
{
	if ($notes[0] == "(insufficient data)") return $notes;
	$result = [];
	foreach ($notes as $key => $note) $result[$key] = "<a href=\"note.php?n=$note\">$note</a>";
	return $result;
}

function all_empirical_pairs_for_receptor($protein, $return_1dim = false, $agonists_only = false)
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
				if ($agonists_only)
				{
					if (isset($acv[$protein]['adjusted_curve_top']) && @$acv[$protein]['adjusted_curve_top'] <= 0) continue;
					if (@$acv[$protein]['type'] == "na" || @$acv[$protein]['type'] == "ia") continue;
				}

				if (!isset($array[$oid]))
				{
					if (@$acv[$protein]['type'] == 'na') $acv[$protein]['adjusted_curve_top'] = 0;
					$array[$oid] = $acv[$protein];
					if (isset($array[$oid]['adjusted_curve_top'])) $array[$oid]['top_ref'] = $ref;
					if (isset($array[$oid]['ec50'])) $array[$oid]['ec50_ref'] = $ref;
				}
				else
				{
					if (@$acv[$protein]['adjusted_curve_top'] > (@$array[$oid]['adjusted_curve_top']?:0))
					{
						$array[$oid]['adjusted_curve_top'] = $acv[$protein]['adjusted_curve_top'];
						$array[$oid]['top_ref'] = $ref;
					}
					if (isset($acv[$protein]['ec50']) && @$acv[$protein]['ec50'] < (@$array[$oid]['ec50'] ?: 0))
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
					$value -= $array[$oid]['ec50']*1.6;
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

function check_isomers($ligname, $randomize=true)
{
    $odor = find_odorant($ligname);
	if (!$odor) die("Odorant not found $ligname.\n");
	if (!isset($odor["isomers"])) return false;

	$result = [];
	foreach (array_keys($odor["isomers"]) as $iso) $result[] = "$iso-{$odor["full_name"]}";
	return $result;
}

function ensure_sdf_exists($ligname)
{
	global $odors;

	$o = find_odorant($ligname);
	if (!$o) die("Odorant not found $ligname.\n");
	$isomers = check_isomers($ligname);
	$lignamei = $isomers ? $isomers[0] : $ligname;

	$pwd = getcwd();
	chdir(__DIR__);
	chdir("..");
	$sdfname = str_replace(" ", "_", "sdf/$lignamei.sdf");
	if (!file_exists($sdfname) || filesize($sdfname) < 20)
	{
		if (isset($o["isomers"]))
		{
			$fullname = str_replace(" ", "_", $o["full_name"]);
			$first = true;
			foreach ($o["isomers"] as $iso => $ismiles)
			{
				$parts = explode("|", $ismiles);
				$ismiles = $parts[0];
				$isofname = escapeshellarg("sdf/$iso-$fullname.sdf");
				exec("obabel -:\"$ismiles\" --gen3D -osdf -O$isofname");
				if (@$parts[1])
				{
					$sub4 = substr($parts[1], 0, 4);
					$rest = substr($parts[1], 4);
					switch ($sub4)
					{
						case "rflp":
						$cmd = "bin/ringflip $isofname $rest";
						// echo "$cmd\n";
						exec($cmd);
						break;

						default:
						die("Unknown modifier $sub4.\n");
					}
				}

				if ($first)
				{
					copy("sdf/$iso-$fullname.sdf", "sdf/$fullname.sdf");
					$first = false;
				}
			}
		}
		else
		{
			$smiles = $o["smiles"];
			$fullname = str_replace(" ", "_", $o["full_name"]);
			$ffname = escapeshellarg("sdf/$fullname.sdf");
			$cmd = "obabel -:\"$smiles\" --gen3D -osdf -O$ffname";
			// echo "$cmd\n";
			exec($cmd);
		}
	}
	chdir($pwd);
	return;	
}

function get_at_wt($element)
{
	static $elements = [];
	if (!count($elements))
	{
		$cwd = getcwd();
		chdir(__DIR__);
		chdir("..");
		$c = file_get_contents("data/elements.dat");
		$c = explode("\n", $c);
		for ($i=1; $i<256; $i++)
		{
			if (!isset($c[$i])) break;
			$ln = explode(" ", preg_replace("/\\s+/", " ", trim($c[$i])));

			$esym = $ln[1];
			$elements[$esym] = $ln;
		}
		chdir($cwd);
	}

	if (!isset($elements[$element])) return 0;
	return $elements[$element][4];
}

function get_mol_wt($odorant)
{
	$o = find_odorant($odorant);
	ensure_sdf_exists($o['full_name']);
	$sdfname = str_replace(" ", "_", "sdf/{$o['full_name']}.sdf");

	$cwd = getcwd();
	chdir(__DIR__);
	chdir("..");
	$c = file_get_contents($sdfname) or die("FAILED to open $sdfname for reading.");

	$retval = 0.0;
	$c = explode("\n", $c);
	$meta = explode(" ", preg_replace("/\\s+/", " ", trim($c[3])));
	$j = intval($meta[0]) + 4;
	for ($i=4; $i<$j; $i++)
	{
		$ln = explode(" ", preg_replace("/\\s+/", " ", trim($c[$i])));
		$atwt = get_at_wt($ln[3]);
		$retval += $atwt;
	}

	chdir($cwd);
	return $retval;
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
		$namematch = false;
		for ($i=1; isset($o["name$i"]); $i++)
		{
			if ($o["name$i"] == $aroma) $namematch = true;
		}

		if ( $o['smiles'] == $aroma
			 || preg_replace( "/[^a-z0-9]/", "", strtolower($o['full_name']) ) == $aroma1
			 || @$o['iupac'] == $aroma
			 || $namematch
			)
		{
			$retval = $o;
			$retval['oid'] = $oid;
			return $retval;
		}
	}
	return false;
}

