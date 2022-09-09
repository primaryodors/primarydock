<?php

global $odors;

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

function best_empirical_pair($protein, $aroma)
{
	global $odors, $sepyt;
	
	$btyp = $sepyt["?"];
	
	foreach ($odors as $oid => $o)
	{
		if ($o['full_name'] == $aroma
			||
			str_replace(" ", "_", $o['full_name']) == $aroma
		)
		{
			foreach ($o['activity'] as $ref => $acv)
			{
				if (isset($acv[$protein]))
				{
					echo "Reference $ref says $aroma is a {$acv[$protein]['type']} for $protein.\n";
					if ($sepyt[$acv[$protein]['type']] > $btyp || $btyp == $sepyt["?"])
					{
						$btyp = $sepyt[trim($acv[$protein]['type'])];
						// echo "Ligand type set to $btyp.\n";
					}
				}
			}
			
			return $btyp;
		}
	}
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

