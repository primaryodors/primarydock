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
];
$sepyt = array_flip($types);

print_r($sepyt);

function best_empirical_pair($protein, $aroma)
{
	global $odors, $sepyt;
	
	$btyp = 0;
	
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
					if ($sepyt[$acv[$protein]['type']] > $btyp || !$btyp)
					{
						$btyp = $sepyt[trim($acv[$protein]['type'])];
						echo "Ligand type set to $btyp.\n";
					}
				}
			}
			
			return $btyp;
		}
	}
}
