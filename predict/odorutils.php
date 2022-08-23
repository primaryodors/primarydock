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

function best_empirical_pair($protein, $aroma)
{
	global $odors, $sepyt;
	
	$btyp = 0;
	
	foreach ($odors as $oid => $o)
	{
		if ($o['full_name'] == $aroma)
		{
			foreach ($o['activity'] as $ref => $acv)
			{
				if (isset($acv[$protein]))
				{
					if ($sepyt[$acv[$protein]['type']] > $btyp || !$btyp) $btyp = $sepyt[$acv[$protein]['type']];
				}
			}
			
			return $btyp;
		}
	}
}
