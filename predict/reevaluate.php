<?php

// Scan all dock_results and rerun their values through the evaluate_result() fuction.

chdir(__DIR__);
require("../data/protutils.php");
chdir(__DIR__);
require("../data/odorutils.php");
chdir(__DIR__);
require("algorithm_affinity.php");
chdir(__DIR__);
chdir("..");

// Load data
$dock_results = [];
$json_file = "predict/dock_results.json";
$date_cutoff = 1723400000;

if (file_exists($json_file))
{
	$dock_results = json_decode(file_get_contents($json_file), true);
}
else die("No results to scan. Please run predictions for several empirical pairs and retry.\n");

$right = 0;
$wrong = 0;
foreach ($dock_results as $protid => $docks)
{
    foreach ($docks as $ligand => $array)
    {
		if (@$array['version'] < $date_cutoff)
		{
			unset($dock_results[$protid][$ligand]);
			if (!count($dock_results[$protid])) unset($dock_results[$protid]);
			continue;
		}

    	$o = find_odorant($ligand);
    	if (!$o)
    	{
    		echo "Warning: $ligand not found in odorants data.\n";
    		continue;
    	}

    	$r = empirical_response($protid, $o);
		$ligname = $o['full_name'];

    	if (!$r)
    		$actual = "(unknown)";
    	else
    	{
			$a = is_agonist($r);
			if ($a > 0) $actual = "Agonist";
			else if ($a < 0) $actual = "Inverse Agonist";
			else $actual = "Non-Agonist";
    	}
    	$array["Actual"] = $actual;

        $array = make_prediction($array);
		$prediction = $array['Predicted'];
        if ($prediction == $actual) $right++;
        else $wrong++;

        $dock_results[$protid][$ligand] = $array;
    }
}

// if (!$right) die("No correct predictions made; exiting.\n");
$percent = round(100.0 * $right / ($right+$wrong), 1);
echo "$percent% success rate.\n";

chdir(__DIR__);
chdir("..");
$bak_file = str_replace(".json", ".bak.json", $json_file);
if (file_exists($bak_file)) unlink($bak_file);
rename($json_file, $bak_file);

$f = fopen($json_file, "wb");
if (!$f) die("File write FAILED. Make sure have access to write $json_file.");
fwrite($f, json_encode_pretty($dock_results));
fclose($f);

