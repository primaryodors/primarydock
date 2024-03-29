<?php

// Scan all duck_results and rerun their values through the evaluate_result() fuction.

require("../data/protutils.php");
require("../data/odorutils.php");
require("dock_eval.php");

// Load data
$dock_results = [];
$json_file = "predict/dock_results.json";

chdir(__DIR__);
chdir("..");

if (file_exists($json_file))
{
	$dock_results = json_decode(file_get_contents($json_file), true);
}
else die("No results to scan. Please run method_combined.php for several empirical pairs and retry.\n");

$right = 0;
$wrong = 0;
foreach ($dock_results as $protid => $docks)
{
    foreach ($docks as $ligand => $array)
    {
    	$o = find_odorant($ligand);
    	if (!$o)
    	{
    		echo "Warning: $ligand not found in odorants data.\n";
    		continue;
    	}

    	$r = empirical_response($protid, $o);

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

        $array = evaluate_result($array);
		$prediction = $array['Prediction'];
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

