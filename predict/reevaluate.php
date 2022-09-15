<?php

// Scan all duck_results and rerun their values through the evaluate_result() fuction.

require("protutils.php");
require("odorutils.php");
require("dock_eval.php");

// Load data
$dock_results = [];
$json_file = "predict/dock_results.json";

$odors = json_decode(file_get_contents("data/odorant.json"), true);

if (file_exists($json_file))
{
	$dock_results = json_decode(file_get_contents($json_file), true);
}
else die("No results to scan. Please run method_combined.php for several empirical pairs and retry.\n");

$right = 0;
$wrong = 0;
foreach ($dock_results as $protid => $docks)
{
    if (substr($protid, 0, 4) != "TAAR") continue;
    foreach ($docks as $ligand => $array)
    {
        if (isset($array['Node 4'])) $pocketnode = 3;
        else $pocketnode = 2;
        $prediction = evaluate_result($array, $pocketnode);
        if ($prediction == $array["Actual"]) $right++;
        else $wrong++;

        $dock_results[$protid][$ligand]['Prediction'] = $prediction;
    }
}

if (!$right) die("No correct predictions made; exiting.\n");
$percent = round(100.0 * $right / ($right+$wrong), 1);
echo "$percent% success rate.\n";

$f = fopen($json_file, "wb");
if (!$f) die("File write FAILED. Make sure have access to write $json_file.");
fwrite($f, json_encode_pretty($dock_results));
fclose($f);

