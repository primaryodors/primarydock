<?php

/*
// Axial tumble
// No mol fullrot iter
$filename = "dock_results_sweep6.twcuog.json";
$threshold = 100;
$leniency  = 500;
// OR51E2: 9 correct predictions out of 11 total (6 agonists / 5 other), 1 false positives / 1 false negatives, 81.82% accuracy.

*/

$filename = "dock_results_sweep6.json";
$threshold = 100;
$leniency  = 500;

chdir(__DIR__);
$data = json_decode(file_get_contents($filename), true);

$analysis = [];

foreach ($data as $rcpid => $ligdata)
{
    foreach ($ligdata as $ligand => $r)
    {
        for ($t = 2; $t <= 7; $t++)
        {
            if (isset($r["TM{$t}_d"]))
            {
                $v = $r["TM{$t}_d"];
                $a = @$r["Actual"];
                if (!$a) continue;

                if (!isset($analysis["TMR$t"][$a]["min"])) $analysis["TMR$t"][$a]["min"] = $v;
                else if ($analysis["TMR$t"][$a]["min"] > $v) $analysis["TMR$t"][$a]["min"] = $v;
                
                if (!isset($analysis["TMR$t"][$a]["max"])) $analysis["TMR$t"][$a]["max"] = $v;
                else if ($analysis["TMR$t"][$a]["max"] < $v) $analysis["TMR$t"][$a]["max"] = $v;   
                
                if (!isset($analysis["TMR$t"][$a]["avg"])) $analysis["TMR$t"][$a]["avg"] = $v;
                else $analysis["TMR$t"][$a]["avg"] += $v;    
                
                if (!isset($analysis["TMR$t"][$a]["count"])) $analysis["TMR$t"][$a]["count"] = 1;
                else $analysis["TMR$t"][$a]["count"]++;
            }
        }
    }
}

foreach ($analysis as $region => $rdat)
{
    foreach ($rdat as $a => $r)
    {
        if (@$r["count"]) $analysis[$region][$a]["avg"] /= $r["count"];
    }
}

print_r($analysis);

$correct = $total = 0;
$agonists = $nonagonists = 0;
$fp = $fn = 0;
foreach ($data as $rcpid => $ligdata)
{
    foreach ($ligdata as $ligand => $r)
    {
        if (!preg_match("/_acid$/", $ligand)) continue;

        $quantified = (-$r["TM6_d"]) * ($leniency + max(-@$r["TM2_d"], -@$r["TM3_d"], -@$r["TM4_d"], -@$r["TM5_d"])) / 10000;

        $data[$rcpid][$ligand]["Predicted"] = ($quantified > $threshold) ? "Agonist" : "Non-Agonist";

        if (!isset($r["Actual"])) continue;
        if ($r["Actual"] == "Inverse Agonist") $r["Actual"] = "Non-Agonist";
        if ($r["Actual"] == $data[$rcpid][$ligand]["Predicted"]) $correct++;
        $total++;
        if ($r["Actual"] == "Agonist")
        {
            $agonists++;
            if ($r["Actual"] != $data[$rcpid][$ligand]["Predicted"]) $fn++;
        }
        else
        {
            $nonagonists++;
            if ($r["Actual"] != $data[$rcpid][$ligand]["Predicted"]) $fp++;
        }
    }
}

echo "\n";

if (!$total) echo "Not enough data processed.\n\n";
else echo "$correct correct predictions out of $total total ($agonists agonists / $nonagonists other), $fp false positives / $fn false negatives, "
  . round(100.0 * $correct / $total, 2) . "% accuracy.\n\n";




















