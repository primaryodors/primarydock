<?php

$filename = "dock_results_sweep6.json";

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

