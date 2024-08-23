<?php

chdir(__DIR__);
chdir("..");
$cwd = getcwd();
require_once("data/protutils.php");
chdir($cwd);
require_once("data/odorutils.php");
chdir($cwd);

$data = json_decode(file_get_contents("predict/dock_results.json"), true);

$a = @$argv[1] ?: "(Z)-3-hexen-1-ol";

$prot = false;
$odor = find_odorant($a);
if (!$odor) $prot = @$prots[$a];
if (!$odor && !$prot) die("Unrecognized input $a. Please specify an odorant or a receptor.\n");
if ($prot && substr($prot['id'], 0, 4) == "MS4A") die("No prediction method exists for MS4A receptors.\n");

if ($odor)
{
    foreach ($prots as $protid => $p)
    {
        if (substr($protid, 0, 4) == "MS4A") break;
        if ($protid == $startat) $started = true;
        if (!$started) continue;
        $cmd = "./run_prediction.sh $protid \"$o\"";
        echo "$cmd\n\n";
        passthru($cmd);
    }
}
else if ($prot)
{
    foreach (all_empirical_pairs_for_receptor($prot['id'], true, false) as $oid => $acv)
    {
        //
    }
}
