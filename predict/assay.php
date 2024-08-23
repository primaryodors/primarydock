<?php

chdir(__DIR__);
chdir("..");
$cwd = getcwd();
require_once("data/protutils.php");
chdir($cwd);
require_once("data/odorutils.php");
chdir($cwd);

$data = json_decode(file_get_contents("predict/dock_results.json"), true);

$o = @$argv[1] ?: "(Z)-3-hexen-1-ol";
$odor = find_odorant($o);
if (!$odor) die("Odorant not found: $o\n");
$o = $odor['full_name'];
if (@$argv[2])
{
    $started = false;
    $startat = $argv[2];
}
else $started = true;

foreach ($prots as $protid => $p)
{
    if (substr($protid, 0, 4) == "MS4A") break;
    if ($protid == $startat) $started = true;
    if (!$started) continue;
    $cmd = "./run_prediction.sh $protid \"$o\"";
    echo "$cmd\n\n";
    passthru($cmd);
}
