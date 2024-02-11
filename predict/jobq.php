<?php

$procs = [];
exec("ps -ef | grep obabel | grep -v grep", $procs);
if (count($procs)) die("Waiting for obabel.");

$procs = [];
exec("ps -ef | grep bin/fyg_activate_or | grep -v grep", $procs);
if (count($procs)) die("Waiting for activation.");

$max_concurrent = 2;
$queue = [];

chdir(__DIR__);
chdir("..");
$dock_results = json_decode(file_get_contents("predict/dock_results.json"), true);
$version = max(filemtime("bin/primarydock"),
    filemtime("predict/methods_common.php"),
    filemtime("predict/method_directmdl.php"),
    filemtime("predict/method_fygactive.php")
    );
$lines = explode("\n", @file_get_contents("jobq"));

foreach ($lines as $ln)
{
    $ln = trim($ln);
    if (!$ln) continue;
    if (substr($ln, 0, 1) == '#') continue;
    $pieces = explode(" ", $ln);

    if ($pieces[0] == "MAX") $max_concurrent = intval($pieces[1]);
    else if ($pieces[0] == "PRDT") $queue[] = [$pieces[1], $pieces[2]];
}
if (!count($queue)) die("No jobs.\n");

$procs = [];
exec("ps -ef | grep \" bin/primarydock \" | grep -v grep", $procs);
if (count($procs) >= $max_concurrent) die("Maximum processes.\n");

// Example:
// primary+ 1342875 1342799 94 18:09 pts/0    01:24:43 bin/primarydock tmp/prediction.OR51E2_i.pelargonic_acid.config
foreach ($queue as $q)
{
    $already = false;
    $prot = $q[0];
    $lig = $q[1];

    if (isset($dock_results[$prot][$lig]))
    {
        if (intval($dock_results[$prot][$lig]['version']) >= $version) continue;
    }

    foreach ($procs as $p)
    {
        if (false===strpos($p, ".{$prot}_")) continue;
        if (false===strpos($p, ".{$lig}.")) continue;
        $already = true;
        break;
    }

    if (!$already)
    {
        passthru("./run_prediction.sh $prot $lig &");
        exit;
    }
}

die("All done.\n");

