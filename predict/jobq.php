<?php

$procs = [];
exec("ps -ef | grep obabel | grep -v grep", $procs);
if (count($procs)) die("Waiting for obabel.");

$dbg = false;
foreach ($argv as $a) if ($a == 'dbg') $dbg = true;

$max_concurrent = 2;
$queue = [];

chdir(__DIR__);
chdir("..");
require_once("data/odorutils.php");
chdir(__DIR__);
chdir("..");
$dock_results = json_decode(file_get_contents("predict/dock_results.json"), true);
$version = max(filemtime("bin/primarydock"),
    filemtime("predict/methods_common.php"),
    filemtime("predict/method_directmdl.php")
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
exec("ps -ef | grep \" bin/primarydock \" | grep -v \"sh -c \" | grep -v grep", $procs);
if (count($procs) >= $max_concurrent) die("Maximum processes.\n");

function uptodate($prot, $lig)
{
    global $dock_results, $version;
    return isset($dock_results[$prot][$lig])
        && intval($dock_results[$prot][$lig]['version']) >= $version
        && isset($dock_results[$prot][$lig]['Predicted']);
}

function running($prot, $lig)
{
    global $procs, $dbg;
    
    foreach ($procs as $p)
    {
        if ($dbg) echo "$p\n$prot | $lig\n\n";
        if (false!==strpos($p, ".{$prot}_") && false!==strpos($p, ".{$lig}.")) return true;
    }

    return false;
}

// Example:
// primary+ 1342875 1342799 94 18:09 pts/0    01:24:43 bin/primarydock tmp/prediction.OR51E2_i.pelargonic_acid.config
foreach ($queue as $q)
{
    $already = false;
    $prot = $q[0];
    $lig = $q[1];

    if ($lig == "*")
    {
        $pairs = all_empirical_pairs_for_receptor($prot, true, false);
        foreach (array_keys($pairs) as $oid)
        {
            $lig = $odors[$oid]['full_name'];
            if (!uptodate($prot, $lig)) if (!running($prot, $lig)) break;
        }
    }

    $lig = @find_odorant($lig)["full_name"];
    if (!$lig) continue;
    $lig = trim(str_replace(" ", "_", $lig));

    if (uptodate($prot, $lig)) continue;

    $already = running($prot, $lig);
    if ($dbg) echo "already? $already\n";

    if (!$already)
    {
        if (!$dbg) passthru("./run_prediction.sh $prot '$lig' &");
        exit;
    }
}

die("All done.\n");

