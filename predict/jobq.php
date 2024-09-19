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
    global $dock_results, $version, $dbg;
    $lig = trim(str_replace(" ", "_", $lig));
    if ($dbg) echo "Checking $prot $lig \n";
    if (!isset($dock_results[$prot][$lig])) return false;
    if ($dbg) echo "$prot $lig {$dock_results[$prot][$lig]['version']} $version\n";
    if (@intval($dock_results[$prot][$lig]['version']) < $version) return false;
    return isset($dock_results[$prot][$lig]['Predicted']);
}

function running($prot, $lig)
{
    global $procs, $dbg;
    $procs = [];
    exec("ps -ef | grep \" bin/primarydock \" | grep -v \"sh -c \" | grep -v grep", $procs);

    foreach ($procs as $p)
    {
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
            if ($dbg) echo uptodate($prot, $lig) ? "$prot $lig is up to date.\n" : "$prot $lig is NOT up to date.\n";
            if ($dbg) echo running($prot, $lig) ? "$prot $lig is running.\n" : "$prot $lig is NOT running.\n";
            if (!uptodate($prot, $lig)) if (!running($prot, $lig)) break;
        }
    }

    $lig = @find_odorant($lig)["full_name"];
    if (!$lig)
    {
        if ($dbg) echo "Ligand not found.\n";
        continue;
    }

    if (uptodate($prot, $lig)) continue;
    $lig = trim(str_replace(" ", "_", $lig));
    $already = running($prot, $lig);

    if (!$already)
    {
        if (!$dbg) passthru("./run_prediction.sh $prot '$lig' &");
        exit;
    }
}

die("All done.\n");

