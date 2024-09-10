<?php

chdir(__DIR__);
chdir("..");

require("predict/statistics.php");
require("data/odorutils.php");
require("data/protutils.php");

$colors =
[
    "PAILV" => "\033[38;5;247m",
    "G" => "\033[38;5;243m",
    "M" => "\033[38;5;185m",
    "C" => "\033[38;5;220m",
    "STNQ" => "\033[38;5;49m",
    "DE" => "\033[38;5;196m",
    "KR" => "\033[38;5;27m",
    "H" => "\033[38;5;99m",
    "FWY" => "\033[38;5;171m"
];
$colorless = "\033[0m";

$orphan = "\033[2m\033[3m";
$deorphan = ""; // "\033[1m";
$reset = "\033[22m\033[23m\033[24m";

echo "Legend: {$deorphan}Deorphaned receptor$reset {$orphan}Orphan receptor$reset\n\n";

if (!@$argv[1] || is_numeric($argv[1])) die("Error no receptors specified.\n");
$rcppatt = $argv[1];
$rcppl = strlen($rcppatt);

$bwnos = [];
for ($i=2; $i<$argc; $i++)
{
    $bwnos[] = $argv[$i];
}

echo "          ".implode(" ", $bwnos)."\n";
foreach ($prots as $rcpid => $p)
{
    if (substr($rcpid, 0, $rcppl) != $rcppatt
        && !preg_match("/$rcppatt/", $rcpid)
        ) continue;

    echo str_pad($rcpid, 10);
    foreach ($bwnos as $bw)
    {
        $bwl = strlen($bw)+1;
        $resno = resno_from_bw($rcpid, $bw);
        if ($resno > 0)
        {
            $letter = substr($p['sequence'], $resno-1, 1);
            foreach ($colors as $l => $c) if (false!==strpos($l, $letter)) { echo $c; break; }
            echo str_pad($letter, $bwl);
            echo $colorless;
        }
        else echo str_pad("-", $bwl);
    }
    echo "\n";
}
