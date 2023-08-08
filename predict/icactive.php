<?php

$rcpid = @$argv[1] or die("Usage: php -f predict/icactive.php receptor_id\n");

$energy_threshold = 0.5;
$distance_threshold = 0.1;

if (substr($rcpid, 0, 2) == "OR")
{
    $active_pdbfn = "pdbs/Gprot/hOR51E2-hGNAS2.pdb";
    $active_chain = "R";

    $inactive_pdbfn = "pdbs/OR51/OR51E2.upright.pdb";
    $inactive_chain = "A";
}
else if (substr($rcpid, 0, 4) == "TAAR")
{
    $active_pdbfn = "pdbs/Gprot/mTAAR9-hGNAS2-pea.pdb";
    $active_chain = "R";

    $inactive_pdbfn = "pdbs/mTAAR9.pdb";
    $inactive_chain = "A";
}
else die("Unsupported protein.\n");

if (!file_exists($active_pdbfn) || !file_exists($inactive_pdbfn)) passthru("bin/pepteditor predict.gprot.pepd");

$active_pairs_raw = [];
exec("bin/ic $active_pdbfn $active_chain $energy_threshold", $active_pairs_raw);

$inactive_pairs_raw = [];
exec("bin/ic $inactive_pdbfn $inactive_chain $energy_threshold", $inactive_pairs_raw);

$active_pairs = [];
$inactive_pairs = [];

foreach ($active_pairs_raw as $pair)
{
    $pieces = explode(" ", $pair);
    $pairid = str_replace(':', '', $pieces[0]);
    $energy = floatval($pieces[1]);
    $radius = floatval(preg_replace("/[^0-9.]/", "", $pieces[2]));

    $active_pairs[$pairid] = [$energy, $radius];
}

foreach ($inactive_pairs_raw as $pair)
{
    $pieces = explode(" ", $pair);
    $pairid = str_replace(':', '', $pieces[0]);
    $energy = floatval($pieces[1]);
    $radius = floatval(preg_replace("/[^0-9.]/", "", $pieces[2]));

    $inactive_pairs[$pairid] = [$energy, $radius];
}


$contacts_made = [];
$contacts_broken = [];
$repacked_pairs = [];

foreach ($active_pairs as $pairid => $pair)
{
    if (!isset($inactive_pairs[$pairid])) continue;

    // Filter out pairs whose energy is positive in both the active and inactive states.
    if ($inactive_pairs[$pairid][0] > 10 && $pair[0] > 10) continue;

    $energy_difference = $pair[0] - $inactive_pairs[$pairid][0];
    $radius_difference = $pair[1] - $inactive_pairs[$pairid][1];

    // echo "$pairid: $energy_difference $radius_difference\n";

    // Energy improvements with decrease in radius means contacts made during activation.
    if ($energy_difference < -$energy_threshold && $radius_difference < $distance_threshold)
        $contacts_made[$pairid] = $radius_difference;

    // Energy closer to zero with increase in radius means contacts broken during activation.
    else if ($energy_difference > $energy_threshold && $radius_difference > $distance_threshold)
        $contacts_broken[$pairid] = $radius_difference;

    // Energy changes without a significant change in radius means repacking.
    else if (abs($energy_difference) > $energy_threshold)
        $repacked_pairs[$pairid] = $radius_difference;
}

print_r($contacts_made);

print_r($contacts_broken);

print_r($repacked_pairs);
