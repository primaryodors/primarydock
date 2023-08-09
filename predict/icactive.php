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
    if ($inactive_pairs[$pairid][0] > 5 && $pair[0] > 5) continue;

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

$get_bw50_resnos = "";
for ($i=1; $i<=7; $i++) $get_bw50_resnos .= "ECHO \"$i: \" %$i.50\n";

$gpcr_contacts = [ "3x40-6x48", "5x51-6x44", "5x55-6x41", "3x43-7x49", "3x43-7x53", "5x58-6x40", "3x46-7x53", "1x53-7x54", "5x62-6x37", "3x50-7x53" ];
$measure_contacts = "";
foreach ($gpcr_contacts as $contacts)
{
    $bwnos = explode('-', $contacts);
    $measure_contacts .= <<<heredoc
LET @distance = @A.{$bwnos[0]} - @A.{$bwnos[1]}
LET &active_distance = @distance

LET @distance = @I.{$bwnos[0]} - @I.{$bwnos[1]}
LET &inactive_distance = @distance

LET &change = &active_distance - &inactive_distance
ECHO "{$bwnos[0]}|{$bwnos[1]}|" &change


heredoc;
}

$pepd_script = <<<heredoc

LOAD $active_pdbfn $active_chain A
BWCENTER

LOAD $inactive_pdbfn $inactive_chain I
BWCENTER

$get_bw50_resnos
ECHO ""

$measure_contacts
ECHO ""


heredoc;

$pepd_fname = "tmp/icactive.pepd";
$fp = fopen($pepd_fname, "w") or die("File access or permissions error.\n");
fwrite($fp, $pepd_script);
fclose($fp);

$result = [];
$cmd = "bin/pepteditor $pepd_fname";
exec($cmd, $result);

print_r($result);
