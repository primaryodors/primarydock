<?php

chdir(__DIR__);
chdir("..");
require_once("data/protutils.php");

function residues_compatible($aa1, $aa2)
{
    // Hydrophobic.
    if (false!==strpos("MAILVGCPHFWY", $aa1) && false!==strpos("MAILVGCPHFWY", $aa2)) return true;

    // Polar uncharged.
    if (false!==strpos("STYNQ", $aa1) && false!==strpos("STYNQ", $aa2)) return true;

    if (false!==strpos("STYNQEDHRK", $aa1) && false!==strpos("STYNQEDHRK", $aa2))
    {
        if (false!==strpos("KR", $aa1) && false!==strpos("KR", $aa2)) return false;
        else if (false!==strpos("DE", $aa1) && false!==strpos("DE", $aa2)) return false;
        else return true;
    }

    return false;
}

function get_distance($xyz1, $xyz2)
{
    return sqrt(
        pow($xyz1["x"]-$xyz2["x"], 2)
        + pow($xyz1["y"]-$xyz2["y"], 2)
        + pow($xyz1["z"]-$xyz2["z"], 2)
    );
}

// Note this script is unfinished.

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

if (!file_exists("pdbs/Gprot")) mkdir("pdbs/Gprot");
if (!file_exists($active_pdbfn) || !file_exists($inactive_pdbfn)) passthru("bin/pepteditor predict/gprot.pepd");

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

$get_bw50_resnos = "";
for ($i=1; $i<=7; $i++) $get_bw50_resnos .= "ECHO \"$i: \" %$i.50\n";

$gpcr_contacts_a = [ "45x51-6x55", "45x53-6x55", "3x40-6x48", "5x47-6x48", "5x51-6x44", "5x55-6x41", "3x43-7x49", "3x43-7x53", "5x58-6x40", "5x58-7x53", "3x46-7x53", "1x53-7x54", "5x62-6x37", "3x50-7x53" ];
$gpcr_contacts_i = [ "3x43-6x40", "3x43-6x41", "6x40-7x49", "3x46-6x37", "1x53-7x53", "7x53-8x50", "7x54-8x51", "3x50-6x37" ];
$skip_compat_chk = [ "45.53-6.55" ];

$measure_contacts = "";
foreach ($gpcr_contacts_a as $contacts)
{
    $bwnos = explode('-', $contacts);
    if ($bwnos[0] == "ligand" || $bwnos[1] == "ligand") continue;

    $measure_contacts .= <<<heredoc
LET @distance = @A.{$bwnos[0]} - @A.{$bwnos[1]}
LET &active_distance = @distance

LET @distance = @I.{$bwnos[0]} - @I.{$bwnos[1]}
LET &inactive_distance = @distance

LET &change = &active_distance - &inactive_distance
ECHO "+{$bwnos[0]}|{$bwnos[1]}|" &change


heredoc;
}
foreach ($gpcr_contacts_i as $contacts)
{
    $bwnos = explode('-', $contacts);
    if ($bwnos[0] == "ligand" || $bwnos[1] == "ligand") continue;

    $measure_contacts .= <<<heredoc
LET @distance = @A.{$bwnos[0]} - @A.{$bwnos[1]}
LET &active_distance = @distance

LET @distance = @I.{$bwnos[0]} - @I.{$bwnos[1]}
LET &inactive_distance = @distance

LET &change = &active_distance - &inactive_distance
ECHO "-{$bwnos[0]}|{$bwnos[1]}|" &change


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

$bw50 = [];
$contacts_made_bw = [];
$contacts_broken_bw = [];

foreach ($result as $line)
{
    $pieces = explode(": ", $line);
    if (count($pieces) > 1)
    {
        $region_no = intval($pieces[0]);
        $resno = intval($pieces[1]);
        $bw50[$region_no] = $resno;
        continue;
    }

    $pieces = explode("|", $line);
    if (count($pieces) >= 3)
    {
        $c = substr($pieces[0], 0, 1);
        if ($c == '+') $contacts_made_bw[substr($pieces[0], 1)."-".$pieces[1]] = floatval($pieces[2]);
        else if ($c == '-') $contacts_broken_bw[substr($pieces[0], 1)."-".$pieces[1]] = floatval($pieces[2]);
    }
}

foreach ($contacts_made as $resnos => $change)
{
    $residues = explode("-", $resnos);
    $resno1 = intval(preg_replace("/[^0-9]/", "", $residues[0]));
    $resno2 = intval(preg_replace("/[^0-9]/", "", $residues[1]));

    $closest_bw50_helixno = $closest_bw50_resno_delta = 0;
    foreach ($bw50 as $regno => $resno50)
    {
        $delta = abs($resno1 - $resno50);
        if (!$closest_bw50_helixno || $delta < $closest_bw50_resno_delta)
        {
            $closest_bw50_resno_delta = $delta;
            $closest_bw50_helixno = $regno;
        }
    }

    $bw1 = $closest_bw50_helixno . '.' . ($resno1 - $bw50[$closest_bw50_helixno] + 50);

    $closest_bw50_helixno = $closest_bw50_resno_delta = 0;
    foreach ($bw50 as $regno => $resno50)
    {
        $delta = abs($resno2 - $resno50);
        if (!$closest_bw50_helixno || $delta < $closest_bw50_resno_delta)
        {
            $closest_bw50_resno_delta = $delta;
            $closest_bw50_helixno = $regno;
        }
    }

    $bw2 = $closest_bw50_helixno . '.' . ($resno2 - $bw50[$closest_bw50_helixno] + 50);

    $contacts_made_bw["$bw1-$bw2"] = floatval($change);
}

foreach ($contacts_broken as $resnos => $change)
{
    $residues = explode("-", $resnos);
    $resno1 = intval(preg_replace("/[^0-9]/", "", $residues[0]));
    $resno2 = intval(preg_replace("/[^0-9]/", "", $residues[1]));
    if ($resno1 == "ligand" || $resno2 == "ligand") continue;

    $closest_bw50_helixno = $closest_bw50_resno_delta = 0;
    foreach ($bw50 as $regno => $resno50)
    {
        $delta = abs($resno1 - $resno50);
        if (!$closest_bw50_helixno || $delta < $closest_bw50_resno_delta)
        {
            $closest_bw50_resno_delta = $delta;
            $closest_bw50_helixno = $regno;
        }
    }

    $bw1 = $closest_bw50_helixno . '.' . ($resno1 - $bw50[$closest_bw50_helixno] + 50);

    $closest_bw50_helixno = $closest_bw50_resno_delta = 0;
    foreach ($bw50 as $regno => $resno50)
    {
        $delta = abs($resno2 - $resno50);
        if (!$closest_bw50_helixno || $delta < $closest_bw50_resno_delta)
        {
            $closest_bw50_resno_delta = $delta;
            $closest_bw50_helixno = $regno;
        }
    }

    $bw2 = $closest_bw50_helixno . '.' . ($resno2 - $bw50[$closest_bw50_helixno] + 50);

    $contacts_broken_bw["$bw1-$bw2"] = floatval($change);
}

// Output results.
/*echo "Contacts made: ";
print_r($contacts_made_bw);
echo "\nContacts broken: ";
print_r($contacts_broken_bw);*/


// Check residues of $rcpid and list which made/broken contacts are compatible.
$working_pdbfn = filename_protid($rcpid);
$bw_lookup_residue = ["6.28", "6.59", "5.33", "5.68", "7.48", "4.57"];
foreach (array_merge($contacts_made_bw, $contacts_broken_bw) as $key => $value)
{
    $residues = explode('-', $key);

    for ($i=0; $i<=1; $i++)
    {
        // Pepteditor doesn't recognize GPCRDB numbering yet, so we have to assume Ballesteros-Weinstein.
        // For ORs, the two are almost always the same.
        $residue_bw = str_replace('x', '.', $residues[$i]);
        if ($residue_bw == "ligand") continue;
        $bw_lookup_residue[] = $residue_bw; // = "ECHO \"$residue_bw|\" \$R.$residue_bw %R.$residue_bw \"|\" @R.$residue_bw \"|\" \$A.$residue_bw";
    }
}

$bw_lookup = [];
foreach ($bw_lookup_residue as $residue_bw)
    $bw_lookup[$residue_bw] = "ECHO \"$residue_bw|\" \$R.$residue_bw %R.$residue_bw \"|\" @R.$residue_bw \"|\" \$A.$residue_bw";

$optimized = [];
foreach ($contacts_made_bw as $key => $value)
{
    $residues = explode('-', $key);
    if ($residues[0] == "ligand" || $residues[1] == "ligand") continue;
    $optimized[$key] = "MOC %R.{$residues[0]} %R.{$residues[1]} &optimized\nECHO \"{$residues[0]}\" \"-\" \"{$residues[1]}\" \"|\" &optimized";
}

$bw_lookup = implode("\n", $bw_lookup);
$optimized = implode("\n", $optimized);

$pepd_script = <<<heredoc

STRAND A
LOAD "$inactive_pdbfn"

STRAND R
LOAD "$working_pdbfn"

$bw_lookup

ECHO "-----"

$optimized

ECHO "-----"

CANMOVE %5.43 %5.54 TO %6.44 %6.55 &TMR5x
ECHO "TMR5x: " &TMR5x

heredoc;


$pepd_fname = "tmp/icactive.pepd";
$fp = fopen($pepd_fname, "w") or die("File access or permissions error.\n");
fwrite($fp, $pepd_script);
fclose($fp);

$result = [];
$cmd = "bin/pepteditor $pepd_fname";
exec($cmd, $result);

print_r($result);


$residue_info = [];
$contact_spacing = [];

$mode = 0;
foreach ($result as $line)
{
    if (trim($line) == "-----")
    {
        $mode++;
        continue;
    }

    switch ($mode)
    {
        case 0:
        $line = explode("|", $line);
        $xyz = explode(',',preg_replace("/[\\[\\]]/", '', $line[2]));
        $line[2] = ["x" => floatval($xyz[0]), "y" => floatval($xyz[1]), "z" => floatval($xyz[2])];
        $residue_info[$line[0]] = $line;
        break;

        case 1:
        $line = explode("|", $line);
        $line[0] = str_replace("x", ".", $line[0]);
        $contact_spacing[$line[0]] = floatval($line[1]);
        break;

        case 2:
        $line = explode(':', $line);
        $varname = $line[0];
        $$varname = is_numeric($line[1]) ? floatval($line[1]) : $line[1];
        break;

        default:
        break 2;
    }
}

ksort($residue_info);
// print_r($contact_spacing);

foreach ($contact_spacing as $contacts => $spacing)
{
    $residues = explode('-', $contacts);
    $delete = false;

    if (!isset($residue_info[$residues[0]])) $delete = true;
    if (!isset($residue_info[$residues[1]])) $delete = true;
    if (!@$residue_info[$residues[0]][3]) $delete = true;
    if (!@$residue_info[$residues[1]][3]) $delete = true;

    $tmr0 = intval(explode('.', explode('x', $residues[0])[0])[0]);
    $tmr1 = intval(explode('.', explode('x', $residues[1])[0])[0]);

    if ($tmr0 == $tmr1) $delete = true;

    if (!$delete)
    {
        $aa0 = substr($residue_info[$residues[0]][1], 0, 1);
        $aa1 = substr($residue_info[$residues[1]][1], 0, 1);

        if (/*($aa0 != $residue_info[$residues[0]][3] || $aa1 != $residue_info[$residues[1]][3]) &&*/ !isset($skip_compat_chk[$contacts]))
        {
            if (!residues_compatible($aa0, $aa1))
            {
                echo "Delete incompatible $aa0{$residues[0]}-$aa1{$residues[1]}\n";
                $delete = true;
            }
        }
    }

    if ($delete) unset($contact_spacing[$contacts]);
}

// print_r($residue_info);
// print_r($contact_spacing);


// Types of TMR6 activation motion:
$tmr6type = "?";
$hasR6x59 = false;
$exrbend = 0;
$dynamics = [];
$tmr6_distance_h = 7.0;

foreach ($contact_spacing as $key => $value)
{
    if ($key == "45.53-6.55" && isset($contact_spacing["45.51-6.55"])) continue;
    
    $residues = explode('-', $key);

    $aa0 = substr($residue_info[$residues[0]][1], 0, 1);
    $aa1 = substr($residue_info[$residues[1]][1], 0, 1);

    if (residues_compatible($aa0, $aa1))
        $dynamics[] = "BRIDGE ".str_replace('-', ' ', $key);
}

// This is important for OR51E2 ligand specificity.
// Motion should be prevented by V5.46, however the docker seems to not be avoiding internal clashes.
/*$aa4x57 = substr($residue_info["4.57"][1], 0, 1);
if ($aa4x57 == 'F' || $aa4x57 == 'W' || $aa4x57 == 'Y')
{
    $dynamics[] = "STCR 4.57";
}*/

// 6.48 Rock: If there is no 45.51-6.55 contact, and there is room for the EXR end of TMR6 to move toward the EXR2 helix, then TMR6 pivots at 6.48,
// creating a rift in the cytoplasmic end and closing around the ligand. If R6.59 is present, it partially uncoils and rotates its side chain
// inward to make contact with the ligand. Else if 45.53 is compatible with 6.55, then that pair becomes the contact.
// Examples of 6.48 rock receptors: OR51E2 (with R6.59), OR1G1 (N45.53).
if (!isset($contact_spacing["45.51-6.55"])
    &&
    (   floatval(@$contact_spacing["45.53-6.55"])
        || substr($residue_info["6.59"][1], 0, 1) == 'R'
    )
)
{
    $tmr6type = "6.48 Rock";
    $m6name = "rock6";
    if (substr($residue_info["6.59"][1], 0, 1) == 'R')
    {
        $hasR6x59 = true;
        $dynamics[] = "ATOMTO 6.59 CZ 4.60";
        $dynamics[] = "STCR 6.59";
        // $dynamics[] = "DYNAMIC WIND wind6 6.56 6.59 -8";
    }

    $distance_v = get_distance($residue_info["6.48"][2], $residue_info["6.28"][2]);
    $distance_h = $tmr6_distance_h;
    $angle = round(45.0 * $distance_h / $distance_v, 2);            // Approximation.
    $dynamics[] = "DYNAMIC BEND $m6name 56.50 6.59 6.48 7.43 -$angle";
}

// 6.48 Bend: If there is a 45.51-6.55 contact and a 3.40-6.48 contact, then TMR6 bends at the 6.48 position to create a cytoplasmic rift, but does
// not move the extracellular end unless necessary to complete the 45.51-6.55 contact.
// Examples of 6.48 bend receptors: OR1A1 (no EXR bend), OR5AN1 (EXR bend).
else if (isset($contact_spacing["45.51-6.55"]) && isset($contact_spacing["3.40-6.48"]))
{
    $tmr6type = "6.48 Bend";
    $m6name = "bend6";
    $exrbend = $contact_spacing["45.51-6.55"];

    $distance_v = get_distance($residue_info["6.48"][2], $residue_info["6.28"][2]);
    $distance_h = $tmr6_distance_h;
    $angle = round(45.0 * $distance_h / $distance_v, 2);
    $dynamics[] = "DYNAMIC BEND $m6name 56.50 6.48 6.48 7.43 -$angle";
}

// 6.55 Bend: If there is a 45.51-6.55 contact but no strong 3.40-6.48 contact, then TMR6 bends at the 6.55 position and most of its length moves
// as a unit to create the cytoplasmic rift.
// Examples of possible 6.55 bend receptors: OR8D1, OR14A2.
else if (isset($contact_spacing["45.51-6.55"]) && !isset($contact_spacing["3.40-6.48"]))
{
    $tmr6type = "6.55 Bend";
    $m6name = "bend6";

    $distance_v = get_distance($residue_info["6.55"][2], $residue_info["6.28"][2]);
    $distance_h = $tmr6_distance_h;
    $angle = round(45.0 * $distance_h / $distance_v, 2);
    $dynamics[] = "DYNAMIC BEND $m6name 56.50 6.55 6.55 2.64 -$angle";
}

// The motion of TMR6 shall be sufficient to move the side chain of 6.40 out of the way for 5.58 and 7.53 to make contact. The side chain of 6.40
// moves to point to 7.52. In OR51E2, the CA of 6.27 moves about 7A away from the CYT2 loop.
echo "TMR6 type: $tmr6type" . ($hasR6x59 ? " with R6.59" : "") . ($exrbend ? " with EXR bend $exrbend" : "") . ".\n";


// TMR5 activation motion:
// TMR5 generally bends at the extracellular end, to keep up with the motion of the cytoplasmic end of TMR6. TMR5 also moves slightly toward TMR6,
// about 1A at the extracellular end and about 3.5A at the cytoplasmic end, the difference being due to the bend.
$distance_v = get_distance($residue_info["5.33"][2], $residue_info["5.68"][2]);

// TODO: The horizontal motion of TMR5+TMR6 is driven by the 5.58-7.53 contact.
// 3.43 probably limits the inward bend of TMR7.
// The overall horizontal motion of TMR5 is limited by multiple residues between 5.43-5.54 and 6.44-6.55.
// The remaining portion of the gap is closed by rotating all of TMR5 to bring 5.58 near 7.53 and rotating
// TMR6 about the Y axis so its CYT end aligns with TMR5's CYT end.
// These motions are important for e.g. OR51E2 in which TMR5 closes onto the ligand binding site, limiting the pocket volume from larger SCFAs.
$distance5h = $TMR5x;
$distance_h = $tmr6_distance_h - $distance5h;
$angle = round(45.0 * $distance_h / $distance_v, 2);
$dynamics[] = "DYNAMIC BEND bend5 5.33 56.49 5.33 2.64 -$angle SYNC $m6name";
$dynamics[] = "DYNAMIC MOVE move5 5.33 56.49 5.43 6.55 $distance5h SYNC bend5";


// TMR7 activation motion:
// The near universal 5.58-7.53 contact is improved by a bend at 7.48. In receptors that lack this contact, the TMR7 bend likely would not occur.
if (@$contact_spacing["5.58-7.53"])
{
    $distance_v = get_distance($residue_info["7.48"][2], $residue_info["7.53"][2]);
    $distance_h = floatval($contact_spacing["5.58-7.53"]);
    $angle = round(45.0 * $distance_h / $distance_v, 2);
    $dynamics[] = "DYNAMIC BEND bend7 7.48 7.56 7.48 2.50 $angle MAX $m6name";
    $dynamics[] = "DYNAMIC WIND wind7 7.51 7.56 -10 SYNC bend7";
}


// EXR2 activation motion:
// If R6.59 exists, the inward motion of its side chain causes 45.53 to point "downward" in the cytoplasmic direction, with a bend of the helix
// in the range of the 45.52-45.54 positions.
if ($hasR6x59)
{
    $dynamics[] = "DYNAMIC BEND bend45 45.53 45.54 45.53 45.54 30 MIN $m6name";
}


echo "\nPRIMARYDOCK_CONFIG_PARAMS:\n".implode("\n", $dynamics)."\n\n";
