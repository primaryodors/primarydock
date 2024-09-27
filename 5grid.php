<?php

chdir(__DIR__);
require_once("data/protutils.php");
chdir(__DIR__);
require_once("data/odorutils.php");
chdir(__DIR__);

function do_5wind($pdbname, $amount)
{
  $cmd = "bin/pepteditor hm/fivewinder.pepd \"$pdbname\" \"$amount\"";
  echo "$cmd\n";
  passthru($cmd);
}

$outjson = "hm/5grid.json";
$results = [];
if (file_exists($outjson)) $results = json_decode(file_get_contents($outjson), true);

foreach ($prots as $rcpid => $p)
{
  $type = substr($rcpid, 0, 4);
  if ($type == "OR51" || $type == "OR52" || $type == "OR56") continue;      // class I
  if ($type == "TAAR") continue;
  if ($type == "VN1R") continue;
  if ($type == "MS4A") continue;

  // Find strongest ligand. If orphan receptor or only weak ligands, skip.
  $pairs = all_empirical_pairs_for_receptor($rcpid, false, true);
  // print_r($pairs); exit;
  $oid = array_keys($pairs)[0];
  if (isset($pairs[$oid]['adjusted_curve_top']) && floatval($pairs[$oid]['adjusted_curve_top']) < 4) continue;
  else if (isset($pairs[$oid]['ec50']) && floatval($pairs[$oid]['ec50']) > -4) continue;
  $lig = $odors[$oid]['full_name'];

  $initw5 =  -5;
  $w5step = 0.5;
  $endw5  =   5;

  if (isset($results[$rcpid])) $initw5 = 0.1*max(array_keys($results[$rcpid])) + $w5step;
  if ($initw5 > $endw5) continue;

  // Copy the active model to a temporary file.
  $fam = family_from_protid($rcpid);
  $acvfn = "pdbs/$fam/$rcpid.active.pdb";
  $bkpfn = "pdbs/$fam/$rcpid.active.bak.pdb";
  copy($acvfn, $bkpfn);

  do_5wind($acvfn, $initw5);
  for ($w5 = $initw5; $w5 <= $endw5; $w5 += $w5step)
  {
    $cmd = "./run_prediction.sh $rcpid \"$lig\"";
    passthru($cmd);

    // Get dock score from JSON.
    $data = json_decode(file_get_contents("predict/dock_results.json"), true);
    if (!isset($data[$rcpid][$lig]['DockScore'])) $dockscore = "ERROR";
    else $dockscore = floatval($data[$rcpid][$lig]['DockScore']);

    $results[$rcpid][10*$w5] = $dockscore;
    $f = fopen($outjson, "wb");
    if (!$f) die("File write FAILED. Make sure have access to write $outjson.");
    fwrite($f, json_encode_pretty($results));
    fclose($f);    

    unlink($acvfn);
    copy($bkpfn, $acvfn);
    exec("git commit -am \"Update.\"");
    exec("git push");
    do_5wind($acvfn, $w5);
  }

  unlink($acvfn);
  rename($bkpfn, $acvfn);
}













