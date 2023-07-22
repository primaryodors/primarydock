<?php

chdir(__DIR__);
chdir("..");

// Includes
require("cputemp.php");
require("data/protutils.php");

foreach (@$argv as $a)
{
    $a = explode('=',$a,2);
    $_REQUEST[$a[0]] = (count($a)>1) ? $a[1] : true;
}

$max_simultaneous_couples = 4;
$cpl_dir = "pdbs/coupled";

if (@$_REQUEST['simul']) $max_simultaneous_couples = intval($_REQUEST['simul']);

if (!file_exists($cpl_dir)) mkdir($cpl_dir);

if (@$_REQUEST['next'])
{
    $cmd = "ps -ef | grep ':[0-9][0-9] bin/pepteditor' | grep -v grep";
    exec($cmd, $results);
    if (!@$_REQUEST['force'] && trim(@$results[$max_simultaneous_couples-1])) die("Already running.\n".print_r($results, 1));
    $already = implode("\n", $results);

    $gpcrid = @$_REQUEST['rcp'] ?: false;
    $gpid = @$_REQUEST['gprot'] ?: false;

    foreach (array_keys($prots) as $rcpid)
    {
        if (false!==strpos($already, "/$rcpid"))
        {
            continue;
        }

        if ($gpcrid && $gpcrid != $rcpid && !preg_match("/^$gpcrid$/", $rcpid) ) continue;
        foreach ($gprots as $gprotid => $gpdata)
        {
            if ($gpid && $gpid != $gprotid) continue;

            $fam = family_from_protid($rcpid);
            if (!file_exists("$cpl_dir/$fam")) mkdir("$cpl_dir/$fam");

            $fname = "$cpl_dir/$fam/{$rcpid}_{$gprotid}.pdb";
            if (!file_exists($fname))
            {
                $gpcrid = $rcpid;
                $gpid = $gprotid;
                goto _found_next;
            }
        }
    }
	  
    die("All done!");
}
else
{
    $gpcrid = @$_REQUEST['rcp'] ?: "OR1A1";
    $gpid = @$_REQUEST['gprot'] ?: "hGNAL";
}
	
_found_next:
;

$fam = family_from_protid($gpcrid);
$orfn = "pdbs/$fam/{$gpcrid}.upright.pdb";
$fname = "$cpl_dir/$fam/{$gpcrid}_{$gpid}.pdb";

if (!file_exists("$cpl_dir/$fam")) mkdir("$cpl_dir/$fam");

die_if_too_hot();

$cmd = "bin/pepteditor predict/couple.pepd $orfn $gpid | tee tmp/cplout";
echo "Running: $cmd\n";
passthru($cmd);

$fp = fopen($fname, "wb");
if ($fp)
{
    fwrite($fp, file_get_contents("tmp/coupled.{$gpcrid}.pdb"));
    fclose($fp);
    unlink("tmp/coupled.{$gpcrid}.pdb");
}








