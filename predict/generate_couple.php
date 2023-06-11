<?php

chdir(__DIR__);
chdir("..");

// Includes
require("predict/protutils.php");

$Gprots = ["hGNAL", "hGNAS2"];

foreach (@$argv as $a)
{
	$a = explode('=',$a,2);
	$_REQUEST[$a[0]] = (count($a)>1) ? $a[1] : true;
}

$max_simultaneous_couples = 4;
$cpl_dir = "pdbs/coupled";

if (!file_exists($cpl_dir)) mkdir($cpl_dir);

if (@$_REQUEST['next'])
{
	$cmd = "ps -ef | grep ':[0-9][0-9] bin/couple' | grep -v grep";
	exec($cmd, $results);
	if (!@$_REQUEST['force'] && trim(@$results[$max_simultaneous_couples-1])) die("Already running.\n".print_r($results, 1));
	$already = implode("\n", $results);

	$gpcrid = @$_REQUEST['rcp'] ?: false;
	$gpid = @$_REQUEST['gprot'] ?: false;
	
	foreach (array_keys($prots) as $rcpid)
	{
	    if (false!==strpos($already, "/$rcpid."))
		{
			continue;
		}

		if ($gpcrid && $gpcrid != $rcpid && !preg_match("/^$gpcrid$/", $rcpid) ) continue;
		foreach ($Gprots as $gprotid)
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
$cfgf = "tmp/{$gpcrid}_{$gpid}.cplcfg";
$fname = "$cpl_dir/$fam/{$gpcrid}_{$gpid}.pdb";

if (!file_exists("$cpl_dir/$fam")) mkdir("$cpl_dir/$fam");

$cfg = <<<heredoc

PROT1 pdbs/$fam/$gpcrid.upright.pdb
PROT2 pdbs/Gprot/$gpid.pdb
# TEMPLATE pdbs/OR51/OR51E2.8f76.pdb

CONTACT 1:RK6.58~1 1:DE45.51
CONTACT N4.38~2 DEQRK(LLLLGAGESGKSTIVKQMRILH-8)
CONTACT RK3.59~1 DE(GIXETXF+9)~1
CONTACT KR7.57~2 DE(HFTCATDTXNXXFVF+30)~3
CONTACT DE6.29 RK(HFTCATDTXNXXFVF+23)~2
CONTACT H3.59~1 Y(HFTCATDTXNXXFVF+29)~2
CONTACT STNQ56.50 Y(HFTCATDTXNXXFVF-4)~2
CONTACT STNQ7.59 STNQ(LKQYELL+2)
CONTACT MAILV5.65 MAILV(LKQYELL+0)
CONTACT MAILV5.61 MAILV(LKQYELL+5)
CONTACT MAILV5.68 MAILV(LKQYELL+6)
CONTACT MAILV3.54 MAILV(LKQYELL+6)

SEGMENT 1 1.27 0 1.58
SEGMENT 1.59 2.37 1.28 2.65
SEGMENT 2.66 3.20 2.38 3.56
SEGMENT 3.57 4.38 3.21 4.64
SEGMENT 4.65 5.31 4.39 5.67
SEGMENT 5.68 6.27 5.32 6.59 5.32 6.48
SEGMENT 6.60 7.29 6.28 7.56
SEGMENT 7.57 end 7.30

ITER 50

OUT $fname



heredoc;

$fp = fopen($cfgf, "wb") or die("FAILED to open $cfgf for writing.");
fwrite($fp, $cfg);
fclose($fp);

$cmd = "bin/couple $cfgf";
echo "Running: $cmd\n";
passthru($cmd);
unlink($cfgf);







