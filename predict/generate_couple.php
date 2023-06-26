<?php

chdir(__DIR__);
chdir("..");

// Includes
require("cputemp.php");
require("predict/protutils.php");

// $Gprots = ["hGNAL", "hGNAS2"];
$Gprots = [];

$c = file_get_contents("data/gprots_aligned.txt");
foreach (explode("\n", $c) as $ln)
{
    if (substr($ln, 0, 1) == '#') continue;
    $gpid = trim(substr($ln, 0, 10));
    if ($gpid)
    {
        if (file_exists("pdbs/Gprot/$gpid.pdb")) $Gprots[] = $gpid;
    }
}

// print_r($Gprots);

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
    $cmd = "ps -ef | grep ':[0-9][0-9] bin/couple' | grep -v grep";
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
TEMPLATE pdbs/OR51/OR51E2.8f76.pdb pdbs/OR51/OR51E2.upright.pdb

# CONTACT 1:RK6.58~1 1:DE45.51                        # EXR salt bridge for fish-like ORs ***CAUSES CLASHES***
CONTACT 1:YHNQ6.55 1:DE45.51                        # EXR hbond for tetrapod ORs
CONTACT N4.38~2 DEQRK(LLLLGAGESGKSTIVKQMRILH-8)     # CYT2 and N-terminus "index finger" helix
CONTACT KR7.57~2 DE(HFTCATDTXNXXFVF+30)~3           # TMR7/tail joint and tip of C-terminus "thumb" helix
CONTACT DE6.29 RK(HFTCATDTXNXXFVF+23)~2             # TMR6 and middle of "thumb" helix
CONTACT H3.59~1 Y(HFTCATDTXNXXFVF+29)~2             # TMR3 and "thumb" tip
# CONTACT STNQED56.50 K(IIQRMHLKQYELL+7)              # CYT3 and "thumb" tip (olfactory G proteins)
# CONTACT STNQED56.50 K(DAVTDIIIKENLKDCGLF+8)         # CYT3 and "thumb" tip (non-olfactory G proteins)
CONTACT STNQED7.59 STNQ(LKQYELL+2)                  # TMR7/tail joint and tip of C-terminus "thumb" helix
CONTACT MAILV5.65 MAILV(LKQYELL+0)                  # TMR2/3/5 and "thumb" tip hydrophobic pocket
CONTACT MAILV5.61 MAILV(LKQYELL+5)                  # TMR2/3/5 and "thumb" tip hydrophobic pocket
CONTACT MAILV5.68 MAILV(LKQYELL+6)                  # TMR2/3/5 and "thumb" tip hydrophobic pocket
CONTACT MAILV3.54 MAILV(LKQYELL+6)                  # TMR2/3/5 and "thumb" tip hydrophobic pocket
CONTACT MAILV2.39 MAILV(LKQYELL+5)                  # TMR2/3/5 and "thumb" tip hydrophobic pocket
CONTACT MAILV3.46 MAILV(LKQYELL+5)                  # TMR2/3/5 and "thumb" tip hydrophobic pocket
CONTACT MAILV2.42 MAILV(LKQYELL+5)                  # "thumb" tip deeper into receptor
CONTACT MAILV2.43 MAILV(LKQYELL+5)                  # "thumb" tip deeper into receptor
CONTACT MAILV2.46 MAILV(LKQYELL+5)                  # "thumb" tip deeper into receptor
CONTACT MAILV2.49 MAILV(LKQYELL+5)                  # "thumb" tip deeper into receptor
CONTACT X3.58 X(LLLLGAGESGKSTIVKQMRILH-8)           # Dummy contact to align "index finger" with CYT2.
CONTACT X3.62 X(LLLLGAGESGKSTIVKQMRILH-4)           # Dummy contact to align "index finger" with CYT2.
CONTACT ERKSTYNQD56.50 Y(YCYPHFTCAVDTENIRR+0)		# CYT3 region to strand beside "thumb" (olfactory G proteins)
CONTACT ERKSTYNQD56.50 Y(YSHFTCATDTQNIQFVF+0)		# CYT3 region to strand beside "thumb" (non-olfactory G proteins)
CONTACT MAILV3.62~2 FHMAILVP(HFTCAXDTXNIXFVF+14)~1	# CYT2 region to middle of "thumb"

SEGMENT 1 1.27 0 1.58
SEGMENT 1.59 2.37 1.28 2.65
# SEGMENT 2.66 3.20 2.38 3.56
SEGMENT 3.57 4.38 3.21 4.64
# SEGMENT 4.65 5.31 4.39 5.67
SEGMENT 4.65 4.66 4.39 45.51
SEGMENT 5.30 5.31 45.58 5.67
SEGMENT 5.68 6.27 5.32 6.59 5.32 6.48
SEGMENT 6.60 7.29 6.28 7.56
SEGMENT 7.57 end 7.30

ITER 50

OUT $fname



heredoc;

$fp = fopen($cfgf, "wb") or die("FAILED to open $cfgf for writing.");
fwrite($fp, $cfg);
fclose($fp);

die_if_too_hot();

$cmd = "bin/couple $cfgf | tee tmp/cplout";
echo "Running: $cmd\n";
passthru($cmd);
unlink($cfgf);

if (!file_exists($fname))
{
    $fp = fopen($fname, "wb");
    if ($fp)
    {
        fwrite($fp, file_get_contents("tmp/cplout"));
        fclose($fp);
    }
}







