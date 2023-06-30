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

$iters = 50;
if ($gpid == "hGNAO1") $iters = 60;

$fam = family_from_protid($gpcrid);
$cfgf = "tmp/{$gpcrid}_{$gpid}.cplcfg";
$fname = "$cpl_dir/$fam/{$gpcrid}_{$gpid}.pdb";

if (!file_exists("$cpl_dir/$fam")) mkdir("$cpl_dir/$fam");

/*
$cfg = <<<heredoc

PROT1 pdbs/$fam/$gpcrid.upright.pdb
PROT2 pdbs/Gprot/$gpid.pdb
TEMPLATE pdbs/OR51/OR51E2.8f76.pdb pdbs/OR51/OR51E2.upright.pdb

CONTACT 1:YHNQ6.55 1:DE45.51                                    # EXR hbond for tetrapod ORs
CONTACT STYNQ4.38~2 ERKSTYNQD(LLLLGAGESGKSTIVKQMRILH-8)~2       # CYT2 and N-terminus "index finger" helix
CONTACT KR7.57~2 DE(HFTCATDTXNXXFVF+30)~3                       # TMR7/tail joint and tip of C-terminus "thumb" helix
CONTACT DE7.59~2 KR(NKEIYCHXTCATDTXNXXXVF+33)~2                 # TMR7/tail joint and tip of C-terminus "thumb" helix
CONTACT DE6.29 KR(HFTCATDTXNXXFVF+24)~4                         # TMR6 and middle of "thumb" helix

CONTACT NQRHFWY3.59~1 FWY(HFTCATDTXNXXFVF+29)~3                 # TMR3 and "thumb" tip
CONTACT NQRHFWY3.59~1 FWY(DIIIANNLR+10)~3                       # TMR3 and "thumb" tip

CONTACT STNQED7.59 STNQED(LKQYELL+2)                            # TMR7/tail joint and tip of C-terminus "thumb" helix
CONTACT STNQED7.59 STNQED(NLKXCGLF+3)                           # TMR7/tail joint and tip of C-terminus "thumb" helix

CONTACT MAILV5.65 MAILV(LKQYELL+0)                              # TMR2/3/5 and "thumb" tip hydrophobic pocket
CONTACT MAILV5.61 MAILV(LKQYELL+5)                              # TMR2/3/5 and "thumb" tip hydrophobic pocket
CONTACT MAILV5.68 MAILV(LKQYELL+6)                              # TMR2/3/5 and "thumb" tip hydrophobic pocket
CONTACT MAILV3.54 MAILV(LKQYELL+6)                              # TMR2/3/5 and "thumb" tip hydrophobic pocket
CONTACT MAILV2.39 MAILV(LKQYELL+5)                              # TMR2/3/5 and "thumb" tip hydrophobic pocket
CONTACT MAILV3.46 MAILV(LKQYELL+5)                              # TMR2/3/5 and "thumb" tip hydrophobic pocket

CONTACT MAILV3.54 MAILV(DAVTDVIIKNNLK+2)                        # TMR2/3/5 and "thumb" tip hydrophobic pocket
CONTACT MAILV5.65 MAILV(DAVTDVIIKNNLK+5)                        # TMR2/3/5 and "thumb" tip hydrophobic pocket
CONTACT MAILV5.65 MAILV(DAVTDVIIKNNLK+6)                        # TMR2/3/5 and "thumb" tip hydrophobic pocket
CONTACT MAILV5.68 MAILV(DAVTDVIIKNNLK+7)                        # TMR2/3/5 and "thumb" tip hydrophobic pocket
CONTACT MAILV2.39 MAILV(DAVTDVIIKNNLK+11)                       # TMR2/3/5 and "thumb" tip hydrophobic pocket

CONTACT MAILV2.42 MAILV(LKQYELL+5)                              # "thumb" tip deeper into receptor
CONTACT MAILV2.43 MAILV(LKQYELL+5)                              # "thumb" tip deeper into receptor

CONTACT X3.58 X(LLLLGAGESGKSTIVKQMRILH-8)                       # Dummy contact to align "index finger" with CYT2.
CONTACT X3.62 X(LLLLGAGESGKSTIVKQMRILH-4)                       # Dummy contact to align "index finger" with CYT2.

CONTACT ERKSTYNQD56.50 Y(YCYPHFTCAVDTENIRR+0)		            # CYT3 region to strand beside "thumb" (olfactory G proteins)
CONTACT ERKSTYNQD56.50 Y(YSHFTCATDTQNIQFVF+0)		            # CYT3 region to strand beside "thumb" (non-olfactory G proteins)
CONTACT MAILV3.62~2 FHMAILVP(HFTCAXDTXNIXFVF+14)~1            	# CYT2 region to middle of "thumb"

SEGMENT 1 1.27 0 1.58
SEGMENT 1.59 2.37 1.28 2.65
# SEGMENT 2.66 3.20 2.38 3.56
SEGMENT 3.57 4.38 3.21 4.64
# SEGMENT 4.65 5.31 4.39 5.67
# SEGMENT 4.65 4.66 4.39 45.51
# SEGMENT 5.30 5.31 45.58 5.67
# SEGMENT 5.68 6.27 5.32 6.59 5.32 6.48
SEGMENT 5.68 6.27 5.32 6.59 5.32 6.59
SEGMENT 6.60 7.29 6.28 7.56
SEGMENT 7.57 end 7.30

MAKESURE Y6.55 DE45.51 6.27 6.59 6.27

ITER $iters

OUT $fname



heredoc;
*/

$cfg = <<<heredoc

PROT1 pdbs/$fam/$gpcrid.upright.pdb
PROT2 pdbs/Gprot/$gpid.pdb
TEMPLATE pdbs/OR51/OR51E2.8f76.pdb pdbs/OR51/OR51E2.upright.pdb

CONTACT X7.57 STYNQRKH(XXXXXX$+3)~2
CONTACT X7.59 X(XXXXXX$+0)
CONTACT MAILVCP3.46 MAILVCPFWY(XXXX$+2)~1
CONTACT X3.60 X(KLLLLGAGESGKST-4)
CONTACT X3.54 MAILVCP(XXXXXXX$+0)~1
CONTACT X5.66 MAILVCP(XXXXXXX$+0)~1
CONTACT X56.50 H(YXHFTCATDTXNI+2)~1
CONTACT X3.58 X(LLLLGAGESGKST-8)                                # Dummy contact to align "index finger" with CYT2.
CONTACT X3.62 X(LLLLGAGESGKST-4)                                # Dummy contact to align "index finger" with CYT2.
CONTACT X3.58 X(LLLLGAGESGKST-3)                                # Dummy contact to align "index finger" with CYT2.
CONTACT X3.62 X(LLLLGAGESGKST-3)                                # Dummy contact to align "index finger" with CYT2.
CONTACT 1:YHNQ6.55 1:DE45.51                                    # EXR hbond for tetrapod ORs

SEGMENT 1 1.27 0 1.58
SEGMENT 1.59 2.37 1.28 2.65
# SEGMENT 2.66 3.20 2.38 3.56
SEGMENT 3.57 4.38 3.21 4.64
# SEGMENT 4.65 5.31 4.39 5.67
# SEGMENT 4.65 4.66 4.39 45.51
# SEGMENT 5.30 5.31 45.58 5.67
# SEGMENT 5.68 6.27 5.32 6.59 5.32 6.48
SEGMENT 5.68 6.27 5.32 6.59 5.32 6.59
SEGMENT 6.60 7.29 6.28 7.56
SEGMENT 7.57 end 7.30

MAKESURE Y6.55 DE45.51 6.27 6.59 6.27

ITER $iters

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







