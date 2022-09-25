<?php

// sequence_update.php
//
// Reads the sequences_aligned.txt file and updates all the PDBs with region info and binding site info.
//

chdir(__DIR__);
$c = file_get_contents('sequences_aligned.txt', 'rb');
if (!$c) die ('NO INPUT FILE');

$lines = explode("\n", $c);
$ells = array();
$sites = array();
$fam = isset($_REQUEST['family']) ? $_REQUEST['family'] : false;
$dosql = isset($_REQUEST['dosql']) ? intval($_REQUEST['dosql']) : false;
$ballwein = [];     // https://www.sciencedirect.com/science/article/abs/pii/S1043947105800497

$lln1 = '';
foreach ($lines as $i => $ln)
{	
	$ln = strtoupper(str_replace("\t", '    ', $ln));
	
	if (substr($ln,0,1) == '#') continue;
	
	if (substr($ln,0,1) == '%') continue;               // TODO: Sequence alignment check.
	
	if (preg_match("/^[\sLSM]+$/", $ln))
	{
        $ells = array(); $j = -1;
		$jarr = [];
		do
		{
			preg_match('/[LSM]/', $ln, $jarr, PREG_OFFSET_CAPTURE, $j+1);
			foreach ($jarr as $je)
			{	
                $j = $je[1];
				$ells[$j] = substr($ln, $j, 1); // $j;
			}
		} while (count($jarr));
		//print_r($ells);
		
		continue;
	}
	
	if (preg_match('/TMR[0-9][-]+/', $ln))
	{   
        $ln .= '     ';
	    $rgnames = ['HEAD', 'TMR1', 'CYT1', 'TMR2', 'EXR1', 'TMR3', 'CYT2', 'TMR4', 'EXR2', 'TMR5', 'CYT3', 'TMR6', 'EXR3', 'TMR7', 'TAIL'];
	    $rgidx = 0; $k = 0;
		$rgns = [];
	    do
	    {	
            $lk = $k;
		    $j  = strpos($ln, 'TMR', $k);
		    if (!$j)
		    {	
                $j = strlen($ln)+100;
			    $rgns[$rgidx]  = [$lk, $j-1, $rgnames[$rgidx]];
			    break;
		    }
		    $k = strpos($ln, ' ', $j);
		    
		    $l = strpos($ln, '|', $lk);
		    if ($l > 0 && $l < $j) $ballwein[$rgnames[$rgidx]] = $l;
		    
		    $l = strpos($ln, '|', $j);
		    if ($l > 0) $ballwein[$rgnames[$rgidx+1]] = $l;
		    
		    $rgns[$rgidx]  = [$lk, $j-1, $rgnames[$rgidx]];
		    $rgidx++;
		    
		    $rgns[$rgidx]  = [$j, $k-1, $rgnames[$rgidx]];
		    $rgidx++;
	    } while ($k);
	
	    // print_r($rgns); exit;
	    // print_r($ballwein); exit;
		
		continue;
	}
	
	
	$match = "/^(OR[0-9]{1,2}[A-Z]{1,2}[0-9]{1,2}\s+[.AC-IK-NP-TVWY]+)|(TAAR[0-9]\s+[.AC-IK-NP-TVWY]+)|(VN1R[0-9]\s+[.AC-IK-NP-TVWY]+)/";
	
	if (preg_match($match, $ln))
    {
		$regionse = [];
		
		$r = 0; $pos1 = $posL = $posM = $posS = 0; $emmed = $spaced = $lastrgn = false; $lsites = array();
		$or = trim(substr($ln, 0, 9));
		if ($j = strpos($or, ' ')) $or = substr($or, 0, $j);
		$or = trim($or);

        $rcpbw = [];
        $rcpbs = [];
		
		foreach(str_split($ln) as $pos => $chr)
		{	
			if ($chr == ' ') $spaced = true;
			if ($spaced && ($chr == 'M' || $chr == '.')) $emmed = true;
			if ($emmed)
			{	if ($chr >= 'A' && $chr <= 'Y') $r++;
				$lsites[] = $chr.$r;
				$strgn = '????';
				foreach ($rgns as $rgi => $rg)
				{	if ($pos >= $rg[0] && $pos <= $rg[1]) 
					{	$strgn = $rg[2];
						$regionse[$rgi][0] = $strgn;
						if (!isset($regionse[$rgi][1])) $regionse[$rgi][1] = $r; 
						$regionse[$rgi][2] = $r;							
						break;
					}
				}
				
				if ($strgn != $lastrgn) $posL = $posM = $posS = 0;
				
				if (substr($strgn,0,3) == 'TMR' && $pos == @$ballwein[$strgn])
				{
					$tmrno = intval(substr($strgn, 3));
                    $rcpbw[$tmrno] = $r;
				}
				
				if (substr($strgn,0,3) == 'EXR' && $pos == @$ballwein[$strgn])
				{
					$exrno = intval(substr($strgn, 3));
					$a = $exrno * 2;
					$b = $a + 1;
                    $rcpbw[$a.$b] = $r;
				}
				
				if (isset($ells[$pos]) && $ells[$pos] != '') 
				{				
                    $rcpbs[] = $r;
				}

            }
        }

        // TODO: Read the PDB and then rewrite it with updated contents.
        $rcpid = explode(" ", trim(substr($ln, 0, 11)))[0];
        $subdir = substr($rcpid, 0, 4);
        if (substr($rcpid, 0, 2) == "OR") $subdir = "OR".intval(substr($rcpid,2,2));

        chdir(__DIR__);
        $pdbname = "../pdbs/$subdir/$rcpid.upright.pdb";
        if (file_exists("../pdbs/$subdir/$rcpid.rotated.pdb"))
            rename("../pdbs/$subdir/$rcpid.rotated.pdb", $pdbname);

        if (!file_exists($pdbname)) continue;
        
        $pdbdat = explode("\n",file_get_contents($pdbname));
        $remarks = [];
        $not_remarks = [];
        foreach ($pdbdat as $pdbln)
        {
            if (substr($pdbln, 0, 6) == "REMARK")
            {
                $remno = intval(substr($pdbln, 7));
                if ($remno < 225) $remarks[] = $pdbln;
            }
            else $not_remarks[] = $pdbln;
        }

        $remarks[] = "REMARK 650";
        foreach ($regionse as $se)
        {
            $region = $se[0];
            if (substr($region, 0, 3) == "TMR")
                $remarks[] = "REMARK 650 HELIX $region {$se[1]} {$se[2]}";
        }

        $remarks[] = "REMARK 800";
        foreach ($rcpbw as $region => $bw50)
        {
            $remarks[] = "REMARK 800 SITE BW $region.50 $bw50";
        }

        $remarks[] = "REMARK 800";
        foreach ($rcpbs as $bsr)
        {
            $remarks[] = "REMARK 800 SITE LIGAND_BINDING $bsr";
        }

        $pdbdat = implode("\n", $remarks) . "\n" . implode("\n", $not_remarks) . "\n\n";

        $f = fopen($pdbname, "wb");
        if (!$f) die("FAILED to open $pdbname for writing; ensure have access.\n");
        fwrite($f, $pdbdat);
        fclose($f);

        // exit;       // debug step.
    }
}