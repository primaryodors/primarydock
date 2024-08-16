<?php

chdir(__DIR__);
require_once("../data/protutils.php");
chdir(__DIR__);

if (!file_exists('/usr/lib/libmodeller.so')) die("This feature requires MODELLER.\nPlease see: https://salilab.org/modeller/ to obtain this third party software.\n\n");

$rcpid = @$argv[1];
$p = $prots[$rcpid];
if (!$rcpid) die("Usage:\nphp -f hm/dohm.php PROTID\n\n");

$helices = "";
foreach ($p["region"] as $rgname => $rgnse)
{
    $nmsub3 = substr($rgname, 0, 3);
    if ($nmsub3 != "TMR" && $nmsub3 != "HXR") continue;
    $rgs = $rgnse['start'];
    $rge = $rgnse['end'];
    $helices .= "        rsr.add(secondary_structure.Alpha(self.residue_range('$rgs:A', '$rge:A')))\n";
}

$disulfs = "";
$xlinx = [ ["3.25", "45.62"] ];
foreach ($xlinx as $xl)
{
    try
    {
        $rno1 = resno_from_bw($rcpid, $xl[0]);
        $rno2 = resno_from_bw($rcpid, $xl[1]);
    }
    catch (Exception $e)
    {
        continue;
    }

    $raa1 = substr($p['sequence'], $rno1-1, 1);
    $raa2 = substr($p['sequence'], $rno2-1, 1);
    if ($raa1 == 'C' && $raa2 == 'C') $disulfs .=
         "        self.patch(residue_type='DISU', residues=(self.residues['$rno1:A'],\n"
        ."                                                  self.residues['$rno2:A']))\n";
}

if ($disulfs) $disulfs = "    def special_patches(self, aln):\n$disulfs";

$mdlcls = "DOPEHRLoopModel";

$CLASSI = "'8f76', '8hti'";
$TAAR1 = "'8jln', '8jlo', '8jlp', '8jlq', '8jlr', '8jso'";
$MTAAR9 = "'8iwe', '8iwm', '8itf', '8iw4', '8iw9'";
// $ADORA2A = "'6gdg'";
// $ADRB2 = "'7dhr', '8gej'";
$LPAR1 = "'7td0', '7yu3'";
$TAS2R = "'7xp6'";
$CB = "'5xr8', '5xra'";

$fam = family_from_protid($rcpid);
switch ($fam)
{
    case 'TAAR':
    if ($rcpid == "TAAR1")
    {
        $knowns = $TAAR1;
        $mdlcls = "AutoModel";
    }
    else if ($rcpid == "TAAR9") $knowns = $MTAAR9;
    else $knowns = "$MTAAR9, $TAAR1";
    break;

    case 'VN1R':
    $knowns = "$TAS2R, $CB";
    break;

    case 'MS4A':
    die("No HM templates known for MS4A receptors.\n");
    break;

    case 'OR51':
    case 'OR52':
    case 'OR56':
    if ($rcpid == "OR51E2") $knowns = "'8f76'";
    else $knowns = "$CLASSI";
    break;

    default:        // Class II ORs
    $knowns = "$CLASSI, $MTAAR9, $CB, $LPAR1";
}


// Generate Python script
$py = <<<natrixs

from modeller import *
from modeller.automodel import *

# log.verbose()
env = Environ()

class MyModel($mdlcls):
    def special_restraints(self, aln):
        rsr = self.restraints
        at = self.atoms
$helices
$disulfs

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

a = $mdlcls(env,
              alnfile  = 'allgpcr.ali',
              knowns   = ($knowns),
              sequence = '$rcpid')
a.starting_model = 0
a.ending_model   = 9
a.library_schedule = autosched.slow
a.max_var_iterations = 300
# a.md_level = refine.very_slow
a.make()

natrixs;

$fp = fopen("$rcpid.hm.py", "w");
fwrite($fp, $py);
fclose($fp);

@unlink("hm.out");
passthru("python3 $rcpid.hm.py | tee $rcpid.hm.out");
$c = file_get_contents("$rcpid.hm.out");
$best_energy = 1e9;
$pyoutfn = false;
$mode = false;
foreach (explode("\n", $c) as $ln)
{
    if (false !== strpos($ln, "Summary of successfully produced models:")) $mode = true;
    if (false !== strpos($ln, "Summary of successfully produced loop models:")) $mode = true;
    if (preg_match("/Filename\\s+molpdf/i", $ln)) $mode = true;
    if (!$mode) continue;

    $ln = preg_replace("/\\s+/", " ", $ln);
    $pieces = explode(" ", $ln);
    if (count($pieces) < 2) continue;
    if (!is_numeric($pieces[1])) continue;
    $e = floatval($pieces[1]);
    if ($e < $best_energy)
    {
        $best_energy = $e;
        $pyoutfn = $pieces[0];
    }
}

if (!$pyoutfn) die("FAIL.\n");

$famno = intval(preg_replace("/[^0-9]/", "", $fam));
$adjustments = "";
if ($famno < 50) $adjustments .= "IF $3.37 != \"G\" THEN ATOMTO %3.37 EXTENT @6.48\n";
else if ($famno == 51 || $famno == 52) $adjustments .= "ATOMTO %6.59 EXTENT @4.57\n";
else if ($rcpid == "OR56B2") $adjustments .= "ATOMTO %6.58 EXTENT @4.57\n";

$pepd = <<<blixtos

LET \$rcpid = "$rcpid"

LET \$inpf = "pdbs/$fam/$rcpid.upright.pdb"

LET \$mdld = "hm/$pyoutfn"

LOAD \$inpf A I
LOAD \$mdld A A

BWCOPY I A

STRAND I
UPRIGHT I
BWCENTER

STRAND A

IF $5.58 != "Y" GOTO _not_57_hbond
IF $7.53 != "Y" GOTO _not_57_hbond
MEASURE %5.58 "OH" %7.53 "OH" &d57
ECHO "Y5.58 and Y7.53 are " &d57 "A apart."
IF &d57 > 5.3 DIE "FAIL"
_not_57_hbond:

A100 &a100
ECHO "A100 Score: " &a100
IF &a100 < 10 DIE "FAIL"

DELETE 1 %1.20
HYDRO
UPRIGHT A
BWCENTER
$adjustments
MINC
$adjustments

# LET \$mdlf = "pdbs/" + \$rcpid
# LET \$mdlf += ".models.pdb"
# SAVE \$mdlf

UNCHAIN I
UNCHAIN O
STRAND A

LET \$outf = "pdbs/$fam/$rcpid.active.pdb"
SAVE \$outf

blixtos;

$fp = fopen("$rcpid.hm.pepd", "w");
fwrite($fp, $pepd);
fclose($fp);

chdir(__DIR__);
chdir("..");
passthru("./bin/pepteditor hm/$rcpid.hm.pepd");

chdir(__DIR__);
foreach (glob("$rcpid.*") as $doomed) unlink($doomed);
