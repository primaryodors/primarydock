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

$CLASSI = "'8f76', '8hti'";
$TAAR1 = "'8jln', '8jlo', '8jlp', '8jlq', '8jlr', '8jso'";
$MTAAR9 = "'8iwe', '8iwm', '8itf', '8iw4', '8iw9'";
$ADRB2 = "'7dhr', '8gej'";
$ADORA2A = "'6gdg'";
$LPAR1 = "'7td0', '7yu3'";

$fam = family_from_protid($rcpid);
switch ($fam)
{
    case 'TAAR':
    $knowns = "$MTAAR9, $TAAR1";
    break;

    case 'VN1R':
    $knowns = "$LPAR1, $TAAR1";
    break;

    case 'MS4A':
    die("No HM templates known for MS4A receptors.\n");
    break;

    case 'OR51':
    case 'OR52':
    case 'OR56':
    $knowns = "$CLASSI";
    break;

    default:        // Class II ORs
    $knowns = "$ADORA2A, $MTAAR9, $CLASSI";                 // This combination seems to produce the most energetically favorable models.
}


// Generate Python script
$py = <<<natrixs

from modeller import *
from modeller.automodel import *

log.verbose()
env = Environ()

class MyModel(DOPEHRLoopModel):
    def special_restraints(self, aln):
        rsr = self.restraints
        at = self.atoms
$helices
$disulfs

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

a = DOPEHRLoopModel(env,
              alnfile  = 'allgpcr.ali',
              knowns   = ($knowns),
              sequence = '$rcpid')
a.starting_model= 1
a.ending_model  = 1
a.library_schedule = autosched.slow
a.max_var_iterations = 300
# a.md_level = refine.very_slow
a.make()

natrixs;

$fp = fopen("hm.py", "w");
fwrite($fp, $py);
fclose($fp);

passthru("python3 hm.py");

$pyoutfn = false;
$dir = dir(__DIR__);
while (false !== ($entry = $dir->read()))
{
    if (preg_match("/^{$rcpid}[.]BL[0-9]+[.]pdb$/", $entry))
    {
        $pyoutfn = $entry;
        break;
    }
}
$dir->close();

if (!$pyoutfn) die("FAIL.\n");

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
DELETE 1 %1.20
HYDRO
UPRIGHT A
BWCENTER
MINC

# LET \$mdlf = "pdbs/" + \$rcpid
# LET \$mdlf += ".models.pdb"
# SAVE \$mdlf

UNCHAIN I
UNCHAIN O

LET \$outf = "pdbs/$fam/$rcpid.active.pdb"
SAVE \$outf

blixtos;

$fp = fopen("hm.pepd", "w");
fwrite($fp, $pepd);
fclose($fp);

chdir(__DIR__);
chdir("..");
passthru("./bin/pepteditor hm/hm.pepd");

chdir(__DIR__);
foreach (glob("$rcpid.*") as $doomed) unlink($doomed);
