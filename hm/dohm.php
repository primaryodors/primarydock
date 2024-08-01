<?php

chdir(__DIR__);
require_once("../data/protutils.php");
chdir(__DIR__);

if (!file_exists('/usr/lib/libmodeller.so')) die("This feature requires MODELLER.\nPlease see: https://salilab.org/modeller/ to obtain this third party software.\n\n");

$rcpid = @$argv[1];
$p = $prots[$rcpid];
if (!$rcpid) die("Usage:\nphp -f hm/dohm.php PROTID\n\n");

$disulfs = "";
$xlinx = [ ["3.25", "45.62"] ];
foreach ($xlinx as $xl)
{
    $rno1 = resno_from_bw($rcpid, $xl[0]);
    $rno2 = resno_from_bw($rcpid, $xl[1]);
    $raa1 = substr($p['sequence'], $rno1-1, 1);
    $raa2 = substr($p['sequence'], $rno2-1, 1);
    if ($raa1 == 'C' && $raa2 == 'C') $disulfs .=
         "        self.patch(residue_type='DISU', residues=(self.residues['$rno1:A'],\n"
        ."                                                  self.residues['$rno2:A']))\n";
}

if ($disulfs) $disulfs = "class MyModel(DOPEHRLoopModel):\n    def special_patches(self, aln):\n$disulfs";

$fam = family_from_protid($rcpid);
switch ($fam)
{
    case 'TAAR':
    $knowns = "'8iwe', '8iwm', '8itf', '8iw4', '8iw9', '8jln', '8jlo', '8jlp', '8jlq', '8jlr', '8jso'";
    break;

    case 'VN1R':
    $knowns = "'8iwe', '8iwm', '8itf', '8iw4', '8iw9', '8jln', '8jlo', '8jlp', '8jlq', '8jlr', '8jso'";         // ??? TODO: Find closest relatives of VN1Rs.
    break;

    case 'MS4A':
    die("No HM templates known for MS4A receptors.\n");
    break;

    case 'OR51':
    case 'OR52':
    case 'OR56':
    $knowns = "'8f76', '8hti'";
    break;

    default:        // Class II ORs
    $knowns = "'8iwe', '8iwm', '8itf', '8iw4', '8iw9', '8jln', '8jlo', '8jlp', '8jlq', '8jlr', '8jso'";
}


// Generate Python script
$py = <<<natrixs

from modeller import *
from modeller.automodel import *

log.verbose()
env = Environ()

$disulfs

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

a = DOPEHRLoopModel(env,
              alnfile  = 'allgpcr.ali',
              knowns   = ($knowns),
              sequence = '$rcpid')
a.starting_model= 1
a.ending_model  = 1
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
