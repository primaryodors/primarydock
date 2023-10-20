<?php
chdir(__DIR__);

$c = file_get_contents("../viewer.htm");

if (@$_REQUEST['view'] == "pred")
{
    require("../data/protutils.php");
    $protid = $_REQUEST["prot"];
    $fam = family_from_protid($protid);
    $odor = $_REQUEST["odor"];
    $mode = $_REQUEST["mode"];      // active or inactive.
    $n = @$_REQUEST["n"] ?: 1;

    $path = "../output/$fam/$protid/$protid.$odor.$mode.model$n.pdb";
    if (!file_exists($path)) $path = "../output/$fam/$protid/$protid.$odor.model$n.pdb";
    if (!file_exists($path)) die("Something went wrong.");
    $pdb = file_get_contents($path);

    $dock = "../output/$fam/$protid/$protid.$odor.$mode.dock";

    $ligbs = "	var lligbs = [";
    $pdblines = explode("\n", $pdb);
    $atom_lines = [];
    $hetatm_xyz = [];
    foreach ($pdblines as $line)
    {
        if (substr($line, 0, 7) == "ATOM   ") $atom_lines[] = $line;
        if (substr($line, 0, 7) == "HETATM ")
        {
            $x = floatval(substr($line, 31, 7));
            $y = floatval(substr($line, 39, 7));
            $z = floatval(substr($line, 47, 7));
            $hetatm_xyz[] = [$x, $y, $z];
        }
    }

    $ligand_residues = [];
    foreach ($atom_lines as $line)
    {
        $resno = intval(substr($line, 23, 3));
        if (isset($ligand_residues[$resno])) continue;
        $x = floatval(substr($line, 31, 7));
        $y = floatval(substr($line, 39, 7));
        $z = floatval(substr($line, 47, 7));

        foreach ($hetatm_xyz as $xyz)
        {
            $r = sqrt(pow($x-$xyz[0], 2) + pow($y-$xyz[1], 2) + pow($z-$xyz[2], 2));
            if ($r <= 5.0)
            {
                if (count($ligand_residues)) $ligbs .= ", ";
                $ligbs .= "$resno";
                $ligand_residues[$resno] = $resno;
                break;
            }
        }
    }
    $ligbs .= "];\n";

    $c = str_replace("	var lligbs = get_ligbs_from_orid();\n", $ligbs, $c);
    $c = str_replace("var literal_pdb = false;\n", "var literal_pdb = `$pdb`;\n", $c);
    $c = str_replace("var literal_fname = \"\";\n", "var literal_fname = \"$protid.$odor.$mode.model$n.pdb\";\n", $c);

    $d = explode("\n", file_get_contents($dock));
    $dockdisp = [];
    $ddidx = 0;
    $dockdisp[$ddidx] = "";
    foreach ($d as $lineno => $line)
    {
        $line2 = trim($d[$lineno+2]);
        if ($line2 == "PDBDAT:") break;

        $dockdisp[$ddidx] .= "$line\n";
    }

    $c .= <<<dockdata

<script>
$('#dockfloat span')[0].innerText = `{$dockdisp[0]}`;
</script>
dockdata;

}

echo $c;
