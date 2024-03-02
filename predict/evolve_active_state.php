<?php
chdir(__DIR__);
require_once("../data/protutils.php");
require_once("../data/odorutils.php");

$mutation_rate = 0.25;

if (!file_exists("../../OR5K1_binding_site/OR5K1_IFD3_models/MM_IFD3/1015_MM_1_IFD3_cmpd1_1.pdb"))
{
    if (!file_exists("../../OR5K1_binding_site/OR5K1_IFD3_models.zip"))
    {
        chdir("../..");
        exec("git clone https://github.com/dipizio/OR5K1_binding_site.git");
        chdir("OR5K1_binding_site");
        exec("unzip OR5K1_IFD3_models.zip");
        chdir(__DIR__);
    }
}

if (!file_exists("../../OR5K1_binding_site/OR5K1_IFD3_models/MM_IFD3/1015_MM_1_IFD3_cmpd1_1.pdb"))
{
    echo "To use this feature, please clone the following repository alongside primarydock:\n";
    echo "git clone https://github.com/dipizio/OR5K1_binding_site.git\n";
    echo "Make sure to uncompress the two .zip files.\n";
    exit;
}

function frand($min, $max)
{
    return 1e-6 * rand($min*1e6, $max*1e6);
}

function runpepd($values, $save = false)
{
    global $tne, $argv, $prots, $odors, $pairs, $orid, $usrmesgd;
    $pepd = [];
    foreach ($tne as $ln)
    {
        if (substr($ln, 0, 5) == "SAVE ") continue;
        else if (substr($ln, 0, 5) == "LOAD ")
        {
            $orid = (@$argv[1] && isset($prots[$argv[1]])) ? $argv[1] : explode(".upright.",explode("/", $ln)[2])[0];
            $fam = family_from_protid($orid);
            $ln = "LOAD \"pdbs/$fam/$orid.upright.pdb\"";
        }
        else if (substr($ln, 0, 7) == "SCREPL ")
        {
            $orid = (@$argv[1] && isset($prots[$argv[1]])) ? $argv[1] : explode(".upright.",explode("/", $ln)[2])[0];
            $fam = family_from_protid($orid);
            chdir(__DIR__);
            if (file_exists("../tmp/$orid.screpl.pdb"))
                $ln = "LOAD tmp/$orid.screpl.pdb";
            else
                $ln = "SCREPL \"pdbs/$fam/$orid.upright.pdb\"\nSAVE tmp/$orid.screpl.pdb";
        }
        else if (substr($ln, 0, 16) == "LET \$usrmesg = \"")
        {
            $quot = strpos($ln, '"');
            $usrmesg = str_replace('"', '', substr($ln, $quot));
            if (!$usrmesgd) echo "Message to user: $usrmesg\n";
            $usrmesgd = true;
        }
        else if (substr($ln, 0, 14) == "LET \$ligand = ")
        {
            if (!$pairs) $pairs = array_values(all_empirical_pairs_for_receptor($orid));
            if (!$pairs || !count($pairs)
                || @$pairs[0]["adjusted_curve_top"] < 0
                || (!@$pairs[0]["adjusted_curve_top"] && !$pairs[0]["ec50"]))
            {
                if (@$argv[2]) $ligand = $argv[2];
                else die("No agonists for $orid.\n");
            }
            else
            {
                $ligand = $odors[$pairs[0]['oid']]['full_name'];
            }

            $ligfn = ensure_sdf_exists($ligand);
            $ln = "LET \$ligand = \"$ligfn\"";
        }
        else if (substr($ln, 0, 4) == "LET ")
        {
            $words = explode(" ", $ln);
            if ($words[2] == "=")
            {
                $varname = $words[1];
                $vartyp = substr($varname, 0, 1);
                $varname = substr($varname, 1);
                if (isset($values[$varname]))
                {
                    $words[3] = $values[$varname];
                    $ln = implode(" ", $words);
                }
                else if ($vartyp == "@")
                {
                    if (isset($values["$varname.x"]) && isset($values["$varname.y"]) && isset($values["$varname.z"]))
                    {
                        $words[3] = "[{$values["$varname.x"]},{$values["$varname.y"]},{$values["$varname.z"]}]";
                        $ln = implode(" ", $words);
                    }
                }
            }
        }
        $pepd[] = $ln;
    }

    chdir(__DIR__);
    chdir("..");
    if ($save) $pepd[] = "SAVE pdbs/$fam/$orid.evolved.pdb";
    $pepdname = $save ? "pdbs/$fam/$orid.evolved.pepd" : "tmp/evolving.pepd";
    $fp = fopen($pepdname, "wb");
    if (!$fp) die("FAILED to open $pepdname for writing.\n");
    fwrite($fp, implode("\n", $pepd));
    fclose($fp);

    $output = [];
    exec("bin/pepteditor $pepdname", $output);
    if (false!==strpos(implode("\n", $output), " splice.php "))
    {
        exec("php -f predict/splice.php");
        $output = [];
        exec("bin/pepteditor $pepdname", $output);
    }
    return $output;
}

function score_result($result)
{
    $score = 0.0;
    foreach ($result as $ln)
    {
        if (false!==strpos($ln, " should be "))
        {
            $words = explode(" ", $ln);
            $shouldbe0 = 0;
            $shouldbe1 = 0;
            $is = 0;

            $mult = 1.0;
            if (false!==strpos($ln, "clashes")) $mult = 0.001;

            foreach ($words as $i => $w)
            {
                if (substr($w, -1) == ':') $is = floatval($words[$i+1]);
                else if ($w == "-" || $w == " to " || $w == "and")
                {
                    $shouldbe0 = floatval($words[$i-1]);
                    $shouldbe1 = floatval($words[$i+1]);
                }
                else if ($w == "less" && $words[$i+1] == "than")
                {
                    $shouldbe0 = 0;
                    $shouldbe1 = floatval($words[$i+2]);
                }
                else if ($w == "<")
                {
                    $shouldbe0 = 0;
                    $shouldbe1 = floatval($words[$i+1]);
                }
            }

            if ($shouldbe1 && $is)
            {
                if ($shouldbe0)
                {
                    if ($score >= $shouldbe0 && $score <= $shouldbe1)
                        $score += $mult * (0.5+0.5*cos((($is-$shouldbe0) / ($shouldbe1-$shouldbe0) - 0.5)*6.28));
                    else $score -= $mult * (min(abs($shouldbe0 - $score), abs($score - $shouldbe1)) / ($shouldbe1-$shouldbe0));
                }
                else $score += $mult * (($shouldbe1-$is) / $shouldbe1);
            }
        }
    }
    return $score;
}

chdir(__DIR__);
chdir("..");
$tne = explode("\n", file_get_contents((@$argv[1] && file_exists($argv[1])) ? $argv[1] : "predict/evolvable_activation.pepd"));

$evparams = [];
$evtypes = [];
$varsearch = false;
$pairs = false;
$usrmesgd = false;
foreach ($tne as $ln)
{
    if (trim($ln) == "# Configurable variables.") $varsearch = true;
    else if (trim($ln) == "# Internal variables.") $varsearch = false;
    else if ($varsearch && (substr($ln, 0, 4) == "LET "))
    {
        $words = explode(" ", $ln);
        if ($words[2] != "=") continue;
        $varname = $words[1];
        $vartyp = substr($varname, 0, 1);
        $varname = substr($varname, 1);

        switch ($vartyp)
        {
            case "$":
            continue 2;

            case "&":
            $evparams[$varname] = floatval($words[3]);
            $evtypes[$varname] = "f";
            break;

            case "%":
            $evparams[$varname] = intval($words[3]);
            $evtypes[$varname] = "i";
            break;

            case "@":
            $xyz = explode(",",str_replace("[", "", str_replace("]", "", $words[3])));
            $evparams["$varname.x"] = floatval($xyz[0]);
            $evtypes["$varname.x"] = "x";
            $evparams["$varname.y"] = floatval($xyz[1]);
            $evtypes["$varname.y"] = "y";
            $evparams["$varname.z"] = floatval($xyz[2]);
            $evtypes["$varname.z"] = "z";
            break;

            default:
            continue 2;
        }
    }
}

if (!count($evparams)) die("Error: parameters not delineated.\n");

for ($generation=1; $generation<=1000000; $generation++)
{
    echo "Beginning generation $generation...\n";

    // Create 20 hybrids. If there are no parents yet, create 20 mutants instead.
    $population = [];
    if (isset($best) && isset($secondbest))
    {
        for ($i=0; $i<20; $i++)
        {
            foreach (array_keys($evparams) as $param)
            {
                $value = rand(0,1) ? $secondbest[$param] : $best[$param];
                if (frand(0,1) <= $mutation_rate) $value += frand(-1, 1);
                $population[$i][$param] = $value;
            }
        }
    }
    else
    {
        for ($i=0; $i<20; $i++)
        {
            foreach ($evparams as $param => $value)
            {
                if (frand(0,1) <= $mutation_rate) $value += frand(-1, 1);
                $population[$i][$param] = $value;
            }
        }
    }

    echo "Spawned population.\n";

    // Run each individual through pepteditor.
    foreach ($population as $i => $individual)
    {
        $result = runpepd($individual);
        echo ".";

        // Score the results.
        $score = score_result($result);

        // Choose the best two individuals for the next generation.
        if (!isset($best) || $score > $best_score)
        {
            if (isset($best))
            {
                $secondbest = $best;
                $secondbest_score = $best_score;
            }

            $best_score = $score;
            $best = $individual;
        }
        else if (!isset($secondbest) || $score > $secondbest_score)
        {
            $secondbest_score = $score;
            $secondbest = $individual;
        }
    }
    echo "\n";

    echo "Best score: $best_score\n";
    runpepd($best, true);
}

