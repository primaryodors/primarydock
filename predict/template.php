<?php
global $template, $protid, $cryoem, $has_rock6, $has_fyg, $args, $accuracy_receptors;

$accuracy_receptors = ["OR51E2", "TAAR1"];
$template = [];
chdir(__DIR__);
chdir("..");
$cryoem = json_decode(file_get_contents("data/cryoem_motions.json"), true);

foreach (@$argv as $a)
{
	$a = explode('=',$a,2);
	$_REQUEST[$a[0]] = (count($a)>1) ? $a[1] : true;
}

$max_simul = 4;
if (@$_REQUEST["simul"]) $max_simul = intval($_REQUEST["simul"]);

function build_template()
{
    global $template, $protid, $cryoem, $has_rock6, $has_fyg, $args;

    if (substr($protid, 0, 4) == "TAAR")
    {
        $template = $cryoem["mTAAR9"];
        // TODO: Blend mTAAR9 with TAAR1 depending on sequence similarity.
    }
    else if (substr($protid, 0, 2) == "OR")
    {
        $fam = intval(substr($protid, 2, 2));
        if ($fam == 51)
        {
            $template = $cryoem["OR51E2"];
        }
        else if ($fam == 52)
        {
            foreach ($cryoem["OR52"] as $hxno => $metrics)
            {
                foreach ($metrics as $metric => $dimensions)
                {
                    $template[$hxno][$metric] = $cryoem["OR52"][$hxno][$metric];

                    $template[$hxno][$metric]["sigma"]
                        + 0.4 * $cryoem["TAAR1"][$hxno][$metric]["sigma"]
                        + 0.6 * $cryoem["mTAAR9"][$hxno][$metric]["sigma"];
                }
            }
        }
        else
        {
            foreach ($cryoem["OR51E2"] as $hxno => $metrics)
            {
                foreach ($metrics as $metric => $dimensions)
                {
                    foreach (array_keys($dimensions) as $dimension)
                    {
                        if ($dimension == "sigma")
                        {
                            $template[$hxno][$metric][$dimension]
                                = 0.4 * $cryoem["TAAR1"][$hxno][$metric][$dimension]
                                + 0.6 * $cryoem["mTAAR9"][$hxno][$metric][$dimension];
                        }
                        else
                        {
                            // TODO: Make the proportions dependent on sequence similarity.
                            $template[$hxno][$metric][$dimension]
                                = 0.40 * $cryoem["OR51E2"][$hxno][$metric][$dimension]
                                + 0.30 * $cryoem["OR52"][$hxno][$metric][$dimension]
                                + 0.10 * $cryoem["TAAR1"][$hxno][$metric][$dimension]
                                + 0.20 * $cryoem["mTAAR9"][$hxno][$metric][$dimension];
                        }
                    }
                }
            }
        }
    }

    // Get typology. Do a no-save run of fyg_activate_or to determine whether protein has FYG or rock6 capabilities.
    $output = [];
    exec("bin/fyg_activate_or --nosave $protid", $output);
    $has_rock6 = $has_fyg = false;
    foreach ($output as $line)
    {
        if (false!==strpos($line, "Performing rock6")) $has_rock6 = true;
        if (false!==strpos($line, "Performing FYG activation")) $has_fyg = true;
    }
}

function check_already_fyg_activating($protid)
{
    global $max_simul;
    $output = [];
    exec("ps -ef | grep -e \"[0-9] bin/fyg_\" | grep -v grep", $output);
    foreach ($output as $o)
    {
        if (false!== strpos($o, $protid)) return true;                // This will give a false positive for e.g. OR2T2 if OR2T29 is running, but this is a minor bug.
    }

    if (count($output) > $max_simul) die("Too many simultaneous activations.");

    if (file_exists("tmp/$protid.fyg"))
    {
        $pid = file_get_contents("tmp/$protid.fyg");
        if (file_exists("/proc/$pid")) return true;
        else unlink("tmp/$protid.fyg");
    }
    else return false;
}

function do_templated_activation()
{
    global $template, $args, $has_fyg, $has_rock6, $protid, $accuracy_receptors;

    $mypid = getmypid();
    exec("echo \"$mypid\" > tmp/$protid.fyg");
    $args = "$protid";

    build_template();

    foreach ($template as $hxno => $metrics)
    {
        foreach ($metrics as $metric => $dimensions)
        {
            if ($hxno == 6)
            {
                if ($metric == "cyt" && ($has_fyg || $has_rock6)) continue;
                else if ($metric == "exr" && $has_rock6) continue;
            }

            $cmdarg = "--" . substr($metric, 0, 1) . $hxno;
            $args .= " $cmdarg {$dimensions['x']} {$dimensions['y']} {$dimensions['z']}";
        }
    }
    $cmd = "bin/fyg_activate_or $args";
    echo "$cmd\n";
    passthru($cmd);
    unlink("tmp/$protid.fyg");

    if (in_array($protid, $accuracy_receptors) && file_exists("devenv"))
    {
        $c = file_get_contents("predict/model_accuracy.sh");
        $lines = explode("\n", $c);
        foreach ($lines as $lno => $ln)
        {
            if (false!==strpos($ln, "bin/fyg_activate_or $protid "))
            {
                $lines[$lno] = "make bin/fyg_activate_or && $cmd && bin/pepteditor predict/model_accuracy.pepd";
            }
        }
        $c = implode("\n", $lines);
        $fp = fopen("predict/model_accuracy.sh", "w");
        if ($fp)
        {
            fwrite($fp, $c);
            fclose($fp);
        }
    }
}
