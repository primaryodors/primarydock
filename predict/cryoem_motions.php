<?php

chdir(__DIR__);
require_once("../data/protutils.php");
chdir("..");

$data_keys = ["x.50", "exr", "cyt", "vec", "rot"];

$output = [];
echo "Executing pepd script...\n";
exec("bin/pepteditor predict/bnav2.pepd", $output);

$prot = false;
$data = [];
foreach ($output as $line)
{
    if (substr($line, 0, 6) == "Wrote ")
    {
        $pieces = explode("/", $line);
        if (count($pieces) > 1)
        {
            $pieces = explode(".", $pieces[1]);
            if (count($pieces) >= 2) $prot = trim($pieces[0]);
        }
    }
    else if (trim($line) == "") $prot = false;
    else if ($prot && substr($line, 1, 3) == ": [")
    {
        $hxno = intval(substr($line, 0, 2));
        $pieces = explode("|", substr($line, 3));

        if ($hxno == 1) $hxdata = [];
        foreach($pieces as $k => $v)
        {
            $v = trim($v);
            if (substr($v, 0, 1) == '[')
            {
                $v = str_replace("[", "", trim($v));
                $v = str_replace("]", "", trim($v));

                $xyz = explode(",",$v);
                $hxdata[$hxno][$data_keys[$k]] = ['x' => $xyz[0], 'y' => $xyz[1], 'z' => $xyz[2]];
            }
            else
                $hxdata[$hxno][$data_keys[$k]] = floatval($v);
        }
        if ($hxno == 7) $data[$prot][] = $hxdata;
    }
}

$average = [];
foreach ($data as $prot => $values)
{
    foreach ($values as $helices)
    {
        foreach ($helices as $helix => $metrics)
        {
            foreach ($metrics as $metric => $xyz)
            {
                if (is_array($xyz) || is_object($xyz))
                {
                    foreach ($xyz as $dimension => $cartesian)
                    {
                        if (!isset($average[$prot][$helix][$metric][$dimension])) $average[$prot][$helix][$metric][$dimension] = 0.0;
                        $average[$prot][$helix][$metric][$dimension] += $cartesian;
                    }
                }
                else
                {
                    if (!isset($average[$prot][$helix][$metric])) $average[$prot][$helix][$metric] = 0.0;
                    $average[$prot][$helix][$metric] += $xyz;
                    // echo "$metric = $xyz\n";
                }
            }
        }
    }

    $count = count($values);
    if ($count)
    {
        foreach ($helices as $helix => $metrics)
        {
            foreach ($metrics as $metric => $xyz)
            {
                if (is_array($xyz))
                    foreach ($xyz as $dimension => $cartesian)
                    {
                        $average[$prot][$helix][$metric][$dimension] /= $count;
                    }
                else
                    $average[$prot][$helix][$metric] /= $count;
            }
        }
    }
}

$deviation = [];
foreach ($data as $prot => $values)
{
    foreach ($values as $helices)
    {
        foreach ($helices as $helix => $metrics)
        {
            foreach ($metrics as $metric => $xyz)
            {
                if (is_array($xyz))
                {
                    if (!isset($deviation[$prot][$helix][$metric])) $deviation[$prot][$helix][$metric] = 0.0;
                    $radius_sqrd = 0.0;
                    foreach ($xyz as $dimension => $cartesian)
                    {
                        $delta = $cartesian - $average[$prot][$helix][$metric][$dimension];
                        $delta *= $delta;
                        
                        $radius_sqrd += $delta;
                    }
                    $deviation[$prot][$helix][$metric] += $radius_sqrd;
                }
            }
        }
    }
}

$output = $average;

foreach ($deviation as $prot => $helices)
{
    $count = count($data[$prot]);
    foreach ($helices as $helix => $metrics)
    {
        foreach ($metrics as $metric => $xyz)
        {
            if (is_array($output[$prot][$helix][$metric])) foreach ($output[$prot][$helix][$metric] as $dimension => $cartesian) $output[$prot][$helix][$metric][$dimension] = round($cartesian, 4);

            if ($count) $deviation[$prot][$helix][$metric] /= $count;
            $output[$prot][$helix][$metric]["sigma"] = round(sqrt($deviation[$prot][$helix][$metric]), 5);
        }
    }
}

$ofname = "data/cryoem_motions.json";
$fp = fopen($ofname, "wb");
if (!$fp) die("Failed to open $ofname for writing.\n");
fwrite($fp, json_encode_pretty($output));
fclose($fp);
echo "Wrote $ofname.\n";

