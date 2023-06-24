<?php

function die_if_too_hot()
{
    $results = [];
    exec("which sensors", $results);            // sudo apt-get install lm-sensors
    if (count($results))
    {
        $results = [];
        exec("sensors", $results);
        foreach ($results as $ln)
        {
            if (preg_match("/:\\s+[+-][0-9]{1,3}[.][0-9]*[^(]*[(]high =/", $ln))
            {
                $ln = explode("(high = ", $ln);
                $ln[0] = preg_replace("/^.*:\\s*/", "", $ln[0]);
                $limit = floatval($ln[1]);
                if ($limit)
                {
                    $limit -= 6;
                    if (floatval($ln[0]) >= $limit) die("Computer too hot.\n");
                }
            }
            else if (preg_match("/:\\s+[+-][0-9]{1,3}[.][0-9]*[^(]*[(]crit =/", $ln))
            {
                $ln = explode("(crit = ", $ln);
                $ln[0] = preg_replace("/^.*:\\s*/", "", $ln[0]);
                $limit = floatval($ln[1]);
                if ($limit)
                {
                    $limit -= 20;
                    if (floatval($ln[0]) >= $limit) die("Computer too hot.\n");
                }
            }
        }
    }
}
