<?php

chdir(__DIR__);
chdir("..");
require_once("predict/odorutils.php");

foreach (@$argv as $k => $a)
{
    if ($k > 0) ensure_sdf_exists($a);
}

