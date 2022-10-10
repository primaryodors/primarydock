<?php
chdir(__DIR__);
require_once("../predict/protutils.php");
require_once("../predict/odorutils.php");

include("header.php");

?>
<h1>Odorants List</h1>

<div class="scrollh" style="height: 780px;">
<?php

foreach ($odors as $odor)
{
    echo "<a href=\"odorant.php?o={$odor['full_name']}\">{$odor['full_name']}</a><br>\n";
}