<?php
chdir(__DIR__);
require_once("../data/protutils.php");
require_once("../data/odorutils.php");

include("header.php");

?>
<h1>References</h1>

<div class="scrollh" style="height: 750px;">
<?php
foreach ($refs as $url => $r)
{
    if (substr($r['short_name'], 0, 1) == '(') continue;
    echo "<p class=\"ref\">\n";
    echo $r['citation'];
    echo " <a href=\"$url\">&#127891;</a>\n";
    echo "</p>\n\n";
}

