<?php 
$sidelinks = [];

$cwd = getcwd();
chdir(__DIR__);
if (file_exists('sidelinks.json')) $sidelinks = json_decode(file_get_contents('sidelinks.json'), true);
chdir($cwd);

?>
<ul>
    <li><a href="index.php">Home</a></li>
    <li><a href="receptors.php">Receptors</a></li>
    <li><a href="odorants.php">Odorants</a></li>
    <li><a href="viewer.php" target="_VIEWER">3D Viewer</a></li>

    <?php
    foreach ($sidelinks as $sl)
    {
        echo "<li><a href=\"{$sl['url']}\"";
        if (@$sl['target']) echo " target=\"{$sl['target']}\"";
        echo ">{$sl['name']}</a></li>\n";
    }
    ?>
</ul>
