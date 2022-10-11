<?php
chdir(__DIR__);
require_once("../predict/protutils.php");
require_once("../predict/odorutils.php");

include("header.php");

?>
<h1>Home</h1>

<p>This web application provides visual representations of the JSON files located in the data subfolder.
</p>

<center>
<table class="idxbuttons">
    <tr>
        <td>
            <a href="receptors.php"><img src="assets/receptors.png"></a>
        </td>
        <td>
            <a href="odorants.php"><img src="assets/odorants.png"></a>
        </td>
    </tr>
</table>
</center>
