<?php
chdir(__DIR__);
require_once("../data/protutils.php");
require_once("../data/odorutils.php");

include("header.php");

?>
<h1>Home</h1>

<?php if (@$customizations['welcome']) echo $customizations['welcome'];
else
{ ?>
<p>This web application provides visual representations of the JSON files located in the data subfolder.
</p>
<?php } ?>

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
