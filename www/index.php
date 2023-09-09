<?php
chdir(__DIR__);
require_once("../data/protutils.php");
require_once("../data/odorutils.php");

if (file_exists("tags.php")) include("tags.php");

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
            <a href="receptors.php"><img src="assets/receptors.png" alt="List of Receptors"></a>
        </td>
        <td>
            <a href="odorants.php"><img src="assets/odorants.png" alt="List of Odorants"></a>
        </td>
    </tr>
</table>
</center>
