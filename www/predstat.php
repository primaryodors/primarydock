<?php
chdir(__DIR__);
require_once("../data/protutils.php");
require_once("../data/odorutils.php");
require_once("../predict/statistics.php");

$dr = json_decode(file_get_contents("../predict/dock_results.json"), true);

$corr = $fail = 0;

foreach ($dr as $rcpid => $ligands)
{
    foreach ($ligands as $ligand => $result)
    {
        $actual = strtolower(@$result["Actual"]);
        $predict = strtolower(@$result["Predicted"]);

        if (!$actual || $actual == "(unknown)")
        {
            continue;
        }
        else if ($actual == "agonist")
        {
            if ($predict == "agonist") $corr++;
            else $fail++;
        }
        else
        {
            if ($predict == "agonist") $fail++;
            else $corr++;
        }
    }
}

$percent = round(100.0 * $corr / ($corr+$fail), 2);

$page_title = "$percent% Prediction Status";
include("header.php");

?><h1>Prediction Status</h1>

<label for="pa">Prediction Accuracy:</label><br>
<progress id="pa" value="<?php echo $percent; ?>" max="100" style="width: 90%;"></progress> <?php echo $percent; ?>%

<script language="JavaScript">
window.setTimeout(function()
{
    window.location.reload();
}, 5*60*1000);
</script>
