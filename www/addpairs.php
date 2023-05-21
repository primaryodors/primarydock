<?php

chdir(__DIR__);
require_once("../predict/protutils.php");
require_once("../predict/odorutils.php");
$odors = json_decode(file_get_contents("../data/odorant.json"), true);

$page_title = "Add Empirical Pairs";
include("header.php");

if (!file_exists("../www/security.json"))
{
    if (isset($_POST['apw']))
    {
        $salt = substr(md5(random_bytes(53)), 13, 2);
        $security = ['apw' => crypt($_POST['apw'], $salt), 'salt' => $salt];
        $fp = fopen("../www/security.json", "wb");
        if (!$fp) die("Unable to write security file; please check permissions.");
        else
        {
            fwrite($fp, json_encode($security));
            fclose($fp);
        }
    }
    else
    {
        ?>
        <form action="" method="POST">
        Please set an administrator password:
        <br>
        <input type="password" name="apw">
        <br>
        <input type="submit">
        </form>
        <?php
        exit;
    }
}
else
{
    $security = json_decode(file_get_contents("../www/security.json"), true);
}

if (!isset($security)) exit;

if (!isset($_REQUEST["token"]))
{
    if (isset($_POST["pw"]))
    {
        if (hash_equals($security["apw"], crypt($_POST["pw"], $security["salt"])))
        {
            $token = bin2hex(random_bytes(128));
            $security['token'] = $token;
            $security['ip'] = $_SERVER["REMOTE_ADDR"];
            $fp = fopen("../www/security.json", "wb");
            if (!$fp) die("Unable to update security file; please check permissions.");
            else
            {
                fwrite($fp, json_encode($security));
                fclose($fp);
            }
        }
        else die("Password incorrect.");
    }
    else
    {
        _relogin:
        ?>
        <form action="" method="POST">
        Please enter the administrator password:
        <br>
        <input type="password" name="pw">
        <br>
        <?php
        foreach ($_POST as $k => $v)
        {
            if ($k != "token") echo "<input type=\"hidden\" name=\"$k\" value=\"$v\">\n";
        }
        ?>
        <input type="submit">
        </form>
        <?php
        exit;
    }
}
else $token = $_REQUEST["token"];

if ($security['token'] != $token || $security['ip'] != $_SERVER["REMOTE_ADDR"]) goto _relogin;

if (isset($_POST["refurl"]))
{
    $new_ref = [ "short_name" => $_POST["sname"], "citation" => $_POST["cite"] ];
    if (!isset($refs[$_POST["refurl"]]))
    {
        if (!$_POST["sname"] || !$_POST["cite"]) die("Short name and citation are required.");
        else $refs[$_POST["refurl"]] = $new_ref;
    }

    $tallest = 0.0;
    for ($i=1; isset($_POST["prot$i"]); $i++)
    {
        $f = floatval(@$_POST["top$i"]);
        if (abs($f) > $tallest) $tallest = $f;
    }
    $scale = $tallest ? (10.0 / $tallest) : 1;

    for ($i=1; isset($_POST["prot$i"]); $i++)
    {
        if ($_POST["prot$i"] && $_POST["odor$i"])
        {
            $new_odor =
            [
                "full_name" => $_POST["odorfn$i"],
                "smiles" => $_POST["smiles$i"],
            ];

            if (@$_POST["olept$i"])
            {
                $olept = explode(" ", $_POST["olept$i"]);
                $new_odor["aroma"]["http://www.thegoodscentscompany.com"] = $olept;
            }

            if (@$_POST["top$i"]) $new_odor["activity"][$_POST["refurl"]][$_POST["prot$i"]]["adjusted_curve_top"]
                = floatval($_POST["top$i"]) * $scale;
            
            else if (@$_POST["top$i"] === "0") $new_odor["activity"][$_POST["refurl"]][$_POST["prot$i"]]["adjusted_curve_top"] = "0";

            if (@$_POST["ec$i"])
            {
                $ec50 = floatval($_POST["ec$i"]);

                // TODO: type of ec50, e.g. log10, mM, μM, etc.
                $ec50 = pow(10.0, $ec50);           // log10.

                if ($ec50 > 0 && $ec50 <= 1) $new_odor["activity"][$_POST["refurl"]][$_POST["prot$i"]]["ec50"] = $ec50;
            }

            if (@$_POST["ant$i"]) $new_odor["activity"][$_POST["refurl"]][$_POST["prot$i"]]["antagonist"] = 1;

            if ($_POST["odor$i"] == "new")
            {
                $oid = md5($new_odor["smiles"]);
            }
            else
            {
                $oid = $_POST["odor$i"];
            }

            if (isset($odors[$oid]))
            {
                $odors[$oid]["activity"][$_POST["refurl"]] = $new_odor["activity"][$_POST["refurl"]];
            }
            else
            {
                $odors[$oid] = $new_odor;
            }

            // echo "<pre>".json_encode_pretty($new_odor)."</pre>";
        }
    }

    $fp = fopen("../data/refs.json", "wb");
    if (!$fp) die("Unable to open data/refs.json; check permissions.");
    else
    {
        fwrite($fp, json_encode_pretty($refs));
        fclose($fp);
    }

    $fp = fopen("../data/odorant.json", "wb");
    if (!$fp) die("Unable to open data/odorant.json; check permissions.");
    else
    {
        fwrite($fp, json_encode_pretty($odors));
        fclose($fp);
    }

    echo "<h2>Success!</h2><p>You may optionally add more pairs if desired.</p>";
}

?>
<style>
body
{
    overflow: auto!important;
}

textarea
{
    border-radius: 4px;
    background-color: #aabbb7;
    color: #08101b;
    border: 1px solid #214577;
    margin: 2px;
}

div.pair
{
    margin: 15px;
    border: 1px solid #069;
    background-color: #235;
    padding: 15px;
}

div.pair label
{
    width: 7%;
    display: inline-block;
}
</style>

<script language="Javascript">
var lastselor = false;
function new_pair()
{
    var pairno = $("div.pair").length + 1;

    var ctnr = $("#pairs")[0];
    var div = document.createElement("div");
    div.className = "pair";
    div.id = "pair" + pairno;

    var lbl = document.createElement("label");
    lbl.innerText = "Prot:";
    div.appendChild(lbl);

    var ctl = document.createElement("select");
    ctl.name = ctl.id = "prot" + pairno;
    var opt = document.createElement("option");
    opt.value = "";
    opt.innerText = "-----";
    ctl.appendChild(opt);
    <?php
    foreach (array_keys($prots) as $rcpid)
    {
        ?>
        opt = document.createElement("option");
        opt.value = "<?php echo $rcpid; ?>";
        opt.innerText = opt.value;
        ctl.appendChild(opt);

        <?php
    }
    ?>
    if (lastselor) ctl.value = lastselor;
    ctl.setAttribute("onchange", "lastselor = this.value;");
    div.appendChild(ctl);
    div.appendChild(document.createElement("br"));

    var lbl = document.createElement("label");
    lbl.innerText = "Odor:";
    div.appendChild(lbl);

    ctl = document.createElement("select");
    ctl.name = ctl.id = "odor" + pairno;
    opt = document.createElement("option");
    opt.value = "";
    opt.innerText = "-----";
    ctl.appendChild(opt);
    <?php
    foreach ($odors as $oid => $odor)
    {
        ?>
        opt = document.createElement("option");
        opt.value = "<?php echo $oid; ?>";
        opt.innerText = "<?php echo $odor["full_name"]; ?>";
        ctl.appendChild(opt);

        <?php
    }
    ?>
    opt = document.createElement("option");
    opt.value = "new";
    opt.innerText = "(add new)";
    ctl.appendChild(opt);

    ctl.setAttribute("onchange", "(this.value == 'new') ? new_odorant(this.parentElement) : new_pair();");
    div.appendChild(ctl);
    div.appendChild(document.createElement("br"));

    lbl = document.createElement("label");
    lbl.innerText = "Top:";
    div.appendChild(lbl);
    ctl = document.createElement("input");
    ctl.name = ctl.id = "top" + pairno;
    div.appendChild(ctl);
    div.appendChild(document.createElement("br"));

    lbl = document.createElement("label");
    lbl.innerHTML = "EC<sub>50</sub>:";
    div.appendChild(lbl);
    ctl = document.createElement("input");
    ctl.name = ctl.id = "ec" + pairno;
    div.appendChild(ctl);
    div.appendChild(document.createElement("br"));

    lbl = document.createElement("label");
    lbl.innerHTML = "Antagonist:";
    div.appendChild(lbl);
    ctl = document.createElement("input");
    ctl.type = "checkbox";
    ctl.name = ctl.id = "ant" + pairno;
    div.appendChild(ctl);
    div.appendChild(document.createElement("br"));

    ctnr.appendChild(div);
}

function new_odorant(div)
{
    var pairno = parseInt(div.id.substr(4));

    lbl = document.createElement("label");
    lbl.innerHTML = "Odor:";
    div.appendChild(lbl);
    ctl = document.createElement("input");
    ctl.name = ctl.id = "odorfn" + pairno;
    div.appendChild(ctl);
    div.appendChild(document.createElement("br"));

    lbl = document.createElement("label");
    lbl.innerHTML = "SMILES:";
    div.appendChild(lbl);
    ctl = document.createElement("input");
    ctl.name = ctl.id = "smiles" + pairno;
    ctl.setAttribute("onchange", "new_pair();");
    div.appendChild(ctl);
    div.appendChild(document.createElement("br"));

    lbl = document.createElement("label");
    lbl.innerHTML = "Aroma:";
    div.appendChild(lbl);
    ctl = document.createElement("input");
    ctl.name = ctl.id = "olept" + pairno;
    ctl.setAttribute("onchange", "new_pair();");
    div.appendChild(ctl);
    div.appendChild(document.createElement("br"));
}
</script>
<h2>Add Empirical Receptor-Ligand Pairs</h2>

<form action="" method="POST">
<input type="hidden" name="token" value="<?php echo $token; ?>">
<label for="refurl">Reference URL:</label>
<input type="text" name="refurl" style="width: 81%;">
<br>
<label for="sname">Short Name:</label>
<input type="text" name="sname" style="width: 31%;">
<br>
<label for="cite">Citation:</label>
<textarea name="cite" style="width: 100%; height: 111px;"></textarea>
<br>

<div id="pairs">
</div>
<br>

<input type="submit" value="Save">

<script>
new_pair();
</script>
