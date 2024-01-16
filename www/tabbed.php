<?php

require_once("../data/protutils.php");

$method = "icactive";

switch (@$_REQUEST["obj"])
{
    case "model":
        http_response_code(401);
    
    case "dock":
        $prot  = @$_REQUEST["prot"];
        $lig   = @$_REQUEST["odor"];
        $mode  = @$_REQUEST["mode"];

        $fam = family_from_protid($prot);

        $fn = "$prot.$lig.$mode.dock";
        $path = "../output/$fam/$prot/$fn";
        break;
    
    case "json":
        http_response_code(401);
    
    default:
        http_response_code(401);
        exit;
}

chdir(__DIR__);
if (!file_exists($path))
{
    http_response_code(404);
    exit;
}

$c = file_get_contents($path);
$c = trim(explode("Original PDB:", $c)[0]);
$poses = explode("Pose: ", $c);

$html = "<pre>".trim($poses[0])."</pre>\n";
$htmltop = "<div class=\"tab\" style=\"display: inline-block; width: 100%;\">\n";
$htmltop .= "<button class=\"tabstatic\" id=\"tabPoses\"><b>Poses:</b>";
$htmlbdy = "";

$first = true;
foreach ($poses as $p)
{
    $pno = intval($p);
    if (!$pno) continue;

    $tabclass = $first ? "tablinks active" : "tablinks";
    $htmltop .= "<button class=\"$tabclass\" id=\"tabPose$pno\"";
    if ($first) $htmltop .= " default";
    $htmltop .= " onclick=\"openTab(this, 'Pose$pno');\">$pno</button>\n";

    $p = explode("\n", $p, 2)[1];
    $htmlbdy .= "<div id=\"Pose$pno\" class=\"tabcontent\" style=\"height: 546px; overflow-y: scroll;";
    if ($first) $htmlbdy .= " display: block;";
    $htmlbdy .= "\">\n";
    $htmlbdy .= "<textarea style=\"width: 100%; height: 100%;\">$p</textarea>\n";
    $htmlbdy .= "</div>\n";

    $first = false;
}
$htmltop .= "</div>";

$html .= $htmltop . $htmlbdy;
echo $html;
