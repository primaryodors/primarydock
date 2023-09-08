<?php

require_once("../data/protutils.php");

$method = "icactive";

switch (@_REQUEST["obj"])
{
    case "model":
        $prot = @$_REQUEST["prot"];
        $lig  = @$_REQUEST["lig" ];
        $mode = @$_REQUEST["mode"];

        $fam = family_from_protid($prot);

        $fn = "$prot.$lig.$mode.model$mdlno.pdb";
        $path = "../output/$fam/$prot/$fn";
        break;
    
    case "dock":
        $prot  = @$_REQUEST["prot"];
        $lig   = @$_REQUEST["lig" ];
        $mode  = @$_REQUEST["mode"];
        $mdlno = @$_REQUEST["mdl"] ?: 1;

        $fam = family_from_protid($prot);

        $fn = "$prot.$lig.$mode.dock";
        $path = "../output/$fam/$prot/$fn";
        break;
    
    case "json":
        $prot  = @$_REQUEST["prot"];
        $lig   = @$_REQUEST["lig" ];

        $json_file = "predict/dock_results_$method.json";
        $dock_results = json_decode(file_get_contents($json_file), true);

        if (!isset($dock_results[$prot][$lig]))
        {
            http_response_code(404);
            exit;
        }

        $result = [];
        $result[$prot][$lig] = $dock_results[$prot][$lig];
        header("Content-Disposition: attachment; filename=\"$fn\"");
        echo json_encode_pretty($result);

        exit;
    
    default:
        http_response_code(401);
        exit;
}

chdir(__DIR__);
if (!file_exists($fn))
{
    http_response_code(404);
    exit;
}

$c = file_get_contents($fn);
header("Content-Disposition: attachment; filename=\"$fn\"");
echo $c;
