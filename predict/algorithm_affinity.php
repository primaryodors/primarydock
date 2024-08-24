<?php

function make_prediction($data)
{
    global $protid, $ligname, $pose;

    if (isset($data["a_Pose1"]) || isset($data["a_BindingEnergy"]))
    {
        $ascore = max(0, -floatval(@$data['a_Pose1']));
        $iscore = max(0, -floatval(@$data['i_Pose1']));

        $aavg = max(0, -floatval(@$data['a_BindingEnergy']));
        $iavg = max(0, -floatval(@$data['i_BindingEnergy']));

        $aa100 = floatval(@$data['a_A100']);
        $ia100 = floatval(@$data['i_A100']);
        $attns = floatval(@$data['a_occlusion'] ?: 1);
        $ittns = floatval(@$data['i_occlusion'] ?: 1);

        if ($ascore > 0 && $aa100 > 0 && $ascore > $iscore)
        {
            $data['Predicted'] = 'Agonist';
            $data['Affinity'] = round($ascore, 4);
            $data['A100'] = round($aa100, 4);
            $data['DockScore'] = round(min($ascore, 150) * $attns * $aa100 / 100, 4);
        }
        else if ($iscore > 0)
        {
            $data['Affinity'] = round($iscore, 4);
            $data['A100'] = round($ia100, 4);
            $data['Predicted'] = 'Inverse Agonist';
            $data['DockScore'] = round(min($iscore-$ascore, 150) * $ittns * ($ia100 - 20) / 100, 4);
        }
        else
        {
            $data['Affinity'] = 0;
            $data['A100'] = round($ia100, 4);
            $data['Predicted'] = 'Non-agonist';
            $data['DockScore'] = 0.0;
        }

        echo "\nProtein: $protid\nLigand: $ligname";
        echo "\nResult: " . print_r($data, true) . "\n";
    }
    else if (isset($data["a_POSES"]) && !$data["a_POSES"])
    {
        $data['Affinity'] = 0;
        $data['A100'] = "(unknown)";
        $data['Predicted'] = 'Non-agonist';
        $data['DockScore'] = 0.0;

        echo "\nProtein: $protid\nLigand: $ligname";
        echo "\nResult: " . print_r($data, true) . "\n";
    }

    $data['CalculateDate'] = date('Y-m-d H:i:s');

    return $data;
}
