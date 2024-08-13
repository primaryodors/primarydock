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

        if ($ascore > 0 && $aa100 > 0 && ($ascore > $iscore || $aavg > $iavg))
        {
            $data['Predicted'] = 'Agonist';
            $data['Affinity'] = $ascore;
            $data['A100'] = $aa100;
            $data['DockScore'] = min($ascore, 35) * $aa100 / 10;
        }
        else if ($iscore > 0)
        {
            $data['Affinity'] = $iscore;
            $data['A100'] = $ia100;
            $data['Predicted'] = 'Inverse Agonist';
            $data['DockScore'] = ($ia100 - 20) / 10;
        }
        else
        {
            $data['Affinity'] = 0;
            $data['A100'] = $ia100;
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

    return $data;
}
