<?php

function make_prediction($data)
{
    global $protid, $ligname, $pose;

    if (isset($data["a_Pose1"]) || isset($data["a_BindingEnergy"]))
    {
        $ascore = max(0, -floatval(@$data['a_Pose1']));
        $iscore = max(0, -floatval(@$data['i_Pose1']));

        $aa100 = floatval(@$data['a_A100']);
        $ia100 = floatval(@$data['i_A100']);

        $dock_score = $ascore - $iscore;
        if ($aa100 && $ia100) $dock_score *= max(0, min(1, ($aa100 - $ia100) / 50));

        if ($ascore > 0 && $dock_score > 0)
        {
            $data['Predicted'] = 'Agonist';
            $data['DockScore'] = $dock_score;
        }
        else if ($iscore > 0 && $dock_score < 0)
        {
            $data['Predicted'] = 'Inverse Agonist';
            $data['DockScore'] = $dock_score;
        }
        else
        {
            $data['Predicted'] = 'Non-agonist';
            $data['DockScore'] = 0.0;
        }

        echo "\nProtein: $protid\nLigand: $ligname";
        echo "\nResult: " . print_r($data, true) . "\n";
    }
    else if (isset($data["a_POSES"]) && !$data["a_POSES"])
    {
        $data['Predicted'] = 'Non-agonist';
        $data['DockScore'] = 0.0;

        echo "\nProtein: $protid\nLigand: $ligname";
        echo "\nResult: " . print_r($data, true) . "\n";
    }

    return $data;
}
