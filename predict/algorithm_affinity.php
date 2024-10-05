<?php

// https://www.masterorganicchemistry.com/2023/08/02/equilibrium-constant-delta-g-calculations-organic-chemistry/
define("R", 0.00831446261815324);
define("body_temperature", 310.2);

function DeltaG($K)
{
    return (-R * body_temperature)*log($K);
}

function K($DeltaG)
{
    return exp($DeltaG / (-R * body_temperature));
}

function equilibrium($kJmol1, $kJmol2)
{
    $DeltaG = $kJmol2 - $kJmol1;
    $K = K($DeltaG);
    return $K/($K+1);
}

function make_prediction($data)
{
    global $protid, $ligname, $pose;

    if (isset($data["a_Pose1"]) || isset($data["a_BindingEnergy"]))
    {
        $baseline = 0.05;                       // this is a guess.
        $baseDeltaG = DeltaG($baseline);

        $aenergy = floatval(@$data['a_Pose1']);
        $ienergy = floatval(@$data['i_Pose1']);
        $afound = floatval(@$data['a_POSES']); // / intval($pose);
        $ifound = floatval(@$data['i_POSES']); // / intval($pose);

        $aeq = equilibrium(0, $aenergy);                    // how much of the active-state population will be ligand-bound
        $ieq = equilibrium(0, $ienergy);                    // how much of the inactive-state population will be ligand-bound
        $bound_acv = equilibrium($ienergy, $aenergy-$baseDeltaG);       // how much of the ligand-bound population will be active
        $activation = $aeq*$bound_acv + (1.0-$baseline)*(1.0-$aeq) - $baseline;

        $aa100 = floatval(@$data['a_A100']);
        $ia100 = floatval(@$data['i_A100']);
        $attns = floatval(@$data['a_occlusion'] ?: 1);
        $ittns = floatval(@$data['i_occlusion'] ?: 1);

        $data['DockScore'] = 0;
        if ($activation >= $baseline) $data['DockScore'] = round(-$aenergy * $activation * min($aa100, 20)/20 * $attns * $afound, 4);

        if ($activation >= $baseline && $data['DockScore'] > 0)
        {
            $data['Predicted'] = 'Agonist';
            $data['Affinity'] = round(-$aenergy, 4);
            $data['A100'] = round($aa100, 4);
            $data['DockScore'] = round(-$aenergy * $activation * min($aa100, 20)/20 * $attns * $afound, 4);
        }
        else if ($ienergy > 0)
        {
            $data['Affinity'] = round(-$ienergy, 4);
            $data['A100'] = round($ia100, 4);
            $data['Predicted'] = 'Inverse Agonist';
            $data['DockScore'] = round(-$ienergy * ($baseline-$activation) * $ia100/20 * $ittns * $ifound, 4);
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
