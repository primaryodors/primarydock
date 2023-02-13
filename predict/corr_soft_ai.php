<?php

/******************************************************************************/
/*                                                                            */
/* To use this PHP script, you must have Composer installed, as well as the   */
/* php-ai/php-ml library:                                                     */
/*                                                                            */
/* sudo apt-get install composer                                              */
/* composer require php-ai/php-ml                                             */
/*                                                                            */
/******************************************************************************/

chdir(__DIR__);
require_once '../vendor/autoload.php';
chdir(__DIR__);

use Phpml\Classification\MLPClassifier;
use Phpml\ModelManager;
use Phpml\NeuralNetwork\ActivationFunction\PReLU;
use Phpml\NeuralNetwork\ActivationFunction\Sigmoid;

require_once("odorutils.php");

echo "Process begin ".date("Y-m-d H:i:s").".\n";

$training_set   = [];
$evaluation_set = [];
$prediction_set = [];

$labels = ['Non-Agonist',
          'Agonist',
          'Inverse Agonist'
         ];

$json_file = "dock_results_soft.json";
$dock_data = [];
if (file_exists($json_file))
{
    $dock_data = json_decode(file_get_contents($json_file), true);
}

$cp = 0;
$metrics = [];
$dock_data_by_oid = [];
foreach ($dock_data as $prot => $ligdat)
{
    foreach ($ligdat as $lig => $data)
    {
        $odor = find_odorant($lig);
        $oid = $odor['oid'];
        $dock_data_by_oid[$prot][$oid] = $data;

        if (isset($data['Actual']))
        {
            $cp++;
            if ((intval(substr($oid, 12, 1)) & 1)
                &&
                (intval(substr($oid, 24, 1)) & 1)
               )
            {
                $arrvar = "training_set";
            }
            else $arrvar = "evaluation_set";

            $$arrvar["$prot.$oid"] = $data['Actual'];
        }

        foreach ($data as $key => $value)
        {
            if (substr($key, 0, 3) == "TMR")
            {
                $metrics[$key] = $key;
            }
        }
    }
}

// print_r($training_set); exit;
// print_r($evaluation_set); exit;

$ctd = count($training_set);
$minctd = 20;
if ($ctd < $minctd) die("Not enough training data: only $ctd pairs, require at least $minctd.\n");

$cm = count($metrics);

echo "Total $cm metrics for $cp pairs.\n";

$aifile = __DIR__."/corr_soft.ai";
$modelManager = new ModelManager();
if (file_exists($aifile)) 
{
    $mlp = $modelManager->restoreFromFile($aifile);
}
else
{ 
    $mlp = new MLPClassifier($cm, [[100, new PReLU]], $labels);
}

// Train the AI.
$iters = 5;                // How many training iterations.

for ($i=1; $i<=$iters; $i++)
{
    set_time_limit(1800);
    echo "Training iteration $i...";
    
    foreach ($training_set as $pair => $actual)
    {
        $pessia = explode(".", $pair);
        $prot = $pessia[0];
        $oid = $pessia[1];
        $sample = [];
        foreach ($metrics as $m)
        {
            if (isset($dock_data_by_oid[$prot][$oid][$m])) $sample[$m] = floatval($dock_data_by_oid[$prot][$oid][$m]);
            else $sample[$m] = 0;
        }
        $mlp->partialTrain( [ array_values($sample) ], [ $actual ]);
        echo '.';
    }
    
    echo "\n";
}


$modelManager->saveToFile($mlp, $aifile);


$correct_pred = 0;
foreach ($evaluation_set as $pair => $actual)
{
    $pessia = explode(".", $pair);
    $prot = $pessia[0];
    $oid = $pessia[1];
    $sample = [];
    foreach ($metrics as $m)
    {
        if (isset($dock_data_by_oid[$prot][$oid][$m])) $sample[$m] = floatval($dock_data_by_oid[$prot][$oid][$m]);
        else $sample[$m] = 0;
    }

    $pred = $mlp->predict([array_values($sample)]);
    // print_r($pred); exit;
    
    if ($pred[0] == $actual) $correct_pred++;
}

$percent = round(100.0 * $correct_pred / count($evaluation_set), 2);

echo "Got $percent% of evaulation set right.\n";

















