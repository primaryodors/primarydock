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

use Phpml\Classification\MLPClassifier;
use Phpml\ModelManager;
use Phpml\NeuralNetwork\ActivationFunction\PReLU;
use Phpml\NeuralNetwork\ActivationFunction\Sigmoid;

require_once("odorutils.php");

$training_set   = [];
$evaluation_set = [];
$prediction_set = [];

$clsss = ['Non-Agonist',
          'Agonist',
          'Inverse Agonist'
         ];

$json_file = "dock_results_soft.json";
$dock_data = [];
if (file_exists($json_file))
{
    $dock_data = json_decode(file_get_contents($json_file), true);
}

$metrics = [];
foreach ($dock_data as $prot => $ligdat)
{
    foreach ($ligdat as $lig => $data)
    {
        $odor = find_odorant($lig);
        $oid = $odor[$oid];
        
        if (isset($data['Actual']))
        {
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
            $metrics[$key] = $key;
        }
    }
}

$aifile = __DIR__."/corr_soft.ai";
if (file_exists($file)) 
{
    $mlp = $modelManager->restoreFromFile($file);
}
else
{ 
    $mlp = new MLPClassifier(count($metrics), [[100, new PReLU]], $clsss);
}



















