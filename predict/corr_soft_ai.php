<?php

/******************************************************************************/
/*                                                                            */
/* To use this PHP script, you must have Composer installed, as well as the   */
/* php-ai/php-ml library.                                                     */
/*                                                                            */
/* sudo apt-get install composer                                              */
/* composer require php-ai/php-ml                                             */
/*                                                                            */
/******************************************************************************/

use Phpml\Classification\MLPClassifier;
use Phpml\ModelManager;
use Phpml\NeuralNetwork\ActivationFunction\PReLU;
use Phpml\NeuralNetwork\ActivationFunction\Sigmoid;


$json_file = "dock_results_soft.json";
$dock_data = [];
if (file_exists($json_file))
{
    $dock_data = json_decode(file_get_contents($json_file), true);
}

$training_set = [];
$evaluation_set = [];
$prediction_set = [];
