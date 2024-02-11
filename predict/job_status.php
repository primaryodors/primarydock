<?php

function set_color($r, $g, $b)
{
    $lr = intval(max(0,min(5,$r*6)));
    $lg = intval(max(0,min(5,$g*6)));
    $lb = intval(max(0,min(5,$b*6)));

    $ccode = intval(16 + $lb + 6*$lg + 36*$lr);

    echo "\x1b[38;5;{$ccode}m";
}

function clear_color()
{
    echo "\x1b[39m";
}

function colorful_text($text, $rgb)
{
    $r = intval( $rgb / 0x10000);
    $g = intval(($rgb & 0xff00) / 0x100);
    $b = intval( $rgb & 0xff);
    $m = 1.0 / 256;
    set_color($m*$r,$m*$g,$m*$b);
    echo $text;
    clear_color();
}

$queue = [];

chdir(__DIR__);
chdir("..");
$dock_results = json_decode(file_get_contents("predict/dock_results.json"), true);
$version = max(filemtime("bin/primarydock"),
    filemtime("predict/methods_common.php"),
    filemtime("predict/method_directmdl.php"),
    filemtime("predict/method_fygactive.php")
    );
$lines = explode("\n", @file_get_contents("jobq"));

foreach ($lines as $ln)
{
    $ln = trim($ln);
    if (!$ln) continue;
    if (substr($ln, 0, 1) == '#') continue;
    $pieces = explode(" ", $ln);

    if ($pieces[0] == "MAX") $max_concurrent = intval($pieces[1]);
    else if ($pieces[0] == "PRDT") $queue[] = [$pieces[1], $pieces[2]];
}
if (!count($queue)) die("No jobs.\n");

foreach ($queue as $q)
{
    $prot = $q[0];
    $lig = $q[1];

    $text = "unknown";
    $color = 0x996699;

    if (isset($dock_results[$prot][$lig]))
    {
        $dr = $dock_results[$prot][$lig];
        if (intval($dr['version']) >= $version)
        {
            if (isset($dr['Actual']))
            {
                if (($dr['Actual'] == "Agonist") == ($dr['Predicted'] == "Agonist"))
                {
                    $text = "good!";
                    $color = 0x00ff00;
                }
                else
                {
                    $text = "FAIL";
                    $color = 0xff0000;
                }
            }
            else
            {
                $text = @$dr['Predicted'] ?: "incomplete";
                $color = 0x0066ff;
            }
        }
        else
        {
            $text = "outdated";
            $color = 0x666666;
        }
    }
    else
    {
        $text = "incomplete";
        $color = 0x999999;
    }

    colorful_text(mb_chr(0x2588).mb_chr(0x2588).mb_chr(0x2588).mb_chr(0x2588).mb_chr(0x2588), $color);
    echo " $prot - $lig: $text";
    echo "\n";
}
