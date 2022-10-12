<?php
chdir(__DIR__);
require_once("../predict/protutils.php");
require_once("../predict/odorutils.php");

$rcpid = @$_REQUEST['r'];
if (!$rcpid)
{
    header("Location: receptors.php");
    exit;
}

$receptor = @$prots[$rcpid];
if (!$receptor)
{
    header("Location: receptors.php");
    exit;
}

$page_title = $rcpid;
$extra_js = ['js/tabs.js'];
$extra_css = ['assets/tabs.css'];

include("header.php");

?>
<div class="tab" style="display: inline-block; margin-top: 30px;">
    <button class="tabstatic" id="tabGene"><?php echo $rcpid; ?></button>
	<button class="tablinks" id="tabInfo" onclick="openTab(this, 'Info');">Info</button>
    <button class="tablinks default" id="tabLigands" onclick="openTab(this, 'Ligands');">Ligands</button>
</div>

<div id="Info" class="tabcontent">
    <?php 
    if (substr($rcpid,0,2) == 'OR')
    { 
        $substr = substr($rcpid,2);
        $lettre = preg_replace("/[^A-Z]/", "", $substr);
        $numero = explode($lettre, $substr);

        echo "<h3>Olfactory Receptor family {$numero[0]} subfamily $lettre member {$numero[1]}.</h3>";
    } 
    else if (substr($rcpid,0,1) == 'TAAR')
    { 
        $numero = substr($rcpid,-1);
        echo "<h3>Trace Amine Associated Receptor $numero.</h3>";
    }
    else if (substr($rcpid,0,1) == 'VN1R')
    { 
        $numero = substr($rcpid,-1);
        echo "<h3>Vomeronasal type 1 Receptor $numero.</h3>";
    }
    
    ?>
    <h3>Sequence:</h3>
    <pre id="protseq"><?php 
    $seq = $receptor['sequence']; 
    $sl = strlen($seq);

    $tmrcols = [1 => '#99c', '#09f', '#0c6', '#bc9', '#fd0', '#f60', '#f06'];

    $nxtmr = 1;
    $nums = $lets = "";
    $bold = 0;
    for ($i=0; $i<$sl; $i++)
    {
        if (!($i % 10)) $nums .= str_pad($i+10, 10, ' ', STR_PAD_LEFT).' ';
        if ($nxtmr > 7)
        {
            $lets .= substr($seq,$i,1);
            goto _tail;
        }

        if (($i+1) == $receptor['region']["TMR$nxtmr"]['start']) 
        {
            $col = $tmrcols[$nxtmr];
            $lets .= "<b style=\"color: $col;\">";
            $bold=1;
        }

        if (($i+1) == $receptor['bw']["$nxtmr.50"]) 
        $lets .= "<span style=\"background-color: #ddd; color: #000;\">".substr($seq,$i,1)."</span>";
        else
        $lets .= substr($seq,$i,1);

        if (($i+1) == $receptor['region']["TMR$nxtmr"]['end']) 
        {
            $lets .= "</b>";
            $nxtmr++;
            $bold=0;
        }

        _tail:
        if (($i % 10) == 9) $lets .= ' ';

        if (($i % 100) == 99) 
        {
            if ($bold) $lets .= "</b>";
            echo "$nums\n$lets\n\n";
            $nums = $lets = "";
            if ($bold) $lets .= "<b style=\"color: $col;\">";
        }
    }
    if ($lets) echo "$nums\n$lets\n\n";

    ?>
    </pre>

    <?php
    $rgntext = [];
    foreach ($receptor['region'] as $rgn => $r)
    {
        $st  = $r['start'];
        $end = $r['end'];
        
        if (substr($rgn,0,3) == 'TMR')
        {
            $tmr = intval(substr($rgn,-1));
            for ($i=$st; $i<=$end; $i=$j)
            { 
                $j = $i + 4;
                if ($j > $end) 
                { 
                    $rgntext[$rgn][] = str_pad(substr($seq, $i-1, $end-$i+1), 4, ' ', STR_PAD_RIGHT);
                    break;
                }
                else $rgntext[$rgn][] = substr($seq, $i-1, $j-$i);
            
                $i = $j;
                $j = $i + 3;
                if ($j > $end) 
                { 
                    $rgntext[$rgn][] = str_pad(substr($seq, $i-1, $end-$i+1), 3, ' ', STR_PAD_RIGHT);
                    break;
                }
                else $rgntext[$rgn][] = substr($seq, $i-1, $j-$i);
            }
            
            if (!($tmr & 1))
            { 
                foreach ($rgntext[$rgn] as $k => $v) $rgntext[$rgn][$k] = strrev($v);
                $rgntext[$rgn] = array_values(array_reverse($rgntext[$rgn]));
            }
        }
        else $rgntext[$rgn][] = substr($seq, $st-1, $end-$st+1);
    }

    // Binding site residues from https://doi.org/10.1110%2Fps.03296404
    $bsr = array_flip(
    [
        "2.53",
        "3.29", "3.32", "3.33", "3.36", "3.37", "3.40", "3.41",
        "4.53", "4.57", "4.60",
        "45.49", "45.52",
        "5.39", "5.42", "5.43", "5.46", "5.47",
        "6.48", "6.51",
        "7.38", "7.39", "7.42",
    ]);

    echo "<!-- ".print_r($bsr, true)." -->\n";

    ?>

    <h3>Transmembane Helices:</h3>
    <pre style="backgroumnd-color: #def; padding: 5px;"><?php 
    $mxrt = 0; 
    foreach ($rgntext as $rgn => $text) if (substr($rgn,0,3) == 'TMR') 
    { 
        $col = $tmrcols[intval(substr($rgn,-1))];
        echo "<span style=\"color: $col; font-weight: bold;\">".str_pad("  $rgn", 9, ' ', STR_PAD_RIGHT)."</span>"; 
        $len = count($text);
        if ($len > $mxrt) $mxrt = $len;
    }

    echo "\n";

    $ltl = [];
    $rc = 
    [
        1 =>  intval($receptor['region']['TMR1']['start']),
        2 =>  intval($receptor['region']['TMR2']['end']),
        3 =>  intval($receptor['region']['TMR3']['start']),
        4 =>  intval($receptor['region']['TMR4']['end']),
        5 =>  intval($receptor['region']['TMR5']['start']),
        6 =>  intval($receptor['region']['TMR6']['end']),
        7 =>  intval($receptor['region']['TMR7']['start']),
    ];

    echo "<!-- ".print_r($rc, true)." -->\n";

    for ($i=0; $i<$mxrt; $i++)
    {
        foreach ($rgntext as $rgn => $text) if (substr($rgn,0,3) == 'TMR') 
        {
            $tmr = intval(substr($rgn,-1));
            if (isset($text[$i]))
            {
                $tl = strlen($text[$i]);
                if ($tl==3 || @$ltl[$rgn] == 4) echo " ";
                
                for ($j=0; $j<$tl; $j++) 
                { 
                    $ch = substr($text[$i],$j,1);
                    // $col = aacolor($ch);
                    echo "<span ";
                    
                    $bw = $rc[$tmr] + 50 - $receptor['bw']["$tmr.50"];
                    
                    if ($ch != ' ') echo "title=\"$tmr.$bw {$aminos[$ch]}{$rc[$tmr]}\" ";
                    echo "style=\"";

                    if (isset($bsr["$tmr.$bw"]))
                    {
                        if ($ch == 'F' || $ch == 'W' || $ch == 'Y') echo "background-color: #c9f; ";
                        if ($ch == 'S' || $ch == 'T' || $ch == 'N' || $ch == 'C' || $ch == 'Q' || $ch == 'Y') echo "color: #066; ";
                        if ($ch == 'H' || $ch == 'K' || $ch == 'R') echo "background-color: #69f; color: #039; ";
                        if ($ch == 'D' || $ch == 'E') echo "background-color: #f99; color: #900; ";
                        if ($ch == 'S' || $ch == 'T' || $ch == 'N' || $ch == 'Q' ) echo "background-color: #6fd; ";
                        if ($ch == 'M' || $ch == 'A' || $ch == 'I' || $ch == 'L' || $ch == 'V' || $ch == 'P' || $ch == 'F' || $ch == 'W') echo "color: #234; ";
                        if ($ch == 'A' || $ch == 'I' || $ch == 'L' || $ch == 'V' || $ch == 'P' || $ch == 'G') echo "background-color: #9bd; ";
                        if ($ch == 'C' || $ch == 'M') echo "background-color: #ec0; ";
                        if ($ch == 'G') echo "color: #333; ";

                        echo "font-weight: bold; ";
                    }
                    echo "\">";
                    echo substr($text[$i],$j,1);
                    echo "</span>";
                    echo " ";
                    
                    if ($ch != ' ')
                    { 
                        if ($tmr & 1) $rc[$tmr]++;
                        else          $rc[$tmr]--;
                    }
                }
                if ($tl==3 || @$ltl[$rgn] == 4) echo " ";
                
                $ltl[$rgn] = $tl;
            }
            else echo "        ";
            
            echo " ";
        }
    
        echo "\n";
    }
    ?>
    </table>

</div>

<div id="Ligands" class="tabcontent">

<div class="scrollh" style="height: 780px;">
<table class="liglist">
    <tr>
        <th>Odorant</th>
        <th>EC<sub>50</sub></th>
        <th>Adjusted Top</th>
        <th>Aroma Notes</th>
    </tr>

<?php 

$pairs = all_empirical_pairs_for_receptor($rcpid);
// die("<pre>".print_r($pairs,true));

$refs = [];
foreach ($pairs as $oid => $pair)
{
    $odor = $odors[$oid];

    if (@$pair['ec50_ref'])
    {
        $refno_ec50 = array_search($pair['ec50_ref'], $refs);
        if (false === $refno_ec50)
        {
            $refno_ec50 = count($refs)+1;
            $refs[$refno_ec50] = $pair['ec50_ref'];
        }
    }
    if (@$pair['top_ref'])
    {
        $refno_top = array_search($pair['top_ref'], $refs);
        if (false === $refno_top)
        {
            $refno_top = count($refs)+1;
            $refs[$refno_top] = $pair['top_ref'];
        }
    }

    $pq = [];
    foreach ($odor['aroma'] as $refurl => $notes) $pq = array_merge($pq, $notes);
    $pq = array_unique($pq);

    $ufn = urlencode($odor['full_name']);
    echo "<tr>\n";
    echo "<td><a href=\"odorant.php?o=$oid\">{$odor['full_name']}</a></td>\n";

    echo "<td>" . $dispec50 = (@$pair['ec50'] ? ("{$pair['ec50']} <sup><a href=\"{$pair['ec50_ref']}\">$refno_ec50</a></sup>") : "-") . "</td>\n";
    echo "<td>" . $disptop = (@$pair['adjusted_curve_top'] ? (round(@$pair['adjusted_curve_top'], 4) . " <sup><a href=\"{$pair['top_ref']}\">$refno_top</a>") : "-") . "</sup></td>\n";

    echo "<td>" . implode(", ",$pq) . "</td>\n";
    echo "</tr>\n";
}
