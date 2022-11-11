<?php
chdir(__DIR__);
require_once("../predict/protutils.php");
require_once("../predict/odorutils.php");
require_once("../predict/statistics.php");

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

// Binding site residues from https://doi.org/10.1110%2Fps.03296404
$bsr = array_flip(
    [
        "2.53",
        "3.29", "3.32", "3.33", "3.36", "3.37", "3.40", "3.41",
        "4.53", "4.57", "4.60",
        "45.49", "45.52",
        "5.39", "5.43", "5.46", "5.47",
        "6.48", "6.51",
        "7.38", "7.39", "7.42",
    ]);

$fam = family_from_protid($rcpid);
if ($fam == 'TAAR') $bsr['5.42'] = count($bsr);

// Copper binding sites for e.g. OR2T11
// http://pubs.acs.org/doi/abs/10.1021/jacs.6b06983
$cub =
[
	"2.39" => "MCHDENQR",
	"3.46" => "MCHDENQR",
	"4.37" => "MCHDENQR",
	"4.39" => "R",
	"4.42" => "CHMDENQ",
	"6.37" => "CHMDENQ",
	"6.40" => "HCMDENQ",
];

$page_title = $rcpid;
$extra_js = ['js/tabs.js'];
$extra_css = ['assets/tabs.css'];

$pairs = all_empirical_pairs_for_receptor($rcpid);

include("header.php");

?>
<script>
var viewer_loaded = false;
function load_viewer(obj)
{
    openTab(obj, 'Structure');
    if (!viewer_loaded)
    {
        window.setTimeout( function()
        {
            $('#viewer').on('load', function()
            {
                var embdd = $('#viewer')[0];
                $("[type=file]", embdd.contentDocument).hide();
                var filediv = $("#filediv", embdd.contentDocument)[0];

                filediv.innerText = "<?php echo @$rcpid; ?>";

                window.setTimeout( function()
                {
                    <?php
                    $tmo = 81;
                    foreach (array_keys($bsr) as $bw)
                    {
                        try
                        {
                            $resno = resno_from_bw($rcpid, $bw);
                            echo "window.setTimeout( function()\n";
                            echo "{\n";
                            echo "embdd.contentWindow.showSideChain($resno);\n";
                            echo "}, $tmo);\n";
                            $tmo += 53;
                        }
                        catch (Exception $ex)
                        {
                            ;
                        }
                    }
                    ?>
                }, 1234);
            });
            $('#viewer')[0].src = '<?php echo "viewer.php?url=pdb.php&prot=$rcpid"; ?>'; 
        }, 259); 
        viewer_loaded = true;
    }
}
</script>

<div class="tab" style="display: inline-block; margin-top: 30px;">
    <button class="tabstatic" id="tabGene"><?php echo $rcpid; ?></button>
	<button class="tablinks <?php if (!count($pairs)) echo "default"; ?>" id="tabInfo" onclick="openTab(this, 'Info');">Info</button>
    <?php if (count($pairs)) { ?>
    <button class="tablinks default" id="tabLigands" onclick="openTab(this, 'Ligands');">Ligands</button>
    <?php } ?>
    <button class="tablinks" id="tabRefs" onclick="openTab(this, 'Refs');">References</button>
    <button class="tablinks" id="tabStructure" onclick="load_viewer(this);">3D Structure</button>
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

        $between = ($nxtmr-1)."$nxtmr.50";
        if (($i+1) == resno_from_bw($rcpid, "$nxtmr.50")
            ||
            ( isset($prots[$rcpid]["bw"][$between]) && ($i+1) == resno_from_bw($rcpid, $between) )
            )
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

    $tmrstartoff = [];
    for ($i=1; $i<=7; $i++)
    {
        $x = $y = 0.0;
        foreach (array_keys($bsr) as $l => $bw)
        {
            if (substr($bw, 0, 2) == "$i.")
            {
                $resno = resno_from_bw($rcpid, $bw);
                $offset = $resno - intval($receptor['region']["TMR$i"]['start']);
                $angle = (pi()*4 / 7 * $offset) % (pi()*2);
                $x += sin($angle);
                $y += cos($angle);
            }
        }
        $angle = find_angle($x, $y);
        $tmrstartoff[$i] = intval(round($angle * 7 / (pi()*4)));
    }

    $tmrstartoff[6]++;

    $rgntext = [];
    foreach ($receptor['region'] as $rgn => $r)
    {
        $st  = $r['start'];
        $end = $r['end'];
        
        if (substr($rgn,0,3) == 'TMR')
        {
            $tmr = intval(substr($rgn,-1));
            $tso = -$tmrstartoff[$tmr];
            for ($i=$st; $i<=$end; $i=$j)
            { 
                $j = $i + 4;
                if ($j > $end) 
                { 
                    $rgntext[$rgn][] = str_pad(substr($seq, $i-1+$tso, $end-$i+1), 4, ' ', STR_PAD_RIGHT);
                    break;
                }
                else $rgntext[$rgn][] = substr($seq, $i-1+$tso, $j-$i);

                $i = $j;
                $j = $i + 3;
                if ($j > $end) 
                { 
                    $rgntext[$rgn][] = str_pad(substr($seq, $i-1+$tso, $end-$i+1), 3, ' ', STR_PAD_RIGHT);
                    break;
                }
                else $rgntext[$rgn][] = substr($seq, $i-1+$tso, $j-$i);
            }
            
            if (!($tmr & 1))
            { 
                foreach ($rgntext[$rgn] as $k => $v) $rgntext[$rgn][$k] = strrev($v);
                $rgntext[$rgn] = array_values(array_reverse($rgntext[$rgn]));
            }
        }
        else $rgntext[$rgn][] = substr($seq, $st-1, $end-$st+1);
    }

$matched = 0;
$lcub = [];
foreach ($cub as $bw => $allowed)
{
	$res = "<span title=\"$bw\">";
	$resno = resno_from_bw($rcpid, $bw);
	$c = substr($seq, $resno-1, 1);
	if (false !== strpos($allowed, $c))
	{
		$matched++;
		$res .= "<b>$c</b>";
	}
	else $res .= $c;

	$res .= "</span>";

	$lcub[] = $res;
}

echo "<p>$matched residues match deep copper binding site: " . implode(" ",$lcub) . ". ";
if ($matched == count($cub)) echo "This receptor is likely to respond to thiols that can access deep inside the ligand binding pocket.";
echo "</p>";


    ?>

    <h3>Transmembane Helices:</h3>
    <pre style="padding: 5px;"><?php 
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

    for ($i=1; $i<=7; $i++) $rc[$i] -= $tmrstartoff[$i];

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
                    
                    $bw = bw_from_resno($rcpid, $rc[$tmr]); // $rc[$tmr] + 50 - $receptor['bw']["$tmr.50"];
                    
                    if ($ch != ' ') echo "title=\"$bw {$aminos[$ch]}{$rc[$tmr]}\" ";
                    echo "style=\"";

                    if (isset($bsr["$bw"]))
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

    <?php
    if ($uniprot = @$receptor['uniprot_id'])
    {
        $links[] = "<a href=\"https://www.uniprot.org/uniprot/$uniprot\" target=\"_top\">UniProt</a>";
        if (preg_match("/^(OR[0-9]{1,2}[A-Z]{1,2}[0-9]{1,2})|(TAAR[0-9])|(VN1R[1-5])$/", $rcpid))
            $links[] = "<a href=\"https://zhanggroup.org/GPCR-EXP/pdb/hgmod/$uniprot/{$uniprot}_results/\" target=\"_top\">I-TASSER</a>";
        $links[] = "<a href=\"https://alphafold.ebi.ac.uk/entry/$uniprot\" target=\"_top\">AlphaFold</a>";
    }

    if (preg_match("/^OR[0-9]{1,2}[A-Z]{1,2}[0-9]{1,2}$/", $rcpid))
    {
        $links[] = "<a href=\"https://genome.weizmann.ac.il/horde/card/index/symbol:$rcpid\" target=\"_top\">HORDE</a>";
    }

    if (1)
    {
        $links[] = "<a href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=$rcpid\" target=\"_top\">GeneCards</a>";
    }

    if (count($links)) echo implode(" | ", $links);
    ?>

</div>

<div id="Ligands" class="tabcontent">

<div class="scrollh" style="height: 780px;">
<table class="liglist">
    <tr>
        <th>Odorant</th>
        <th>EC<sub>50</sub></th>
        <th>Adjusted Top</th>
        <th>Antagonist?</th>
        <th>Aroma Notes</th>
    </tr>

<?php 

$lrefs = [];
foreach ($pairs as $oid => $pair)
{
    $odor = $odors[$oid];

    if (@$pair['ec50_ref'])
    {
        $refno_ec50 = array_search($pair['ec50_ref'], $lrefs);
        if (false === $refno_ec50)
        {
            $refno_ec50 = count($lrefs)+1;
            $lrefs[$refno_ec50] = $pair['ec50_ref'];
        }
    }
    if (@$pair['top_ref'])
    {
        $refno_top = array_search($pair['top_ref'], $lrefs);
        if (false === $refno_top)
        {
            $refno_top = count($lrefs)+1;
            $lrefs[$refno_top] = $pair['top_ref'];
        }
    }

    $pq = [];
    foreach ($odor['aroma'] as $refurl => $notes) $pq = array_merge($pq, $notes);
    $pq = array_unique($pq);

    $ufn = urlencode($odor['full_name']);
    echo "<tr>\n";
    echo "<td><a href=\"odorant.php?o=$oid\" style=\"white-space: nowrap;\">{$odor['full_name']}</a></td>\n";

    echo "<td>" . 
        ($dispec50 = (@$pair['ec50']
            ? ("{$pair['ec50']} <sup><a href=\"#\" onclick=\"openTab($('#tabRefs')[0], 'Refs');\">$refno_ec50</a></sup>")
            : "-")
        ) . "</td>\n";
    echo "<td>" . 
        ($disptop = (@$pair['adjusted_curve_top']
            ? (round(@$pair['adjusted_curve_top'], 4) . " <sup><a href=\"#\" onclick=\"openTab($('#tabRefs')[0], 'Refs');\">$refno_top</a>")
            : "-")
        ) . "</sup></td>\n";

    if (@$pair['antagonist']) echo "<td>Y</td>";
    else echo "<td>&nbsp;</td>";

    echo "<td style=\"white-space: nowrap;\">" . implode(", ",$pq) . "</td>\n";
    echo "</tr>\n";
}

?>
</table>
</div>
</div>

<div id="Refs" class="tabcontent">
<?php
foreach ($lrefs as $idx => $refurl)
{
    echo "<a href=\"$refurl\"><p>\n";
    echo "$idx.) ";
    echo $refs[$refurl]['citation'];
    echo "</p></a>\n";
}
?>
</div>

<div id="Structure" class="tabcontent">
    <iframe id="viewer"></iframe>
</div>
