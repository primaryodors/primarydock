<?php
chdir(__DIR__);
require_once("../data/protutils.php");
require_once("../data/odorutils.php");
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

$predictions = [];
$pred_shown = [];
$predname = [];
$predate = [];
$preddock = [];

function predcmp($a, $b)
{
    if ($a < $b) return 1;
    else if ($b < $a) return -1;
    else return 0;
}

chdir(__DIR__);
$dock_results = json_decode(file_get_contents("../predict/dock_results.json"), true);
if (isset($dock_results[$rcpid]))
{
    foreach ($dock_results[$rcpid] as $ligname => $dock)
    {
        $odor = find_odorant($ligname);
        $oid = $odor['oid'];
        $predname[$oid] = $ligname;
        $predate[$oid] = date("Y-m-d H:i:s", $dock['version']);
        $preddock[$oid] = $dock['docker'];
        if (isset($dock['DockScore'])) $predictions[$oid] = floatval($dock['DockScore']);
        else if (isset($dock['a_Pose1']) && isset($dock['i_Pose1']))
            $predictions[$oid] = (floatval($dock['i_Pose1']) - floatval($dock['a_Pose1'])) / 2;
    }
}

uasort($predictions, "predcmp");

$cmd = "ps -ef | grep -E ':[0-9][0-9] (bin/primarydock|bin/pepteditor|obabel)' | grep -v grep";
$result = [];
exec($cmd, $result);

$ligands_inprogress = [];
foreach ($result as $line)
{
    $pos = strpos($line, ".$rcpid.");
    if (false == $pos) $pos = strpos($line, ".{$rcpid}_");
    if (false == $pos) continue;
    $pos2 = strpos($line, ".", $pos+1);
    $ligname = substr($line, $pos2);
    $ligname = preg_replace("/[.]config/", "", $ligname);
    $odor = find_odorant($ligname);
    $ligands_inprogress[$odor['oid']] = true;
}

// Copper binding sites for e.g. OR2T11
// http://pubs.acs.org/doi/abs/10.1021/jacs.6b06983
$cub =
[
	"2.39" => "MCHDENQR",
	"3.46" => "MCHDENQR",
	"3.50" => "R",
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
<style>
#skeletal
{
    margin: 10px;
    width: 300px;
    height: 300px;
    position: absolute;
    box-shadow: 25px 25px 35px rgba(0,0,0,0.5);
    z-index: 10000;
    background: #234;
}

#dlmenu
{
    position: absolute;
    box-shadow: 25px 25px 35px rgba(0,0,0,0.5);
    z-index: 10000;
    background: #234;
    padding: 20px;
    font-size: small;
}

.ctxmenu
{
    margin: 0px;
    padding-inline-start: 0px;
}

.ctxmenu li
{
    display: block;
    margin-bottom: 2px;
}

.ctxmenu a
{
    display: inline-block;
    border-top: 2px solid #68a;
    border-left: 2px solid #468;
    background-color: #246;
    border-right: 2px solid #123;
    border-bottom: 2px solid #012;
    padding: 0px 5px;
    font-size: small;
}

#txtview, #htmview
{
    position: absolute;
    z-index: 10000;
    width: 80%;
    height: 80%;
    margin: 5% 10%;
    border-top: 2px solid #68a;
    border-left: 2px solid #468;
    background-color: #246;
    border-right: 2px solid #123;
    border-bottom: 2px solid #012;
    box-shadow: 0px 0px 50px #000, 0px 0px 100px #000, 0px 0px 200px #000;
}

#htmview
{
    padding: 13px;
    padding-top: 4px;
}

#txttext, #htmhtml
{
    width: 100%;
    height: 100%;
}

.closebtn
{
    display: inline-block;
    border-top: 1px solid #ff0;
    border-left: 1px solid #f90;
    background-color: #f00;
    border-right: 1px solid #900;
    border-bottom: 1px solid #300;
    color: #fff!important;
    text-decoration: none!important;
    cursor: pointer;
    padding: 0px 5px;
}
</style>
<script>

var embdd;
function show_bsr_chains()
{
<?php
    $tmo = 81;
    foreach (array_keys($bsr) as $bw)
    {
        try
        {
            $resno = resno_from_bw($rcpid, $bw);
            echo "\twindow.setTimeout( function() { embdd.contentWindow.showSideChain($resno); }, $tmo);\n";
            $tmo += 53;
        }
        catch (Exception $ex)
        {
            ;
        }
    }
?>
}

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
                embdd = $('#viewer')[0];
                $("[type=file]", embdd.contentDocument).hide();
                var filediv = $("#filediv", embdd.contentDocument)[0];

                <?php chdir(__DIR__);
                if (file_exists("../pdbs/$fam/$rcpid.active.pdb"))
                {
                    echo "filediv.innerHTML = \"$rcpid";
                    echo " &#xb7; <small><a id=\\\"mdl3di\\\" class=\\\"mdlsel\\\" href=\\\"#\\\" onclick=\\\"load_remote_pdb('pdb.php?prot=$rcpid'); $('#mdl3da').removeClass('mdlsel'); this.className = 'mdlsel'; window.parent.setTimeout( window.parent.show_bsr_chains, 1234);\\\">inactive</a></small>";
                    echo " &#xb7; <small><a id=\\\"mdl3da\\\" href=\\\"#\\\" onclick=\\\"load_remote_pdb('pdb.php?prot=$rcpid&mod=a'); $('#mdl3di').removeClass('mdlsel'); this.className = 'mdlsel'; window.parent.setTimeout( window.parent.show_bsr_chains, 1234);\\\">active</a></small>";
                    echo "\";";
                }
                else
                { ?>
                filediv.innerText = "<?php echo @$rcpid; ?>";
                <?php } ?>
                
                window.setTimeout( show_bsr_chains, 1234);
            });
            $('#viewer')[0].src = '<?php echo "viewer.php?url=pdb.php&prot=$rcpid"; ?>'; 
        }, 259); 
        viewer_loaded = true;
    }
}

function showSkeletal(e, smiles)
{
    var skeletal = $("#skeletal")[0];
    var top = e.pageY;
    if (top > document.body.scrollHeight - 300) top -= 300;
    skeletal.style.left = `${e.pageX}px`;
    skeletal.style.top = `${top}px`;
    skeletal.innerHTML = svg_from_smiles(smiles, 300, 300);

    $(skeletal).show();
}

function show_dlmenu(e, prot, lig, v, d)
{
    var dlmenu = $("#dlmenu")[0];

    $("#dl_mnu_prot")[0].innerText = prot;
    $("#dl_mnu_lig")[0].innerText = decodeURIComponent(lig);
    $("#dl_cd_v")[0].innerText = v;
    if (d == "pd") $("#dl_dock")[0].innerText = "PrimaryDock";
    else if (d == "vina") $("#dl_dock")[0].innerText = "AutoDock Vina";
    else $("#dl_dock")[0].innerText = "(unknown)";
    $("#dl_acv_mdl")[0].setAttribute("href", "download.php?obj=model&prot="+prot+"&odor="+lig+"&mode=active");
    $("#vw_acv_mdl_3d")[0].setAttribute("href", "viewer.php?view=pred&prot="+prot+"&odor="+lig+"&mode=active");
    $("#dl_iacv_mdl")[0].setAttribute("href", "download.php?obj=model&prot="+prot+"&odor="+lig+"&mode=inactive");
    $("#vw_iacv_mdl_3d")[0].setAttribute("href", "viewer.php?view=pred&prot="+prot+"&odor="+lig+"&mode=inactive");
    $("#dl_acv_dc")[0].setAttribute("href", "download.php?obj=dock&prot="+prot+"&odor="+lig+"&mode=active");
    $("#vw_acv_dc_3d")[0].setAttribute("href", "viewer.php?view=dock&prot="+prot+"&odor="+lig+"&mode=active");
    $("#dl_iacv_dc")[0].setAttribute("href", "download.php?obj=dock&prot="+prot+"&odor="+lig+"&mode=inactive");
    $("#vw_iacv_dc_3d")[0].setAttribute("href", "viewer.php?view=dock&prot="+prot+"&odor="+lig+"&mode=inactive");
    $("#dl_json")[0].setAttribute("href", "download.php?obj=json&prot="+prot+"&odor="+lig);

    var x = e.pageX + 5, y = e.pageY;
    if (y > window.innerHeight - 400) y = window.innerHeight - 400;

    dlmenu.style.left = `${x}px`;
    dlmenu.style.top = `${y}px`;

    $(dlmenu).show();
}

function view_file(url)
{
    var txt = $('#txttext')[0];
    var vw = $('#txtview');
    txt.value = "Loading...";

	$.ajax(
	{
		url: url,
		cache: false,
		success: function(result)
		{
            txt.value = result;
            vw.show();
        }
    });
    vw.show();
}

function view_html_file(url)
{
    var htm = $('#htmhtml')[0];
    var vw = $('#htmview');
    htm.innerHTML = "Loading...";

	$.ajax(
	{
		url: url,
		cache: false,
		success: function(result)
		{
            htm.innerHTML = result;
            vw.show();
        }
    });
    vw.show();
}
</script>

<div class="tab" style="display: inline-block; margin-top: 30px;">
    <button class="tabstatic" id="tabGene"><b><?php echo $rcpid; ?></b>:


<?php

if (@$receptor["btree"])
{
    foreach ($treenodes as $nodeid => $nodename)
    {
        $nodeid = str_replace('*', '', $nodeid);
        $len = strlen($nodeid);
        if (substr($receptor["btree"], 0, $len) === $nodeid)
        {
            echo "<a href=\"receptors.php?tree=1&bt=$nodeid\">$nodename</a>, ";
        }
    }
    echo $rcpid;
}
else
{
    if (substr($fam, 0, 2) == "OR")
    {
        echo "Olfactory receptor ";
        $fmn = intval(preg_replace("/[^0-9]/", "", $fam));
        echo "<a href=\"receptors.php?f=sOR$fmn\">family $fmn</a>, ";
        $sub = preg_replace("/[^A-Z]/", "", substr($rcpid, strlen($fam)) );
        echo "<a href=\"receptors.php?f=sOR$fmn$sub\">subfamily $sub</a>, ";
        $mbr = intval(preg_replace("/[^0-9]/", "", substr($rcpid, strlen($fam)) ));
        echo "member $mbr";
    }
    else
    {
        switch ($fam)
        {
            case 'TAAR':
            echo "Trace amine-associated receptor ";
            $mbr = intval(preg_replace("/[^0-9]/", "", substr($rcpid, strlen($fam)) ));
            echo "$mbr";
            break;

            case 'VN1R':
            echo "Vomeronasal type 1 receptor ";
            $mbr = intval(preg_replace("/[^0-9]/", "", substr($rcpid, strlen($fam)) ));
            echo "$mbr";
            break;

            case 'MS4A':
            echo "Membrane-spanning 4A receptor ";
            $mbr = intval(preg_replace("/[^0-9]/", "", substr($rcpid, strlen($fam)) ));
            echo "$mbr";
            break;

            default:
            ;
        }
    }
}

?>
    </button>
	<button class="tablinks <?php if (!count($pairs)) echo "default"; ?>" id="tabInfo" onclick="openTab(this, 'Info');">Info</button>
    <?php if (count($pairs) || count($predictions)) { ?>
    <button class="tablinks default" id="tabLigands" onclick="openTab(this, 'Ligands');">Ligands</button>
    <?php } ?>
    <button class="tablinks" id="tabRefs" onclick="openTab(this, 'Refs');">References</button>
    <?php if (count($pairs)) { ?>
    <button class="tablinks <?php if (@$_REQUEST['cmp']) echo "default"; ?>" id="tabComparison" onclick="openTab(this, 'Comparison');">Comparison</button>
    <?php } ?>
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
        $aaletter = substr($seq,$i,1);

        if (@$_REQUEST["abc"]) $lc = chr(65+strpos("ARNDCEQGHILKMFPSTWYV", $aaletter));
        else $lc = $aaletter;

        if (!($i % 10)) $nums .= str_pad($i+10, 10, ' ', STR_PAD_LEFT).' ';
        if ($nxtmr > 7)
        {
            $lets .= $lc;
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
        {
            $lets .= "<span style=\"background-color: #ddd; color: #000;\">$lc</span>";
        }
        else
        {
            $lets .= $lc;
        }

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

/*echo "<p>$matched residues match deep copper binding site: " . implode(" ",$lcub) . ". ";
if ($matched == count($cub)) echo "This receptor is likely to respond to thiols that can access deep inside the ligand binding pocket.";
echo "</p>";*/


    ?>

    <h3>Transmembrane Helices:</h3>
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
                $tl1 = strlen(trim($text[$i]));
                if ($tl==3 || @$ltl[$rgn] == 4) echo " ";
                
                if ($tmr & 1)
                {
                    $text[$i] = strrev($text[$i]);
                    $rc[$tmr] += $tl1 - 1;
                }

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
                    if (@$_REQUEST["abc"]) echo chr(65+strpos("ARNDCEQGHILKMFPSTWYV", substr($text[$i],$j,1)));
                    else echo substr($text[$i],$j,1);
                    echo "</span>";
                    echo " ";
                    
                    if ($ch != ' ')
                    { 
                        /*if ($tmr & 1) $rc[$tmr]++;
                        else*/          $rc[$tmr]--;
                    }
                }
                if ($tl==3 || @$ltl[$rgn] == 4) echo " ";
                
                $ltl[$rgn] = $tl;
            }
            else echo "        ";
            
            echo " ";
        
            if ($tmr & 1) $rc[$tmr] += $tl1+1;
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

    if (isset($receptor['ordb']))
    {
        if (is_array($receptor['ordb']))
        {
            $i=1;
            $str="";
            foreach ($receptor['ordb'] as $ordb)
            {
                $ordb = intval($ordb);
                if ($i==1) $str = "<a href=\"https://ordb.biotech.ttu.edu/ORDB/data/$ordb\" target=\"_top\">ORDB</a>";
                else $str .= " <a href=\"https://ordb.biotech.ttu.edu/ORDB/data/$ordb\" target=\"_top\">$i</a>";
                $i++;
            }
            $links[] = $str;
        }
        else
        {
            $ordb = intval($receptor['ordb']);
            $links[] = "<a href=\"https://ordb.biotech.ttu.edu/ORDB/data/$ordb\" target=\"_top\">ORDB</a>";
        }
    }

    if (count($links)) echo implode(" | ", $links);
    ?>

</div>

<div id="Ligands" class="tabcontent">

<div class="box">
<div class="row content scrollh">

<?php if (count($predictions)) { ?>
<div style="background-color: #fc9; color: 000; font-size: 0.7em; padding: 2px;">
This page features beta versions of predictions, including links to viewable 3D models. 
The numbers and models are not yet guaranteed to be accurate,
and we welcome collaboration on our codebase.
The repository is accepting pull requests at
<a style="color: #00c;" href="https://github.com/primaryodors/primarydock">https://github.com/primaryodors/primarydock</a>.
</div>
<?php } ?>

<table class="liglist">
    <tr>
        <th width="20%">Odorant</th>
        <th width="10%">EC<sub>50</sub></th>
        <th width="10%">Adjusted Top</th>
        <th width="10%">Antagonist?</th>
        <?php
        if (count($predictions)) echo "<th width=\"1%\">Predicted</th>";
        ?>
        <th width="50%">Aroma Notes</th>
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
    echo "<td><a href=\"odorant.php?o=$oid\" style=\"white-space: nowrap;\"";

    $smilesesc = str_replace("\\", "\\\\", $odor['smiles']);
    echo " onmouseenter=\"showSkeletal(event, '$smilesesc');\"";
    echo " onmouseout=\"$('#skeletal').hide();\"";
    echo ">{$odor['full_name']}</a>";
    echo "</td>\n";

    echo "<td>" . 
        ($dispec50 = (@$pair['ec50']
            ? ("{$pair['ec50']} <sup><a href=\"#\" onclick=\"openTab($('#tabRefs')[0], 'Refs');\">$refno_ec50</a></sup>")
            : "-")
        ) . "</td>\n";
    echo "<td>" . 
        ($disptop = ((@$pair['adjusted_curve_top'] || !@$pair['ec50'])
            ? (round(@$pair['adjusted_curve_top'], 4) . " <sup><a href=\"#\" onclick=\"openTab($('#tabRefs')[0], 'Refs');\">$refno_top</a>")
            : "-")
        ) . "</sup></td>\n";

    if (@$pair['antagonist']) echo "<td>Y</td>";
    else echo "<td>&nbsp;</td>";

    if (count($predictions))
    {
        if (isset($predictions[$oid]))
        {
            echo "<td><a href=\"#\" onclick=\"show_dlmenu(event, '$rcpid', '".urlencode($predname[$oid])."', '{$predate[$oid]}', '{$preddock[$oid]}');\">"
                .round($predictions[$oid], 2)."</a></td>";
            $pred_shown[$oid] = true;
        }
        else if (@$ligands_inprogress[$oid]) echo "<td><span title=\"A prediction is currently running for this receptor+ligand pair.\">&#x23F3;</span></td>";
        else echo "<td>&nbsp;</td>";
    }

    echo "<td style=\"white-space: nowrap;\">" . implode(", ",make_clickable_notes($pq)) . "</td>\n";
    echo "</tr>\n";
}   // Each pair.

if (count($predictions))
{
    foreach ($predictions as $oid => $pred)
    {
        if (!@$pred_shown[$oid])
        {
            $lodor = $odors[$oid];
            echo "<tr>\n";
            echo "<td><a href=\"odorant.php?o=$oid\" style=\"white-space: nowrap;\"";
        
            $smilesesc = str_replace("\\", "\\\\", $lodor['smiles']);
            echo " onmouseenter=\"showSkeletal(event, '$smilesesc');\"";
            echo " onmouseout=\"$('#skeletal').hide();\"";
            echo ">{$lodor['full_name']}</td>\n";
            echo "<td>&nbsp;</td>\n";
            echo "<td>&nbsp;</td>\n";
            echo "<td>&nbsp;</td>\n";

            echo "<td><a href=\"#\" onclick=\"show_dlmenu(event, '$rcpid', '".urlencode($predname[$oid])."', '{$predate[$oid]}', '{$preddock[$oid]}');\">"
                .round($predictions[$oid], 2)."</a></td>";
            echo "</td>";
            $pred_shown[$oid] = true;

            $pq = [];
            foreach ($lodor['aroma'] as $refurl => $notes) $pq = array_merge($pq, $notes);
            $pq = array_unique($pq);
            echo "<td style=\"white-space: nowrap;\">" . implode(", ",make_clickable_notes($pq)) . "</td>\n";
            echo "</tr>\n";
        }
    }
}

?>
</table>
</div>
</div>
</div>

<div id="Comparison" class="tabcontent">
<div class="box">
<div class="row content scrollh">

<?php
if (!@$_REQUEST['cmp'])
{
    ?><form action="" method="GET"><?php
    foreach ($_REQUEST as $k => $v)
    {
        if ($k != "cmp") echo "<input type=\"hidden\" name=\"$k\" value=\"$v\">\n";
    }
    ?><select name="cmp">
        <?php
            foreach ($prots as $protid => $p)
            {
                if ($protid == $rcpid) continue;
                echo "<option value=\"$protid\">$protid</option>\n";
            }
        ?>
    </select>
    <input type="submit" value="Compare">
    </form>
    <?php
}
else
{
    $cmp = $_REQUEST['cmp'];
    $tode  = all_empirical_pairs_for_receptor($rcpid, true);
    $touto = all_empirical_pairs_for_receptor($cmp,   true);

    $max = max(max($tode), max($touto));

    $arr = [];
    foreach ($tode as $k => $v) $arr[$k] = [$v, 0];
    foreach ($touto as $k => $v) $arr[$k][1] = $v;

    function cmpcmp($a, $b)
    {
        $aa = @$a[1] - @$a[0];
        $bb = @$b[1] - @$b[0];

        if ($aa > $bb) return 1;
        else if ($bb > $aa) return -1;
        else return 0;
    }

    uasort($arr, "cmpcmp");

    $xvals = [];

    $arrnotes = [];
    $allnotes = [];
    $notefreq = [];
    foreach ($arr as $oid => $vals)
    {
        $xvals[$oid] = @$vals[1] - @$vals[0];
        $notes = [];
        if (@$odors[$oid]['aroma'])
        {
            foreach ($odors[$oid]['aroma'] as $ref => $a)
            {
                if (false!==strpos($ref, "primaryodors")) continue;
                foreach ($a as $n)
                {
                    $notes[$n] = $n;
                    $allnotes[$n] = $n;

                    if (isset($notefreq[$n])) $notefreq[$n]++;
                    else $notefreq[$n] = 1;
                }
            }
        }
        $arrnotes[$oid] = $notes;
    }

    $notecorrs = [];
    foreach ($allnotes as $note)
    {
        if ($notefreq[$note] < 5) continue;
        $yvals = [];
        foreach ($arrnotes as $oid => $an) $yvals[$oid] = in_array($note, $an) ? 1 : 0;
        $corr = correlationCoefficient($xvals, $yvals);
        $notecorrs[$note] = $corr;
    }

    ?><table>
        <tr>
            <th style="text-align: left;">Odorant</th>
            <th style="text-align: right;"><?php echo $rcpid; ?></th>
            <th style="text-align: center;">|</th>
            <th style="text-align: left;"><?php echo $cmp; ?></th>
            <th style="text-align: left;">Notes</th>
        </tr>
    <?php
    foreach ($arr as $oid => $vals)
    {
        if (max($vals) <= 0) continue;

        foreach ($vals as $i => $v)
            if ($v < 0)
            {
                $vals[1-$i] += $v;
                $vals[$i] -= $v;
            }
        ?>
        <tr>
            <td style="text-align: left;">
                <a href="odorant.php?o=<?php echo $oid; ?>">
                    <?php echo @$odors[$oid]['full_name']; ?>
                </a>
            </td>
            <td style="text-align: right; color: #39f;"><?php
            for ($i=0; $i<$vals[0]; $i += 0.5) echo "&block;";
            ?></td>
            <td style="text-align: center;">|</td>
            <td style="text-align: left; color: #983;"><?php
            for ($i=0; $i<$vals[1]; $i += 0.5) echo "&block;";
            ?></td>
            <?php
            $notes = $arrnotes[$oid];

            foreach ($notes as $k => $n)
            {
                if (@$notecorrs[$n] > 0)
                {
                    if ($notecorrs[$n] >= 0.75) $notes[$k] = "<b style=\"color: #f00;\">$n</b>";
                    else if ($notecorrs[$n] >= 0.5) $notes[$k] = "<b style=\"color: #f60;\">$n</b>";
                    else if ($notecorrs[$n] >= 0.25) $notes[$k] = "<b style=\"color: #f90;\">$n</b>";
                    else if ($notecorrs[$n] >= 0.1) $notes[$k] = "<span style=\"color: #f93;\">$n</span>";
                    else if ($notecorrs[$n] >= 0.05) $notes[$k] = "<span style=\"color: #b97;\">$n</span>";
                }
                else if (@$notecorrs[$n] < 0)
                {
                    if ($notecorrs[$n] <= -0.75) $notes[$k] = "<b style=\"color: #00f;\">$n</b>";
                    else if ($notecorrs[$n] <= -0.5) $notes[$k] = "<b style=\"color: #06f;\">$n</b>";
                    else if ($notecorrs[$n] <= -0.25) $notes[$k] = "<b style=\"color: #09f;\">$n</b>";
                    else if ($notecorrs[$n] <= -0.1) $notes[$k] = "<span style=\"color: #39f;\">$n</span>";
                    else if ($notecorrs[$n] <= -0.05) $notes[$k] = "<span style=\"color: #79b;\">$n</span>";
                }
            }

            $notes = implode(", ", $notes);
            echo "<td>$notes</td>\n";
            ?>
        </tr>
        <?php
    }

    ?></table><?php
}
?>

</div>
</div>
</div>

<div id="Refs" class="tabcontent">
<div class="box">
<div class="row content scrollh">
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
</div>
</div>

<div id="Structure" class="tabcontent">
    <iframe id="viewer"></iframe>
</div>


</div>
<div id="skeletal"></div>
<script>
$('#skeletal').hide();
</script>

<div id="dlmenu" onclick="$('#dlmenu').hide();">
    <div align="right" style="margin-bottom: -25px;">
        <a href="#" class="closebtn" onclick="$('#dlmenu').hide();">&#xd7;</a>
    </div>
    <div id="dl_mnu_prot">&nbsp;</div>
    <div id="dl_mnu_lig">&nbsp;</div>
    <br>
    Code version:<div id="dl_cd_v">&nbsp;</div><br>
    Docker used:<div id="dl_dock">&nbsp;</div><br>
    Files:<br>
    <table class="ctxmenu">
        <tr><td>Active model:</td>
            <td><a href="#" onclick="view_file($('#dl_acv_mdl')[0].href); $('#dlmenu').hide();">text</a></td>
            <td><a id="vw_acv_mdl_3d" href="" target="_3d">3D</a></td>
            <td><a id="dl_acv_mdl" href="" target="_dl" onclick="$('#dlmenu').hide();">download</a></td>
        </tr>
        <tr><td>Inactive model:</td>
            <td><a href="#" onclick="view_file($('#dl_iacv_mdl')[0].href); $('#dlmenu').hide();">text</a></td>
            <td><a id="vw_iacv_mdl_3d" href="" target="_3d">3D</a></td>
            <td><a id="dl_iacv_mdl" href="" target="_dl" onclick="$('#dlmenu').hide();">download</a></td>
        </tr>
        <tr><td>Active dock:</td>
            <td><a href="#" onclick="view_html_file($('#dl_acv_dc')[0].href.replace('download.php','tabbed.php')); $('#dlmenu').hide();">text</a></td>
            <td><a id="vw_acv_dc_3d" href="" target="_3d">3D</a></td>
            <td><a id="dl_acv_dc" href="" target="_dl" onclick="$('#dlmenu').hide();">download</a></td>
        </tr>
        <tr><td>Inactive dock:</td>
            <td><a href="#" onclick="view_html_file($('#dl_iacv_dc')[0].href.replace('download.php','tabbed.php')); $('#dlmenu').hide();">text</a></td>
            <td><a id="vw_iacv_dc_3d" href="" target="_3d">3D</a></td>
            <td><a id="dl_iacv_dc" href="" target="_dl" onclick="$('#dlmenu').hide();">download</a></td>
        </tr>
        <tr><td>JSON entry:</td>
            <td><a href="#" onclick="view_file($('#dl_json')[0].href); $('#dlmenu').hide();">text</a></td>
            <td>&nbsp;</td>
            <td><a id="dl_json" href="" target="_dl" onclick="$('#dlmenu').hide();">download</a></td>
        </tr>
    </table>
</div>
<script>
$('#dlmenu').hide();
</script>

<div id="txtview" style="display: none;">
    <div align="right">
        <a href="#" class="closebtn" onclick="$('#txtview').hide();">&#xd7;</a>
    </div>
    <textarea id="txttext" enabled="false">
        Loading...
    </textarea>
</div>

<div id="htmview" style="display: none;">
    <div align="right">
        <a href="#" class="closebtn" onclick="$('#htmview').hide();">&#xd7;</a>
    </div>
    <div id="htmhtml">
        Loading...
    </div>
</div>

