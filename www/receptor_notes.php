<?php 
$cwd = getcwd();
chdir(__DIR__);
require_once("../predict/protutils.php");
require_once("../predict/odorutils.php");
require_once("../predict/statistics.php");
chdir($cwd);

global $prots, $odors;

function correlate_receptors_aromanotes()
{
	global $prots, $odors;

    set_time_limit(3600);

    $xvals = [];            // 2-dimensional per-receptor, per-odorant activity of the given receptor.
    $yvals = [];            // 2 dimensional per-note, per-odorant presence or absence of each perceptual quality.

    $ycount = [];

    foreach (array_keys($prots) as $rcpid)
    {
        $xvals[$rcpid] = [];
        foreach ($odors as $oid => $odor)
        {
            $pair = best_empirical_pair($rcpid, $oid, true);
            $xval = 0.0;
            $samples = 0;
            if (is_array($pair))
            {
                $adjustment = 1;
                if ($odor['full_name'] == 'indole' || $odor['full_name'] == 'skatole') $adjustment = 0.25;        // Two odorants that are dominating lists they oughtn't.
                if (isset($pair['adjusted_curve_top']))
                {
                    $xval += $adjustment * min($pair['adjusted_curve_top'], 10);
                    $samples++;
                }
                if (isset($pair['ec50']))
                {
                    $xval -= $adjustment * $pair['ec50'];
                    $samples++;
                }
                if ($samples) $xval /= $samples;
            }
            $xvals[$rcpid][$oid] = $xval;
        }
    }

    $notes = [];
    foreach ($odors as $oid => $odor)
    {
        foreach ($odor['aroma'] as $refurl => $pqlist)
        {
            $notes = array_merge($notes, $pqlist);
        }
    }
    $notes = array_unique($notes);

    foreach ($notes as $pq)
    {
        $yvals[$pq] = [];
        $ycount[$pq] = 0;
        foreach ($odors as $oid => $odor)
        {
            $yvals[$pq][$oid] = 0;
            foreach ($odor['aroma'] as $refurl => $pqlist)
            {
                if ($refurl == "http://www.primaryodors.org") continue;
                if (in_array($pq, $pqlist))
                {
                    $yvals[$pq][$oid] = 1;
                    $ycount[$pq]++;
                }
            }
        }
    }

    $correlations = [];

    foreach ($xvals as $rcpid => $xv)
    {
        foreach ($yvals as $pq => $yv)
        {
            if ($ycount[$pq] < 3) continue;
            $corr = correlationCoefficient($xv, $yv);
            if ($corr > 0) $correlations[$rcpid][$pq] = $corr;
        }

        if ($correlations[$rcpid]) arsort($correlations[$rcpid]);
    }

    return $correlations;
}

function get_notes_for_receptor($rcpid, $correlations)
{
    if (!isset($correlations[$rcpid])) return "(none)";

    $maxcorr = floatval(max($correlations[$rcpid]));
    $threshold = $maxcorr / 2;
    $retval = [];
    foreach ($correlations[$rcpid] as $pq => $corr)
    {
        if ($corr >= $threshold) $retval[] = $pq;
    }

    $retval = implode(", ", $retval);
    return $retval;
}
