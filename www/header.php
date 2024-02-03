<?php 
global $extra_css, $extra_js, $head_tags;

$customizations = [];

$cwd = getcwd();
chdir(__DIR__);
if (file_exists('custom.json')) $customizations = json_decode(file_get_contents('custom.json'), true);
chdir($cwd);

if ($extra_css && !is_array($extra_css)) $extra_css = [$extra_css];
if ($extra_js  && !is_array($extra_js )) $extra_js  = [$extra_js ];

?>
<html>
    <head>
        <title><?php echo @$page_title ?: @$customizations['title'] ?: "PrimaryDock Web App"; ?></title>
        <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.6.0/jquery.min.js"></script>
        <script src="https://www.lactame.com/lib/openchemlib/5.2.0/openchemlib-minimal.js"></script>
        <script>
        function svg_from_smiles(smiles, w, h)
        {
            var molecule=OCL.Molecule.fromSmiles(smiles);
            return molecule.toSVG(w, h, Math.random.toString(36), {fontWeight: 900})
                .replace(/rgb\(0,0,0\)/g,"rgb(255,255,255)")
                .replace(/fill=\"rgb\(160,0,0\)\">.*<\/text/g, '></text')
                .replace(/rgb\(160,0,0\)/g,"rgb(170,187,204)")
                ;
        }
        </script>
        <link rel="stylesheet" href="assets/style.css?<?php echo time();?>">

        <?php 
        if ($extra_js)
            foreach ($extra_js as $js)
                echo "<script src=\"$js\"></script>\n";
        if ($extra_css)
            foreach ($extra_css as $css)
                echo "<link rel=\"stylesheet\" href=\"$css\">\n";
        if ($head_tags) echo "\n$head_tags\n\n\n";
        ?>
    </head>
    <body>
        <div id="logo">
            <?php if (@$customizations['logo']['href']) echo "<a href=\"{$customizations['logo']['href']}\">"; ?>
            <img src="<?php echo @$customizations['logo']['url'] ?: "assets/logo.png"; ?>"
               alt="<?php echo @$customizations['logo']['alt'] ?: "Main Logo"; ?>"
               style="float: left;">            
            <?php if (@$customizations['logo']['text']) echo "<h2>{$customizations['logo']['text']}</h2>"; ?>
            <?php if (@$customizations['logo']['href']) echo "</a>"; ?>
        </div>
        <div id="side"<?php if (@$customizations['logo']['height'])
        {
            $offset = intval($customizations['logo']['height']) - 37;
            echo " style=\"top: {$offset}px;\"";
        } ?>>
            <?php include("sidebar.php"); ?>
        </div>
        <div id="main"<?php if (@$customizations['logo']['height']) echo " style=\"top: {$customizations['logo']['height']}px;\""; ?>>
