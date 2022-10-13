<?php 
global $extra_css, $extra_js;

if ($extra_css && !is_array($extra_css)) $extra_css = [$extra_css];
if ($extra_js  && !is_array($extra_js )) $extra_js  = [$extra_js ];

?>
<html>
    <head>
        <title><?php echo @$page_title ?: "PrimaryDock Web App"; ?></title>
        <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.6.0/jquery.min.js"></script>
        <link rel="stylesheet" href="assets/style.css">

        <?php 
        if ($extra_js)
            foreach ($extra_js as $js)
                echo "<script src=\"$js\"></script>\n";
        if ($extra_css)
            foreach ($extra_css as $css)
                echo "<link rel=\"stylesheet\" href=\"$css\">\n";
        ?>
    </head>
    <body>
        <div id="logo">
            <img src="assets/logo.png">
        </div>
        <div id="side">
            <?php include("sidebar.php"); ?>
        </div>
        <div id="main">