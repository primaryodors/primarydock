<?php 

?>
<html>
    <head>
        <title><?php echo @$page_title ?: "PrimaryDock Web App"; ?></title>
        <link rel="stylesheet" href="assets/style.css">
    </head>
    <body>
        <div id="logo">
            <img src="assets/logo.png">
        </div>
        <div id="side">
            <?php include("sidebar.php"); ?>
        </div>
        <div id="main">