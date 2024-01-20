<?php

if (isset($_REQUEST['smiles']))
{
    chdir(__DIR__);
    $smilesu = escapeshellarg($_REQUEST['smiles']);
    $cmd = "obabel --gen3d -osdf -O\"test.sdf\" -:$smilesu 2>&1";
    echo ("$cmd\n");
    passthru($cmd);
    exit;
}

?><html>
<head>
<title>SMILES viewer</title>
<style>
body
{
    background-color: #000;
    color: #fff;
}

input
{
    border: 1px solid #0ff;
    background-color: #000;
    color: #fff;
}

textarea
{
    border: 1px solid #0ff;
    background-color: #000;
    height: 80%;
    width: 100%;
    color: #0f0;
    font-family: Pet Me 2Y, Courier New, Monospace;
    overflow-x: auto;
    overflow-y: auto;
}

.flex-outer
{
    display: flex;
    flex-direction: row;
    height: 95%;
}

.flex-inner
{
    display: flex;
    flex-direction: column;
    border: 1px solid #00f;
    height: 100%;
    padding: 15px;
}
</style>
<script src="https://code.jquery.com/jquery-3.6.0.min.js" integrity="sha256-/xUj+3OJU5yExlq6GSYGSHk7tPXikynS7ogEvDej/m4=" crossorigin="anonymous"></script>
<script src="https://unpkg.com/ngl@0.10.4/dist/ngl.js"></script>
<script>
var stage;

function set_stage()
{
    stage = new NGL.Stage("viewport");

    stage.setParameters(
    {
        'backgroundColor': '#020408',
        'ambientColor': '#59f',
        'ambientIntensity': 0.6,
        'lightColor': '#fc8',
        'lightIntensity': 0.8,
    });

    window.onresize = function()
    {
        var w = parseInt((window.visualViewport.width - $("#viewport")[0].getBoundingClientRect().top - 5) * 0.666);
        var h = parseInt(window.visualViewport.height - 5);
        $('#viewport').css("width", w+"px").css("height",h+"px");
        stage.setSize(w, h);
    };
    window.onresize();
}

var color_scheme = NGL.ColormakerRegistry.addScheme(function (params)
{
    this.atomColor = function(atom)
    {
		if (typeof atom.atomType == "undefined") return 0xff99ff;
    	var elem = atom.atomType.element;
        if (atom.resno <= 1)
        {
        	elem = atom.atomname.replace(/[^A-Za-z]/g, '');
        	elem = elem.charAt(0).toUpperCase() + elem.slice(1).toLowerCase();
    	}
        switch (elem)
        {
        	case 'H': return atom.isBackbone() ? 0x79949b : 0xe1faff;
			case 'He': return 0xfcaa93;
			case 'Ne': return 0xff4e00;
            
        	case 'C':
            return 0x666666;

        	case 'N': return 0x0033FF;
        	case 'O': return 0xFF2222;
        	case 'P': return 0xff8add;
        	case 'S': return 0xffe022;
        	case 'Se': return 0xffb022;
        	case 'Cl': return 0x339900;
        	case 'Br': return 0x993300;
        	case 'I': return 0x990033;
        	case 'Cu': return 0xe5a072;
        	case 'Zn': return 0x92b0cc;
        	case 'Na': return 0xffcc00;
        	case 'Mg': return 0x99ffcc;
        	case 'K': return 0x9966ff;
        	case 'Ca': return 0xf7f7f7;
        	case 'Fe': return 0x845f42;

        	default: return 0x99aabb;
        }

        return 0xcc00ff;
    };
});


function go(smiles_string)
{
    stage.removeAllComponents();
    url = "<?php echo $_SERVER["PHP_SELF"]; ?>";
	$.ajax(
	{
		url: url,
		cache: false,
		data:
		{
			smiles: smiles_string
		},
		success: function(result)
		{
            $('#result')[0].value = result;
            url = "test.sdf";
            $.ajax(
            {
                url: url,
                cache: false,
                success: function(response)
                {
			        var stringBlob = new Blob( [ response ], { type: 'text/plain'} );
                    stage.loadFile(stringBlob, { ext: "sdf" }).then( function( comp )
                    {
                        var rparam =
                        {
                            multipleBond: true,
                            colorScheme: color_scheme
                        }
                        comp.addRepresentation("ball+stick", rparam);
                        comp.autoView();
                    });
                }
            });
        }
    });
}
</script>
</head>
<body onload="set_stage();">
<div class="flex-outer">
<div class="flex-inner" style="width:33%;">
<label for="smiles">SMILES string:</label>
<br>
<input type="text" name="smiles" id="smiles" value="CCO">
<br>
<input type="button" value="Go" onclick="go($('#smiles')[0].value);">
<br><br>
<textarea id="result">
</textarea>
</div>
<div class="flex-inner" style="width:66%;">
<div id="viewport" style="width:100%; height:100%; background-color: #020408;"></div>
</div>
</div>

