<?php

if (isset($_REQUEST['smiles']))
{
    chdir(__DIR__);

    $smilesu = escapeshellarg($_REQUEST['smiles']);
    $cmd = "obabel --gen3d -osdf -O\"test.sdf\" -:$smilesu 2>&1";
    echo ("$cmd\n");
    passthru($cmd);

    if (@$_REQUEST["rflp"])
    {
        chdir("..");
        $rflp = escapeshellarg($_REQUEST['rflp']);
        $cmd = "bin/ringflip \"www/test.sdf\" $rflp";
        echo ("$cmd\n");
        passthru($cmd);
    }

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
    font-family: Monospace;
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
            
        	case 'C':
            return 0x666666;

            case 'Li': return 0xff1852;
            case 'Be': return 0x1e4cff;
            case 'B': return 0x774422;
        	case 'N': return 0x0033FF;
        	case 'O': return 0xFF2222;
            case 'F': return 0xddff66;
			case 'Ne': return 0xff4022;
        	case 'Na': return 0xffcc00;
        	case 'Mg': return 0x99ffcc;
            case 'Al': return 0x7aadbd;
            case 'Si': return 0x5d6886;
        	case 'P': return 0xff8add;
        	case 'S': return 0xffe022;
        	case 'Cl': return 0x339900;
            case 'Ar': return 0xd46aff;
        	case 'K': return 0x9966ff;
        	case 'Ca': return 0xf7f7f7;
            case 'Sc': return 0xe0a5ff;
            case 'Ti': return 0x738cad;
            case 'V': return 0x6b6455;
            case 'Cr': return 0x1fad3a;
            case 'Mn': return 0x53264f;
        	case 'Fe': return 0x845f42;
            case 'Co': return 0x3735eb;
            case 'Ni': return 0x05a47c;
        	case 'Cu': return 0xe5a072;
        	case 'Zn': return 0x92b0cc;
            case 'Ga': return 0x656187;
            case 'Ge': return 0x666a59;
            case 'As': return 0x339cff;
        	case 'Se': return 0xffb022;
        	case 'Br': return 0x993300;
            case 'Kr': return 0xffc2ef;
            case 'Rb': return 0x990022;
            case 'Sr': return 0xd60862;
            case 'Y': return 0xb78d6f;
            case 'Zr': return 0x688d87;
            case 'Nb': return 0x749ebc;
            case 'Mo': return 0xb5c24b;
            case 'Tc': return 0x3a78b0;
            case 'Ru': return 0x4c4cb1;
            case 'Rh': return 0x8f48c2;
            case 'Pd': return 0xa5bab9;
            case 'Ag': return 0xdbe5e1;
            case 'Cd': return 0xffcc00;
            case 'In': return 0x2520e1;
            case 'Sn': return 0x908d85;
            case 'Sb': return 0xd38d65;
            case 'Te': return 0xa3142d;
        	case 'I': return 0x990033;
            case 'Xe': return 0x7a95ff;
            case 'Cs': return 0x3467fb;
            case 'Ba': return 0x95e104;
            case 'La': return 0xbcb1ff;
            case 'Ce': return 0xffa429;
            case 'Pr': return 0x96c72e;
            case 'Nd': return 0xd1346c;
            case 'Pm': return 0xcd5a19;
            case 'Sm': return 0xff5308;
            case 'Eu': return 0xff4502;
            case 'Gd': return 0x9d5b6a;
            case 'Tb': return 0x77ff38;
            case 'Dy': return 0xeaff41;
            case 'Ho': return 0x57bf86;
            case 'Er': return 0x18ff2b;
            case 'Tm': return 0x2d46ff;
            case 'Yb': return 0x8cc83b;
            case 'Lu': return 0x3c9fa6;
            case 'Hf': return 0x636b80;
            case 'Ta': return 0x516169;
            case 'W': return 0x2d5d6c;
            case 'Re': return 0x0a7cbd;
            case 'Os': return 0x325565;
            case 'Ir': return 0x5d528a;
            case 'Pt': return 0xa4b3c3;
            case 'Au': return 0xda9d17;
            case 'Hg': return 0x2f6569;
            case 'Tl': return 0x0d8b24;
            case 'Pb': return 0x2d4242;
            case 'Bi': return 0x2960ff;
            case 'Po': return 0x325da3;
            case 'At': return 0x111133;
            case 'Rn': return 0xff19a1;
            case 'Fr': return 0x771117;
            case 'Ra': return 0xff5a64;
            case 'Ac': return 0xbab462;
            case 'Th': return 0x50989f;
            case 'Pa': return 0x794a7e;
            case 'U': return 0x25c24c;
            case 'Np': return 0xffbf29;
            case 'Pu': return 0xfe40d7;

        	default: return 0x99aabb;
        }

        return 0xcc00ff;
    };
});


function go(smiles_string)
{
    $('#gobtn').prop('disabled', true);
    stage.removeAllComponents();
    stage.setSpin(false);
    url = "<?php echo $_SERVER["PHP_SELF"]; ?>";
	$.ajax(
	{
		url: url,
		cache: false,
		data:
		{
			smiles: smiles_string,
            rflp: $('#rflp')[0].value,
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
                        $('#gobtn').prop('disabled', false);
                    });
                }
            });
        }
    });
}
</script>
</head>
<body onload="set_stage(); go($('#smiles')[0].value);">
<div class="flex-outer">
<div class="flex-inner" style="width:33%;">
<label for="smiles">SMILES string:</label>
<br>
<input type="text" name="smiles" id="smiles" value="<?php echo @$_REQUEST['s'] ?: "CCO"; ?>" onchange="$('#rflp')[0].value = '';">
<br>
<label for="rflp">Ring flips (optional):</label>
<br>
<input type="text" name="rflp" id="rflp" value="">
<br>
<input type="button" id="gobtn" value="Go" onclick="go($('#smiles')[0].value);" onchange="go($('#smiles')[0].value);">
<br><br>
<textarea id="result">
</textarea>
</div>
<div class="flex-inner" style="width:60%;">
<div id="viewport" style="width:100%; height:100%; background-color: #020408;"></div>
</div>
</div>

