
NAME="$1"
SMILES="$2"

obabel -:$SMILES --gen3D -osdf -Osdf/$NAME.sdf
