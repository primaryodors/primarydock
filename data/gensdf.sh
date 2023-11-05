
NAME="$1"
SMILES="$2"

echo "Hash for odorants.json:"
echo -n $SMILES | md5sum | sed '/[0-9a-fA-F]/s/[ -]//g'

obabel -:$SMILES --gen3D -osdf -Osdf/$NAME.sdf
