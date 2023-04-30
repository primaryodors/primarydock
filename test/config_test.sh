clear
make code || exit 1

LIG=$(cat $1 | grep 'LIG ')
for (( i=2; i <= "$#"; i++ )); do
    if [ "${!i}" = "--lig" ]
    then
        i=$((i+1))
        LIG="${!i}"
    fi
done
LIG="${LIG/LIG /}"
LIG="${LIG/sdf\//}"
LIG="${LIG/.sdf/}"

WHPHP=$(which php)
if [ $WHPHP ]
then
    RESULT=$(php -f data/fetch_sdf.php $LIG)

    if echo $RESULT | grep -q 'Odorant not found'; then
        echo $RESULT
        exit 1
    fi
fi

bin/primarydock "$@"

WHSOX=$(which play)
if [ $WHSOX ]
then
    play chime.mp3
fi

