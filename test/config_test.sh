clear
make code || exit 1
bin/primarydock "$@"

WHSOX=$(which play)
if [ $WHSOX ]
then
    play chime.mp3
fi

