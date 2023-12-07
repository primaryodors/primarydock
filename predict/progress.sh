watch -n 5 "ps -ef | grep -E ':[0-9][0-9] (bin/primarydock|bin/pepteditor|bin/ic|bin/fyg|obabel)' | grep -v grep; echo; WHSENS=$(which sensors); if [ \$WHSENS ]; then sensors; fi"

