
DIR=$(dirname "$(readlink -f "$0")")
crontab -l > /tmp/mycron
echo "* * * * * php -f \"$DIR/predict/jobq.php\" 2>&1 > /dev/null" >> /tmp/mycron
crontab /tmp/mycron
rm /tmp/mycron
