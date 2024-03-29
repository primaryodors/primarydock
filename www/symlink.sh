
# To use this script, you must have Apache2 or some other web server installed and enabled.
# Make sure to check that the destination path is correct. For Apache2, the default public folder
# is /var/www/html/ but if your configuration is different or you are using another web server
# package then you should modify this file to point to the correct path.
#
# The following packages are also required: php, php-gd, php-curl.
#
# You may have to run this script under root if you are installing it on a dedicated machine, vs.
# for example a web host where you would have access to the public folder by default.
#

# Obtain the correct working directory.
cd "$(dirname "$0")"

# Create the symlink.
ln -s "$(pwd)" /var/www/html/primarydock
