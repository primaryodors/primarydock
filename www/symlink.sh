
# To use this script, you must have Apache2 or some other web server installed and enabled.
# Make sure to check that the destination path is correct. For Apache2, the default public folder
# is /var/www/html/ but if your configuration is different or yo uare using another web server
# package then you should modify this file to point to the correct path.

# Obtain the correct working directory.
cd "$(dirname "$0")"

# Create the symlink.
sudo ln -s "$(pwd)" /var/www/html/primarydock
