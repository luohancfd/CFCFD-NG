#!/bin/sh
# send_to_web_server.sh
#
# Extracted from update_monc_files.sh
# PJ, 21-Oct-05
#

SRC=/data/tunnels
DEST=proba:/home2/moncdata/
FACILITY_LIST="T4 X2 X3 X1 drummond"

# --------------------------------------------------------
echo "Begin sending files to the web server."

for FACILITY in $FACILITY_LIST
do
    echo "Sending $FACILITY to web server;"
    echo "(Enter password at the prompt if necessary.)"
    rsync -av -e ssh --delete $SRC/$FACILITY $DEST
done 

echo "Done sending files to the web server."

