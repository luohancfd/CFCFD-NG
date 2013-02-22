#!/bin/bash
#
# rsync_data_X2.sh
#
# This shell script was written to back up the tunnel data to the triton 
# server
#
# Carolyn Jacobs
# 16-Feb-07

# Disallow the use of undefined variables
shopt -s -o nounset

# Global declarations
declare -rx SCRIPT=${0##*/}           # SCRIPT is the name of this script
declare -rx rsync="/usr/bin/rsync"    # the rsync command - man 1 rsync

# Sanity checks - the redirect is to the standard error file (>&2)
if test -z "$BASH" ; then
    printf "$SCRIPT:$LINENO: please run this script with the BASH shell\n" >&2
    exit 192
fi
if test ! -x "$rsync" ; then
    printf "$SCRIPT:$LINENO: the command $rsync is not available - \
aborting\n" >&2
    exit 192
fi

# Step 1: User defines the shot number to be backed up
read -p "Shot number to be backed up (eg. x2s100):" -r SHOT
printf "Shot number to be backed up is:%s \n" $SHOT

# Step 2: Change the group ownership of the shot directory to tunnels
chgrp -R tunnels /data/X2/$SHOT

# Step 3: Copy the data across to the archive directory on Triton
#  - preserving the group 
rsync -avg /data/X2/$SHOT/ triton:/archive1/tunnels/X2/$SHOT

# Cleanup
exit 0    # all is well
