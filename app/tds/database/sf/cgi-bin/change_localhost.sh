#! /bin/sh
# change_localhost.sh
# Replaces all occurrences of localhost with the mech web server address
# in all the files that are specified on the command line.

for f in $@
do
    echo "processing $f"
    cp $f $f.old
    sed -e 's/localhost/www.mech.uq.edu.au/' $f > xxxx
    mv xxxx $f
done
