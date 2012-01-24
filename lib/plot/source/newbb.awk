# newbb.awk
# Change the Boundingbox comment in a postscript file.
# P.J.
# 18-Apr-99
#
# Invoke with
# awk -v x1=0 -v y1=0 -v x2=500 -v y2=300 -f newbb.awk infile > outfile
# where the numbers for x1..y2 are examples.
#
BEGIN { if (x2 <= x1 || y2 <= y1) {
            print "Invalid bounding box values." > /dev/stderr
            exit; 
        }
      }
$1 == "%%BoundingBox:" { print $1, x1, y1, x2, y2 }
$1 != "%%BoundingBox:" { print $0 }
