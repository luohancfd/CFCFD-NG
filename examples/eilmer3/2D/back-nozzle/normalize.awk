# normalize.awk
# Normalize the surface pressure over the length of the nozzle.
BEGIN {
    p0 = 418.7e3 
    print "# Normalized surface pressure for the Back nozzle (simulation)"
    print "# x(inches) p/pt"
}

$1 != "#" {  # For non-comment lines in the data file do... 
   p = $9
   r = $2
   x = $1
   print x/0.0254, p/p0
}
