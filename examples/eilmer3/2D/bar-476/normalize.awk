# normalize.awk
# Normalize the surface pressure over the centreline static pressure.
BEGIN {
  p_centre = -1.0;
}

$1 != "#" { 
  p = $9;
  r = $2;
  if (p_centre < 0.0) p_centre = p;
  print r/0.005, p/p_centre;
}
