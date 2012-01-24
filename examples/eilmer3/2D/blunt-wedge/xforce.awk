# xforce.awk
# Extract the simulation times and axial force values from the log file.
# 

BEGIN {
    print "# time (microseconds)  x-force-total  only-cylinder only-wedge";
}

/XFORCE/ && $5 == 0 {
    # Select just the simulation time and the pressure forces for block 0.
    t = $3;  # in seconds
    fx_p_0 = $9;  # force on cylinder in Newtons
    # Don't do anything until we pick up the wedge data (block 1).
}

/XFORCE/ && $5 == 1 {
    # Select just the simulation time and the pressure forces for block 1.
    t = $3;  # in seconds
    fx_p_1 = $9; # wedge surface in Newtons
    print t*1.0e6, fx_p_0 + fx_p_1, fx_p_0, fx_p_1;
}
