# make_syn_cmds.awk
# Makes alternative Tcl command names for original commands specified
# in a data file (via the command line).
#
# PJ, 04-Oct-00
# -------------------------------------------------------------------

BEGIN {
    print "# moc_syn_cmds.tcl"
    print "#"
    print "#@doc"
    print "# Contains automatically generated wrapper procedures."
    print "# These procedures provide synonymous commands for"
    print "# a number of the IMOC procedures and functions."
    print "#@end"
    print " "
    print "package provide imoc 0.1"
    print " "
    print "# ----------- begin synonymous commands -------------"
}

$1 != "#" { # for each non-comment line in the input data
    original_name = $1
    for (i = 2; i <= NF; ++i) {
        nick_name = $i

        print "proc ", nick_name, " { args } { "
        print "    #@proc"
        print "    #@doc"
        print "    # Equivalent to the ", original_name, " command."
        print "    #@end"
        print "    eval ", original_name, " $args"
        print "}; #end proc"
        print "#"
    }
}

END {
   print "# --------- end of synonymous commands -----------"
}

# ------------ end of make_syn_cmds.awk -----------------------------
