# moc_kernel.tcl

#@doc
# Tcl part of the unit process module.
# These procedures include compound unit processes and
# some more convenient versions of the C functions.
#@end

package provide imoc 0.1

#
# Some composite processes, built up of sequences of unit processes.
#

proc MarchAlongCMinus { old_first new_first direction } {
    #@proc
    #@doc
    # Purpose: Generate nodes along a new C- characteristic line. <BR>
    # Input  : <BR>
    # old_first : a starting point on an existing C- line <BR>
    # new_first : the starting point on the new C- line 
    #             on which the new nodes are to be generated. <BR>
    # direction : up or down <BR>
    # Output : <BR>
    # Returns a list of node indices on the new C- curve.
    # in the order that they are generated.
    #@end
    global gd
    set nodelist [list $new_first]
    set n1 $new_first
    set n2 $old_first
    set more_to_do 1
    while {$more_to_do == 1} {
        if { $gd(echoCommands) == 1 } {
            puts "SubCommand: InteriorNode $n1 $n2 -1"
        }; # end if
        set n4 [InteriorNode $n1 $n2 -1]
        if {[GetDebugLevel] >= 1} {
            puts "Made node $n4 from nodes $n1 $n2"
        }; # end if
        lappend nodelist $n4
        set n1 $n4
        if {[string compare $direction "up"] == 0} {
            set n2 [GetNodeDataValue $n2 CMinusUp]
        } else {
            set n2 [GetNodeDataValue $n2 CMinusDown]
        }; # end if
        if {[GetDebugLevel] >= 1} {
            puts "Puts next pair $n1 $n2"
        }; # end if
        if {$n2 == -1} { set more_to_do 0 }
    }; # end while
    return $nodelist
}; # end proc MarchAlongCMinus

#---------------------------------------------------------------------

proc MarchAlongCPlus { old_first new_first direction } {
    #@proc
    #@doc
    # Purpose: Generate nodes along a new C+ characteristic line. <BR>
    # Input  : <BR>
    # old_first : a starting point on an existing C+ line <BR>
    # new_first : the starting point on the new C+ line
    #             on which the new nodes are to be generated. <BR>
    # direction : up or down <BR>
    # Output : <BR>
    # Returns a list of node indices on the new C+ curve.
    # in the order that they are generated.
    #@end
    global gd
    set nodelist [list $new_first]
    set n1 $old_first
    set n2 $new_first
    set more_to_do 1
    while {$more_to_do == 1} {
        if { $gd(echoCommands) == 1 } {
            puts "SubCommand: InteriorNode $n1 $n2 -1"
        }; # end if
        set n4 [InteriorNode $n1 $n2 -1]
        if {[GetDebugLevel] >= 1} {
            puts "Made node $n4 from nodes $n1 $n2"
        }; # end if
        lappend nodelist $n4
        set n2 $n4
        if {[string compare $direction "up"] == 0} {
            set n1 [GetNodeDataValue $n1 CPlusUp]
        } else {
            set n1 [GetNodeDataValue $n1 CPlusDown]
        }; # end if
        if {[GetDebugLevel] >= 1} {
            puts "Puts next pair $n1 $n2"
        }; # end if
        if {$n1 == -1} { set more_to_do 0 }
    }; # end while
    return $nodelist
}; # end proc MarchAlongCPlus

#---------------------------------------------------------------------
