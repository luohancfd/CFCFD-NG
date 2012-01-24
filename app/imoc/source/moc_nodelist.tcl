# moc_nodelist.tcl
# Graphical User-Interface part of the MOC program:

#@doc
# Several operations on the characteristic mesh require
# the indices of particular (selected) nodes.
# The following functions handle the book-keeping for
# a list of "selected" nodes.
# Procedures which do the actual computations can then
# get their node identifiers from this list. <BR>
#
# It is intended that these procedures be invoked from
# the menu system.
#
# No error trapping is done; be careful.
#@end

package provide imoc 0.1

#------------------------------------------------------------

proc nodelist_AddNode { nid } {
    #@proc
    global selectedNodeList
    lappend selectedNodeList $nid
}; # end proc

# --------------------------------------------------------

proc nodelist_Clear {} {
    #@proc
    global selectedNodeList
    foreach nid $selectedNodeList {
        recolourNode $nid
    }; # end foreach
    set selectedNodeList {}
}; # end proc

# --------------------------------------------------------

proc nodelist_GetNodeId { {i 0} {list_action keep} } {
    #@proc
    # Returns the i-th node index from the selected-node list.
    # Remember that elements are counted from zero to n-1.
    # If list_action is "delete", the entry is removed from
    # the selectedNodeList; default action is to keep.
    # If the list is empty or if $i is too large, returns -1.
    #
    global selectedNodeList
    set list_length [llength $selectedNodeList]
    if { $list_length == 0 || $i >= $list_length } { 
        return -1 
    } else {
        set nid [lindex $selectedNodeList $i]
        if { [string compare $list_action "delete"] == 0 } {
            set selectedNodeList [lreplace $selectedNodeList $i $i]
        }; # end if
        return $nid
    }; # end if
}; # end proc

# --------------------------------------------------------

proc nodelist_DeleteNode {} {
    #@proc
    # If the nodelist is not empty, delete the first node 
    # in the selection list and call the procedure again.
    global gd
    set nid [nodelist_GetNodeId 0 delete]
    if { $nid < 0 } {
        puts "Node list is empty."
    } else {
        if { $gd(debug) == 1 } {
            puts "Delete node $nid"
        }; # end if
        if { $gd(echoCommands) == 1 } {
            puts "SubCommand: DeleteNode $nid"
        }; # end if
        DeleteNode $nid
        nodelist_DeleteNode
    }; # end if
}; # end proc

# --------------------------------------------------------

proc nodelist_InteriorNode {} {
    #@proc
    # Pick up the first and second nodes from the selection list
    # and create a new internal node.
    # Assumption: first node is on the C- characteristic
    #             second node is on the C+ characteristic
    #
    global gd
    if { $gd(debug) == 1 } {
        puts "Create an interior node."
    }; # end if
    set n1 [nodelist_GetNodeId 0]
    set n2 [nodelist_GetNodeId 1]
    if {$n1 < 0 || $n2 < 0} {
        puts "nodelist_InteriorNode: Invalid nodes in selection list: n1=$n1, n2=$n2"
        return
    }; # end if
    if { $gd(echoCommands) == 1 } {
        puts "SubCommand: InteriorNode $n1 $n2 -1"
    }; # end if
    set n4 [InteriorNode $n1 $n2 -1]
    if { $gd(debug) == 1 } {
        puts "New node is n4=$n4."
    }; # end if
}; # end proc

# --------------------------------------------------------

proc nodelist_InsertNodeHalfway {} {
    #@proc
    # Pick up the first and second nodes from the selection list
    # and create a new internal node halfway between them.
    #
    global gd
    if { $gd(debug) == 1 } {
        puts "Create an interior node halfway between 2 others."
    }; # end if
    set n1 [nodelist_GetNodeId 0]
    set n2 [nodelist_GetNodeId 1]
    if {$n1 < 0 || $n2 < 0} {
        puts "nodelist_InsertNodeHalfway: "
        puts "Invalid nodes in selection list: n1=$n1, n2=$n2"
        return
    }; # end if
    if { $gd(echoCommands) == 1 } {
        puts "SubCommand: InsertNode $n1 $n2 -1 0.5"
    }; # end if
    set n4 [InsertNode $n1 $n2 -1 0.5]
    if { $gd(debug) == 1 } {
        puts "New node is n4=$n4."
    }; # end if
}; # end proc

# --------------------------------------------------------

proc nodelist_WallNode { Ctype iw } {
    #@proc
    # Pick up the first node from the selection list
    # and create a new node on Wall iw.
    # Ctype indicates the type of the incoming characteristic:
    # CPlus or CMinus
    #
    global gd
    if { [WallIsPresent $iw] == 0 } {
        puts "nodelist_WallNode: Wall $iw is not present."
        return
    }; # end if
    set n1 [nodelist_GetNodeId 0]
    if {$n1 < 0} {
        puts "nodelist_WallNode: Invalid node in selection list: n1=$n1"
        return
    }; # end if
    if {[string compare $Ctype CMinus] == 0} {
        if { $gd(echoCommands) == 1 } {
            puts "SubCommand: CMinusWallNode $iw $n1 -1"
        }; # end if
        set n4 [CMinusWallNode $iw $n1 -1]
    } else {
        if { $gd(echoCommands) == 1 } {
            puts "SubCommand: CPlusWallNode $iw $n1 -1"
        }; # end if
        set n4 [CPlusWallNode $iw $n1 -1]
    }; # end if
    if { $gd(debug) == 1 } {
        puts "New node is n4=$n4."
    }; # end if
}; # end proc

# --------------------------------------------------------

proc nodelist_FreeBndyNode { Ctype } {
    #@proc
    # Pick up the first node from the selection list as 
    # being on a free boundary and the second node as being
    # an internal node on a characteristic which intersects
    # the free boundary.
    # Create a new node on the boundary and on the characteristic.
    # Ctype indicates the type of the incoming characteristic:
    # CPlus or CMinus
    #
    global gd
    set n0 [nodelist_GetNodeId 0]
    if {$n0 < 0} {
        puts "nodelist_FreeBndyNode: Invalid node in selection list: n0=$n0"
        return
    }; # end if
    set nChar [nodelist_GetNodeId 1]
    if {$nChar < 0} {
        puts "nodelist_FreeBndyNode: Invalid node in selection list: nChar=$nChar"
        return
    }; # end if
    if {[string compare $Ctype CMinus] == 0} {
        if { $gd(echoCommands) == 1 } {
            puts "SubCommand: CMinusFreeBndyNode $n0 $nChar -1"
        }; # end if
        set n4 [CMinusFreeBndyNode $n0 $nChar -1]
    } else {
        if { $gd(echoCommands) == 1 } {
            puts "SubCommand: CPlusFreeBndyNode $n0 $nChar -1"
        }; # end if
        set n4 [CPlusFreeBndyNode $n0 $nChar -1]
    }; # end if
    if { $gd(debug) == 1 } {
        puts "New node is n4=$n4."
    }; # end if
}; # end proc

# --------------------------------------------------------

proc nodelist_ExtendStreamline {} {
    #@proc
    # Pick up the first (n0), second (n1) and third (n2) nodes from 
    # the selection list and extend a streamline from n0 to 
    # the line joining n1 and n2.
    #
    global gd
    set n0 [nodelist_GetNodeId 0]
    set n1 [nodelist_GetNodeId 1]
    set n2 [nodelist_GetNodeId 2]
    if { $gd(debug) >= 1 } {
        puts "Extend a streamline from node $n0 to nodes $n1, $n2."
    }; # end if
    if {$n1 < 0 || $n2 < 0 || $n0 < 0} {
        puts "nodelist_ExtendStreamline: Invalid nodes in selection list:"
        puts "    n0=$n0, n1=$n1, n2=$n2"
        return
    }; # end if
    if { $gd(echoCommands) == 1 } {
        puts "SubCommand: AddStreamNode $n0 $n1 $n2 -1 0"
    }; # end if
    set n4 [AddStreamNode $n0 $n1 $n2 -1 0]
    if { $gd(debug) >= 1 } {
        puts "New node is n4=$n4."
        puts "Node data: [GetNodeData $n4]"
    }; # end if
}; # end proc

# --------------------------------------------------------

proc nodelist_ExtendStreamlineShepard {} {
    #@proc
    # Pick up the first (n0) from the selection list and 
    # extend a streamline from n0 by a previously set step-size.
    #
    global gd
    set n0 [nodelist_GetNodeId 0]
    set dL $gd(StreamStepSize)
    if { $gd(debug) >= 1 } {
        puts "Extend a streamline from node $n0 by step $dL."
    }; # end if
    if {$n0 < 0} {
        puts "nodelist_ExtendStreamlineShepard: Invalid node in selection list:"
        puts "    n0=$n0"
        return
    }; # end if
    if { $gd(echoCommands) == 1 } {
        puts "SubCommand: StepStreamNode $n0 -1 $dL"
    }; # end if
    set n4 [StepStreamNode $n0 -1 $dL]
    if { $gd(debug) >= 1 } {
        puts "New node is n4=$n4."
        puts "Node data: [GetNodeData $n4]"
    }; # end if
}; # end proc

# --------------------------------------------------------

proc nodelist_MarchAlongCharacteristic { Ctype direction } {
    #@proc
    # Pick up the first node and second node from the selection list
    # Start at n1 and create a new characteristic 
    # marching along the old characteristic marked by n0
    # in the direction specified.
    # Ctype indicates the type of the characteristic: CPlus or CMinus
    #
    global gd
    set n0 [nodelist_GetNodeId 0]
    set n1 [nodelist_GetNodeId 1]
    if { $gd(debug) >= 1 } {
        puts "nodelist_MarchAlongCharacteristic:"
        puts "Starting nodes $n0, $n1, dir=$direction, Ctype=$Ctype"
    }; # end if
    if {$n0 < 0 || $n1 < 0} {
        puts "nodelist_MarchAlongCharacteristic: "
        puts "Invalid nodes in selection list: n0=$n0, n1=$n1"
        return
    }; # end if
    if {[string compare $Ctype CMinus] == 0} {
        if { $gd(echoCommands) == 1 } {
            puts "SubCommand: MarchAlongCMinus $n0 $n1 $direction"
        }; # end if
        set new_nodes [MarchAlongCMinus $n0 $n1 $direction]
    } else {
        if { $gd(echoCommands) == 1 } {
            puts "SubCommand: MarchAlongCPlus $n0 $n1 $direction"
        }; # end if
        set new_nodes [MarchAlongCPlus $n0 $n1 $direction]
    }; # end if
    if { $gd(debug) >= 1 } {
        puts "New nodes generated: $new_nodes."
    }; # end if
}; # end proc

# --------------------------------------------------------
