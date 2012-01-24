# moc_plot.tcl

#@doc
# Routines to initialize and close the display canvas
# and to render the nodes onto the canvas.
#@end

package provide imoc 0.1

#---------------------------------------------------------

proc initDisplay {} {
    #@proc
    #@doc
    # The main window contains a display (canvas),
    # a couple of scroll bars for that canvas and
    # a status/hints line below the display. <BR>
    # This procedure also sets the bindings for the canvas
    # so that left-clicking a node will select it and
    # right-clicking a node will display its properties in
    # a dialogue window.
    #@end

    global gd

    if {$gd(debug) == 1} {
        puts "initializeDisplay: begin."
    }; # end if

    set displayFrame .df
    set display      .df.display
    set xsbar        .df.xsbar
    set ysbar        .df.ysbar
    set statusFrame  .sf
    set statusline   .sf.status
    set coordDisplay .sf.cd

    set canvasXsize $gd(canvasXsize)
    set canvasYsize $gd(canvasYsize)
    set vportXsize  $gd(vportXsize)
    set vportYsize  $gd(vportYsize)

    if {$gd(debug) == 1} {
        puts "initializeDisplay: Set up scrolled canvas."
    }; # end if
    frame $displayFrame
    canvas $display -width $vportXsize -height $vportYsize \
        -xscrollcommand "$xsbar set" \
        -yscrollcommand "$ysbar set" \
        -background white
    set sr [list 0 0 $canvasXsize $canvasYsize]
    $display configure -scrollregion $sr
    set fraction [expr 1.0 - ($vportYsize / $canvasYsize)]
    $display yview moveto $fraction

    scrollbar $xsbar -orient horizontal \
        -command "$display xview"
    scrollbar $ysbar -orient vertical \
        -command "$display yview"

    grid $display -row 0 -column 0 -sticky nsew
    grid $ysbar -row 0 -column 1 -sticky ns
    grid $xsbar -row 1 -column 0 -sticky ew
    grid columnconfigure $displayFrame 0 -weight 1
    grid rowconfigure $displayFrame 0 -weight 1

    if {$gd(debug) == 1} {
        puts "initializeDisplay: set up status and coordinate line."
    }; # end if
    frame $statusFrame
    label $statusline -relief sunken -borderwidth 2 \
        -text "Status line..." -anchor w -width 40
    label $coordDisplay -relief sunken -borderwidth 2 \
        -text "X: ... Y: ..." -anchor w -width 25
    pack $coordDisplay -side right
    pack $statusline -side left -fill x -expand yes

    # Pack the status frame first so that it is not lost as the
    # user resizes the window to be too small.
    pack $statusFrame -expand no -fill x -side bottom
    pack $displayFrame -expand yes -fill both -side top

    # Set up a couple of bindings to echo coordinates
    # and select nodes (later)
    bind $display <Motion> { displayCoordinates %x %y }
    bind $display <ButtonPress-1> { pickSomething %x %y }
    bind $display <ButtonPress-3> { 
        nodelist_Clear
        selectNode %x %y 
        editNodeData
        nodelist_Clear
    }

    # remember the name of the canvas and status line
    set gd(display_canvas)  $display
    set gd(status_line)     $statusline
    set gd(coord_display)   $coordDisplay

    if {$gd(debug) == 1} {
        puts "initDisplay: end."
    }; # end if

}; # end proc initDisplay

#----------------------------------------------------------

proc closeDisplay {} {
    #@proc
    # destroy all of the objects created in initDisplay
    set displayFrame .df
    set display      .df.display
    set xsbar        .df.xsbar
    set ysbar        .df.ysbar
    set statusline   .status
    destroy $statusline
    destroy $ysbar
    destroy $xsbar
    destroy $display
    destroy $displayFrame
}; #end proc closeDisplay

#-------------------------------------------------------------

proc showStatusMsg { msg } {
    #@proc
    global gd
    $gd(status_line) configure -text $msg
    update idletasks
}; # end proc showStatusMsg

#-------------------------------------------------------------

proc displayCoordinates { vportX vportY } {
    #@proc
    # Get the current mouse position, convert the coordinates to 
    # world coordinates and display.
    # If we are zooming in on our view, display a set of cross hairs also.
    global gd

    set canvas $gd(display_canvas)
    set cX [$canvas canvasx $vportX]
    set cY [$canvas canvasy $vportY]
    set x [worldX $cX]
    set y [worldY $cY]
    $gd(coord_display) configure -text [format "X:%.3f Y:%.3f" $x $y]

    set cX_min [canvasX $gd(xmin)]
    set cY_min [canvasY $gd(ymin)]
    set cX_max [canvasX $gd(xmax)]
    set cY_max [canvasY $gd(ymax)]
    if { [ string compare $gd(pickAction) pickCorner1 ] == 0 } {
        $canvas delete cursor1
        $canvas create line $cX $cY_min $cX $cY_max -fill red  -tags cursor1
        $canvas create line $cX_min $cY $cX_max $cY -fill red  -tags cursor1
    }; # end cursor1
    if { [ string compare $gd(pickAction) pickCorner2 ] == 0 } {
        $canvas delete cursor2
        $canvas create line $cX $cY_min $cX $cY_max -fill red  -tags cursor2
        $canvas create line $cX_min $cY $cX_max $cY -fill red  -tags cursor2
    }; # end cursor2

}; # end proc displayCoordinates

proc eraseCursors {} {
    #@proc
    global gd
    set canvas $gd(display_canvas)
    $canvas delete cursor1
    $canvas delete cursor2
}; # end proc eraseCursors

#-------------------------------------------------------------

proc plotAxes {} {
    #@proc
    global gd

    set canvas $gd(display_canvas)

    set xmin $gd(xmin)
    set ymin $gd(ymin)
    set xmax $gd(xmax)
    set ymax $gd(ymax)
    set dx   $gd(dx)
    set dy   $gd(dy)

    set x1 [canvasX $xmin]
    set y1 [canvasY $ymin]

    set x2 [canvasX $xmax]
    set y2 [canvasY $ymin]
    $canvas create line $x1 $y1 $x2 $y2 -fill black -tags axes

    set x2 [canvasX $xmin]
    set y2 [canvasY $ymax]
    $canvas create line $x1 $y1 $x2 $y2 -fill black -tags axes

    for {set x $xmin} {$x <= [expr $xmax + 1.0e-6]} {set x [expr $x + $dx]} {
        set x1 [canvasX $x]
        set x2 $x1
        set y1 [canvasY $ymin]
        set y2 [expr $y1 + 5]
        $canvas create line $x1 $y1 $x2 $y2 -fill black -tags axes
        $canvas create text [expr $x2] [expr $y2+5] -fill black \
            -text [format "%.3g" $x] -anchor n -tags axes
    }; # end for x

    for {set y $ymin} {$y <= [expr $ymax + 1.0e-6]} {set y [expr $y + $dy]} {
        set x1 [canvasX $xmin]
        set x2 [expr $x1 - 5]
        set y1 [canvasY $y]
        set y2 $y1
        $canvas create line $x1 $y1 $x2 $y2 -fill black -tags axes
        $canvas create text [expr $x2-5] [expr $y2] -fill black \
            -text [format "%.3g" $y] -anchor e -tags axes
    }; # end for x

}; # end proc plotAxes

proc eraseAxes {} {
    #@proc
    global gd
    set canvas $gd(display_canvas)
    $canvas delete axes
}; # end proc eraseAxes

#-------------------------------------------------------------

proc plotWalls {} {
    #@proc
    global gd

    set canvas $gd(display_canvas)

    set n_step 100
    set dt [expr 1.0 / $n_step]

    for {set iw 0} {$iw < 2} {incr iw} {
        if [WallIsPresent $iw] {
            set t 0.0
            set x1 [canvasX [WallPosX $iw $t]]
            set y1 [canvasY [WallPosY $iw $t]]
            for {set j 1} {$j <= 100} {incr j} {
                set t [expr $j * $dt]
                set x2 [canvasX [WallPosX $iw $t]]
                set y2 [canvasY [WallPosY $iw $t]]
                $canvas create line $x1 $y1 $x2 $y2 \
                    -fill brown -width 2.0 -tags walls
                set x1 $x2
                set y1 $y2
            }; # end for j
        }; # end if
    }; # end for iw

}; # end proc plotWalls

proc eraseWalls {} {
    #@proc
    global gd
    set canvas $gd(display_canvas)
    $canvas delete walls
}; # end proc eraseWalls

#-------------------------------------------------------------

proc coordinatesAreBad { x y } {
    #@proc
    # Look for the substring ERROR in either value.
    # This indicates that GetNodeDataValue has detected an error.
    set resultx [string match *ERROR* $x]
    set resulty [string match *ERROR* $y]
    set result [expr $resultx || $resulty]
    return $result
}; # end proc coordinatesAreBad

#-------------------------------------------------------------

proc plotMesh {} {
    #@proc
    global gd
    global nodeIndex
    global objectIndex

    set display $gd(display_canvas)

    if [info exist nodeIndex] {
        unset nodeIndex
    }; # end if
    if [info exist objectIndex] {
        unset objectIndex
    }; # end if

    if {$gd(debug) == 1} {
        puts "plotMesh: Begin plotting nodes, mesh and streamlines..."
    }; # end if
    set nid -1;  # start at the beginning
    set nid [GetNextNodeId $nid]
    set node_count 0;  # used to trigger GUI updates
    while {$nid != -1} {
        set x [GetNodeData_C $nid X]
        set y [GetNodeData_C $nid Y]
        if [coordinatesAreBad $x $y] {
            puts "Bad coordinates ($x, $y) for node $nid"
            set nid [GetNextNodeId $nid]
            continue; # Jump straight to next pass
        }; # end if
        if {$gd(debug) == 1} {
            puts "Plot node $nid at ($x, $y)."
        }; # end if

        set x [canvasX $x]
        set y [canvasY $y]
        if {$gd(debug) == 1} {
            puts "Canvas coordinates are ($x, $y)."
        }; # end if

        # Decide on the size of the node to allow room for
        # a node number if requested
        if { $gd(showNodeNumbers) == 1 } {
            set r 15
        } else {
            set r 3
        }; #end if

        # Draw the node as a filled circle.
        set objectid [ $display create oval \
            [expr $x-$r] [expr $y-$r] \
            [expr $x+$r] [expr $y+$r] \
            -outline black -fill gray -tags node ]
        # Remember which node is associated with each oval object
        set nodeIndex($objectid) $nid
        set objectIndex($nid) $objectid
        # Add the node id, maybe.
        if { $gd(showNodeNumbers) == 1 } {
            $display create text $x $y \
                -text $nid -anchor center -tags nodeid
        }; #end if

        if { $gd(showCharMesh) == 1 } {
            # Find the coordinates of any upstream node and
            # plot the segment of characteristic connecting it
            # to the present node
            set nidUp [GetNodeData_C $nid CPlusUp]
            if { $nidUp != -1 } {
                set xUp [GetNodeData_C $nidUp X]
                set yUp [GetNodeData_C $nidUp Y]
                if [coordinatesAreBad $xUp $yUp] {
                    puts "Bad coordinates ($xUp, $yUp) for node $nidUp"
                    set nid [GetNextNodeId $nid]
                    continue; # Jump straight to next pass
                }; # end if
                if {$gd(debug) == 1} {
                    puts "CPlusUp node $nidUp at ($xUp, $yUp)."
                }; # end if
                set xUp [canvasX $xUp]
                set yUp [canvasY $yUp]
                set this_line [$display create line \
                    $xUp $yUp $x $y -fill red -tags lines]
                # Arrange the nodes to sit over the line segments
                # by sending the line to the bottom of display list.
                $display lower $this_line
            }; # end if
            set nidUp [GetNodeData_C $nid CMinusUp]
            if { $nidUp != -1 } {
                set xUp [GetNodeData_C $nidUp X]
                set yUp [GetNodeData_C $nidUp Y]
                if [coordinatesAreBad $xUp $yUp] {
                    puts "Bad coordinates ($xUp, $yUp) for node $nidUp"
                    set nid [GetNextNodeId $nid]
                    continue; # Jump straight to next pass
                }; # end if
                if {$gd(debug) == 1} {
                    puts "CMinusUp node $nidUp at ($xUp, $yUp)."
                }; # end if
                set xUp [canvasX $xUp]
                set yUp [canvasY $yUp]
                set this_line [$display create line \
                    $xUp $yUp $x $y -fill blue -tags lines]
                # Arrange the nodes to sit over the line segments
                # by sending the line to the bottom of display list.
                $display lower $this_line
            }; # end if
        }; # end if

        if { $gd(showStreamlines) == 1 } {
            # Find the coordinates of any upstream node and
            # plot the segment of streamline connecting it
            # to the present node
            set nidUp [GetNodeData_C $nid CZeroUp]
            if { $nidUp != -1 } {
                set xUp [GetNodeData_C $nidUp X]
                set yUp [GetNodeData_C $nidUp Y]
                if [coordinatesAreBad $xUp $yUp] {
                    puts "Bad coordinates ($xUp, $yUp) for node $nidUp"
                    set nid [GetNextNodeId $nid]
                    continue; # Jump straight to next pass
                }; # end if
                if {$gd(debug) == 1} {
                    puts "CZeroUp node $nidUp at ($xUp, $yUp)."
                }; # end if
                set xUp [canvasX $xUp]
                set yUp [canvasY $yUp]
                set this_line [$display create line \
                    $xUp $yUp $x $y -fill green -width 2.0 -tags lines]
                # Arrange the nodes to sit over the line segments
                # by sending the line to the bottom of display list.
                $display lower $this_line
            }; # end if
        }; # end if

        set nid [GetNextNodeId $nid]; # prepare for next pass
        incr node_count
        if { $node_count == 10 } {
            # update idletasks
            # This turns out to be a bad idea because
            # it prevents Tk from accumulating a number of graphics
            # commands and eliminating redundant renderings.
            set node_count 0
        }; # end if
    }; # end while

    if {$gd(debug) == 1} {
        puts "plotMesh: Finished plotting nodes, mesh and streamlines."
    }; # end if
    showStatusMsg "Nodes and Mesh replotted."

}; # end proc plotMesh

#------------------------------------------------------------

proc eraseMesh {} {
    #@proc
    global gd
    global nodeIndex
    global objectIndex

    set display $gd(display_canvas)

    $display delete node
    $display delete nodeid
    $display delete lines

    if [info exist nodeIndex] {
        unset nodeIndex
    }; # end if
    if [info exist objectIndex] {
        unset objectIndex
    }; # end if

}; # end proc eraseMesh

#----------------------------------------------------------

proc refreshDisplay {} {
    #@proc
    global gd

    eraseCursors
    eraseAxes
    eraseWalls
    eraseMesh

    plotAxes
    plotWalls
    plotMesh
}; # end proc refreshDisplay

#----------------------------------------------------------

proc pickSomething { vportX vportY } {
    #@proc
    global gd
    set can $gd(display_canvas)

    if { [ string compare $gd(pickAction) pickNode ] == 0 } {
        selectNode $vportX $vportY
    } elseif { [ string compare $gd(pickAction) pickCorner1 ] == 0 } {
        set cX [$can canvasx $vportX]
        set cY [$can canvasy $vportY]
        set gd(newX1) [worldX $cX]
        set gd(newY1) [worldY $cY]
        puts "Selected corner 1: $gd(newX1), $gd(newY1)"
        set gd(pickAction) pickCorner2
        showStatusMsg "Left mouse button to select second corner."
    } elseif { [ string compare $gd(pickAction) pickCorner2 ] == 0 } {
        set cX [$can canvasx $vportX]
        set cY [$can canvasy $vportY]
        set gd(newX2) [worldX $cX]
        set gd(newY2) [worldY $cY]
        puts "Selected corner 2: $gd(newX2), $gd(newY2)"
        setXYRanges $gd(newX1) $gd(newY1) $gd(newX2) $gd(newY2)
        refreshDisplay
        set gd(pickAction) pickNode
        showStatusMsg "Finished zooming window."
    }; # end of alternative pickActions

}; # end proc pickSomething

#----------------------------------------------------------

proc startZoomWindow {} {
    #@proc
    global gd
    set gd(pickAction) pickCorner1
    showStatusMsg "Left button to select first corner."
}; # end proc startZoomWindow

#----------------------------------------------------------

proc findNodeBoundingBox {} {
    #@proc
    global gd

    set xmin 0.0
    set ymin 0.0
    set xmax 0.0
    set ymax 0.0

    set this_is_first 1
    set nid -1;  # start at the beginning
    set nid [GetNextNodeId $nid]
    while {$nid != -1} {
        set x [GetNodeData_C $nid X]
        set y [GetNodeData_C $nid Y]
        if {$gd(debug) == 1} {
            puts "Found node $nid at ($x, $y)."
        }; # end if

        if { $this_is_first == 1 } {
            set xmin $x
            set ymin $y
            set xmax $x
            set ymax $y
            set this_is_first 0
        } else {
            if { $x < $xmin } { set xmin $x }
            if { $y < $ymin } { set ymin $y }
            if { $x > $xmax } { set xmax $x }
            if { $y > $ymax } { set ymax $y }
        }; # end if

        set nid [GetNextNodeId $nid]
    }; # end while

    if {$gd(debug) == 1} {
        puts "Bounding box ($xmin, $ymin) -- ($xmax, $ymax)."
    }; # end if
    return "$xmin $ymin $xmax $ymax"
}; # end proc findNodeBoundingBox

#----------------------------------------------------------

proc selectNode { vportX vportY } {
    #@proc
    global gd
    global nodeIndex
    global objectIndex

    set can $gd(display_canvas)
    set cX [$can canvasx $vportX]
    set cY [$can canvasy $vportY]
    set gd(selectedNode) -1
    set selectedNode -1

    # Make a list of all items that overlap the specified pixel
    set objectList [$can find overlapping $cX $cY $cX $cY]

    # Search the list in reverse order to find a node
    set last [llength $objectList]
    set foundNodeList {}
    for {set i [expr $last-1]} {$i >= 0} {incr i -1} {
        set object [lindex $objectList $i]
        set objectType [$can type $object]
        if {$objectType == "oval"} {
            lappend foundNodeList $nodeIndex($object)
        }; # end if
    }; # end for

    if { [llength $foundNodeList] > 1 && $gd(displayDialogForCoincidentNodes) == 1 } {
        # Put up a list of the possible nodes so that the user can select one.
        selectNodeDialog $foundNodeList
        vwait gd(selectedNode)
        # Once the user has finished with the dialog, go and
        # collect the selected node id.
        set selectedNode $gd(selectedNode)
    } elseif { [llength $foundNodeList] == 1 } {
        # Just take the first element in the found list
        set selectedNode [lindex $foundNodeList 0]
    }; # end if

    if { $selectedNode >= 0 } {
        $can itemconfigure $objectIndex($selectedNode) -fill red
        nodelist_AddNode $selectedNode
        showStatusMsg "Node $selectedNode selected."
    }; # end if
}; # end proc selectNode
 
# --------------------------------------------------------

proc selectNodeDialog { listOfNodes } {
    #@proc
    global gd

    if { $gd(debug) == 1 } {
        puts "Start up the selectNodeDialog with the following list of nodes:"
        puts $listOfNodes
    }; # end if

    if { [winfo exists .selectNodeDialog] == 0 } {
        # Create the widget
        toplevel .selectNodeDialog
        wm title .selectNodeDialog "Select Node"

        set sF [frame .selectNodeDialog.selectNodeFrame]
        set bF [frame .selectNodeDialog.buttonFrame]

        set lb [listbox $sF.listbox -height 5 -width 10 -selectmode single ]
        eval $lb insert end $listOfNodes
        $lb selection set 0
        set sb [scrollbar $sF.yscroll -command [list $lb yview] -orient vertical]
        $lb configure -yscrollcommand [list $sb set]
        pack $sb -side right -fill y
        pack $lb -side left -fill both
        $lb see 0
        set gd(selectNodeListbox) $lb

        set bAccept [button $bF.accept -text Accept -command { 
            set xx_lb $gd(selectNodeListbox)
            set xx_nl $gd(selectNodeListbox.listOfNodes)
            set gd(selectedNode) [lindex $xx_nl [$xx_lb curselection]]
            destroy .selectNodeDialog
        } ]
        set bDismiss [button $bF.dismiss -text Dismiss -command {
            destroy .selectNodeDialog
            set gd(selectedNode) -1
        } ]
        pack $bAccept $bDismiss -side left -padx 4 -pady 4

        pack $sF $bF
    } else {
        # the dialog already exists; bring it up with the new list
        set lb $gd(selectNodeListbox)
        $lb delete 0 end
        eval $lb insert end $listOfNodes
        $lb selection set 0
        $lb see 0
        raise .selectNodeDialog
    }; # end if

    set gd(selectNodeListbox.listOfNodes) $listOfNodes
}; # end proc

#----------------------------------------------------------

proc recolourNode { nid } {
    #@proc
    global gd
    global objectIndex

    set can $gd(display_canvas)
    set listofnames [array names objectIndex $nid]
    if { [llength $listofnames] >= 1 } {
        set object $objectIndex($nid)
        $can itemconfigure $object -fill gray
    }; # end if
}; # end proc unselectNode
 
#----------------------------------------------------------

proc plotPSFile { filename } {
    #@proc
    global gd

    set display $gd(display_canvas)
    set w       $gd(canvasXsize)
    set h       $gd(canvasYsize)

    $display postscript -file $filename -x 0 -width $w -y 0 -height $h \
        -pagewidth 150m

}; # end proc

#----------------------------------------------------------

