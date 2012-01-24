# moc_menu.tcl
#@doc
# Procedures for posting menus and dialogue windows.
# Since they are all related to invoking GUI elements from the menus, 
# the user should not need to call these procedures directly.
#@end

package provide imoc 0.1

# --------------------------------------------------------

proc initMenu {} {
    #@proc
    #@doc
    # Set up menu bar for the main window.
    #@end
    # m0 is the top-level menubar
    # m1 is a first-level menu
    # m2 is a second-level menu

    global gd

    set m0 [menu .menubar]
    . config -menu $m0
    wm title . "IMOC: Interactive Method of Characteristics."

    set m1 [menu $m0.file -tearoff 0]
    $m0 add cascade -label File -menu $m1
    $m1 add command -label "New" \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: nodelist_Clear; DeleteAll; refreshDisplay"
            }; # end if
            nodelist_Clear; DeleteAll; refreshDisplay
            showStatusMsg "Ready to start over."
        }
    $m1 add command -label "Source Script..." \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: sourceUserScript"
            }; # end if
            showStatusMsg "Begin reading user script..."
            sourceUserScript 
            showStatusMsg "Begin reading user script...Done."
        }
    $m1 add separator
    $m1 add command -label "Load Node Data..." \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: loadNodeData; nodelist_Clear; refreshDisplay"
            }; # end if
            showStatusMsg "Load node data from file..."
            loadNodeData; nodelist_Clear; refreshDisplay 
            showStatusMsg "Load node data from file...Done."
        }
    $m1 add command -label "Load Wall 0..." \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: loadWall_gui 0; nodelist_Clear; refreshDisplay"
            }; # end if
            showStatusMsg "Load wall 0 data from file..."
            loadWall_gui 0; nodelist_Clear; refreshDisplay 
            showStatusMsg "Load wall 0 data from file...Done."
        }
    $m1 add command -label "Load Wall 1..." \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: loadWall_gui 1; nodelist_Clear; refreshDisplay"
            }; # end if
            showStatusMsg "Load wall 1 data from file..."
            loadWall_gui 1; nodelist_Clear; refreshDisplay 
            showStatusMsg "Load wall 1 data from file...Done."
        }
    $m1 add separator
    $m1 add command -label "Save Node Data..." \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: saveNodeData"
            }; # end if
            showStatusMsg "Save node data to file..."
            saveNodeData 
            showStatusMsg "Save node data to file...Done."
        }
    $m1 add command -label "Save Wall 0..." \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: saveWall_gui 0"
            }; # end if
            showStatusMsg "Save wall 0 data to file..."
            saveWall_gui 0 
            showStatusMsg "Save wall 0 data to file...Done."
        }
    $m1 add command -label "Save Wall 1..." \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: saveWall_gui 1"
            }; # end if
            showStatusMsg "Save wall 1 data to file..."
            saveWall_gui 1 
            showStatusMsg "Save wall 1 data to file...Done."
        }
    $m1 add separator
    $m1 add command -label Exit -command "exit"

    #----------------------------------------------------
    set m1 [menu $m0.edit -tearoff 0 -title "Edit Menu"]
    $m0 add cascade -label Edit -menu $m1
    $m1 add command -label "Create Node..." \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: editNodeData [CreateNode -1]; refreshDisplay"
            }; # end if
            showStatusMsg "Creating new node..."
            editNodeData [CreateNode -1]; refreshDisplay 
            showStatusMsg "Creating new node...Done."
        }
    $m1 add command -label "Edit Node Data..." \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: editNodeData; nodelist_Clear"
            }; # end if
            showStatusMsg "Editing node data..."
            editNodeData; nodelist_Clear 
            showStatusMsg "Editing node data...Done."
        }
    $m1 add command -label "Delete Selected Nodes" \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: nodelist_DeleteNode; nodelist_Clear; refreshDisplay"
            }; # end if
            showStatusMsg "Deleting selected nodes..."
            nodelist_DeleteNode; nodelist_Clear; refreshDisplay 
            showStatusMsg "Deleting selected nodes...Done."
        }
    $m1 add separator
    $m1 add command -label "Edit Wall 0..." \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: editWallData 0; refreshDisplay"
            }; # end if
            showStatusMsg "Editing wall 0 data..."
            editWallData 0; refreshDisplay 
            showStatusMsg "Editing wall 0 data...Done."
        }
    $m1 add command -label "Edit Wall 1..." \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: editWallData 1; refreshDisplay"
            }; # end if
            showStatusMsg "Editing wall 1 data..."
            editWallData 1; refreshDisplay 
            showStatusMsg "Editing wall 1 data...Done."
        }
    $m1 add separator
    $m1 add command -label "Clear List of Selected Nodes" \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: nodelist_Clear"
            }; # end if
            showStatusMsg "Clearing list of selected nodes..."
            nodelist_Clear 
            showStatusMsg "Clearing list of selected nodes...Done."
        }
    $m1 add separator
    set m2 [menu $m1.options -tearoff 0 -title "Options Menu"]
    $m1 add cascade -label "Options..." -menu $m2
    $m2 add check -label "Echo Commands to Console" -variable gd(echoCommands)
    $m2 add check -label "Display Dialog for Coincident Nodes" \
        -variable gd(displayDialogForCoincidentNodes)

    #----------------------------------------------------
    set m1 [menu $m0.compute -tearoff 1 -title "Compute Menu"]
    $m0 add cascade -label Compute -menu $m1
    $m1 add command -label "Interior Node" \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: nodelist_InteriorNode; nodelist_Clear; refreshDisplay"
            }; # end if
            showStatusMsg "Compute a new interior node..."
            nodelist_InteriorNode; nodelist_Clear; refreshDisplay 
            showStatusMsg "Compute a new interior node...Done."
        }
    $m1 add command -label "Insert Node Halfway" \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: nodelist_InsertNodeHalfway; nodelist_Clear; refreshDisplay"
            }; # end if
            showStatusMsg "Insert node halfway..."
            nodelist_InsertNodeHalfway; nodelist_Clear; refreshDisplay 
            showStatusMsg "Insert node halfway...Done."
        }

    set m2 [menu $m1.wall -tearoff 0 -title "Compute Wall-Node Menu"]
    $m1 add cascade -label "Wall Node" -menu $m2
    $m2 add command -label "C- to Wall 0" \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: nodelist_WallNode CMinus 0; nodelist_Clear; refreshDisplay"
            }; # end if
            showStatusMsg "New along down C- to wall 0..."
            nodelist_WallNode CMinus 0; nodelist_Clear; refreshDisplay 
            showStatusMsg "New along down C- to wall 0...Done."
        }
    $m2 add command -label "C- to Wall 1" \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: nodelist_WallNode CMinus 1; nodelist_Clear; refreshDisplay"
            }; # end if
            showStatusMsg "New node along C- to wall 1..."
            nodelist_WallNode CMinus 1; nodelist_Clear; refreshDisplay 
            showStatusMsg "New node along C- to wall 1...Done"
        }
    $m2 add command -label "C+ to Wall 0" \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: nodelist_WallNode CPlus 0; nodelist_Clear; refreshDisplay"
            }; # end if
            showStatusMsg "New along along C+ to wall 0..."
            nodelist_WallNode CPlus 0; nodelist_Clear; refreshDisplay 
            showStatusMsg "New along along C+ to wall 0...Done."
        }
    $m2 add command -label "C+ to Wall 1" \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: nodelist_WallNode CPlus 1; nodelist_Clear; refreshDisplay"
            }; # end if
            showStatusMsg "New along along C+ to wall 1..."
            nodelist_WallNode CPlus 1; nodelist_Clear; refreshDisplay 
            showStatusMsg "New along along C+ to wall 1...Done."
        }

    set m2 [menu $m1.free -tearoff 0 -title "Compute Free Boundary Menu"]
    $m1 add cascade -label "Free Boundary Node" -menu $m2
    $m2 add command -label "From C-" \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: nodelist_FreeBndyNode CMinus; nodelist_Clear; refreshDisplay"
            }; # end if
            showStatusMsg "New along along C- to free boundary..."
            nodelist_FreeBndyNode CMinus; nodelist_Clear; refreshDisplay 
            showStatusMsg "New along along C- to free boundary...Done."
        }
    $m2 add command -label "From C+" \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: nodelist_FreeBndyNode CPlus; nodelist_Clear; refreshDisplay"
            }; # end if
            showStatusMsg "New along along C+ to free boundary..."
            nodelist_FreeBndyNode CPlus; nodelist_Clear; refreshDisplay 
            showStatusMsg "New along along C+ to free boundary...Done."
        }

    set m2 [menu $m1.march -tearoff 0 -title "March-Along-Charactaristic Menu"]
    $m1 add cascade -label "March Along Charcteristic" -menu $m2
    $m2 add command -label "Down C-" \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: nodelist_MarchAlongCharacteristic CMinus down;"
                puts "Command: nodelist_Clear; refreshDisplay"
            }; # end if
            showStatusMsg "March downstream along C- ..."
            nodelist_MarchAlongCharacteristic CMinus down; 
            nodelist_Clear; refreshDisplay 
            showStatusMsg "March downstream along C- ...Done."
        }
    $m2 add command -label "Up C-" \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: nodelist_MarchAlongCharacteristic CMinus up;"
                puts "Command: nodelist_Clear; refreshDisplay"
            }; # end if
            showStatusMsg "March upstream along C- ..."
            nodelist_MarchAlongCharacteristic CMinus up; 
            nodelist_Clear; refreshDisplay 
            showStatusMsg "March upstream along C- ...Done."
        }
    $m2 add command -label "Down C+" \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: nodelist_MarchAlongCharacteristic CPlus down;"
                puts "Command: nodelist_Clear; refreshDisplay"
            }; # end if
            showStatusMsg "March downstream along C+ ..."
            nodelist_MarchAlongCharacteristic CPlus down; 
            nodelist_Clear; refreshDisplay 
            showStatusMsg "March downstream along C+ ...Done."
        }
    $m2 add command -label "Up C+" \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: nodelist_MarchAlongCharacteristic CPlus up;"
                puts "Command: nodelist_Clear; refreshDisplay"
            }; # end if
            showStatusMsg "March upstream along C+ ..."
            nodelist_MarchAlongCharacteristic CPlus up; 
            nodelist_Clear; refreshDisplay 
            showStatusMsg "March upstream along C+ ...Done."
        }

    set m2 [menu $m1.stream -tearoff 0 -title "Streamline Menu"]
    $m1 add cascade -label "Streamline Node" -menu $m2
    $m2 add command -label "Use Manually Selected Nodes" \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: nodelist_ExtendStreamline; nodelist_Clear; refreshDisplay"
            }; # end if
            showStatusMsg "Extend streamline to selected nodes..."
            nodelist_ExtendStreamline; nodelist_Clear; refreshDisplay 
            showStatusMsg "Extend streamline to selected nodes...Done."
        }
    $m2 add command -label "Extend by Fixed Step Size" \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: nodelist_ExtendStreamlineShepard;"
                puts "Command: nodelist_Clear; refreshDisplay"
            }; # end if
            showStatusMsg "Extend streamline by fixed step..."
            nodelist_ExtendStreamlineShepard; nodelist_Clear; refreshDisplay 
            showStatusMsg "Extend streamline by fixed step...Done."
        }
    $m2 add command -label "Set Step Size..." \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: displayStreamlineStepSize"
            }; # end if
            showStatusMsg "Start dialog for streamline stepsize..."
            displayStreamlineStepSize 
        }

    $m1 add separator
    $m1 add check -label "Axisymmetric Geometry" \
        -variable gd(axiFlag) \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: SetAxiFlag $gd(axiFlag)"
            }; # end if
            SetAxiFlag $gd(axiFlag) 
            showStatusMsg "Axisymmetric flag altered."
        }
    $m1 add command -label "Set Gamma..." \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: displayGamma"
            }; # end if
            showStatusMsg "Start Gamma dialog."
            displayGamma 
        }

    #----------------------------------------------------
    set m1 [menu $m0.plot -tearoff 0 -title "Plot Menu"]
    $m0 add cascade -label Plot -menu $m1
    $m1 add command -label "Refresh Display" \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: nodelist_Clear; refreshDisplay"
            }; # end if
            nodelist_Clear; refreshDisplay 
        }
    $m1 add command -label "Generate PS File..." \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: generatePSFile"
            }; # end if
            showStatusMsg "Generating a postscript file..."
            generatePSFile 
            showStatusMsg "Generating a postscript file...Done."
        }

    $m1 add separator
    $m1 add command -label "Set New X,Y Ranges..." \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: showRangeDialog"
            }; # end if
            showStatusMsg "Adjust display ranges..."
            showRangeDialog 
            showStatusMsg "Adjust display ranges...Done."
        }
    $m1 add command -label "Set Default X,Y Ranges" \
        -command { 
            if { $gd(echoCommands) == 1 } {
                puts "Command: setXYRanges; setXYTics; refreshDisplay"
            }; # end if
            setXYRanges; setXYTics; refreshDisplay 
            showStatusMsg "Have set default display ranges."
        }
    $m1 add command -label "Zoom to Selected Window" \
        -command { startZoomWindow }
    $m1 add command -label "Zoom to Include All" \
        -command { 
            set bb [findNodeBoundingBox] 
            puts "Bounding box found to be: $bb"
            eval setXYRanges $bb
            refreshDisplay
        } 
    $m1 add check -label "Same X,Y Scales" \
        -variable gd(sameScales) -command { setXYScales; refreshDisplay }

    $m1 add separator
    $m1 add check -label "Show Node Numbers" \
        -variable gd(showNodeNumbers) -command { refreshDisplay }
    $m1 add check -label "Show Characteristic Mesh" \
        -variable gd(showCharMesh) -command { refreshDisplay }
    $m1 add check -label "Show Streamlines" \
        -variable gd(showStreamlines) -command { refreshDisplay }

    #----------------------------------------------------
    set m1 [menu $m0.help -tearoff 0]
    $m0 add cascade -label Help -menu $m1
    $m1 add command -label "About..." \
        -command "showAboutBox"
    $m1 add command -label "General Help" \
        -command {
            showStatusMsg "Starting browser with help files..."
            showGeneralHelp
            showStatusMsg "Starting browser with help files...Done."
        }

}; # end proc initMenu

#---------------------------------------------------------

proc closeMenu {} {
    #@proc
    destroy .menubar
}; # end proc closeMenu

# --------------------------------------------------------

# Some dialog windows...

proc showAboutBox {} {
    #@proc
    set    msg "IMOC 0.1.4\n"
    append msg "Interactive Method of Characteristics\n"
    append msg "Copyright (C) 2000, 2001\n"
    append msg "P.Jacobs\n"
    append msg "Department of Mechanical Engineering\n"
    append msg "The University of Queensland\n"
    append msg "http://www.mech.uq.edu.au/staff/jacobs/cfcfd/\n"
    tk_messageBox -type ok -title "About IMOC" \
        -message $msg -icon info
}; # end proc showAboutBox

# --------------------------------------------------------

proc showGeneralHelp {} {
    #@proc
    #@doc
    # Start up the HTML viewer.
    #@end
    global gd
    set helpFile [file join $gd(IMOC_HOME) doc index.html]
    set helpFile [file nativename $helpFile]
    exec $gd(htmlViewer) $helpFile &
}; # end proc showGeneralHelp

# --------------------------------------------------------

proc showRangeDialog {} {
    #@proc
    if [winfo exists .rangeDialog] {
        wm deiconify .rangeDialog
        raise .rangeDialog
    } else {
        makeRangeDialog 
    }; # end if
}; # end proc

# --------------------------------------------------------

proc makeRangeDialog {} {
    #@proc
    # Set up the separate window for entering the range values
    global gd 

    # Get current ranges for display
    set gd(display_xmin) $gd(xmin)
    set gd(display_xmax) $gd(xmax)
    set gd(display_dx)   $gd(dx)
    set gd(display_ymin) $gd(ymin)
    set gd(display_ymax) $gd(ymax)
    set gd(display_dy)   $gd(dy)

    toplevel .rangeDialog
    wm title .rangeDialog "Set Plotting Ranges"
    set rF [frame .rangeDialog.rangeFrame]
    set bF [frame .rangeDialog.buttonFrame]

    set lXmin [label $rF.labelXmin -text "Xmin"]
    set lXmax [label $rF.labelXmax -text "Xmax"]
    set ldX   [label $rF.labeldX -text "dX"]
    set lYmin [label $rF.labelYmin -text "Ymin"]
    set lYmax [label $rF.labelYmax -text "Ymax"]
    set ldY   [label $rF.labeldY -text "dY"]
    set eXmin [entry $rF.entryXmin -textvariable gd(display_xmin) \
        -relief sunken -bg white -width 10]
    set eXmax [entry $rF.entryXmax -textvariable gd(display_xmax) \
        -relief sunken -bg white -width 10]
    set edX   [entry $rF.entrydX -textvariable gd(display_dx) \
        -relief sunken -bg white -width 10]
    set eYmin [entry $rF.entryYmin -textvariable gd(display_ymin) \
        -relief sunken -bg white -width 10]
    set eYmax [entry $rF.entryYmax -textvariable gd(display_ymax) \
        -relief sunken -bg white -width 10]
    set edY   [entry $rF.entrydY -textvariable gd(display_dy) \
        -relief sunken -bg white -width 10]
    grid $lXmin $eXmin $lXmax $eXmax $ldX $edX
    grid $lYmin $eYmin $lYmax $eYmax $ldY $edY

    set bAccept [button $bF.accept -text Apply \
        -command {
            setXYRanges $gd(display_xmin) $gd(display_ymin) \
                $gd(display_xmax) $gd(display_ymax) 
            setXYTics $gd(display_dx) $gd(display_dy) 
            refreshDisplay; } ]
    set bClear [button $bF.clear -text "Clear Edits" \
        -command {
            set gd(display_xmin) $gd(xmin); set gd(display_xmax) $gd(xmax); 
            set gd(display_ymin) $gd(ymin); set gd(display_ymax) $gd(ymax); 
            set gd(display_dx)   $gd(dx);   set gd(display_dy)   $gd(dy); } ]
    set bDismiss [button $bF.dismiss -text Minimize \
        -command {wm iconify .rangeDialog} ]
    pack $bAccept $bClear $bDismiss -side left -padx 4 -pady 4

    pack $rF $bF
}; # end proc rangeDialog

# --------------------------------------------------------

proc displayGamma {} {
    #@proc
    # Set up the separate window for 
    # displaying ratio of specific heats
    global gd

    # Get current ranges for display
    set gd(display_gam) $gd(gamma)
    if { $gd(debug) == 1 } {
        puts "displayGamma: on entry, gamma = $gd(display_gam)"
    }; # end if

    toplevel .gammaDialog 
    wm title .gammaDialog "Set Ratio of Specific Heats"
    set gF [frame .gammaDialog.gammaFrame]
    set bF [frame .gammaDialog.buttonFrame]

    set lG [label $gF.labelG -text "Gamma"]
    set eG [entry $gF.entryG -textvariable gd(display_gam) \
        -relief sunken -bg white]
    grid $lG $eG

    set bAccept [button $bF.accept -text Accept \
        -command { set gd(gamma) $gd(display_gam);
                   SetGamma $gd(gamma);
                   destroy .gammaDialog } ]
    set bClear [button $bF.clear -text "Clear Edits" \
        -command { set gd(display_gam) $gd(gamma) } ]
    set bDismiss [button $bF.dismiss -text Dismiss \
        -command {destroy .gammaDialog} ]
    pack $bAccept $bClear $bDismiss -side left -padx 4 -pady 4

    pack $gF $bF
}; # end proc displayGamma

# --------------------------------------------------------

proc displayStreamlineStepSize {} {
    #@proc
    # Set up the separate window for displaying the size of the step
    # to be used for extending the streamline via Shepard interpolation
    global gd
    global stepsize

    # Get current ranges for display
    set stepsize $gd(StreamStepSize)
    if { $gd(debug) == 1 } {
        puts "displayStreamlineStepSize: on entry, stepsize = $stepsize"
    }; # end if

    toplevel .sssDialog 
    wm title .sssDialog "Set Streamline Step Size"
    set sF [frame .sssDialog.sssFrame]
    set bF [frame .sssDialog.buttonFrame]

    set lS [label $sF.labelS -text "Step Size"]
    set eS [entry $sF.entryS -textvariable stepsize \
        -relief sunken -bg white]
    grid $lS $eS

    set bAccept [button $bF.accept -text Accept \
        -command { set gd(StreamStepSize) $stepsize;
                   destroy .sssDialog } ]
    set bClear [button $bF.clear -text "Clear Edits" \
        -command { set stepsize $gd(StreamStepSize) } ]
    set bDismiss [button $bF.dismiss -text Dismiss \
        -command {destroy .sssDialog} ]
    pack $bAccept $bClear $bDismiss -side left -padx 4 -pady 4

    pack $sF $bF
}; # end proc displayStreamlineStepSize

# --------------------------------------------------------

proc loadNodeData {} {
    #@proc
    global gd

    set inputFileName \
        [tk_getOpenFile -title "Load Node Data"]
    set gd(nodeFile) [file tail $inputFileName]
    if { [string compare $gd(nodeFile) ""] == 0 } { return }
    set gd(workDir) [file dirname $inputFileName]
    cd $gd(workDir)
    if { $gd(debug) == 1 } {
        puts "loadNodeData from file: $gd(nodeFile)"
    }; # end if
    LoadNodes $gd(nodeFile)

}; # end proc loadNodeData

# --------------------------------------------------------

proc loadWall_gui { iw } {
    #@proc
    global gd

    set inputFileName \
        [tk_getOpenFile -title "Load Wall $iw"]
    set gd(WallFile) [file tail $inputFileName]
    if { [string compare $gd(WallFile) ""] == 0 } { return }
    set gd(workDir) [file dirname $inputFileName]
    cd $gd(workDir)

    if { $gd(debug) == 1 } {
        puts "Load Wall $iw data from file: $gd(Wall0File)"
    }; # end if
    LoadWall $iw $gd(WallFile)

}; # end proc loadWall_gui

# --------------------------------------------------------

proc saveNodeData {} {
    #@proc
    global gd

    set inputFileName \
        [tk_getSaveFile -title "Save Node Data" \
            -initialfile $gd(nodeFile) ]
    set gd(nodeFile) [file tail $inputFileName]
    if { [string compare $gd(nodeFile) ""] == 0 } { return }
    set gd(workDir) [file dirname $inputFileName]
    cd $gd(workDir)
    puts "Selected node file: $gd(nodeFile)"

    if { $gd(debug) == 1 } {
        puts "saveNodeData to file: $gd(nodeFile)"
    }; # end if
    SaveNodes $gd(nodeFile)

}; # end proc saveNodeData

# --------------------------------------------------------

proc saveWall_gui { iw } {
    #@proc
    global gd

    set inputFileName \
        [tk_getSaveFile -title "Save Wall $iw" \
            -initialfile $gd(WallFile)]
    set gd(WallFile) [file tail $inputFileName]
    if { [string compare $gd(WallFile) ""] == 0 } { return }
    set gd(workDir) [file dirname $inputFileName]
    cd $gd(workDir)
    puts "Selected file: $gd(WallFile)"

    if { $gd(debug) == 1 } {
        puts "Save Wall $iw data to file: $gd(WallFile)"
    }; # end if
    SaveWall $iw $gd(WallFile)

}; # end proc saveWall_gui

# --------------------------------------------------------

proc editNodeData { {id -1} } {
    #@proc
    global gd
    global new_property_value old_property_value
    global property_list edit_node_id
    global mach_prandtl_option

    set mach_prandtl_option setMach

    if { $id < 0 } {
        set edit_node_id [nodelist_GetNodeId 0]
    } else {
        set edit_node_id $id
    }; # end if
    if { $edit_node_id < 0 } {
        # We still don't have a valid node id.
        return
    }; # end if
    if { $gd(debug) == 1 } {
        puts "Edit the properties of node $edit_node_id"
        puts "Node data: [GetNodeData $edit_node_id]"
    }; # end if

    set node_data [split [GetNodeData $edit_node_id]]
    set property_list {}
    foreach {property value} $node_data {
        lappend property_list $property
        set new_property_value($property) $value
        set old_property_value($property) $value
    }; # end foreach

    if { [winfo exists .propertyDialog] == 0 } {
        # Create the widget
        toplevel .propertyDialog
        set pF [frame .propertyDialog.propertyFrame]
        set oF [frame .propertyDialog.optionFrame]
        set bF [frame .propertyDialog.buttonFrame]

        foreach property $property_list {
            set label [label $pF.label$property -text "$property"]
            set entry [entry $pF.entry$property \
                -textvariable new_property_value($property) \
                -relief sunken -bg white -width 15]
            grid $label $entry -sticky w
        }; # end foreach

        set MPMLabel [label $oF.label -text Set: ]
        set MFromNu [radiobutton $oF.mFromNuButton \
            -text "Mach Number" \
            -variable mach_prandtl_option -value setMach ]
        set NuFromM [radiobutton $oF.nuFromMButton \
            -text "PM Function" \
            -variable mach_prandtl_option -value setNu ]
        pack $MPMLabel $MFromNu $NuFromM -side left

        set bAccept [button $bF.accept -text Apply \
            -command {
                # Only one of Mach number and Prandtl-Meyer function can be set.
                # Remove the other from the property_list that is to be used
                # in copying back the new values.
                if {[string compare $mach_prandtl_option setMach] == 0} {
                    set pl_i [lsearch $property_list Nu]
                } else {
                    set pl_i [lsearch $property_list Mach]
                }; # end if
                # puts "Omit item $pl_i which is [lindex $property_list $pl_i]"
                set reduced_property_list [lreplace $property_list $pl_i $pl_i]
                foreach property $reduced_property_list {
                    SetNodeData $edit_node_id $property $new_property_value($property)
                }; # end foreach
                # Now, make the text entries for Nu and Mach consistent.
                if {[string compare $mach_prandtl_option setMach] == 0} {
                    set new_property_value(Nu) [GetNodeDataValue $edit_node_id Nu]
                } else {
                    set new_property_value(Mach) [GetNodeDataValue $edit_node_id Mach]
                }; # end if
                refreshDisplay
            } ]
        set bClear [button $bF.clear -text "Clear Edits" \
            -command {
                foreach property $property_list {
                    set new_property_value($property) $old_property_value($property)
                }; # end foreach
            } ]
        set bDismiss [button $bF.dismiss -text Dismiss \
            -command {destroy .propertyDialog} ]
        pack $bAccept $bClear $bDismiss -side left -padx 4 -pady 4

        pack $pF $oF $bF
    } else {
        # Update the widget for the newly selected node
        raise .propertyDialog
    }; # end if
    wm title .propertyDialog "Properties for Node $edit_node_id"
}; # end proc

# --------------------------------------------------------

proc editWallData { iw } {
    #@proc
    global gd
    global ew_iw MAX_N new_N new_X new_Y old_N old_X old_Y
    set MAX_N 10;  # Maximum number of points in a wall

    if { $gd(debug) == 0 } {
        puts "Edit the data points for wall $iw"
    }; # end if

    set ew_iw $iw; # put it into the global namespace for the Accept command

    for {set ip 0} {$ip < $MAX_N} {incr ip} {
        set old_X($ip) ""
        set new_X($ip) ""
        set old_Y($ip) ""
        set new_Y($ip) ""
    }; # end for

    set old_N [WallGetNumberOfPoints $iw]
    set new_N $old_N
    for {set ip 0} {$ip < $old_N} {incr ip} {
        set old_X($ip) [WallGetPointX $iw $ip]
        set new_X($ip) $old_X($ip)
        set old_Y($ip) [WallGetPointY $iw $ip]
        set new_Y($ip) $old_Y($ip)
    }; # end for

    if { [winfo exists .wallDialog] == 0 } {
        # Create the widget
        toplevel .wallDialog
        set nF [frame .wallDialog.numberFrame]
        set pF [frame .wallDialog.pointFrame]
        set bF [frame .wallDialog.buttonFrame]

        set labelN [label $nF.labelN -text "Number of points"]
        set entryN [entry $nF.entryN \
            -textvariable new_N \
            -relief sunken -bg white -width 8]
        pack $labelN $entryN -side left

        for {set ip 0} {$ip < $MAX_N} {incr ip} {
            set labelX [label $pF.labelX$ip -text "X($ip)"]
            set entryX [entry $pF.entryX$ip \
                -textvariable new_X($ip) \
                -relief sunken -bg white -width 8]
            set labelY [label $pF.labelY$ip -text "Y($ip)"]
            set entryY [entry $pF.entryY$ip \
                -textvariable new_Y($ip) \
                -relief sunken -bg white -width 8]
            grid $labelX $entryX $labelY $entryY -sticky w
        }; # end for

        set bAccept [button $bF.accept -text Apply \
            -command {
                WallDeletePoints $ew_iw
                for {set ip 0} {$ip < $new_N} {incr ip} {
                    WallAddPoint $ew_iw $new_X($ip) $new_Y($ip)
                }; # end for
                refreshDisplay
            } ]
        set bClear [button $bF.clear -text "Clear Edits" \
            -command {
                set new_N $old_N
                for {set ip 0} {$ip < $MAX_N} {incr ip} {
                    set new_X($ip) $old_X($ip)
                    set new_Y($ip) $old_Y($ip)
                }; # end for
            } ]
        set bDismiss [button $bF.dismiss -text Dismiss \
            -command {destroy .wallDialog} ]
        pack $bAccept $bClear $bDismiss -side left -padx 4 -pady 4

        pack $nF $pF $bF
    } else {
        # Update the widget for the newly selected wall
        raise .wallDialog
    }; # end if
    wm title .wallDialog "Points Defining Wall $iw"
}; # end proc

# --------------------------------------------------------

proc sourceUserScript {} {
    #@proc
    global gd

    set inputFileName \
        [tk_getOpenFile -title "Select Script File"]
    set gd(scriptFile) [file tail $inputFileName]
    if { [string compare $gd(scriptFile) ""] == 0 } {
        puts "sourceUserScript: did not source file."
        return
    }; # end if
    set gd(workDir) [file dirname $inputFileName]
    cd $gd(workDir)

    if { $gd(debug) == 1 } {
        puts "Source user script from file: $gd(scriptFile)"
    }; # end if
    uplevel #0 "source $gd(scriptFile)"

}; # end proc sourceUserScript

# --------------------------------------------------------

proc generatePSFile {} {
    #@proc
    global gd

    set inputFileName \
        [tk_getSaveFile -title "Generate Postscript File" \
            -initialfile $gd(plotFile) ]
    set gd(plotFile) [file tail $inputFileName]
    if { [string compare $gd(plotFile) ""] == 0 } {
        puts "generatePSFile: did not generate file."
        return
    }; # end if
    set gd(workDir) [file dirname $inputFileName]
    cd $gd(workDir)
    puts "Postscript plot file: $gd(plotFile)"

    if { $gd(debug) == 1 } {
        puts "Generate postscript plot into file: $gd(plotFile)"
    }; # end if
    plotPSFile $gd(plotFile)

}; # end proc generatePSFile

# --------------------------------------------------------
