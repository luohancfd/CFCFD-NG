# td_br_gui_procs.tcl
# Procedures for GUI interaction

puts "Start of td_br_gui_procs.tcl"


proc myBusy { action widget } {
    # Use either the blt::busy widget or cover over the toplevel manually.
    #     action is expected to be either "hold" or "forget"
    #     widget is usually the toplevel "."
    if { $::td(useBltBusy) == 1 } {
        rbc::busy $action $widget
    } else {
	if { [string equal $action hold] } {
	    # Put up a transparent toplevel and give it focus.
	    if { [winfo exists .busyDialog] } {
		raise .busyDialog
	    } else {
		toplevel .busyDialog -background {} -cursor watch
	    }
	    wm geometry .busyDialog [winfo geometry $widget]
	    focus .busyDialog
	} else {
	    # Assume that we want to release the busyDialog
	    wm withdraw .busyDialog
	    focus $widget
	}
    }
    return
}; # end proc myBusy


proc setMouseBindings { {option pickTimeRange} } {
    # option can be one of pickTimeRange (default)
    #                      pickPoint
    #                      pickPointWindowCoord

    bind $::mygraph <Motion> {
        displayCoordinates [%W axis invtransform x %x] \
                           [%W axis invtransform y %y] \
                           [%W axis invtransform y2 %y]
    }

    if { [string equal $option pickTimeRange] } {
        bind $::mygraph <ButtonPress-1> {
            setLeftMarker [%W axis invtransform x %x]
        }
        bind $::mygraph <ButtonPress-3> {
            setRightMarker [%W axis invtransform x %x]
        }
        bind $::mygraph <Control-ButtonPress-1> {
            setTopMarker [%W axis invtransform y %y] [%W axis invtransform y2 %y]
        }
        bind $::mygraph <Control-ButtonPress-3> {
            setBotMarker [%W axis invtransform y %y] [%W axis invtransform y2 %y]
        }
    }; # end if

    if { [string equal $option pickPoint] } {
        bind $::mygraph <ButtonPress-1> {
            recordXYPoint [%W axis invtransform x %x] \
                          [%W axis invtransform y %y] \
                          [%W axis invtransform y2 %y]
        }
        bind $::mygraph <ButtonPress-3> {
            recordXYPoint [%W axis invtransform x %x] \
                          [%W axis invtransform y %y] \
                          [%W axis invtransform y2 %y]
        }
    }; # end if

    if { [string equal $option pickPointWindowCoord] } {
        # We want the raw window coordinates.
        bind $::mygraph <ButtonPress-1> {
            recordXYPoint %x %y %y
        }
        bind $::mygraph <ButtonPress-3> {
            recordXYPoint %x %y %y
        }
    }; # end if

}; # end proc
puts "End proc setMouseBindings"


proc recordXYPoint { x y y2 } {
    global recordXYPointStatus
    global recordXYPointX recordXYPointY
    global annotationText annotationAnchor

    if { [string equal $recordXYPointStatus waitingForTextLocation] } {
        # puts "Write string at ($x, $y)"
        getAnnotationText
        if { ![string equal $annotationText ""] } {
            $::mygraph marker create text -text "$annotationText" \
                -coords [list $x $y] -anchor $annotationAnchor \
                -justify left -font $::td(plotFont)
            # Use the text only once
            set annotationText ""
            set annotationAnchor w
        }; # end if
        set recordXYPointStatus none
    } elseif { [string equal $recordXYPointStatus waitingForFirstPoint] } {
        # puts "First point at ($x, $y)"
        set recordXYPointX $x
        set recordXYPointY $y
        set recordXYPointStatus waitingForSecondPoint
        updateStatusMessage "Waiting for second point."
    } elseif { [string equal $recordXYPointStatus waitingForSecondPoint] } {
        # puts "Second point at ($x, $y)"
        $::mygraph marker create line \
            -coords [list $recordXYPointX $recordXYPointY $x $y]
        set recordXYPointX ""
        set recordXYPointY ""
        set recordXYPointStatus none
    } elseif { [string equal $recordXYPointStatus waitingForSampleTime] } {
        if { [string equal ::whichYAxis($::current) y] } {
            set sampleLevel y
        } else {
            set sampleLevel $y2
        }; # end if
        set sampleTime $x
        snapToReferenceTime $sampleTime
        set recordXYPointX ""
        set recordXYPointY ""
        set recordXYPointStatus none
    } elseif { [string equal $recordXYPointStatus waitingForSampleLevel] } {
        if { [string equal [set ::whichYAxis($::current)] y] } {
            set sampleLevel $y
        } else {
            set sampleLevel $y2
        }; # end if
        snapToReferenceValue $sampleLevel
        set recordXYPointX ""
        set recordXYPointY ""
        set recordXYPointStatus none
    } elseif { [string equal $recordXYPointStatus waitingToPickElement] } {
        # The coordinates are going to be window coordinates.
        # puts "Picking element:"
        if { [$::mygraph element closest $x $y result] == 1 } {
            parray result
            set elem [set result(name)]
            set ::current [set ::traceNumber($elem)]
            set ::td(nhalf) [set ::nhalf($::current)]
            computeAverageValue 
        };
        set recordXYPointX ""
        set recordXYPointY ""
        set recordXYPointStatus none
    } elseif { [string equal $recordXYPointStatus waitingForMarkerLevel] } {
        markTestTime $y
        set recordXYPointX ""
        set recordXYPointY ""
        set recordXYPointStatus none
    } else {
        set recordPointStatus none
    }; # end if

    if { [string equal $recordXYPointStatus none] } {
        # puts "Return to default mouse bindings"
        setMouseBindings pickTimeRange
        updateStatusMessage "Ready."
    }; # end if
}; # end proc
puts "End proc recordXYPoint"


proc getAnnotationText {} {
    global annotationText

    if { [winfo exists .annotationDialog] } {
        # Dialog widget already exists; just show it
        raise .annotationDialog
    } else {
        # Create the dialog widget and show it.
        toplevel .annotationDialog

        set hintLabel [label .annotationDialog.hintLabel \
            -text "Type note using \\n for new line within the text." ]

        set pF [frame .annotationDialog.labelFrame]
        set aF [frame .annotationDialog.anchorFrame]
        set bF [frame .annotationDialog.buttonFrame]

        set xlabel [label $pF.xlabel -text "Text:"]
        set xentry [entry $pF.xentry \
            -textvariable annotationText \
            -relief sunken -bg white -width 50]
        grid $xlabel $xentry -sticky w

        set al [label $aF.label -text "Anchor:"]
        set abw [radiobutton $aF.abw -text West -value w \
                 -variable annotationAnchor]
        set abc [radiobutton $aF.abc -text Center -value center \
                 -variable annotationAnchor]
        set abe [radiobutton $aF.abe -text East -value e \
                 -variable annotationAnchor]
        grid $al $abw $abc $abe -sticky w

        set bApply [button $bF.apply -text Apply \
            -command {
                destroy .annotationDialog
                # replace the backslash-n combinations
                # with actual newline characters
                regsub -all {\\n} $annotationText "\n" annotationText
                set annotationDialogFinished 1
            } ]
        set bDismiss [button $bF.dismiss -text "Close Window" \
            -command {
                set annotationText ""
                set annotationAnchor w
                destroy .annotationDialog
                set annotationDialogFinished 1
            } ]
        pack $bApply $bDismiss -side left -padx 4 -pady 4

        pack $hintLabel $pF $aF $bF
    }; # end if
    wm title .annotationDialog "Add note to graph."
    bind $xentry <KeyPress-Return> { 
        destroy .annotationDialog
        regsub -all {\\n} $annotationText "\n" annotationText
        set annotationDialogFinished 1
    }; # end bind
    focus $pF.xentry
    # Now that the dialog is up, wait for the interaction to be complete.
    vwait annotationDialogFinished
}; # end proc


proc updateTraceNumberList {} {
    global currentEntry traceNumberList
    global numberOfTraces traceIsPlotted

    # Rebuild list from scratch
    set traceNumberList {}
    for {set i 1} {$i <= $numberOfTraces} {incr i} {
        if { $traceIsPlotted($i) } {
            lappend traceNumberList $i
        }
    }; # end for

    $currentEntry configure -values $traceNumberList
}; # end proc


proc getTraceLabelText {} {
    global mygraph traceLabelText current traceName
    global elementName traceLabelDialogFinished

    # get the element name of the current trace
    set elementName [set traceName($current)]
    set traceLabelText [$mygraph element cget $elementName -label]

    if { [winfo exists .traceLabelDialog] } {
        # Dialog widget already exists; just show it
        raise .traceLabelDialog
    } else {
        # Create the dialog widget and show it.
        toplevel .traceLabelDialog

        set hintLabel [label .traceLabelDialog.hintLabel \
            -text "Type label text using \\n for new line within the text." ]

        set pF [frame .traceLabelDialog.labelFrame]
        set bF [frame .traceLabelDialog.buttonFrame]

        set xlabel [label $pF.xlabel -text "Text:"]
        set xentry [entry $pF.xentry \
            -textvariable traceLabelText \
            -relief sunken -bg white -width 50]
        grid $xlabel $xentry -sticky w

        set bApply [button $bF.apply -text Apply \
            -command {
                destroy .traceLabelDialog
                # replace the backslash-n combinations
                # with actual newline characters
                regsub -all {\\n} $traceLabelText "\n" traceLabelText
                if { ![string equal $elementName ""] } {
                    $mygraph element configure $elementName \
                        -label $traceLabelText
                }; # end if
                set traceLabelDialogFinished 1
            } ]
        set bDismiss [button $bF.dismiss -text "Close Window" \
            -command {
                set traceLabelText ""
                destroy .traceLabelDialog
                set traceLabelDialogFinished 1
            } ]
        pack $bApply $bDismiss -side left -padx 4 -pady 4

        pack $hintLabel $pF $bF
    }; # end if
    wm title .traceLabelDialog "Edit label for current trace."
    bind $xentry <KeyPress-Return> { 
        destroy .traceLabelDialog
        regsub -all {\\n} $traceLabelText "\n" traceLabelText
        if { ![string equal $elementName ""] } {
            $mygraph element configure $elementName \
                -label $traceLabelText
        }; # end if
        set traceLabelDialogFinished 1
    }; # end bind
    focus $pF.xentry
    # Now that the dialog is up, wait for the interaction to be complete.
    vwait traceLabelDialogFinished
}; # end proc


proc displayChannelNames {} {
    global td channelNameToNumber
    # This could take some time; put up a busy indicator
    myBusy hold .
    update

    # get the data and put it into the text widget
    if { [string equal $td(dataSource) directRead] } {
        set channel_id_text [readChannelNames]
    } else {
        set channel_id_text [fetchChannelNames]
    }; # end if
    $::headerText delete 1.0 end
    $::headerText insert end $channel_id_text
    
    # Now, extract the channel id numbers and names
    set channelList {}
    array set channelNameToNumber {}
    foreach line_of_text [split $channel_id_text "\n"] {
        # extract the first, white-space delimited token
        # ignoring white-space at the start of the line
	set lineItems [split [string trim $line_of_text] " "]
        set ch_id [string trim [lindex $lineItems 0]]
        set ch_name [string trim [lindex $lineItems 1]]
        if { [string length $ch_id] > 0 } {
            lappend channelList "$ch_id $ch_name"
	    if { [string length $ch_name] > 0 } {
		set channelNameToNumber($ch_name) $ch_id
	    }
        }
    }
    set ::td(channelList) [lsort $channelList]
    lappend ::td(channelList) "list"
    $::channelEntry configure -values $::td(channelList)
    $::normalizeEntry configure -values $::td(channelList)

    # Take down the busy indicator
    myBusy forget .
    update
    updateStatusMessage "Ready."
}; # end proc


proc displayChannelHeader {} {
    global td
    # get the saved text and put it into the text widget
    $::headerText delete 1.0 end
    $::headerText insert end $td(headerText)
}; # end proc


proc displayRunDescription {} {
    global td
    # This could take some time; put up a busy indicator
    myBusy hold .
    update

    # get the saved text and put it into the text widget
    if { [string equal $td(dataSource) directRead] } {
        set td(runDescriptionText) [readRunDescription]
    } else {
        set td(runDescriptionText) [fetchRunDescription]
    }; # end if
    $::headerText delete 1.0 end
    $::headerText insert end $td(runDescriptionText)

    # Take down the busy indicator
    myBusy forget .
    update
    updateStatusMessage "Ready."
}; # end proc


proc displayShotList {} {
    global td
    # This could take some time; put up a busy indicator
    myBusy hold .
    update

    # get the saved text and put it into the text widget
    if { [string equal $td(dataSource) directRead] } {
        updateStatusMessage "Collect shot list from local files."
        set shotListText [buildShotList]
    } else {
        updateStatusMessage "Fetch shot list from http server."
        set shotListText [fetchShotList]
    }; # end if
    $::headerText delete 1.0 end
    $::headerText insert end $shotListText

    # Now, extract the shot ids
    set shotList {}
    foreach lineOfText [split $shotListText "\n"] {
        if { [regexp "facility" $lineOfText] } {
            # ignore the comment line (at the start)
        } else {
            foreach id [split $lineOfText] { 
                set id [string trim $id]
                if { [string length $id] > 0 } {
                    lappend shotList $id
                }
            }; # end foreach id
        }; # end if
    }; # end foreach lineOfText
    set ::td(shotList) [lsort $shotList]
    lappend ::td(shotList) "list"
    $::shotEntry configure -values $::td(shotList)

    # Take down the busy indicator
    myBusy forget .
    update
    updateStatusMessage "Ready."
}; # end proc


proc showAboutBox {} {
    #@proc
    global td
    set    msg "td_browser.tcl $td(versionAndDate)\n"
    append msg "Tunnel Data Browser\n"
    append msg "P.Jacobs, 2002..2009\n"
    tk_messageBox -type ok -title "About td_browser" \
        -message $msg -icon info
}; # end proc showAboutBox


proc updateStatusMessage { {newText ""} } {
    set ::statusText $newText
    update idletasks
}; # end proc


proc displayHelpFile {} {
    global td

    set fp [open $td(helpFileName) "r"]
    set helpText [read $fp]
    close $fp

    append helpText "\nCurrent global data settings are:\n"
    foreach e [array names td] {
        append helpText "td($e) [set td($e)]\n"
    }; # end foreach

    $::headerText delete 1.0 end
    $::headerText insert end $helpText
}; # end proc


proc displayCoordinates { x y y2 } {
    $::xPosLabel configure -text "time: [format "%12.4g" $x]"
    $::yPosLabel configure -text "y-value: [format "%12.4g" $y]"
    $::y2PosLabel configure -text "y2-value: [format "%12.4g" $y2]"
}; # end proc


proc setLeftMarker { t } {
    global current tv
    if { [string equal $t ""] } return
    if { ![info exists ::plotXVec$current] } return
    ::plotXVec$current dup tv
    set dt [expr $tv(1) - $tv(0)]
    set tmin $tv(0)
    set tmax $tv(end)
    if { $t < $tmin } { set t $tmin }
    if { $t > $tmax } { set t $tmax }
    # Do not depend on fixed time increment between samples!
    # Search for samples in a small range and return the first.
    set tA [expr $t + 2.0 * $dt]
    set tIndex [lindex [tv search $t $tA] 0]
    set ::td(t1Index) $tIndex
    set ::td(t1) [set tv($tIndex)]
    if [catch {set ::tdiff [expr $::td(t2)-$::td(t1) ]} result ]  { set ::tdiff 0 }
    $::mygraph marker create line -name t1Mark \
        -coords "$::td(t1) -Inf $::td(t1) Inf"
    computeAverageValue
}; # end proc


proc setRightMarker { t } {
    global current tv
    if { [string equal $t ""] } return
    if { ![info exists ::plotXVec$current] } return
    ::plotXVec$current dup tv
    set dt [expr $tv(1) - $tv(0)]
    set tmin $tv(0)
    set tmax $tv(end)
    if { $t < $tmin } { set t $tmin }
    if { $t > $tmax } { set t $tmax }
    # Do not depend on fixed time increment between samples!
    # Search for samples in a small range and return the first.
    set tA [expr $t + 2.0 * $dt]
    set tIndex [lindex [tv search $t $tA] 0]
    set ::td(t2Index) $tIndex
    set ::td(t2) [set tv($tIndex)]
    if [catch {set ::tdiff [expr $::td(t2)-$::td(t1) ]} result ]  { set ::tdiff 0 }
    $::mygraph marker create line -name t2Mark \
        -coords "$::td(t2) -Inf $::td(t2) Inf"
    computeAverageValue
}; # end proc


proc setTopMarker { y y2 } {
    set ::td(yTop) $y
    set ::td(y2Top) $y2
    $::mygraph marker create line -name yTopMark -coords "-Inf $y Inf $y"
}; # end proc setTopMarker


proc setBotMarker { y y2 } {
    set ::td(yBot) $y
    set ::td(y2Bot) $y2
    $::mygraph marker create line -name yBotMark -coords "-Inf $y Inf $y"
}; # end proc setBotMarker


proc markTestTime { y } {
    # For Ivy, put in a test time marker at the cursor-specified
    # y-position.  
    # The value of y is in current plotting units for the y-axis.
    global td mygraph
    if { [string equal $td(t1) ""] || 
         [string equal $td(t2) ""] ||
	 [string equal $y ""] } {
        return; # do nothing
    } else {
        # Work out the size of the end bars in graph units.
        # Make them effectively 10 pixels high.
        set ypixels [$mygraph axis transform y $y]
        set y2pixels [expr $ypixels - $td(tickSizeInPixels)]
        set y2 [$mygraph axis invtransform y $y2pixels]
        $mygraph marker create line -name testTimeMark1 \
            -coords "$td(t1) $y $td(t2) $y"
        $mygraph marker create line -name testTimeMark2 \
            -coords "$td(t1) $y $td(t1) $y2"
        $mygraph marker create line -name testTimeMark3 \
            -coords "$td(t2) $y $td(t2) $y2"
    }; # end if
    return
}; # end proc markTestTime


proc deleteMarkers { {option all} } {
    # option is one of all
    #                  range
    global td
    if { [string equal $option all] } {
        foreach name [$::mygraph marker names] {
            $::mygraph marker delete $name
        }
    }; # end if
    if { [string equal $option range] } {
        $::mygraph marker delete t1Mark
        $::mygraph marker delete t2Mark
        $::mygraph marker delete yTopMark
        $::mygraph marker delete yBotMark
    }; # end if
    set td(t1Index) ""
    set td(t2Index) ""
    set td(t1) ""
    set td(t2) ""
    set td(yTop) ""
    set td(y2Top) ""
    set td(yBot) ""
    set td(y2Bot) ""
    set td(averageValue) ""
    set td(stddevValue) ""
}; # end proc


proc clearGraph {} {
    global td
    global current numberOfTraces
    deleteMarkers all
    foreach name [$::mygraph element names] {
        $::mygraph element delete $name
    }; # end foreach
    for {set i 1} {$i <= $numberOfTraces} {incr i} {
        catch { rbc::vector destroy ::plotXVec$i }
        catch { rbc::vector destroy ::plotYVec$i }
        set ::traceIsPlotted($i) 0
    };
    set numberOfTraces 0
    set current $numberOfTraces
    $::headerText delete 1.0 end
    set td(xAxisLabel) ""
    set td(yAxisLabel) ""
    set td(y2AxisLabel) ""
    $::mygraph axis configure y2 -hide true
    relabelGraphAxes
    set td(nhalf) 0
    updateTraceNumberList
}; # end if


proc incrementSetOffsetAndScaleFactor {} {
    global td current cumOS incOS
    global g_cum_toffset g_cum_tscale g_cum_yoffset g_cum_yscale
    global g_inc_toffset g_inc_tscale g_inc_yoffset g_inc_yscale

    if { $current == 0 } { return }

    set whichTrace $current
    set g_inc_toffset [lindex [set incOS($current)] 0]
    set g_inc_tscale  [lindex [set incOS($current)] 1]
    set g_inc_yoffset [lindex [set incOS($current)] 2]
    set g_inc_yscale  [lindex [set incOS($current)] 3]
    set g_cum_toffset [lindex [set cumOS($current)] 0]
    set g_cum_tscale  [lindex [set cumOS($current)] 1]
    set g_cum_yoffset [lindex [set cumOS($current)] 2]
    set g_cum_yscale  [lindex [set cumOS($current)] 3]

    if { [winfo exists .offsetDialog] } {
        # Dialog widget already exists; just show it
        raise .offsetDialog
    } else {
        # Create the dialog widget and show it.
        toplevel .offsetDialog
        set eF [frame .offsetDialog.entryFrame]
        set bF [frame .offsetDialog.buttonFrame]

        set ctoffsetlabel [label $eF.ctoffsetlabel -text "cumulative T offset"]
        set ctoffsetentry [entry $eF.ctoffsetentry \
            -textvariable g_cum_toffset \
            -relief sunken -bg gray -width 15 -state disabled]
        set ctscalelabel [label $eF.ctscalelabel -text "cumulative T scale"]
        set ctscaleentry [entry $eF.ctscaleentry \
            -textvariable g_cum_tscale \
            -relief sunken -bg gray -width 15 -state disabled]
        grid $ctoffsetlabel $ctoffsetentry $ctscalelabel $ctscaleentry -sticky w

        set toffsetlabel [label $eF.toffsetlabel -text "incremental T offset"]
        set toffsetentry [entry $eF.toffsetentry \
            -textvariable g_inc_toffset \
            -relief sunken -bg white -width 15]
        set tscalelabel [label $eF.tscalelabel -text "incremental T scale"]
        set tscaleentry [entry $eF.tscaleentry \
            -textvariable g_inc_tscale \
            -relief sunken -bg white -width 15]
        grid $toffsetlabel $toffsetentry $tscalelabel $tscaleentry -sticky w

        set cyoffsetlabel [label $eF.cyoffsetlabel -text "cumulative Y offset"]
        set cyoffsetentry [entry $eF.cyoffsetentry \
            -textvariable g_cum_yoffset \
            -relief sunken -bg gray -width 15 -state disabled]
        set cyscalelabel [label $eF.cyscalelabel -text "cumulative Y scale"]
        set cyscaleentry [entry $eF.cyscaleentry \
            -textvariable g_cum_yscale \
            -relief sunken -bg gray -width 15 -state disabled] 
        grid $cyoffsetlabel $cyoffsetentry $cyscalelabel $cyscaleentry -sticky w

        set yoffsetlabel [label $eF.yoffsetlabel -text "incremental Y offset"]
        set yoffsetentry [entry $eF.yoffsetentry \
            -textvariable g_inc_yoffset \
            -relief sunken -bg white -width 15]
        set yscalelabel [label $eF.yscalelabel -text "incremental Y scale"]
        set yscaleentry [entry $eF.yscaleentry \
            -textvariable g_inc_yscale \
            -relief sunken -bg white -width 15]
        grid $yoffsetlabel $yoffsetentry $yscalelabel $yscaleentry -sticky w

        set bApply [button $bF.apply -text Apply \
            -command {
                set incOS($current) \
                    [list $g_inc_toffset $g_inc_tscale \
                          $g_inc_yoffset $g_inc_yscale] 
                applyOffsetAndScale incremental $current 
                set g_cum_toffset [lindex [set cumOS($current)] 0]
                set g_cum_tscale  [lindex [set cumOS($current)] 1]
                set g_cum_yoffset [lindex [set cumOS($current)] 2]
                set g_cum_yscale  [lindex [set cumOS($current)] 3]
            } ]
        set bDefault [button $bF.default -text "Default Values" \
            -command { 
                resetOffsetAndScale incremental $current
                set g_inc_toffset [lindex [set incOS($current)] 0]
                set g_inc_tscale  [lindex [set incOS($current)] 1]
                set g_inc_yoffset [lindex [set incOS($current)] 2]
                set g_inc_yscale  [lindex [set incOS($current)] 3]
            } ]
        set bDismiss [button $bF.dismiss -text "Close Window" \
            -command {destroy .offsetDialog} ]
        pack $bApply $bDefault $bDismiss -side left -padx 4 -pady 4

        pack $eF $bF
    }; # end if
    wm title .offsetDialog "Increment offset and scale for trace $current."
}; # end proc


proc getDashPattern { i } {
    global td
    if { $i == 1 || $td(useDashedLines) == 0 } {
        set dashlist ""
    } else {
        set dashlist [list 4 [expr ($i - 1) * 2]]
    }; # end if
    return $dashlist
}; # end proc


proc getTraceColour { i } {
    global td
    if { $td(useColours) } {
        set colourList [list darkblue magenta darkgreen red blue brown]
        set nc [llength $colourList]
        set j [expr ($i - 1) % $nc]
        if { $j < 0 } { set j 0 }
        set colour [lindex $colourList $j]
    } else {
        set colour black
    }; # end if
    return $colour
}; # end proc


proc addTraceToGraph {} {
    global td traceName traceNumber current whichYAxis traceIsPlotted

    ::tVec$current dup ::plotXVec$current
    ::vVec$current dup ::plotYVec$current
    filterData

    if { [info exists td(header.transducerName)] } {
        set tname $td(header.transducerName)
        set td(xAxisLabel) "time, $td(currentTimeUnits)"
        set yAxisLabel "$td(header.transducerName), $td(header.dataUnits)"
    } else {
        set tname ""
        set td(xAxisLabel) ""
        set yAxisLabel ""
    }; # end if
    if { [string equal $td(yAxis) y] } {
        set td(yAxisLabel) $yAxisLabel
        set whichYAxis($current) y
        set axisHint L
    } else {
        set td(y2AxisLabel) $yAxisLabel
        # By default, at atartup, the y2 axis is not displayed.
        $::mygraph axis configure y2 -hide false
        set whichYAxis($current) y2
        set axisHint R
    }; # end if
    # Have gone back to NOT having newline characters in the element name
    # so that the postscript file is generated OK.  The element names appear
    # as one-line comments in the postscript file.
    if { $td(addTraceLabels) == 1 } {
        set elementLabel \
            "Shot $td(shot), Ch $td(channel),\n$tname, ($current,$axisHint)"
    } else {
        set elementLabel ""
    }; # end if
    set elem [$::mygraph element create "$td(shot),$td(channel),$current" \
        -xdata ::plotXVec$current -ydata ::plotYVec$current -symbol "" \
        -dashes [getDashPattern $current] -color [getTraceColour $current] \
        -mapy $td(yAxis) -label $elementLabel]
    # Set up cross-reference tables
    set traceName($current) $elem
    set traceNumber($elem) $current
    set traceIsPlotted($current) 1
    puts "Trace $current added to graph as element $elem."
    updateTraceNumberList
}; # end proc


proc deleteTrace { traceIndex } {
    # Delete the spoecified trace from the graph.
    puts -nonewline "Delete trace $traceIndex"
    catch { $::mygraph element delete [set ::traceName($traceIndex)] } junk
    puts " Result: $junk"
    set ::traceIsPlotted($traceIndex) 0
    updateTraceNumberList
}; # end proc


proc refreshGraph { {axisOption both} } {
    global td
    zoomToIncludeAll $axisOption
    displayChannelHeader
    # deleteMarkers range;  # By default, we want to keep user annotations.
    setLeftMarker $td(t1)
    setRightMarker $td(t2)
    relabelGraphAxes
    updateTraceNumberList
}; # end proc


proc relabelGraphAxes {} {
    global td
    $::mygraph axis configure x -title $td(xAxisLabel)
    $::mygraph axis configure y -title $td(yAxisLabel)
    $::mygraph axis configure y2 -title $td(y2AxisLabel)
}; # end proc


proc manuallySetAxisLabels {} {
    global td
    global xAxisLabel_old yAxisLabel_old y2AxisLabel_old

    set xAxisLabel_old $td(xAxisLabel)
    set yAxisLabel_old $td(yAxisLabel)
    set y2AxisLabel_old $td(y2AxisLabel)

    if { [winfo exists .labelDialog] } {
        # Dialog widget already exists; just show it
        raise .labelDialog
    } else {
        # Create the dialog widget and show it.
        toplevel .labelDialog
        set pF [frame .labelDialog.labelFrame]
        set bF [frame .labelDialog.buttonFrame]

        set xlabel [label $pF.xlabel -text "x label"]
        set xentry [entry $pF.xentry \
            -textvariable td(xAxisLabel) \
            -relief sunken -bg white -width 35]
        grid $xlabel $xentry -sticky w
        set ylabel [label $pF.ylabel -text "y label"]
        set yentry [entry $pF.yentry \
            -textvariable td(yAxisLabel) \
            -relief sunken -bg white -width 35]
        grid $ylabel $yentry -sticky w
        set y2label [label $pF.y2label -text "y2 label"]
        set y2entry [entry $pF.y2entry \
            -textvariable td(y2AxisLabel) \
            -relief sunken -bg white -width 35]
        grid $y2label $y2entry -sticky w

        set bApply [button $bF.apply -text Apply \
		-command { relabelGraphAxes } ]
        set bClear [button $bF.clear -text "Clear Edits" \
            -command {
                set td(xAxisLabel) $xAxisLabel_old
                set td(yAxisLabel) $yAxisLabel_old
                set td(y2AxisLabel) $y2AxisLabel_old
            } ]
        set bDismiss [button $bF.dismiss -text "Close Window" \
            -command {destroy .labelDialog} ]
        pack $bApply $bClear $bDismiss -side left -padx 4 -pady 4

        pack $pF $bF
    }; # end if
    wm title .labelDialog "Set the axis labels."
}; # end proc


proc zoomToIncludeAll { {axisOption both} } {
    # Input values: both
    #               onlyY
    #               none
    if { [string equal $axisOption both] } {
        $::mygraph axis configure x -min {} -max {}
        $::mygraph axis configure y -min {} -max {}
        $::mygraph axis configure y2 -min {} -max {}
    }; # end if
    if { [string equal $axisOption onlyY] } {
        $::mygraph axis configure y -min {} -max {}
        $::mygraph axis configure y2 -min {} -max {}
    }; # end if
    # Don't do anything for "none".
}; # end proc


proc zoomToSelectedRange {} {
    global td numberOfTraces whichYAxis

    set i1 $td(t1Index)
    set i2 $td(t2Index)
    if { [string length $i1] > 0 && \
         [string length $i2] > 0 } {
        if { $i2 < $i1 } {
            set ii $i1; set i1 $i2; set i2 $ii
        }; # end if

        # We are going to assume that the most-recently obtained
        # data set is usable.
        # Set the time range with this data.
        rbc::vector create ::subsetT
        ::plotXVec$numberOfTraces dup ::tmpV
        ::subsetT set $::tmpV($i1:$i2)
        set tmin $::subsetT(min)
        set tmax $::subsetT(max)

	set vmin {}
        set vmax {}
        set v2min {}
        set v2max {}
        rbc::vector create ::subsetV
	for {set j 1} {$j <= $numberOfTraces} {incr j} {
            # Check the limits for each trace, skipping those
            # that are not plotted.
            if { [set ::traceIsPlotted($j)] == 0 } { continue }
            ::plotYVec$j dup ::tmpV
            # We shall also assume that the same time-to-index
            # relation also holds.
            ::subsetV set $::tmpV($i1:$i2)
            set vmin_local $::subsetV(min)
            set vmax_local $::subsetV(max)
            if { [string equal [set whichYAxis($j)] y] } {
                if { $vmin == {} } { set vmin $vmin_local }
                if { $vmax == {} } { set vmax $vmax_local }
                if { $vmin_local < $vmin } { set vmin $vmin_local }
                if { $vmax_local > $vmax } { set vmax $vmax_local }
	    } else {
                if { $v2min == {} } { set v2min $vmin_local }
                if { $v2max == {}} { set v2max $vmax_local }
                if { $vmin_local < $v2min } { set v2min $vmin_local }
                if { $vmax_local > $v2max } { set v2max $vmax_local }
	    }; # end if
	}; # end for
        $::mygraph axis configure x -min $tmin -max $tmax
        $::mygraph axis configure y -min $vmin -max $vmax
        $::mygraph axis configure y2 -min $v2min -max $v2max
    }; # end if
    if { [string length $::td(yTop)] > 0 } {
	$::mygraph axis configure y -max $::td(yTop)
	$::mygraph axis configure y2 -max $::td(y2Top)
    }
    if { [string length $::td(yBot)] > 0 } {
	$::mygraph axis configure y -min $::td(yBot)
	$::mygraph axis configure y2 -min $::td(y2Bot)
    }
    deleteMarkers
}; # end if


proc manuallySetAxisRanges {} {
    global td
    global xmin_manual ymin_manual xmax_manual ymax_manual
    global y2min_manual y2max_manual
    global xmin_old ymin_old xmax_old ymax_old y2min_old y2max_old

    set xmin_manual [$::mygraph axis cget x -min]
    set xmax_manual [$::mygraph axis cget x -max]
    set ymin_manual [$::mygraph axis cget y -min]
    set ymax_manual [$::mygraph axis cget y -max]
    set y2min_manual [$::mygraph axis cget y2 -min]
    set y2max_manual [$::mygraph axis cget y2 -max]

    set xmin_old $xmin_manual
    set ymin_old $ymin_manual
    set y2min_old $y2min_manual
    set xmax_old $xmax_manual
    set ymax_old $ymax_manual
    set y2max_old $y2max_manual

    if { [winfo exists .rangeDialog] } {
        # Dialog widget already exists; just show it
        raise .rangeDialog
    } else {
        # Create the dialog widget and show it.
        toplevel .rangeDialog
        set eF [frame .rangeDialog.entryFrame]
        set bF [frame .rangeDialog.buttonFrame]

        set xminlabel [label $eF.xminlabel -text "xmin"]
        set xminentry [entry $eF.xminentry \
            -textvariable xmin_manual \
            -relief sunken -bg white -width 15]
        set xmaxlabel [label $eF.xmaxlabel -text "xmax"]
        set xmaxentry [entry $eF.xmaxentry \
            -textvariable xmax_manual \
            -relief sunken -bg white -width 15]
        grid $xminlabel $xminentry $xmaxlabel $xmaxentry -sticky w

        set yminlabel [label $eF.yminlabel -text "ymin"]
        set yminentry [entry $eF.yminentry \
            -textvariable ymin_manual \
            -relief sunken -bg white -width 15]
        set ymaxlabel [label $eF.ymaxlabel -text "ymax"]
        set ymaxentry [entry $eF.ymaxentry \
            -textvariable ymax_manual \
            -relief sunken -bg white -width 15]
        grid $yminlabel $yminentry $ymaxlabel $ymaxentry -sticky w

        set y2minlabel [label $eF.y2minlabel -text "y2min"]
        set y2minentry [entry $eF.y2minentry \
            -textvariable y2min_manual \
            -relief sunken -bg white -width 15]
        set y2maxlabel [label $eF.y2maxlabel -text "y2max"]
        set y2maxentry [entry $eF.y2maxentry \
            -textvariable y2max_manual \
            -relief sunken -bg white -width 15]
        grid $y2minlabel $y2minentry $y2maxlabel $y2maxentry -sticky w

        set bApply [button $bF.apply -text Apply \
            -command {
                $::mygraph axis configure x -min $xmin_manual -max $xmax_manual
                $::mygraph axis configure y -min $ymin_manual -max $ymax_manual
                $::mygraph axis configure y2 \
                    -min $y2min_manual -max $y2max_manual
            } ]
        set bClear [button $bF.clear -text "Clear Edits" \
            -command {
                set xmin_manual $xmin_old
                set xmax_manual $xmax_old
                set ymin_manual $ymin_old
                set ymax_manual $ymax_old
                set y2min_manual $y2min_old
                set y2max_manual $y2max_old
            } ]
        set bDismiss [button $bF.dismiss -text "Close Window" \
            -command {destroy .rangeDialog} ]
        pack $bApply $bClear $bDismiss -side left -padx 4 -pady 4

        pack $eF $bF
    }; # end if
    wm title .rangeDialog "Manually set the axis ranges."
}; # end proc


proc displayOptionsDialog {} {
    global td

    if { [winfo exists .optionsDialog] } {
        # Dialog widget already exists; just show it
        raise .optionsDialog
    } else {
        # Create the dialog widget and show it.
        toplevel .optionsDialog
        set nb [NoteBook .optionsDialog.notebook -height 200 -width 400]
        set bF [frame .optionsDialog.buttonFrame]

        set page1 [$nb insert end 1 -text "Data Source"]
        set rb1 [radiobutton $page1.rb1 -text "Use HTTP Server" \
            -variable ::td(dataSource) -value httpServer]
        set rb2 [radiobutton $page1.rb2 -text "Use Local Server" \
            -variable ::td(dataSource) -value localServer]
        set rb3 [radiobutton $page1.rb3 \
            -text "Read Local Data Files Directly" \
            -variable ::td(dataSource) -value directRead]
        pack $rb1 $rb2 $rb3 -anchor w

        set page2 [$nb insert end 2 -text "Local Files"]
        set label1 [label $page2.label1 -text "Root Directory:"]
        set entry1 [entry $page2.entry1 \
            -textvariable td(rootDir) \
            -relief sunken -bg white -width 25]
        grid $label1 $entry1 -sticky w
        set cb1 [checkbutton $page2.cb1 \
            -text "Include MetaData When Saving Selected Data" \
            -variable ::td(includeMetaDataWhenSaving)]
        grid $cb1 -sticky w -columnspan 2

        set page3 [$nb insert end 3 -text "Web Access"]
        set label1 [label $page3.label1 -text "Server Name/IP:"]
        set entry1 [entry $page3.entry1 \
            -textvariable td(serverIP) \
            -relief sunken -bg white -width 25]
        grid $label1 $entry1 -sticky w
        set label3 [label $page3.label3 -text "User Name:"]
        set entry3 [entry $page3.entry3 \
            -textvariable td(username) \
            -relief sunken -bg white -width 12]
        grid $label3 $entry3 -sticky w
        set label4 [label $page3.label4 -text "Password:"]
        set entry4 [entry $page3.entry4 \
            -textvariable td(password) \
            -relief sunken -bg white -width 12]
        grid $label4 $entry4 -sticky w
        set cb1 [checkbutton $page3.cb1 -text "Use Proxy with HTTP server" \
            -variable ::td(useProxy)]
        grid $cb1 -sticky w -columnspan 2
        set label2 [label $page3.label2 -text "Proxy Name/IP:"]
        set entry2 [entry $page3.entry2 \
            -textvariable td(proxyhost) \
            -relief sunken -bg white -width 25]
        grid $label2 $entry2 -sticky w

        set page4 [$nb insert end 4 -text "Graph"]
        set cb1 [checkbutton $page4.cb1 \
            -text "Use Colours to Distinguish Traces" \
            -variable ::td(useColours)]
        grid $cb1 -sticky w -columnspan 2
        set cb2 [checkbutton $page4.cb2 \
            -text "Use Dashed Lines to Distinguish Traces" \
            -variable ::td(useDashedLines)]
        grid $cb2 -sticky w -columnspan 2
        set cb3 [checkbutton $page4.cb3 \
            -text "Add traces with auto-generated labels" \
            -variable ::td(addTraceLabels)]
        grid $cb3 -sticky w -columnspan 2
        set label1 [label $page4.label1 -text "Font:"]
        set entry1 [entry $page4.entry1 \
            -textvariable td(plotFont) \
            -relief sunken -bg white -width 25]
        grid $label1 $entry1 -sticky w
        set label2 [label $page4.label2 -text "Tick Size for Test Time:"]
        set entry2 [entry $page4.entry2 \
            -textvariable td(tickSizeInPixels) \
            -relief sunken -bg white -width 12]
        grid $label2 $entry2 -sticky w

        set rb1 [radiobutton $page4.rb1 \
            -text "Rescale ONLY Y-axis when adding new trace" \
            -variable ::td(rescaleAxesOption) -value onlyY]
        set rb2 [radiobutton $page4.rb2 \
            -text "Rescale BOTH axes when adding new trace" \
            -variable ::td(rescaleAxesOption) -value both]
        set rb3 [radiobutton $page4.rb3 \
            -text "Do NOT scale axes when adding new trace" \
            -variable ::td(rescaleAxesOption) -value none -state disabled]
        grid $rb1 -sticky w -columnspan 2
        grid $rb2 -sticky w -columnspan 2
        grid $rb3 -sticky w -columnspan 2

        set page5 [$nb insert end 5 -text "Data"]
        set label1 [label $page5.label1 -text "Default nhalf:"]
        set entry1 [entry $page5.entry1 \
            -textvariable td(default.nhalf) \
            -relief sunken -bg white -width 10]
        grid $label1 $entry1 -sticky w
        set label2 [label $page5.label2 -text "Reference Time:"]
        set entry2 [entry $page5.entry2 \
            -textvariable td(referenceTime) \
            -relief sunken -bg white -width 15]
        grid $label2 $entry2 -sticky w
        set label3 [label $page5.label3 -text "Reference Value:"]
        set entry3 [entry $page5.entry3 \
            -textvariable td(referenceValue) \
            -relief sunken -bg white -width 15]
        grid $label3 $entry3 -sticky w

        NoteBook::compute_size $nb
        $nb raise 1

        set bDismiss [button $bF.dismiss -text "Close Window" \
            -command {destroy .optionsDialog} ]
        pack $bDismiss -side left -padx 4 -pady 4

        pack $nb $bF
    }; # end if
    wm title .optionsDialog "td_browser configuration options."
}; # end proc

proc displayheatDialog {} {
    global calcfrom
    global calcto
    global thermprod
    global transens
    global ambvolt
    set calcfrom 3
    set calcto 1
    set thermprod 1540
    set transens 1e-3
    set ambvolt 1


    if { [winfo exists .heatDialog] } {
        # Dialog widget already exists; just show it
        raise .heatDialog
    } else {
        # Create the dialog widget and show it.
        toplevel .heatDialog
        
        #Create first set of Radiobuttons
        set fra1 [frame .heatDialog.fra1 -borderwidth "5" -relief "groove"]

        set lab1 [label $fra1.lab1 -text "Select to calculate from:"]
        grid $lab1
        set rb1 [radiobutton $fra1.rb1 -text "Heat-transfer Rate"\
               -value 1 -variable calcfrom]
        set rb2 [radiobutton $fra1.rb2 -text "Surface Temperature"\
               -value 2 -variable calcfrom]
        set rb3 [radiobutton $fra1.rb3 -text "Raw EMF"\
               -value 3 -variable calcfrom]
        grid $rb1 -sticky w -columnspan 2
        grid $rb2 -sticky w -columnspan 2
        grid $rb3 -sticky w -columnspan 2

        pack $fra1 -anchor s -padx 4 -pady 4


        #Create second set of Radiobuttons
        set fra2 [frame .heatDialog.fra2 -borderwidth "5" -relief "groove"]

        set lab1 [label $fra2.lab1 -text "Select to calculate to:"]
        grid $lab1
        set rb1 [radiobutton $fra2.rb1 -text "Heat-transfer Rate"\
               -value 1 -variable calcto]
        set rb2 [radiobutton $fra2.rb2 -text "Surface Temperature"\
               -value 2 -variable calcto]
        grid $rb1 -sticky w -columnspan 2
        grid $rb2 -sticky w -columnspan 2
        
        pack $fra2 -anchor s -padx 4 -pady 4

        set lab2 [label .heatDialog.lab2 -text "Thermal Product:"]
        pack $lab2 -padx 4 -pady 4 -anchor s


        #Creat first variable entry section thermal product
        set fra3 [frame .heatDialog.fra3]        
         
        set ent1 [entry $fra3.ent1 -textvariable thermprod \
               -relief sunken -bg white -width 10] 
        set lab3 [label $fra3.lab3 -text "J/m^2/K/s^(1/2)"]       
        grid $ent1 $lab3
      
        pack $fra3 -padx 50 -pady 4 -anchor w


        #Creat second variable entry section tranducer sensitivity
  	  set lab4 [label .heatDialog.lab4 -text "Transducer Sensitivity:"]
        pack $lab4 -padx 4 -pady 4 -anchor s

        set fra4 [frame .heatDialog.fra4]        
         
        set ent2 [entry $fra4.ent2 -textvariable transens \
               -relief sunken -bg white -width 10] 
        set lab5 [label $fra4.lab5 -text ""]       
        grid $ent2 $lab5 

        pack $fra4 -padx 50 -pady 4 -anchor w

  
        #Creat third variable entry section ambient voltage
        set lab6 [label .heatDialog.lab6 -text "Transducer Gain:"]
        pack $lab6 -padx 4 -pady 4 -anchor s

        set fra5 [frame .heatDialog.fra5]        
         
        set ent3 [entry $fra5.ent3 -textvariable ambvolt \
               -relief sunken -bg white -width 10] 
        set lab7 [label $fra5.lab6 -text "V"]       
        grid $ent3 $lab7 
       
        pack $fra5 -padx 50 -pady 4 -anchor w


        #Create calculate and close buttons
        set fra6 [frame .heatDialog.fra6]

        set but1 [button $fra6.but1 -text "CALCULATE"\
                -command {compute_heat_transfer} ]
        set but2 [button $fra6.but2 -text "CLOSE" \
            -command {destroy .heatDialog} ]
        grid $but1 $but2 -sticky w -padx 4 -pady 4
        
        pack $fra6 -anchor s 
    

    }; # end if
    wm title .heatDialog "Heat-transfer Calculations."
}; # end proc


