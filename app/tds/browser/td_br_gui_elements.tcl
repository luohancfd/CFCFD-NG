# td_br_gui_elements.tcl
# Set up frames and widgets in the main part of the window.
#
# A frame to contain statusFrame1 and the graph widget.
# The other items (text window and status lines) will be stacked 
# vertically below this frame.
#
set upperFrame [frame .fu]

#
# A frame containing labels and entry widgets that indicate
# the current shot, channel and averaging.
#
set statusFrame1 [frame $upperFrame.f1]

set facilityLabel [label $statusFrame1.fl -text "Facility"]
set facilityEntry [ComboBox $statusFrame1.fe -textvariable td(facility) \
		       -width 10 -entrybg white -values $td(facilityList) ]
set shotLabel [label $statusFrame1.sl -text "Shot"]
set shotEntry [ComboBox $statusFrame1.se -textvariable td(shot) \
		   -width 10 -entrybg white -values $td(shotList) ]
$shotEntry bind <KeyPress-Return> { getFreshData }
set channelLabel [label $statusFrame1.cl -text "Channel"]
set channelEntry [ComboBox $statusFrame1.ce -textvariable td(channel) \
		      -width 10 -entrybg white -values $td(channelList) ]
$channelEntry bind <KeyPress-Return> { set td(yAxis) y; getFreshData }
set fetchButton [button $statusFrame1.fb -text "Fetch Data" \
    -command { set td(yAxis) y; getFreshData }]

set normalizeLabel [label $statusFrame1.norml -text "Norm. Chan."]
set normalizeEntry [ComboBox $statusFrame1.norme \
			-textvariable ::td(normalizeChannel) \
			-width 10 -entrybg white \
			-values $td(channelList) ]
$normalizeEntry bind <KeyPress-Return> { getNormalizingData }
set tshiftLabel [label $statusFrame1.delayl -text "Time Shift"]
set tshiftEntry [entry $statusFrame1.delaye -textvariable ::td(timeShift) \
    -width 10 -background white]
set normalizeButton [button $statusFrame1.normb -text "Set Norm Chan" \
    -command { getNormalizingData }]

set nhalfLabel [label $statusFrame1.nl -text "nhalf"]
set nhalfEntry [entry $statusFrame1.ne -textvariable ::td(nhalf) \
    -width 10 -background white]
bind $nhalfEntry <KeyPress-Return> { filterData; computeAverageValue }
set filterButton [button $statusFrame1.filtb -text "Apply Filter" \
    -command { filterData; computeAverageValue } ]

set currentLabel [label $statusFrame1.currl -text "Current Trace"]
set currentEntry [ComboBox $statusFrame1.curre -textvariable ::current \
		      -width 10 -entrybg white -values $traceNumberList ]
$currentEntry bind <KeyPress-Return> {
    if { $current < 0 } { set current 0 }
    if { $current > $numberOfTraces } { set current $numberOfTraces }
    set td(nhalf) [set nhalf($current)] 
    computeAverageValue 
}; # end bind

set heatButton [button $statusFrame1.heatb -text "Qdot Function" \
		    -command {displayheatDialog}  ]

set myrow 0
grid $facilityLabel -row $myrow -column 0 -sticky e
grid $facilityEntry -row $myrow -column 1 -sticky ew 
incr myrow
grid $shotLabel -row $myrow -column 0 -sticky e
grid $shotEntry -row $myrow -column 1 -sticky ew
incr myrow
grid $channelLabel -row $myrow -column 0 -sticky e
grid $channelEntry -row $myrow -column 1 -sticky ew
incr myrow
grid $fetchButton -row $myrow -column 1 -sticky ew -pady 5
incr myrow
grid $normalizeLabel -row $myrow -column 0 -sticky e
grid $normalizeEntry -row $myrow -column 1 -sticky ew
incr myrow
grid $tshiftLabel -row $myrow -column 0 -sticky e
grid $tshiftEntry -row $myrow -column 1 -sticky ew
incr myrow
grid $normalizeButton -row $myrow -column 1 -sticky ew -pady 5
incr myrow
grid $nhalfLabel -row $myrow -column 0 -sticky e
grid $nhalfEntry -row $myrow -column 1 -sticky ew
incr myrow
grid $filterButton -row $myrow -column 1 -sticky ew -pady 5
incr myrow
grid $currentLabel -row $myrow -column 0 -sticky e
grid $currentEntry -row $myrow -column 1 -sticky ew
incr myrow
grid $heatButton -row $myrow -column 1 -sticky ew -pady 5

#
# A status line indicating the current cursor position on
# the graph window.
#
set statusFrame2 [frame .f2 -borderwidth 3 -relief groove]
set cursorLabel [label $statusFrame2.c -text "Current Cursor Position: "]
set ::xPosLabel [label $statusFrame2.x -text "time:" -width 20]
set ::yPosLabel [label $statusFrame2.y -text "y-value:" -width 20]
set ::y2PosLabel [label $statusFrame2.y2 -text "y2-value:" -width 20]
pack $cursorLabel $::xPosLabel $::yPosLabel $::y2PosLabel -side left

#
# A status line indicating the marked range,
# together with its average and standard deviation.
#
set statusFrame3 [frame .f3 -borderwidth 3 -relief groove]
set t1Label [label $statusFrame3.t1l -text "Marked Range: t1:"]
set ::t1Entry [entry $statusFrame3.t1e -width 10 \
    -textvariable ::td(t1) -background white]
bind $::t1Entry <KeyPress-Return> { setLeftMarker [$::t1Entry get]}
set t2Label [label $statusFrame3.t2l -text "t2:"]
set ::t2Entry [entry $statusFrame3.t2e -width 10 \
    -textvariable ::td(t2) -background white]
bind $::t2Entry <KeyPress-Return> {setRightMarker [$::t2Entry get]}
if [catch {set ::tdiff [expr $::td(t2)-$::td(t1) ]} result ]  { set ::tdiff 0 } 
set tdiffLabel [label $statusFrame3.tdiffl -text "t2-t1:"]
set ::tdiffEntry [entry $statusFrame3.tdiffe -width 10 \
    -textvariable ::tdiff -background white]
set averageValueLabel [label $statusFrame3.avl -text "Average:"]
set averageValueEntry [entry $statusFrame3.ave -width 20 \
    -textvariable ::td(averageValue)]
set stddevValueLabel [label $statusFrame3.sdvl -text "StdDev:"]
set stddevValueEntry [entry $statusFrame3.sdve -width 20 \
    -textvariable ::td(stddevValue)]
pack $t1Label $t1Entry $t2Label $t2Entry \
    $tdiffLabel $tdiffEntry \
    $averageValueLabel $averageValueEntry \
    $stddevValueLabel $stddevValueEntry \
    -side left

#
# A status line giving some hint as to the current state.
#
set statusFrame4 [frame .f4 -borderwidth 3 -relief flat]
set statusLabel [label $statusFrame4.statl -text "Status :"]
set statusEntry [entry $statusFrame4.state -textvariable ::statusText ]
pack $statusLabel -side left
pack $statusEntry -side left -expand 1 -fill x

#
# A scrolling text window
#
set textFrame [frame .tf]
set ::headerText [text .tf.t -height $td(textHeight) -width 60 \
    -font $td(textFont) -wrap none \
    -yscrollcommand [list $textFrame.vsb set] ]
set textScrollBar [scrollbar .tf.vsb -orient vertical \
    -command {$::headerText yview} ]
pack $::headerText -side left -expand 1 -fill both
pack $textScrollBar -side left -fill y

#
# The plot window featuring the BLT graph widget
#
set graphFrame [frame $upperFrame.gf -borderwidth 3 -relief groove]
set ::mygraph [blt::graph $graphFrame.g -height $td(height) -width $td(width)]
$::mygraph axis configure x -titlefont $td(plotFont)
$::mygraph axis configure y -titlefont $td(plotFont)
$::mygraph axis configure y2 -titlefont $td(plotFont)
$::mygraph axis configure x -tickfont $td(plotFont)
$::mygraph axis configure y -tickfont $td(plotFont)
$::mygraph axis configure y2 -tickfont $td(plotFont)
$::mygraph legend configure -font $td(plotFont)
pack $::mygraph
setMouseBindings pickTimeRange

#
# Pack the major widgets into main window
#
pack $statusFrame1 -side left -fill y
pack $graphFrame -side left
pack $upperFrame 
pack $statusFrame3 -expand 1 -fill x
pack $textFrame -fill both -expand 1
pack $statusFrame2 -fill x -expand 1
pack $statusFrame4 -fill x -expand 1
wm title . "td_browser"

#
# General key bindings
#
event add <<FetchData1>> <KeyPress-F1>
bind . <<FetchData1>> { set td(yAxis) y; getFreshData }
event add <<FetchData2>> <Shift-KeyPress-F1>
event add <<FetchData2>> <Control-KeyPress-F1>
bind . <<FetchData2>> { set td(yAxis) y2; getFreshData }

event add <<FilterData>> <KeyPress-F2>
bind . <<FilterData>> { 
    set td(nhalf) $td(default.nhalf)
    filterData
    computeAverageValue 
}; # end bind
event add <<UnFilterData>> <Shift-KeyPress-F2>
event add <<UnFilterData>> <Control-KeyPress-F2>
bind . <<UnFilterData>> { 
    restoreOriginalTrace
    computeAverageValue 
}

event add <<ZoomIn>> <KeyPress-F3>
bind . <<ZoomIn>> { zoomToSelectedRange }
event add <<ZoomOut>> <Shift-KeyPress-F3>
event add <<ZoomOut>> <Control-KeyPress-F3>
bind . <<ZoomOut>> { zoomToIncludeAll }

bind . <KeyPress-F4> { deleteTrace $current }
event add <<ClearGraph>> <Shift-KeyPress-F4>
event add <<ClearGraph>> <Control-KeyPress-F4>
bind . <<ClearGraph>> { clearGraph }

bind . <KeyPress-F5> { displayChannelNames }

bind . <KeyPress-F6> { displayChannelHeader }

bind . <KeyPress-F7> { displayRunDescription }

bind . <KeyPress-F8> { displayShotList }

bind . <KeyPress-F10> { generateASCIIFiles }

