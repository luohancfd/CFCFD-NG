# td_br_menu.tcl
# Set up menu bar for td_browser.tcl.

menu .mb
. configure -menu .mb

puts "Start of File Menu"
menu .mb.file -tearoff 0
.mb.file add command -label "Save Selected Data to GNUPlot format..." \
    -command { saveSelectedData gnuplot }
.mb.file add command -label "Save Selected Data to CSV format..." \
    -command { saveSelectedData csv }
.mb.file add command -label "Save Displayed Text..." \
    -command { saveDisplayedText }
.mb.file add command -label "Save Postscript Graph..." \
    -command { printPostscriptGraph }
.mb.file add separator
.mb.file add command -label "Generate ASCII Files (defrock)" \
    -accelerator "<F10>" -command { generateASCIIFiles }
.mb.file add separator
.mb.file add command -label "Configure Options..." \
    -command { displayOptionsDialog }
.mb.file add separator
.mb.file add command -label "Exit" -command { tidyUp; exit }
.mb add cascade -label File -menu .mb.file
puts "End of File Menu"

menu .mb.data -tearoff 0
.mb.data add command -label "Pick Trace Using Cursor" \
    -command { 
        set recordXYPointStatus waitingToPickElement
        setMouseBindings pickPointWindowCoord
        updateStatusMessage "Use mouse to select trace."
    }
.mb.data add separator
.mb.data add command -label "Fetch Data (map to y-axis)" \
    -accelerator "<F1>" \
    -command { set td(yAxis) y; getFreshData }
.mb.data add command -label "Fetch Data (map to y2-axis)" \
    -accelerator "<Shift-F1>/<Control-F1>" \
    -command { set td(yAxis) y2; getFreshData }
.mb.data add command -label "Fetch Plain-Old-Data (map to y-axis)" \
    -command { set td(yAxis) y; readPODFile }
.mb.data add command -label "Fetch Plain-Old-Data (map to y2-axis)" \
    -command { set td(yAxis) y2; readPODFile }
.mb.data add separator
.mb.data add command -label "Filter Data with Default nhalf" \
    -accelerator "<F2>" \
    -command { 
        set td(nhalf) $td(default.nhalf);
        filterData; 
        computeAverageValue 
    }
.mb.data add command -label "Restore Original Data" \
    -accelerator "<Shift-F2>/<Control-F2>" \
    -command { restoreOriginalTrace; computeAverageValue }
.mb.data add separator
.mb.data add command -label "Shift Sample Point to Reference Time" \
    -command { 
        set recordXYPointStatus waitingForSampleTime
        setMouseBindings pickPoint
        updateStatusMessage "Use mouse to select sample time."
    }
.mb.data add command -label "Shift Sample Point to Reference Level" \
    -command { 
        set recordXYPointStatus waitingForSampleLevel
        setMouseBindings pickPoint
        updateStatusMessage "Use mouse to select sample level."
    }
.mb.data add command -label "Increment Offset and Scale Factor..." \
    -command { incrementSetOffsetAndScaleFactor }
.mb add cascade -label Data -menu .mb.data
puts "End of Data Menu"

menu .mb.graph -tearoff 0
.mb.graph add command -label "Zoom To Selected Range" \
    -accelerator "<F3>" \
    -command { zoomToSelectedRange }
.mb.graph add command -label "Zoom To Include All" \
    -accelerator "<Shift-F3>/<Control-F3>" \
    -command { zoomToIncludeAll }
.mb.graph add command -label "Manually Set Axis Ranges..." \
    -command { manuallySetAxisRanges }
.mb.graph add command -label "Manually Set Axis Labels..." \
    -command { manuallySetAxisLabels }
.mb.graph add command -label "Edit Current Trace Label..." \
    -command { getTraceLabelText }
.mb.graph add separator
.mb.graph add command -label "Add Text (Note) to Graph..." \
    -command { 
        set recordXYPointStatus waitingForTextLocation
        setMouseBindings pickPoint
        updateStatusMessage "Use mouse to select text location."
    }
.mb.graph add command -label "Add Line Segment to Graph" \
    -command { 
        set recordXYPointStatus waitingForFirstPoint
        setMouseBindings pickPoint
        updateStatusMessage "Use mouse to select first point."
    }
.mb.graph add command -label "Add Test Time Marker to Graph" \
    -command { 
        set recordXYPointStatus waitingForMarkerLevel
        setMouseBindings pickPoint
        updateStatusMessage "Use mouse to select level for test time marker."
    }
.mb.graph add separator
.mb.graph add command -label "Delete Range Markers" \
    -command { deleteMarkers range }
.mb.graph add command -label "Delete All Markers and Notes" \
    -command { deleteMarkers all }
.mb.graph add command -label "Delete Current Trace" \
    -accelerator "<F4>" \
    -command { deleteTrace $current }
.mb.graph add command -label "Clear Graph" \
    -accelerator "<Shift-F4>/<Control-F4>" \
    -command { clearGraph }
.mb add cascade -label Graph -menu .mb.graph
puts "End of Graph Menu"

menu .mb.text -tearoff 0
.mb.text add command -label "Display Channel Names" \
    -accelerator "<F5>" \
    -command { displayChannelNames }
.mb.text add command -label "Display Channel Header" \
    -accelerator "<F6>" \
    -command { displayChannelHeader }
.mb.text add command -label "Display Run Description" \
    -accelerator "<F7>" \
    -command { displayRunDescription }
.mb.text add separator
.mb.text add command -label "Display Shot List" \
    -accelerator "<F8>" \
    -command { displayShotList }
.mb add cascade -label Text-Window -menu .mb.text
puts "End of Text Menu"

menu .mb.help -tearoff 0
.mb.help add command -label "About..." -command { showAboutBox }
.mb.help add command -label "Display Help File" -command { displayHelpFile }
.mb add cascade -label Help -menu .mb.help
puts "End of Help Menu"

