# td_br_config.tcl
# Standard onfiguration settings and initialization of some data.
#
# This file should not normally be edited by the user.
# Put customizations into the file td_browser.ini instead.

set td(versionAndDate) "0.33 05-Apr-09"

# options for the dataSource: httpServer localServer directRead
set td(dataSource) httpServer
set td(serverIP) "www.mech.uq.edu.au"
set td(proxyhost) "http-proxy.uq.edu.au"
set td(proxyport) "80"
set td(useProxy) 0
set td(username) "tdguest"
set td(password) "tdpasswd"
set td(dataDir) ""
set td(dataFile) ""

set td(facility) T4
set td(facilityList) [list T4 X2 X3 drummond]
set td(shot) 7319
set td(shotList) [list "list"]
set td(channel) "spa"
set td(channelList) [list "list"]
set td(normalizeChannel) ""
set td(timeShift) ""
# People using a mix of the newer LabView files with the older
# monc data files might like to set td(currentTimeUnits) "seconds"
# or "milliseconds" or "microseconds".
set td(currentTimeUnits) ""
array set channelNameToNumber {}

set td(headerText) "Dummy header text"
set td(runDescriptionText) "Dummy description text."
set td(textFont) {Courier 10}
set td(textHeight) 12

set td(default.nhalf) 20
set td(nhalf) 0
set td(referenceTime) 0.0
set td(referenceValue) 0.0
# People using a mix of newer LabView files (having t0 corresponding
# to primary diaphragm rupture) and older files with nominally zero t0
# may want to set this following vzlue to 0.
set td(add_t0_to_times) 1

set td(t1) ""
set td(t2) ""
set td(yTop) ""
set td(y2Top) ""
set td(yBot) ""
set td(y2Bot) ""
set td(t1Index) ""
set td(t2Index) ""
set td(averageValue) ""
set td(stddevValue) ""
set td(width) 8i
set td(height) 4i
set td(xAxisLabel) ""
set td(yAxisLabel) ""
set td(y2AxisLabel) ""
set td(yAxis) y
set td(rescaleAxesOption) onlyY; # options: onlyY both none
set td(plotFont) {Helvetica 12}
set td(useColours) 1
set td(useDashedLines) 0
set td(addTraceLabels) 1
set td(tickSizeInPixels) 10

set td(includeMetaDataWhenSaving) 1

# The global variable counts the number of traces that we have loaded.
# It needs to be incremented just as each new data set
# is read or fetched.
set numberOfTraces 0

# The global variable current indicates the "current" trace
# on which we will be working.
# Thus data sets should be indexed 1..$numberOfTraces.
set current $numberOfTraces

if { [string equal $::tcl_platform(platform) windows] } {
    # Home directory for the MONC data on MS-windows boxes
    set td(rootDir) [file join c: /moncdata]
    # A temporary file for decompressing the MONC data files, etc.
    set td(tempDataFile) [file join c: /tdsTempData]
} else {
    # Home directory for the MONC data on Unix boxes
    set td(rootDir) [file join / home moncdata]
    # A temporary file for decompressing the MONC data files, etc.
    set td(tempDataFile) [file join / tmp tdsTempData]
};

set td(workDir) $td(scriptHome)
set td(helpFileName) [file join $td(scriptHome) td_browser.help]

set td(iniFileName) [file join $td(scriptHome) td_browser.ini]
# Locate the server script if possible
set td(serverScriptHome) [file join $td(scriptHome) .. server]
if { ![file exists $td(serverScriptHome)] } {
    set td(serverScriptHome) ""
}; # end if

# Try to make the temporary data file name unique.
append td(tempDataFile) [expr int(rand() * 100000000)]

# We used to use the blt::busy procedure,
# now it sometimes causes trouble.
set td(useBltBusy) 0

# some state data for the procedure that adds annotations
# to the graph
set recordXYPointStatus none
set annotationText ""
set annotationAnchor w

# arrays to hold information on each plotted data set
array set incOS {} ;          # incremental offsets and scales
array set cumOS {} ;          # cumulative offsets and scales
array set whichYAxis {} ;     # holds the value y or y2
array set traceName {} ;      # look up trace number to get graph element name
array set traceNumber {} ;    # look up trace name to get number
array set traceIsPlotted {} ; # boolean to indicate trace status

set traceNumberList {};       # list of plotted traces for the ComboBox

# some other global data
set statusText ""

# Here is a list of the metadata expected at the top of the
# tds-style data files.
set metadataList [list \
    shotName \
    channelId \
    withTimeColumn \
    dataPoints \
    dataType \
    dataUnits \
    timeStart \
    timeAverageWindow \
    timeInterval \
    timeUnits \
    transducerSensitivity \
    transducerSensitivityUnits \
    transducerName \
    transducerSerialNumber \
    gain \
    qfluxgain \
]; # end list
