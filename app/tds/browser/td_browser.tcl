#!/bin/sh
#\
exec wish "$0" ${1+"$@"}

# Use the following on a unix machine if you want to run
# under ActiveTcl which is not your default Tcl.
# exec /usr/local/ActiveTcl/bin/wish "$0" ${1+"$@"}

# File:
# td_browser.tcl
#
# Purpose:
# Tunnel Data Browser
# Ask for some data from the td_server.tcl CGI script and then
# display it using a BLT graph -- proof of principle.
#
# Copyright (C) 2002
# 
# Author: 
# Peter Jacobs
# Department of Mechanical Engineering
# The University of Queensland
# http://www.mech.uq.edu.au/staff/jacobs/cfcfd/
#
# Revisions:
# 31-Dec-01 original code
# 01-Jan-02 Uses Http/CGI server or local script to fetch the data.
#           Most of the user interface functions but it is ugly.
# 03-Jan-02 Clean up the GUI.
# 05-Jan-02 Read the ASCII files directly from the local disk.
# 08-Jan-02 Speed improvements obtained by using lists.
# 10-Jan-02 Initialization file is now sourced.
# 12-Feb-02 Added (hopefully) decent error dialogs.
# 20-Feb-02 Adjust font size for plot elements; 
#           add basic http authorization: username and password.
# For later revisions, see the log messages in the RCS file.

package require http
package require base64
package require rbc
package require BWidget

# --------------------------------------------------------------
puts "Begin preamble..."

if { [string length [info script]] } {
    set td(scriptHome) [file dirname [info script]]
} else {
    set td(scriptHome) [file dirname [info nameofexecutable]]
}; # end if

set suffix [info sharedlibextension]
set sharedLibName [file join $td(scriptHome) td_br_filter$suffix]
puts -nonewline "Load shared module $sharedLibName"
if { [catch {load $sharedLibName Td_br_filter} loadResult] } {
    set td(fastFilter) 0
    puts " Failed: $loadResult"
} else {
    set td(fastFilter) 1
    puts " OK."
}; # end if

# Locate the defrock directory assuming that we are starting
# in the tds directory.
set td(defrockDir) [file join [pwd] .. defrock]
# puts "td(defrockDir) is $td(defrockDir)"
if { [string equal $::tcl_platform(platform) windows] } {
    set td(defrockProgram) defrock
} else {
    set td(defrockProgram) ./defrock
}; # end if

# Remember where we started.
set td(workDir) [pwd]

proc timer_reset {} { set ::timer0 [clock clicks] }
proc timer_now {} { return [expr 1.0e-6 * ([clock clicks] - $::timer0)] }

puts "End preamble."

# ---------------------------------------------------------------
# Most helper procedures are in other files...

set moduleList [list td_br_read.tcl td_br_fetch.tcl td_br_gui_procs.tcl \
    td_br_data.tcl]
foreach f $moduleList {
    puts "Begin source $f..."
    source [file join $td(scriptHome) $f]
    puts "End source $f."
}; 

# ---------------------------------------------------------------
# A couple of miscellaneous procedures...

proc saveDisplayedText {} {
    global td
    set fileName [tk_getSaveFile -title "Save Displayed Text to..." \
        -initialdir $td(workDir) ]
    if {[string length $fileName] > 0} {
        set fp [open $fileName "w"]
        puts $fp [$::headerText get 1.0 end-1c]
        close $fp
        puts "Displayed text written to $fileName."
        set td(workDir) [file dirname $fileName]; # save for next time
    }; # end if    
}; # end proc


proc printPostscriptGraph {} {
    global td
    set fileName [tk_getSaveFile -title "Print postscript to file..." \
        -initialdir $td(workDir) ]
    if {[string length $fileName] > 0} {
        $::mygraph postscript output $fileName -decorations 0 -center 1
        set td(workDir) [file dirname $fileName]; # save for next time
    }; # end if    
}; # end proc


proc tidyUp {} {
    global td
    # do a bit of cleaning up as we exit
    catch { file delete $td(tempDataFile) }
}; # end proc

# ------------------ start of main script --------------------
# This is pretty much just setting up the GUI for later interaction.
# At the end of the script, the Tcl event loop is active and the user
# may then do things via the GUI widgets or via the command line.

if { [string equal $tcl_platform(platform) windows] } {
    console show
}; # end if

source [file join $td(scriptHome) td_br_config.tcl]

if { [file exists $td(iniFileName)] } {
    puts "Read customization file: $td(iniFileName)"
    source $td(iniFileName)
}; # end if

set moduleList [list td_br_menu.tcl td_br_gui_elements.tcl]
foreach f $moduleList {
    puts "Begin source $f..."
    source [file join $td(scriptHome) $f]
    puts "End source $f."
}; 

# Now, the Tcl event loop takes over.
updateStatusMessage "Ready."
puts "**** READY ****"
puts ""
