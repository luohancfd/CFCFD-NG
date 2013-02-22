#!/bin/sh
# getCurrentData.tcl \
exec tclsh "$0" ${1+"$@"}

# Shock-Tunnel Data Management
# Department of Mechanical Engineering
# The University of Queensland
#
# This script will find the shots the are already in the database
# and writes a list to the file current_data.txt
#
# Peter Jacobs, Apr 2004

# -------------------------------------------------------------------------

set ::debugLevel 0
set scriptDir [file dirname [info script]]

# -------------------------------------------------------------------------

proc printUsage {} {
    puts "-------------------------------------------------------------------------------"
    puts "Usage: ./getCurrentData.tcl facilityName archiveDir"
    puts ""
    puts "  facilityName is the name of the shock tunnel / expansion tube."
    puts "               (e.g. T4, X3, X2)"
    puts "  archiveDir   is the directory where the shot information is stored"
    puts ""
    puts "Example: The following command "
    puts "./getCurrentData.tcl T4 /home/moncdata"
    puts "looks for T4 shots in the database and results in"
    puts "the file current_data.txt being written into /home/moncdata/T4/"
    puts "-------------------------------------------------------------------------------"
}

puts "Start processing..."

#
# Start by pulling apart the command line.
#
if {$argc < 2} {
    printUsage
    exit
}

set facilityName [lindex $argv 0]
set archiveDir [lindex $argv 1]
cd $archiveDir; set archiveDir [pwd];  # to ensure that we have the full path
set workDir [file join $archiveDir $facilityName]

set textFileName [file join $workDir current_data_$facilityName.txt]
if {$::debugLevel == 1} {
    puts "scriptDir     : $scriptDir"
    puts "facilityName  : $facilityName"
    puts "workDir       : $workDir"
    puts "textFileName  : $textFileName"
}

# ----------------------------------------------------------------------------

cd $workDir
set queryString "select basic_file_name from shot_descriptions"
append queryString " where facility_name = '$facilityName'"
if {$::debugLevel > 0} {
    puts "queryString: $queryString"
}
set result [exec psql -c $queryString -d moncdata]
if {$::debugLevel > 0} {
    puts "result: $result"
}
set fp [open "current_data.txt" "w"]
foreach line [split $result "\n"] {
    if {[string first "basic_file_name" $line] >= 0} {
	continue
    } elseif {[string first "----------" $line] >= 0} {
	continue
    } elseif {[string first "rows)" $line] >= 0} {
	continue
    } else {
	set basic_file_name [string trim $line]
	if {[string length $basic_file_name] > 0} {
	    puts $fp $basic_file_name
	}
    }
}
close $fp

# -------------------------------------------------------------------------------
puts "End job."














