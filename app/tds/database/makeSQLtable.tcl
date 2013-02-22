#!/bin/sh
# makeSQLtable.tcl \
exec tclsh "$0" ${1+"$@"}

# Shock-Tunnel Data Management
# Department of Mechanical Engineering
# The University of Queensland
#
# This script will find the shots the aren't already in the database
# (according to the file current_data.txt) and then writes them
# into a text file (shot_descriptions.sql) ready for insertion into the database.
#
# Original Effort: Aaron Brandis, Dec 2001 -- Feb 2002
# Rewritten      : Peter Jacobs, Apr 2004

# -------------------------------------------------------------------------

package require textutil

set ::debugLevel 0
set scriptDir [file dirname [info script]]

# -------------------------------------------------------------------------

source [file join $scriptDir fieldsForRunDescription.tcl]
source [file join $scriptDir processRunDescriptionFile.tcl]

proc printUsage {} {
    puts "-------------------------------------------------------------------------------"
    puts "Usage: ./makeSQLtable.tcl facilityName archiveDir shotName ?-createcmd?"
    puts ""
    puts "  facilityName is the name of the shock tunnel / expansion tube."
    puts "               (e.g. T4, X3, X2)"
    puts "  archiveDir   is the directory where the shot information is stored"
    puts "  shotName     can be given the special value all (add in all new shots)"
    puts "               or set to a specific shot name."
    puts "               This is just the name of the shot as it appears in its "
    puts "               specific directory.  It may include letters and digits."
    puts "  -createcmd   indicates that we want to start the SQL file with "
    puts "               a CREATE command."
    puts ""
    puts "Example 1: Do all new shots for X3 and include the CREATE table command."
    puts "./makeSQLtable.tcl X3 /home/moncdata all -createcmd"
    puts "Example 2: Do one specific shot (in this case 7319) for T4."
    puts "./makeSQLtable.tcl T4 /home/moncdata 7319"
    puts "-------------------------------------------------------------------------------"
}

puts "Start processing..."

#
# Start by pulling apart the command line.
#
if {$argc < 3} {
    printUsage
    exit
}

set facilityName [lindex $argv 0]
set archiveDir [lindex $argv 1]
cd $archiveDir; set archiveDir [pwd];  # to ensure that we have the full path
set workDir [file join $archiveDir $facilityName]
set processOption [lindex $argv 2]
if {$processOption == "all"} {
    set shotName ""
} else {
    set shotName $processOption
}
if {$argc > 3} {
    set addCreateCmdFlag [string equal [lindex $argv 3] "-createcmd"]
} else {
    set addCreateCmdFlag 0
}

set mysqlFileName [file join $workDir shot_descriptions.sql]
set problemsFileName [file join $workDir maintenance_problems_$facilityName.txt]
if {$::debugLevel == 1} {
    puts "scriptDir     : $scriptDir"
    puts "facilityName  : $facilityName"
    puts "workDir       : $workDir"
    puts "processOption : $processOption"
    puts "shotName      : $shotName"
    puts "mysqlFile     : $mysqlFileName"
    puts "problemsFile  : $problemsFileName"
}

# ----------------------------------------------------------------------------
# Before getting started, make a list of already-existing entries.
# We assume that the file current_data.txt has been obtained from the
# SQL database already and it has basicFileName as the first column. 
# If this is not so, the list will remain empty.
#
set fileName [file join $workDir current_data.txt]
if { [catch {open $fileName "r"} result] } {
    puts "Warning: $result"
    set existingShotList {}
} else {
    set existingShotList {}
    set existingDataFile $result
    foreach line [split [read $existingDataFile] \n] {
	# Pick up the basicFileName as the shot identifier
	lappend existingShotList [lindex [split $line " \t"] 0]
    }
}    
if { $::debugLevel >= 1 } {
   puts "existingShotList: $existingShotList"
}

# ---------------------------------------------------------------------------
# Now, we are ready to scan the run-description files and collect
# the data into text-file describing the table.
#
set mysqlFile [open $mysqlFileName w]
if { $::debugLevel >= 1 } {
    foreach name $fieldNames {
	puts -nonewline $mysqlFile "$name\t"
    }
    puts $mysqlFile ""
}
if { $addCreateCmdFlag } {
    puts $mysqlFile [makeSQLCreateCmd]
}

cd [file join $workDir rundesc]
if {$processOption == "all"} {
    set files_to_be_searched [glob -nocomplain *.{txt,TXT}]
} else {
    set files_to_be_searched [glob -nocomplain $shotName*.{txt,TXT}]
}
if {$::debugLevel >= 0} {
    puts "files_to_be_searched: $files_to_be_searched"
}

foreach thisFile $files_to_be_searched {
    set basicFileName [file rootname [file tail $thisFile]]
    regexp {[[:digit:]]+} $thisFile shotNumber
    if { [lsearch $existingShotList $basicFileName] == -1 } {
	puts -nonewline "basicFileName: $basicFileName "
	if [processRunDescriptionFile [file join $workDir rundesc $thisFile]] {
	    set ::dataValue(facility_name) $facilityName
	    set ::dataValue(basic_file_name) $basicFileName
	    puts "Write row for SQL table."
	    puts $mysqlFile [makeSQLInsertCmd]
	    writeXMLFileForShotDescription \
		[file join $workDir rundesc $basicFileName.xml]
	} else {
	    puts "Skipped."
	}
    }
}

close $mysqlFile

# -------------------------------------------------------------------------------
puts "End job."














