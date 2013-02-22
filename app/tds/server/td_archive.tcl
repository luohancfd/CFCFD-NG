#!/bin/sh
#\
exec wish "$0" ${1+"$@"}

# File:
# td_archive.tcl
#
# Purpose:
# Scan the shot directories and archive the shot data files
# into one TAR file for each shot.
#
# Usage:
# See below for procedure echoUsage.
# 
# Author: 
# Peter Jacobs
# Department of Mechanical Engineering
# The University of Queensland
#
# Revisions:
# 06-mar-02 original code
#
# ---------------- preamble ----------------------

if { [string length [info script]] } {
    set td(scriptHome) [file dirname [info script]]
} else {
    set td(scriptHome) [file dirname [info nameofexecutable]]
}; # end if
# puts "td(scriptHome) is $td(scriptHome)"
set td(logFileName) td_archive.log
# Remember where we started.
set td(workDir) [pwd]

# ---------------- procedures --------------------

proc echoUsage {} {
    puts "Usage: td_archive.tcl topLevelDir facilityName tarDir ?shot=shotId?"
    puts "where:"
    puts "topLevelDir is the name of the directory which contains"
    puts "            the facilityName directory."
    puts "facilityName is simply the name of the facility."
    puts "            (i.e. topLevelDir/facilityName/ contains"
    puts "            the individual shot directories.)"
    puts "tarDir      is the name of the directory in which we"
    puts "            will create the TAR files."
    puts "By default, all shot directories are archived."
    puts "A single shot can be archived as specified with the"
    puts "shot=shotId argument."
    puts "--------------------------------------------------------"
}; # end proc


proc looksLikeAShotDir { facilityName shotName } {
    # From the top level directory, look for the MONC header file
    # within the shot directory.
    # Try various upper- and lower-case combinations.
    set shotDir [join [list $facilityName / $shotName] ""]
    set shotName1 $shotName
    set pattern1 [join [list $shotDir / $shotName1 . {{hed,HED}}] ""]
    set shotName2 [string toupper $shotDir]
    set pattern2 [join [list $shotDir / $shotName2 . {{hed,HED}}] ""]
    set shotName3 [string tolower $shotDir]
    set pattern3 [join [list $shotDir / $shotName3 . {{hed,HED}}] ""]
    set hedFileList [glob -nocomplain $pattern1 $pattern2 $pattern3]
    return [expr [llength $hedFileList] >= 1]
}; # end proc


proc uniqueElements { {inputList {}} } {
    # Accepts a list which may have multiple entries
    # with the same value.
    # Returns a list containing only unique entries.
    foreach e $inputList {
        set myarray($e) 1
    }; # end foreach
    set uniqueList [array names myarray]
    if { [array exists myarray] } { unset myarray }
    return $uniqueList
}; # end if


# ---------------- start main script -----------------

# Command-line arguments
if { $argc < 3 } {
    echoUsage
    exit
};

set td(topLevelDir) [lindex $argv 0]
if { ![file isdirectory $td(topLevelDir)] } {
    puts "$td(topLevelDir) is not a directory."
    exit
}; # end if

set td(facilityName) [lindex $argv 1]
if { ![file isdirectory $td(topLevelDir)/$td(facilityName)] } {
    puts "$td(topLevelDir)/$td(facilityName) is not a directory."
    exit
}; # end if

set td(tarDir) [lindex $argv 2]
if { ![file isdirectory $td(tarDir)] } {
    puts "$td(tarDir) is not a directory."
    exit
}; # end if

if { $argc >= 4 } {
    # Assume that the fourth command line argument 
    # is there to specify one shot.
    set td(specificShotId) [lindex [split [lindex $argv 3] "="] 1]
} else {
    set td(specificShotId) ""
}; # end if

set logFile [open $td(logFileName) w]
puts $logFile "topLevelDir is $td(topLevelDir)"
puts $logFile "tarDir is $td(tarDir)"
puts $logFile "facilityName is $td(facilityName)"

# Now, sift through the directories within td(topLevelDir)/td(facilityName)
# and process those subdirectories that contain shot data.

cd $td(topLevelDir)/$td(facilityName)
if { $td(specificShotId) == "" } {
    set directoryList [glob -nocomplain *]
} else {
    set directoryList $td(specificShotId)
}; # end if
# puts "directoryList: $directoryList"
cd .. ; # back to topLevelDir

# filter out directories that don't contain shot data
foreach dirName $directoryList {
    if { ![looksLikeAShotDir $td(facilityName) $dirName] } {
        set i [lsearch $directoryList $dirName]
        if { $i >= 0 } {
            set directoryList [lreplace $directoryList $i $i]
        }; # end if
    }; # end if
}; # end foreach

puts $logFile "------------- Start of Report ---------------"
foreach shotName $directoryList {
    puts -nonewline $logFile "shot $shotName "
    puts -nonewline "$shotName "
    flush stdout
    set totalFiles 0

    set shotDir [join [list $td(facilityName) / $shotName] ""]
    incr totalFiles [llength [glob -nocomplain $shotDir/*]]

    # The following patterns should pick up both T4 and X2 style 
    # MOD file names. 
    if { [string compare -nocase -length 1 $shotName s] == 0 } {
        # take off the first "s" or "S"
        set shortName [string range $shotName 1 end]
    } else {
        set shortName $shotName
    }; # end if
    set p1a [join \
        [list $td(facilityName) /descript/ \
            $shotName . {{mod,MOD}}] ""]
    set p1b [join \
        [list $td(facilityName) /descript/ \
            $shortName . {{mod,MOD}}] ""]
    set p1c [join \
        [list $td(facilityName) /descript/ \
            [string toupper $shotName] . {{mod,MOD}}] ""]
    set p1d [join \
        [list $td(facilityName) /descript/ \
            [string toupper $shortName] . {{mod,MOD}}] ""]
    set p1e [join \
        [list $td(facilityName) /descript/ \
            [string tolower $shotName] . {{mod,MOD}}] ""]
    set p1f [join \
        [list $td(facilityName) /descript/ \
            [string tolower $shortName] . {{mod,MOD}}] ""]
    set files1 [uniqueElements \
        [glob -nocomplain $p1a $p1b $p1c $p1d $p1e $p1f]]
    incr totalFiles [llength $files1]

    # The following patterns should get upper- or lower-case txt extensions
    # for both T4 and X2-style names.
    set p2a [join \
        [list $td(facilityName) /rundesc/ \
            $shotName . {{txt,TXT}} ] ""]
    set p2b [join \
        [list $td(facilityName) /rundesc/ \
            {{run,RUN}} $shortName . {{txt,TXT}}] ""]
    set p2c [join \
        [list $td(facilityName) /rundesc/ \
            [string toupper $shotName] . {{txt,TXT}} ] ""]
    set p2d [join \
        [list $td(facilityName) /rundesc/ \
            {{run,RUN}} [string toupper $shortName] . {{txt,TXT}}] ""]
    set p2e [join \
        [list $td(facilityName) /rundesc/ \
            [string tolower $shotName] . {{txt,TXT}} ] ""]
    set p2f [join \
        [list $td(facilityName) /rundesc/ \
            {{run,RUN}} [string tolower $shortName] . {{txt,TXT}}] ""]
    set files2 [uniqueElements \
        [glob -nocomplain $p2a $p2b $p2c $p2d $p2e $p2f]]
    incr totalFiles [llength $files2]

    if { $totalFiles > 0 } {
        # Build a unique tar file name using both the facility and shot names.
        set tarFileName [join \
            [list $td(tarDir) / $td(facilityName) _ $shotName .tar] ""]
        # Need the eval in case the files1 and files2 lists contain 
        # more than one element each.
        if { [catch { eval exec tar cf $tarFileName $shotDir $files1 $files2 } \
             cmdOutput] } {
            puts -nonewline $logFile "; error; command output is $cmdOutput"
        } else {
            puts -nonewline $logFile "; $totalFiles files --> $tarFileName"
        }; # end if
    } else {
        puts -nonewline $logFile "; no files to archive"
    }; # end if

    puts $logFile "" ; # to terminate the line
}; # end foreach

puts "End of archive run."
exit
# ---------------- end main script -----------------
