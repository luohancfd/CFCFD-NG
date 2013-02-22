#!/bin/sh
#\
exec tclsh "$0" ${1+"$@"}

# File:
# td_scan.tcl
#
# Purpose:
# Scan the shot directories and convert the binary data file
# that have been written by MONC into equivalent text files
# and write a shot.index file to aid the server script.
# We'll make use of the defrock program to do the actual conversion.
# 
# Note that the original MONC files should remain intact.
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
# 27-feb-02 original code
# 06-mar-02 added option to process a single shot,
# 13-mar-02 cope with X2 style MOD files
# 23-jan-05 be more careful with new-style shot data that have
#           been collected by son-of-monc
#
# ---------------- preamble ----------------------

if { [string length [info script]] } {
    set td(scriptHome) [file dirname [info script]]
} else {
    set td(scriptHome) [file dirname [info nameofexecutable]]
}; # end if
# puts "td(scriptHome) is $td(scriptHome)"
set td(logFileName) td_scan.log
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

# ---------------- procedures --------------------

proc echoUsage {} {
    puts "Usage: td_scan.tcl dataDir action ?shot=shotId?"
    puts "where:"
    puts "dataDir  is the name of the directory which contains"
    puts "         the individual shot directories."
    puts "action   is one of new|refresh|index"
    puts "         new option creates the ascii files if they are"
    puts "             not already present."
    puts "         refresh option forces the ASCII files to" 
    puts "             be regenerated even if they are already present."
    puts "         index option writes an index file in dataDir"
    puts "By default, all shot directories are scanned."
    puts "A single shot can be scanned as specified with the"
    puts "shot=shotId argument."
    puts "--------------------------------------------------------"
}; # end proc


proc looksLikeAShotDir { shotDir } {
    # From the main data directory, look for either MONC or ASCII files.
    return [expr [MONCFilesArePresent $shotDir] \
	    || [ASCIIFilesArePresent $shotDir]]
}; # end proc


proc MONCFilesArePresent { shotName } {
    # From the main data directory, look for the MONC header file
    # within the shot directory.
    # Try various upper- and lower-case combinations.
    set shotDir $shotName
    set shotName1 $shotName
    set pattern1 [join [list $shotDir / $shotName1 . {{hed,HED}}] ""]
    set shotName2 [string toupper $shotName]
    set pattern2 [join [list $shotDir / $shotName2 . {{hed,HED}}] ""]
    set shotName3 [string tolower $shotName]
    set pattern3 [join [list $shotDir / $shotName3 . {{hed,HED}}] ""]
    set hedFileList [glob -nocomplain $pattern1 $pattern2 $pattern3]
    return [expr [llength $hedFileList] >= 1]
}; # end proc


proc ASCIIFilesArePresent { shotName } {
    # From the main data directory, look into the shot directory
    # to see if the ASCII files are already present in compressed form.
    # Test 1: look for the LST file
    # Test 2: there should be at least one other compressed data file.
    # Both tests need to be true to return a true result.
    set asciiName [join [list $shotName A] ""]
    set i1 [file exists $shotName/$asciiName.LST.gz]
    set fileList [glob -nocomplain $shotName/$asciiName*gz]
    set i2 [expr [llength $fileList] >= 2]
    return [expr $i1 && $i2]
}; # end proc


proc fixFileNameCase { shotDir } {
    # From the main data directory, look to see if the MONC data
    # files have the same case as the directory name.
    # We assume that the shotName and the directory name are the same.
    set shotName $shotDir
    set saveDir [pwd]
    cd $shotDir
    set fileList [glob -nocomplain *]
    foreach f $fileList {
        set oldRootName [file rootname $f]
        set extn [file extension $f]
        if { ![string equal $oldRootName $shotName] } {
            # The names don't match exactly
            if { [string compare -nocase $oldRootName $shotName] == 0 } {
                # ...but they do match except for case of the letters
                catch { file rename $f [join [list $shotName $extn] ""] } junk
            } else {
                # This file name is not of the form shot.nnn; leave it.
            }; # end if
        }; # end if
    }; # end foreach
    cd $saveDir
}; # end proc


proc setDataFilePermissions { shotName } {
    # From the main data directroy, look into the MONC data
    # directory for the shot and make sure that the data files
    # have reasonable permissions.
    # For some reason, defrock seems to need the files to be
    # writeable as well as readable.
    set fileList [glob -nocomplain $shotName/*]
    foreach f $fileList {
        # r=1 w=1 x=0  --> 6
        file attributes $f -permissions 0664
    }; # end foreach
}; # end proc


proc checkForModFile { shotName } {
    # Assume that we are starting from the main data directory.
    # Look for one of the forms of the MOD file.
    # Returns a value of 1 if a standard MOD file already exists
    # or can be copied from a similarly-named MOD file.
    cd descript
    set foundModFile 0
    # Initially, look for a standard (T4) style MOD file name.
    set pattern1 [join [list $shotName . {{mod,MOD}}] ""]
    set fileList [glob -nocomplain $pattern1]
    if { [llength $fileList] == 0 } {
        # Didn't find a T4 standard MOD file.
        # Allow for upper- and lower-case variations.
        set pattern2 [join [list [string toupper $shotName] . {{mod,MOD}}] ""]
        set pattern3 [join [list [string tolower $shotName] . {{mod,MOD}}] ""]

        # Secondly, look for an X2 standard
        # MOD file that has the shot Id starting with "s" or "S"
        # while the MOD filename drops this first character.
        if { [string compare -nocase -length 1 $shotName s] == 0 } {
            # Try removing that first "s".
            set shotName2 [string range $shotName 1 end]
            set pattern1 [join [list $shotName2 . {{mod,MOD}}] ""]
        } else {
            # Try adding an "s".
            set pattern1 [join [list {[sS]} $shotName . {{mod,MOD}}] ""]
        }; # end if

        set fileList [glob -nocomplain $pattern1 $pattern2 $pattern3]
        if { [llength $fileList] >= 1 } {
            # Found a MOD file with a nonstandard name. Copy it.
            set sourceFile [lindex $fileList 0]
            set destinationFile [join [list $shotName .MOD] ""]
            file copy $sourceFile $destinationFile
            set foundModFile 1
        } else {
            # Didn't find an X2 style MOD file, cant' do much else.
            set foundModFile 0
        }; # end if
    } else {
        set foundModFile 1
    }; # end if
    cd ..
    return $foundModFile
}; # end proc


# ---------------- start main script -----------------

# Command-line arguments
if { $argc < 2 } {
    echoUsage
    exit
};

set td(dataDir) [lindex $argv 0]
if { ![file isdirectory $td(dataDir)] } {
    puts "$td(dataDir) is not a directory."
    exit
}; # end if

set td(commandOption) [lindex $argv 1]
if { [lsearch {new refresh index} $td(commandOption)] == -1 } {
    puts "Invalid action: $td(commandOption)"
    exit
}; # end if

if { $argc >= 3 } {
    # Assume that the third command line argument 
    # is there to specify one shot.
    set td(specificShotId) [lindex [split [lindex $argv 2] "="] 1]
} else {
    set td(specificShotId) ""
}; # end if

set logFile [open $td(logFileName) w]
puts $logFile "dataDir is $td(dataDir)"
if { [string equal $td(commandOption) refresh] } {
    puts $logFile "Will refresh ASCII files in shot directories."
} elseif { [string equal $td(commandOption) new] } {
    puts $logFile "Will generate ASCII files in shot directories "
    puts $logFile "only where they don't yet exist."
} elseif { [string equal $td(commandOption) index] } {
    puts $logFile "Will generate a shot.index file."
} else {
    puts "Should not have reached this point."
    puts "td(commandOption) is $td(commandOption)"
    close $logFile
    exit
}; # end if

# Now, sift through the directories within td(dataDir) and
# process those subdirectories that contain shot data.

cd $td(dataDir)
if { $td(specificShotId) == "" } {
    set directoryList [glob -nocomplain *]
} else {
    set directoryList $td(specificShotId)
}; # end if
# puts "directoryList: $directoryList"

# filter out directories that don't contain shot data
foreach dirName $directoryList {
    if { ![looksLikeAShotDir $dirName] } {
        set i [lsearch $directoryList $dirName]
        if { $i >= 0 } {
            set directoryList [lreplace $directoryList $i $i]
        }; # end if
    }; # end if
}; # end foreach

if { [string equal $td(commandOption) index] } {
    # Assming that we are already in $td(dataDir),
    # write the shot.index file with one shot id per line.
    set indexFile [open "shot.index" w]
    puts $indexFile [join [lsort $directoryList] "\n"]
    close $indexFile
    puts $logFile "Index file written."
    close $logFile
    exit
}; # end if


puts $logFile "------------- Start of Report ---------------"
foreach shotName $directoryList {
    puts -nonewline $logFile "shot $shotName "
    puts -nonewline "$shotName "
    flush stdout

    # Assume that the files do not already exist and
    # that we want to make them.
    set doEraseASCIIFiles 0
    set doMakeASCIIFiles  1

    if { [string equal $td(commandOption) new] && \
         [ASCIIFilesArePresent $shotName] } {
        puts -nonewline $logFile \
             "; ASCII files are already present; do nothing"
        set doEraseASCIIFiles 0
        set doMakeASCIIFiles  0
    }; # end if

    if { [string equal $td(commandOption) refresh] } {
        if { [MONCFilesArePresent $shotName] } {
            set doEraseASCIIFiles 1
            set doMakeASCIIFiles  1
	} else {
	    puts -nonewline $logFile \
		"; only ASCII files are present; do nothing"
            set doEraseASCIIFiles 0
            set doMakeASCIIFiles  0
	}
    }; # end if

    if { $doEraseASCIIFiles } {
        puts -nonewline $logFile "; erase ASCII files"
        cd $shotName 
        set targetPattern [join [list $shotName A * gz] ""]
        foreach targetFile [glob -nocomplain $targetPattern] {
            # puts "erase file: $targetFile"
            file attributes $targetFile -permissions 0666
            file delete -force $targetFile
        }; # end foreach
        cd ..
    }; # end if

    if { $doMakeASCIIFiles } {
        if { [checkForModFile $shotName] } {
            # Fix up the file names when shot directory name
            # doesn't have the same case as the files within it.
            fixFileNameCase $shotName

            # Make sure that we have appropriate permissions on the
            # MONC data files.
            # For some reason, defrock needs them to be writeable.
            if { [string equal $::tcl_platform(platform) windows] } {
                # do nothing on windows
            } else {
                setDataFilePermissions $shotName
            }; # end if

            # Build the ASCII files from the MONC files.
            # Be careful with where we start running defrock.
            # It needs the monc.val file to be in the same
            # directory as it is started.
            set saveDir [pwd]
            cd $td(defrockDir)
            catch { exec $td(defrockProgram) $td(dataDir) $shotName } \
                  defrockOutput
            # puts $defrockOutput
            cd $saveDir

            if { [regexp "End of defrock run." $defrockOutput] } {
                puts -nonewline $logFile "; defrock successful"
            } else {
                puts -nonewline $logFile "; defrock unsuccessful"
	    }; # end if

            puts -nonewline $logFile "; gzip ASCII files (if any)"
            cd $shotName
            set targetPattern [join [list $shotName A *] ""]
            foreach targetFile [glob -nocomplain $targetPattern] {
                catch { exec gzip -f $targetFile } gzipOutput
                # puts "; file $targetFile, gzipOutput is $gzipOutput"
            }; # end for
            # Tidy up, just in case we have dumped core here.
            file delete -force core
            cd $saveDir
	} else {
            puts -nonewline $logFile "; cannot find MOD file."
        }; # end if
    }; # end if

    puts $logFile "." ; # to terminate the line
}; # end foreach

puts "End of scan."
exit
# ---------------- end main script -----------------
