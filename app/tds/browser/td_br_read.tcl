# td_br_read.tcl
# Procedures for reading the ASCII MONC data files directly.


proc looksLikeAShotDir { shotName } {
    # From the main data directory, look for the gzipped channel list file
    # within the shot directory.
    set pattern [join [list $shotName / $shotName A . {{LST,lst}} .gz] ""]
    set listFileList [glob -nocomplain $pattern]
    if { [llength $listFileList] > 0 } {
        set result 1
    } else {
        set result 0
    }; # end if
    return $result 
}; # end proc


proc buildShotList {} {
    # Look up the facility directory and extract a list of shot directories.
    # Return the list as a string, with several shot names per line.
    global td

    set saveDir [pwd]
    set dataDir [file join $td(rootDir) $td(facility)]
    if { [file isdirectory $dataDir] } {
        cd $dataDir
        set shotList [glob -nocomplain *]
        # puts "In directory [pwd], shotList is $shotList"
        foreach dirName $shotList {
            if { ![looksLikeAShotDir $dirName] } {
                set i [lsearch $shotList $dirName]
                if { $i >= 0 } {
                    set shotList [lreplace $shotList $i $i]
                }; # end if
            }; # end if
        }; # end foreach
        # puts "After filter, shotList is $shotList"
        cd $saveDir
        # We actually want the end result with one shot per line.
        set resultString "For facility $td(facility), there are "
        append resultString "[llength $shotList] shots with ASCII data.\n"
        foreach {s0 s1 s2 s3 s4 s5 s6 s7 s8 s9} $shotList {
            append resultString "$s0 $s1 $s2 $s3 $s4 $s5 $s6 $s7 $s8 $s9 \n" 
        }; # end foreach
    } else {
        set resultString "td_browser_error : $dataDir is not a directory."
    }; # end if
    return $resultString
}; # end proc


proc readRunDescription {} {
    # Look for the Run Description text file in the 
    # RUNDESC (or rundesc) directory.
    # Take care to look for both upper- and lower-case names.

    global td

    # First, look for run description in the shot directory.
    # Failing that, try looking in the old-monc-style rundesc directory.
    set filePattern0 [join \
        [list $td(rootDir) / $td(facility) / $td(shot) / $td(shot) ".txt"] ""]
    set fileList [glob -nocomplain $filePattern0]
    if { [llength $fileList] == 0 } {
        # Try the old style...
	set dirPattern [join \
	    [list $td(rootDir) / $td(facility) / {{rundesc,RUNDESC}}] ""]
	set dirList [glob -nocomplain $dirPattern]
	if { [llength $dirList] == 0 } {
	    puts "Could not find Run Descriptions directory."
	}; # end if
	set dirName [lindex $dirList 0]

	# T4-style name
	set filePattern1 [join \
	    [list $dirName / $td(shot) . {{txt,TXT}}] ""]
	# X2-style name
	if { [string compare -nocase -length 1 $td(shot) s] == 0 } {
	    # take off the first "s" or "S"
	    set shortName [string range $td(shot) 1 end]
	} else {
	    set shortName $td(shot)
	}; # end if
	set filePattern2 [join \
            [list $dirName / {{run,RUN}} $shortName . {{txt,TXT}}] ""]
	set fileList [glob -nocomplain $filePattern1 $filePattern2]
    }; # end if

    if { [llength $fileList] == 0 } {
        return "td_browser_error : No Run Description file for $td(shot)."
    }; # end if
    set fileName [lindex $fileList 0]

    if { [file exists $fileName] } {
        set fp [open $fileName "r"]
        set content [read $fp]
        # puts $content
        close $fp
        return $content
    } else {
        return "td_browser_error : Run description file $fileName doesn't exist."
    }; # end if
}; # end proc


proc findShotDataDir {} {
    global td
    set dirName [file join $td(rootDir) $td(facility) $td(shot)]
    if { [file exists $dirName] } {
        set td(dataDir) $dirName
        return 1
    } else {
        puts "td_browser_error : Could not find shot directory $dirName"
        set td(dataDir) ""
        return 0
    }; # end if
}; # end proc


proc readChannelNames {} {
    global td
    if { [findShotDataDir] == 0 } { return 0 }
    set fname $td(shot)
    append fname A.LST.gz
    set fileName [file join $td(dataDir) $fname]
    if { [file exists $fileName] } {
        set td(dataFile) $fileName
    } else {
        puts "td_browser_error : Could not find file $fileName"
        set td(dataFile) ""
        return 0
    }; # end if
    if { [catch {exec gzip -d -c $td(dataFile) > $td(tempDataFile)} result] } {
        puts "td_browser_error : gunzip command failed"
        puts $result
        return 0
    }; # end if

    # we have successfully uncompressed the data file
    set fp [open $td(tempDataFile) "r"]
    set listOfNames [read $fp]
    close $fp

    return $listOfNames
}; # end proc


proc readChannelData { channel_id newIndex } {
    # Open the file, read the metadata and 
    # reconstruct the time values if they are not already there.
    # If everything is successful, a value of 1 is returned and
    # two new global vectors ::tVec$newIndex and ::vVec$newIndex 
    # will contain the channel's data.
    # Input: channel_id : channel identity, typically the file extension
    #        newIndex   : index value that will make the vector names unique
    #                     Typically this will be one more than the current
    #                     number of traces.

    global td

    if { [findShotDataDir] == 0 } { return 0 }

    set fname $td(shot)
    append fname A . $channel_id .gz
    set fileName [file join $td(dataDir) $fname]
    if { [file exists $fileName] } {
        set td(dataFile) $fileName
    } else {
        puts "td_browser_error : Could not find file $fileName"
        set td(dataFile) ""
        return 0
    }; # end if

    # puts "Before gzip : [timer_now]"
    if { [catch {exec gzip -d -c $td(dataFile) > $td(tempDataFile)} result] } {
        puts "td_browser_error : gunzip command failed"
        puts $result
        return 0
    }; # end if
    # puts "After gzip : [timer_now]"

    # we have successfully uncompressed the data file
    # puts "Read data from file $td(dataFile)"
    set fp [open $td(tempDataFile) "r"]

    set commentLine [gets $fp]

    # Extract the header information and save it.
    set td(headerText) ""
    while { [eof $fp] == 0 } {
        set headerLine [gets $fp]
        # discard comment character
        regsub {#} $headerLine {} headerLine
        # remove surplus whitespace
        set headerLine [string trim $headerLine]
        regsub -all {\s\s+} $headerLine { } headerLine
        append td(headerText) "$headerLine\n"
        if { [string length $headerLine] == 0 } {
            # we have hit the blank line
            break; 
        } else {
            foreach {name value} [split $headerLine] {
                set td(header.$name) $value
            }; # end for
        }; # end if
    }; # end while
    # puts "After reading header: [timer_now]"

    # Read the data itself.
    set N $td(header.dataPoints)
    set dt $td(header.timeInterval)
    set t0 $td(header.timeStart)
    if {[string equal $td(header.withTimeColumn) "no"]} {
        set generateTime 1
    } else {
        set generateTime 0
    }; # end if

    # Put the new data into fresh vectors
    if { [info exists ::tVec$newIndex] } {
        blt::vector destroy ::tVec$newIndex
    }; # end if
    if { [info exists ::vVec$newIndex] } {
        blt::vector destroy ::vVec$newIndex
    }; # end if
    blt::vector create ::tVec$newIndex
    blt::vector create ::vVec$newIndex
    # puts "After creating vectors: [timer_now]"

    if {[string length $td(currentTimeUnits)] == 0} {
	set td(currentTimeUnits) $td(header.timeUnits)
    }
    set multiplyFactor [getTimeFactor $td(header.timeUnits) $td(currentTimeUnits)]
    set tlist {}
    set vlist {}
    if {$td(add_t0_to_times) == 1} {
	set tstart $t0
    } else {
	set tstart 0.0
    }
    for {set j 0} {$j < $N} {incr j} {
        set dataLine [gets $fp]
        if { [eof $fp] } break
        if { $generateTime } {
            set t [expr $tstart + $j * $dt]
            set v [string trim $dataLine]
        } else {
            # remove surplus whitespace
            set dataLine [string trim $dataLine]
            regsub -all {\s\s+} $dataLine { } dataLine
            foreach {t v} [split $dataLine " \t"] {}
        }; # end if
        lappend tlist [expr $t * $multiplyFactor]
        lappend vlist $v
        # puts "$t $v"
    }; # end for
    ::tVec$newIndex set $tlist
    ::vVec$newIndex set $vlist
    # puts "After saving data in vectors: [timer_now]"
    close $fp

    return 1
}; # end proc readChannelData


proc getTimeFactor {from to} {
    # Multiplying factors to convert between different time units.
    array set factor {}
    set factor(microseconds,microseconds) 1.0
    set factor(microseconds,milliseconds) 0.001
    set factor(microseconds,seconds) 1.0e-6
    set factor(milliseconds,microseconds) 1000.0
    set factor(milliseconds,milliseconds) 1.0
    set factor(milliseconds,seconds) 0.001
    set factor(seconds,microseconds) 1.0e6
    set factor(seconds,milliseconds) 1.0e3
    set factor(seconds,seconds) 1.0
    if { [catch { set multiply_factor $factor($from,$to) } result] } {
        set multiply_factor 1.0
    } 
    return $multiply_factor
}; # end proc getTimeFactor


proc readPlainOldData { newIndex fileName } {
    # Open the file, read the channel data, time and value from each line. 
    # If everything is successful, a value of 1 is returned and
    # two new global vectors ::tVec$newIndex and ::vVec$newIndex 
    # will contain the channel's data.
    # Input: newIndex   : index value that will make the vector names unique
    #                     Typically this will be one more than the current
    #                     number of traces.
    #        filename   : the data file name

    global td metadataList

    if { [file exists $fileName] == 0 } {
        puts "td_browser_error : Could not find file $fileName"
        set td(dataFile) ""
        return 0
    }; # end if

    # We presume that there is no metadata in the file.
    foreach metadataName $metadataList {
        set td(header.$metadataName) ""
    }; # end foreach

    set fp [open $fileName "r"]

    set N 0
    set tlist {}
    set vlist {}
    set td(headerText) ""
    while { [eof $fp] == 0 } {
        set lineText [gets $fp]

        # remove surplus whitespace
        set lineText [string trim $lineText]
        regsub -all {\s\s+} $lineText { } lineText
        if { [string length $lineText] == 0 } {
            # we have hit a blank line
            continue; 
        }; # end if

        if { [string first {#} $lineText] >= 0 } {
            # A sharp character indicates a header/comment line.
            # append td(headerText) "$lineText\n"

            # discard comment character
            regsub {#} $lineText {} lineText
            # remove surplus whitespace
            set lineText [string trim $lineText]
            regsub -all {\s\s+} $lineText { } lineText
            append td(headerText) "$lineText\n"
            if { [string length $lineText] == 0 } {
                # we have hit the blank line, forget it
            } else {
                foreach {name value} [split $lineText] {
                    set indx [lsearch $metadataList $name]
                    if { $indx >= 0 } {
                         # force the name to be precise
                        set name [lindex $metadataList $indx]
                        set td(header.$name) $value
                    }; # end if
                }; # end for
            }; # end if
        } else {
            # Assume two values per line
            foreach {t v} [split $lineText] {}
            lappend tlist $t
            lappend vlist $v
            # puts "$t $v"
            incr N
        }; # end if
    }; # end while

    close $fp

    puts "Read $N data points from $fileName"
    if { $N == 0 } {return 0}; # We managed to load no data.

    # Put the new data into fresh vectors
    if { [info exists ::tVec$newIndex] } {
        blt::vector destroy ::tVec$newIndex
    }; # end if
    if { [info exists ::vVec$newIndex] } {
        blt::vector destroy ::vVec$newIndex
    }; # end if
    blt::vector create ::tVec$newIndex
    blt::vector create ::vVec$newIndex
    ::tVec$newIndex set $tlist
    ::vVec$newIndex set $vlist
    puts "New vectors set up."

    # Fill in some of the metadata that would normally be
    # the head of the tds style data files.
    set td(facility) ""
    set td(shot) $td(header.shotName)
    set td(channel) $td(header.channelId)

    return 1
}; # end proc readPlainOldData


proc generateASCIIFiles {} {
    global td
    set startDir [pwd]
    cd [file join $td(rootDir) $td(facility)]
    set dataDir [pwd]
    set shotName $td(shot)
    puts -nonewline "Generate ASCII files, shot $shotName :"
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
        catch { exec $td(defrockProgram) $dataDir $shotName } \
              defrockOutput
        puts $defrockOutput
        cd $saveDir

        if { [regexp "End of defrock run." $defrockOutput] } {
             puts -nonewline "; defrock successful"
        } else {
            puts -nonewline "; defrock unsuccessful"
        }; # end if

        puts -nonewline "; gzip ASCII files (if any)"
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
        puts -nonewline "; cannot find MOD file."
    }; # end if
    cd $startDir
    puts "...done."
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
                file rename $f [join [list $shotName $extn] ""]
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

