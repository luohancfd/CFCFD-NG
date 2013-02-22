#!/bin/sh
#\
exec tclsh "$0" ${1+"$@"}

# File:
# td_server.tcl
#
# Purpose:
# Send some shock tunnel data via the CGI interface.
#
# Command-line Usage:
# td_server.tcl facility=<facility>+shot=<shot>+channel=<channel>+part=<part>
# where:
# <facility> is one of (T4, X1, X2, X3)
# <shot> is the part of the file name that identifies the shot.
#        MONC data (and the equivalent compressed ASCII data) is stored in 
#        a directory of this name.
#        If a special value of "list" is supplied, a list of shots for
#        the specified facility is returned.
# <channel> is the channel identifier, also used by MONC as the extension 
#        of the channel's data file.
#        A value of "list" will cause the server to return a list of channel
#        numbers and transducer names, one channel per line.
#        A value of "rundescription" will cause the server to return the
#        content of the run description file.
# <part> is one of (info, data, raw)
#        A value of "info" causes the server to return the header information
#        for the channel.
#        A value of "data" will cause the server to return just the history
#        data as a number of "time value" pairs.
#        A value of "raw" will cause the server to return the uncompressed
#        data file unprocessed.
#        This parameter is optional as it has no meaning if a particular
#        channel is not specified and it will default to a value of "data"
#        if a valid channel is specified.
#
# CGI Usage:
# Instead of getting the query string from the command-line,
# the script looks in the environment variable QUERY_STRING.
#
# Author:
# Peter Jacobs
# Department of Mechanical Engineering
# The University of Queensland
#
# Revisions:
# 26-Dec-01 original code
# 29-Dec-01 CGI working, times generation
# 30-Dec-01 documentation, part selector, send list of channels
# 31-Dec-01 options to return shot list and run descriptions
# 01-Jan-02 generalize for Win32 use
# 05-Jan-02 include keyword "td_server_error" in returned error messages
#           return the raw/original data file if requested
# 29-Mar-02 put all results into a content string so that the size
#           may be sent as part of the HTTP header.
# 01-Jun-03 Added allowAccess procedure to filter CGI requests.
# 23-Jan-05 Look for run-description file in shot directory first. 

# ---------------------- procedures ----------------------------

proc initializeData {} {
    global td contentString allowedShotList
    set contentString ""

    if { [string equal $::tcl_platform(platform) windows] } {
        # Home directory for the MONC data
        set td(rootDir) [file join d: /home moncdata]
        # A temporary file for decompressing the MONC data files, etc.
        set td(tempDataFile) [file join d: /tdsTempData]
    } else {
        # Home directory for the MONC data; there is more than one possibility.
        set td(rootDir) [file join / home moncdata]
	if { ![file exists $td(rootDir)] } {
	    set td(rootDir) [file join / home2 moncdata]
	}
        # A temporary file for decompressing the MONC data files, etc.
        set td(tempDataFile) [file join / tmp tdsTempData]
    };
    # Try to make the temporary data file name unique.
    append td(tempDataFile) [expr int(rand() * 100000000)]

    set td(dataDir) ""
    set td(dataFile) ""

    # check if we are in a CGI script environment
    set td(cgiScript) [info exists ::env(QUERY_STRING)]

    # Data for the allowAccess procedure: for each restricted user,
    # there is a list of facility+shotId pairs.
    # We may put in a few sample sets but we'll need to think
    # of something "scalable" it it gets out of hand.
    set allowedShotList(tdguest) [list T4+7319 T4+7320 T4+list X3+list X2+list]
}; # end proc


proc sendContentString {} {
    global td
    global contentString

    if { $td(cgiScript) } {
        puts "Content-type: text/plain"
        puts "Content-length: [string length $contentString]"
        puts ""
    }; # end if
    puts -nonewline $contentString
}; # end proc


proc parseQueryString {} {
    # We may have started this script from the normal shell and
    # have passed the query string via the command line.
    # First, check the QUERY_STRING environment variable and,
    # if there is nothing there, try the command-line.

    global td

    # default values for the expected command-line parameters
    set td(facility) ""
    set td(shot) ""
    set td(channel) ""
    set td(part) "data"

    if { [info exists ::env(QUERY_STRING)] } {
        set queryString $::env(QUERY_STRING)
    } else {
        set queryString [join $::argv "+"]
    }; # end if

    foreach {a} [split $queryString "+"] {
        foreach {name value} [split $a "="] {
            set td($name) $value
            # puts "name: $name, value: [set td($name)]"
        }; # end foreach
    }; # end foreach
}; # end proc


proc sendTestData {} {
    # Generate some dummy data
    global contentString
    for {set t 0.0} {$t < 10.0} {set t [expr $t + 0.1]} {
        set v [expr sin($t)]
        append contentString "$t $v\n"
    }; # end for
}; # end proc


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

        if { [file exists "shot.index"] } {
            # Read the shot list from the file.
            set indexFile [open "shot.index" r]
            set shotList [read -nonewline $indexFile]
            close $indexFile
            # puts "shot.index content: $shotList ."
            set infoSource "(used shot.index file)"
        } else {
            # Build the shot list by searching individual directories.
            set shotList [glob -nocomplain *]
            # puts "In directory [pwd], before filter, shotList is $shotList"
            # Filter out those directories that appear
            # to NOT contain shot data.
            foreach dirName $shotList {
                if { ![looksLikeAShotDir $dirName] } {
                    set i [lsearch $shotList $dirName]
                    if { $i >= 0 } {
                        set shotList [lreplace $shotList $i $i]
                    }; # end if
                }; # end if
            }; # end foreach
            # puts "After filter, shotList is $shotList"
            set infoSource "(searched directories)"
        }; # end if
        cd $saveDir
        # We actually want the end result with one shot per line.
        set resultString "For facility $td(facility), there are "
        append resultString "[llength $shotList] shots with ASCII data."
        append resultString " $infoSource \n"
        foreach {s0 s1 s2 s3 s4 s5 s6 s7 s8 s9} $shotList {
            append resultString "$s0 $s1 $s2 $s3 $s4 $s5 $s6 $s7 $s8 $s9 \n" 
        }; # end foreach
    } else {
        set resultString "td_server_error : $dataDir is not a directory."
    }; # end if
    return $resultString
}; # end proc


proc sendListOfShots {} {
    # This new version used the same procedure as td_browser.
    global contentString
    append contentString [buildShotList]
}; # end proc


proc sendListOfShots_old_version {} {
    # Go to the facility directory and send a list of subdirectories
    # that contain the actual shot data.
    global td
    global contentString
    set dirName [file join $td(rootDir) $td(facility)]
    if { [file exists $dirName] } {
        cd $dirName
        foreach shotDir [glob *] {
            if { ![string equal $shotDir "rundesc"] && \
	         ![string equal $shotDir "descript"] } {
                # send only real shot directory names
                append contentString "$shotDir\n"
            }; # end if
        }; # end foreach
    } else {
        append contentString \
        "td_server_error : Cannot find facility $facility directory.\n"
    }; # end if
}; # end proc


proc sendRunDescription {} {
    # Look for the Run Description text file in the 
    # RUNDESC (or rundesc) directory.
    # Take care to look for both upper- and lower-case names.

    global td
    global contentString

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
	    append contentString \
		"td_server_error : Could not find Run Descriptions directory.\n"
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
        append contentString \
        "td_server_error : No Run Description file for $td(shot).\n"
    }; # end if
    set fileName [lindex $fileList 0]

    if { [file exists $fileName] } {
        set fp [open $fileName "r"]
        append contentString [read $fp]
        close $fp
    } else {
        append contentString \
        "td_server_error : Run description file $fileName doesn't exist.\n"
    }; # end if
}; # end proc


proc sendChannelList {} {
    global td
    global contentString

    if { [catch {exec gzip -d -c $td(dataFile) > $td(tempDataFile)} result] } {
        append contentString "td_server_error : gunzip command failed\n"
        append contentString "$result\n"
    } else {
        # we have successfully uncompressed the data file
        set fp [open $td(tempDataFile) "r"]
        append contentString [read $fp]
        close $fp
    }; # end if
}; # end proc


proc sendChannelData {} {
    global td
    global contentString
    # Open the file, read the metadata and 
    # reconstruct the time values if they are not already there.

    # puts "send data from file $td(dataFile)"
    if { [catch {exec gzip -d -c $td(dataFile) > $td(tempDataFile)} result] } {
        append contentString "td_server_error : gunzip command failed\n"
        append contentString "$result\n"
    } else {
        # we have successfully uncompressed the data file
        set fp [open $td(tempDataFile) "r"]

        if { [string equal $td(part) raw] } {
            # Send the whole data file as it is.
            append contentString [read $fp]
        } else {
            # Process the data file a little more and 
            # send only the requested piece.
            set commentLine [gets $fp]

            # Extract the header information and save it.
            while { [eof $fp] == 0 } {
                set headerLine [gets $fp]
                # discard comment character
                regsub {#} $headerLine {} headerLine
                # remove surplus whitespace
                set headerLine [string trim $headerLine]
                regsub -all {\s\s+} $headerLine { } headerLine
                if { [string length $headerLine] == 0 } {
                    # we have hit the blank line
                    break; 
                } else {
                    foreach {name value} [split $headerLine] {
                        set td(header.$name) $value
                    }; # end for
                }; # end if
            }; # end while

            if { [string equal $td(part) info] } {
                # Send back just the metadata/info 
                foreach {item} [array names td "header.*" ] {
                    regsub {header\.} $item {} itemLabel
                    append contentString "$itemLabel [set td($item)]\n"
                }; # end foreach
            } else {
                # Send back the data itself.
                set N $td(header.dataPoints)
                set dt $td(header.timeInterval)
                set t0 $td(header.timeStart)
                if {[string equal $td(header.withTimeColumn) "no"]} {
                    set generateTime 1
                } else {
                    set generateTime 0
                }; # end if
                for {set j 0} {$j < $N} {incr j} {
                    set dataLine [gets $fp]
                    if { [eof $fp] } break
                    if { $generateTime } {
                        set t [expr $j * $dt]
                        set v [string trim $dataLine]
                    } else {
                        foreach {t v} [split $dataLine] {}
                    }; # end if
                    append contentString "$t $v\n"
                }; # end for
            }; # end if
        }; # end if

        close $fp
    }; # end if
}; # end proc


proc findShotDataDir {} {
    global td
    set dirName [file join $td(rootDir) $td(facility) $td(shot)]
    if { [file exists $dirName] } {
        set td(dataDir) $dirName
        return 1
    } else {
        set td(dataDir) ""
        return 0
    }; # end if
}; # end proc


proc findChannelDataFile {} {
    global td
    set fname $td(shot)
    if { [string equal $td(channel) list] } {
        append fname A.LST.gz
    } else {
        append fname A . $td(channel) .gz
    }; # end if
    set fileName [file join $td(dataDir) $fname]
    if { [file exists $fileName] } {
        set td(dataFile) $fileName
        return 1
    } else {
        set td(dataFile) ""
        return 0
    }; # end if
}; # end proc


proc processRequest {} {
    global td
    global contentString
    if { [string equal $td(shot) "list"] } {
	# We want a list of the shots for a particular facility.
	sendListOfShots
    } else {    
	# We want the data for a particular shot.
	if { [findShotDataDir] } {
	    if { [string equal $td(channel) "rundescription"] } {
		sendRunDescription
	    } elseif { [findChannelDataFile] } {
		if { [string equal $td(channel) "list"] } {
		    sendChannelList
		} else {
		    sendChannelData
		}; # end if
	    } else {
		append contentString "td_server_error : "
		append contentString "Channel $td(channel) not found.\n"
	    }; # end if
	} else {
	    append contentString "td_server_error : "
	    append contentString "Facility: $td(facility) Shot: $td(shot) "
	    append contentString "not found.\n"
	}; # end if
    }; # end if
};  # end proc


proc allowAccess {} {
    # Provides a basic access filter so that the generic guest
    # login that is embedded in the distributed client files
    # can look at a few example data but not more.
    # We can also use this procedure to apply an access policy
    # more generally.
    global td contentString allowedShotList
    if { $td(cgiScript) } {
	# Impose restrictions on the tdguest only.
	# All others let through by the web-server should be valid users.
	set user $::env(REMOTE_USER)
	if [info exists allowedShotList($user)] {
	    set combinedName $td(facility)+$td(shot)
	    set combinedList [set allowedShotList($user)]
	    if { [lsearch -exact $combinedList $combinedName] >= 0 } {
		# restricted user but selection allowed
		return 1
	    } else {
		# restricted user and invalid selection
		return 0
	    }
	} else {
	    # not a restricted user
	    return 1
	}
    } else {
	# Don't impose any restrictions in a command-line environment.
	# The data must already be on the local disk.
	return 1
    }
}


# -------------- main script starts here -----------------

package require ncgi
package require base64

initializeData
parseQueryString
if { [allowAccess] } {
    processRequest
} else {
    append contentString "td_server_error : "
    append contentString "Access to that data is not allowed.\n"
}
sendContentString

# tidy up
catch { file delete $td(tempDataFile) }

