# td_br_data.tcl
# Data manipulation procedures for td_browser.tcl

puts "Start of td_br_data.tcl"


proc readPODFile {} {
    global td
    set fileName [tk_getOpenFile -title "Load Plain-Old-Data File..." \
        -initialdir $td(workDir) ]
    puts "readPODFile: fileName=$fileName"
    if {[string length $fileName] > 0} {
        getFreshData pod $fileName
        set td(workDir) [file dirname $fileName]; # save for next time
    }; # end if    
}; # end proc


proc getFreshData { {typeOfData tds} {fileName {}} } {
    # top level procedure to read the data from a local file
    # or fetch it from the http server.
    # Input: typeOfData : tds  data file is a tds data file with header
    #                     pod  data file is plain old data (two column)

    global td current numberOfTraces nhalf channelNameToNumber

    # timer_reset

    # Intercept the special case where the user is asking for 
    # the lists rather than a valid shot or channel.
    if { [string equal $td(shot) "list"] ||
	 [string equal $td(shot) ""] } {
        displayShotList
        return
    }
    if { [string equal $td(channel) "list"] ||
	 [string equal $td(channel) ""] } {
        displayChannelNames
        return
    }

    if { [string equal $typeOfData tds] } {
        # We are fetching from the tds style files,
        # either locally or from the httpd server.

        # The text in the td(channel) variable may consist of more than one word.
        # Look at just the first.
        set td(channel) [lindex [split [string trim $td(channel)]] 0]

	# The user may have specified the channel by its name rather than
	# its number.  If so try to look up the channel number.
        if {![string is integer $td(channel)]} {
	    updateStatusMessage "Determining channel number."
	    displayChannelNames
            if { [catch { 
		set ch_name $td(channel)
		set ch_id [set channelNameToNumber($ch_name)]
	    } result ] } {
		set msg "Couldn't convert channel name \"$ch_name\" to number."
		puts $msg
		updateStatusMessage $msg
		return
	    }
	    set td(channel) $ch_id
  	}

        if { [string length $td(channel)] == 2 } {
            # Assume that we have card and module numbers only
            # and that we want data stream 0 even if the data
            # was multiplexed. 
            append td(channel) "0"
        }; # end if

        if { [string equal $td(dataSource) directRead] } {
            # puts "Before readChannelData call: [timer_now]"
            updateStatusMessage "Reading data directly."
            set found_data [readChannelData $td(channel) [expr $numberOfTraces + 1]]
            # puts "After readChannelData returns: [timer_now]"
        } else {
	    # This may take some time; put up a busy indicator
	    myBusy hold .
	    update
            updateStatusMessage "Fetching channel header."
            fetchChannelHeader $td(channel)
            updateStatusMessage "Fetching data."
            set found_data [fetchChannelData $td(channel) [expr $numberOfTraces + 1]]
	    # Take down the busy indicator
	    myBusy forget .
	    update
        }; # end if
    } else {
        # Lets try to read plain old data from a file.
        updateStatusMessage "Reading data from plain data file."
        set found_data [readPlainOldData [expr $numberOfTraces + 1] $fileName]
    }; # end if

    updateStatusMessage "Finished getting data."

    if { $found_data == 1 } {
        incr numberOfTraces
        set current $numberOfTraces
        set td(nhalf) 0
        set nhalf($current) $td(nhalf)
        resetOffsetAndScale cumulative $current
        resetOffsetAndScale incremental $current
        if { [string length $td(normalizeChannel)] > 0 } {
            # Assume that we have normalizing data present.
	    normalizeData
        }; # end if
        addTraceToGraph
        refreshGraph $td(rescaleAxesOption) 
    }; # end if

    updateStatusMessage "Ready."
    # puts "After refreshGraph: [timer_now]"
}; # end proc


proc getNormalizingData { {typeOfData tds} {fileName {}} } {
    # top level procedure to read the normalizing data from a local file
    # or fetch it from the http server.
    # Input: typeOfData : tds  data file is a tds data file with header
    #                     pod  data file is plain old data (two column)

    # This function is specifically for the T4 fellows who want to divide
    # the test-section pressure signals by the stagnation pressure signal.
    # We assume that the signals have units of kPa and that the interesting
    # values for the stagnation pressure are in the megaPascal range 
    # (so that the addition of a small amount to avoid divide-by-zero 
    # is not a problem).

    global td current numberOfTraces nhalf channelNameToNumber

    puts "Attempting to get normalizing data channel $td(normalizeChannel)"
    if { [string equal $td(normalizeChannel) "list"] ||
	 [string equal $td(normalizeChannel) ""] } {
        displayChannelNames
        return
    }

    if { [string equal $typeOfData tds] } {
        # We are fetching from the tds style files,
        # either locally or from the httpd server.

        # The text in the td(normalizeChannel) variable may consist of more than one word.
        # Look at just the first.
        set td(normalizeChannel) [lindex [split [string trim $td(normalizeChannel)]] 0]

	# The user may have specified the channel by its name rather than
	# its number.  If so try to look up the channel number.
        if {![string is integer $td(normalizeChannel)]} {
	    updateStatusMessage "Determining normalizing channel number."
	    displayChannelNames
            if { [catch { 
		set ch_name $td(normalizeChannel)
		set ch_id [set channelNameToNumber($ch_name)]
		puts "Using channel $ch_id as normalising channel."
	    } result ] } {
		set msg "Couldn't convert channel name \"$ch_name\" to number."
		puts $msg
		updateStatusMessage $msg
		return
	    }
	    puts "Using channel $ch_id as normalising channel."
	    set td(normalizeChannel) $ch_id
  	}

        if { [string length $td(normalizeChannel)] == 2 } {
            # Assume that we have card and module numbers only
            # and that we want data stream 0 even if the data
            # was multiplexed. 
            append td(normalizeChannel) "0"
        }; # end if

        if { [string equal $td(dataSource) directRead] } {
            updateStatusMessage "Reading data directly."
            set found_data [readChannelData $td(normalizeChannel) "normalizeChannel"]
        } else {
	    # This may take some time; put up a busy indicator
	    myBusy hold .
	    update
            updateStatusMessage "Fetching normalizing channel header."
            fetchChannelHeader $td(normalizeChannel)
            updateStatusMessage "Fetching normalizing data."
            set found_data [fetchChannelData $td(normalizeChannel) "normalizeChannel"]
	    # Take down the busy indicator
	    myBusy forget .
	    update
        }; # end if
    } else {
        # Lets try to read plain old data from a file.
        updateStatusMessage "Reading normalizing data from plain data file."
        set found_data [readPlainOldData "normalizeChannel" $fileName]
    }; # end if

    updateStatusMessage "Finished getting normalizing data."

    # Clean up the normalizing data so that it is all positive.
    ::vVecnormalizeChannel expr "abs(::vVecnormalizeChannel) + 1.0e-6"

    updateStatusMessage "Ready."
}; # end proc


proc computeAverageValue {} {
    global td
    global current
    set i1 $td(t1Index)
    set i2 $td(t2Index)
    if { [string length $i1] > 0 && \
         [string length $i2] > 0 } {
        if { $i2 < $i1 } { set ii $i1; set i1 $i2; set i2 $ii }
        ::plotYVec$current dup ::tmpV
        set td(averageValue) [blt::vector expr mean(::tmpV($i1:$i2))]
        set td(stddevValue) [blt::vector expr sdev(::tmpV($i1:$i2))]
    } else {
        set td(averageValue) ""
    }; # end if
}; # end proc


proc normalizeData {} {
    global td
    global current
    if { [string length $td(normalizeChannel)] > 0 &&
	 [info exists ::vVecnormalizeChannel] } {
	# Proceed to normalize the data channel.
	::vVec$current dup ::vdata
	::tVec$current dup ::tdata
	set dt_data [expr $::tdata(1) - $::tdata(0)]
	::vVecnormalizeChannel dup ::vnorm
	::tVecnormalizeChannel dup ::tnorm
	set dt_norm [expr $::tnorm(1) - $::tnorm(0)]
	puts "dt_data= $dt_data  dt_norm= $dt_norm"

	# If necessary, interpolate the less-frequently-sampled trace.
	set dt_ratio [expr $dt_data / $dt_norm]
	if { $dt_ratio < 0.9 } {
	    puts " Data channel is sampled faster than normalizing channel."
	    set density [expr round($dt_norm / $dt_data)]
	    ::vnorm populate ::vnorm_new $density
	    ::tnorm populate ::tnorm_new $density
	    ::vnorm_new dup ::vnorm
	    ::tnorm_new dup ::tnorm
	    set dt_norm [expr $::tnorm(1) - $::tnorm(0)]
	} elseif { $dt_ratio > 1.1 } {
	    puts "Data channel is sampled slower than normalizing channel."
	    set density [expr round($dt_data / $dt_norm)]
	    ::vdata populate ::vdata_new $density
	    ::tdata populate ::tdata_new $density
	    ::vdata_new dup ::vdata
	    ::tdata_new dup ::tdata
	    set dt_data [expr $::tdata(1) - $::tdata(0)]
	}

	# Time-shifting the normalizing data, if necessary.
	set tshift $td(timeShift)
	if { [string length [string trim $tshift]] == 0 } {
	    set tshift 0.0
	}
	if { [string is integer $tshift] || \
		 [string is double $tshift)] } {
	    puts "Numeric time shift= $tshift"
	    if { $tshift >= 0 } {
		# We will work with later portions of the normalizing trace.
		set nshift [expr round($tshift / $dt_norm)]
		::vnorm expr "::vnorm($nshift:end)"
		::tnorm expr "::tnorm($nshift:end)"
	    } else {
		# We will work with earlier times of the normalizing trace.
		set nshift [expr round(abs($tshift) / $dt_data)]
		::vdata expr "::vdata($nshift:end)"
		::tdata expr "::tdata($nshift:end)"
	    }
	} else {
	    set nshift 0
	}
	# We may end up with a shortened data vector so
	# we need to make the vector lengths consistent.
	set min_length [::vdata length]
	if { [::vnorm length] < $min_length } {
	    set min_length [::vnorm length]
	}
	if { [::vdata length] > $min_length } {
	    ::vdata length $min_length
	    ::tdata length $min_length
	}
	if { [::vnorm length] > $min_length } {
	    ::vnorm length $min_length
	    ::tnorm length $min_length
	}

	# Replace original data with normalized data.
	::vVec$current expr "::vdata / ::vnorm"
	::tVec$current expr "::tdata"

	blt::vector destroy ::vdata ::tdata ::vnorm ::tnorm
    }; # end if
}; # end proc


proc filterData {} {
    global td current

    if { $current == 0 } { return }

    # Apply the moving-average filter.
    if { $td(nhalf) > 0 } {
        # This can be a slow process; put up a busy indicator
        myBusy hold .
        update

        ::vVec$current dup ::vVec

	if { $td(fastFilter) == 1 } {
            cfilter ::vVec $td(nhalf)
            ::vVec dup ::plotYVec$current
        } else {
            if { [info exists ::fVec] } {
                blt::vector destroy ::fVec
            }; # end if
            blt::vector create ::fVec
            set i0 0
            set ilast [expr [::vVec length] - 1]
            for { set i $i0 } { $i <= $ilast } { incr i } {
                set i1 [expr $i - $td(nhalf)]
                if { $i1 < $i0 } { set i1 $i0 }
                set i2 [expr $i + $td(nhalf)]
                if { $i2 > ($ilast - 1) } { set i2 [expr $ilast - 1] }
                # Yes, this is the slow way to do it -- but, it is easy.
                # One day, we should try to do it via the C API so that
                # we can implement a much faster filter.
                ::fVec append [blt::vector expr mean(::vVec($i1:$i2))]
            }; # end for
            ::fVec dup ::plotYVec$current
        }; # end if
    
        # Remove the busy indicator
        myBusy forget .
        update
    } else {
        # For a zero point/one point average, 
        # just copy back the original data.
        ::vVec$current dup ::plotYVec$current
    }; # end if
    puts "Filtered data with $::td(nhalf) points."
    set ::nhalf($current) $::td(nhalf)

    # Get a clean copy of the sample-times vector
    # so that the offsets and scales can be consistently applied.
    ::tVec$current dup ::plotXVec$current
    applyOffsetAndScale cumulative $current
}; # end proc 


proc compute_heat_transfer {} {
    global td current
    global calcfrom
    global calcto
    global thermprod
    global transens
    global ambvolt

    if { $current == 0 } { return }

    # Apply the heat-transfer calculation -- Josh Corbett, 2004.

    # This can be a slow process; put up a busy indicator
    myBusy hold .
    update

    # Get a clean copy of the sample-times vector
    # and the (presumably) temperature data
    ::tVec$current dup ::ttVec
    ::vVec$current dup ::vVec
    
    set dt [expr $::ttVec(1) - $::ttVec(0)]

    if {$calcfrom == 3 && $calcto == 1} {
         c_heat_transfer_emf ::vVec $dt $thermprod $ambvolt $transens
         # Put the qdot data up on the graph (as y-coordinates)
         ::vVec dup ::vVec$current
         ::vVec dup ::plotYVec$current
         set td(yAxisLabel) "W/m*"
         relabelGraphAxes
         set calcfrom 1
         set calcto 2
    } elseif {$calcfrom == 1 && $calcto == 2} {
    	   c_surface_temperature ::vVec $dt $thermprod $ambvolt $transens
         # Put the qdot data up on the graph (as y-coordinates)
         ::vVec dup ::vVec$current
         ::vVec dup ::plotYVec$current
         set td(yAxisLabel) "Degree K"
         relabelGraphAxes
         set calcfrom 2
         set calcto 1
    } elseif {$calcfrom == 2 && $calcto == 1} {
    	   c_heat_transfer_temp ::vVec $dt $thermprod $ambvolt $transens
         # Put the qdot data up on the graph (as y-coordinates)
         ::vVec dup ::vVec$current
         ::vVec dup ::plotYVec$current
         set td(yAxisLabel) "W/m*"
         relabelGraphAxes
         set calcfrom 1
         set calcto 2
    } elseif {$calcfrom == 3 && $calcto == 2} {
    	   c_emf_temp ::vVec $dt $thermprod $ambvolt $transens
         # Put the qdot data up on the graph (as y-coordinates)
         ::vVec dup ::vVec$current
         ::vVec dup ::plotYVec$current
         set td(yAxisLabel) "Degree K"
         relabelGraphAxes
         set calcfrom 2
         set calcto 1
    }

  
    # Remove the busy indicator
    myBusy forget .
    update

    # Get a clean copy of the sample-times vector
    # so that the offsets and scales can be consistently applied.
    ::tVec$current dup ::plotXVec$current
    applyOffsetAndScale cumulative $current
}; # end proc 


proc saveSelectedData { {requestedFormat gnuplot} } {
    # requestedFormat is one of gnuplot, csv
    global td current metadataList
    global current
    set fileName [tk_getSaveFile -title "Save Selected Data to..." \
        -initialdir $td(workDir) ]
    if {[string length $fileName] > 0} {
        set i1 $td(t1Index)
        if { [string length $i1] == 0 } { set i1 0 }
        set i2 $td(t2Index)
        if { [string length $i2] == 0 } { 
            set i2 [expr [::plotYVec$current length] - 1]
        }; # end if
        if { $i1 > $i2 } { set i $i1; set i1 $i2; set i2 $i }
        ::plotXVec$current dup ::tmpX
        ::plotYVec$current dup ::tmpY
        set fp [open $fileName "w"]
        if { $td(includeMetaDataWhenSaving) } {
            foreach metadataName $metadataList {
		if {[info exists td(header.$metadataName)]} {
		    set value [set td(header.$metadataName)]
		    if {[string equal $requestedFormat gnuplot]} {
			puts $fp "\# $metadataName $value"
		    } elseif {[string equal $requestedFormat csv]} {
			if {[string is integer $value] || [string is double $value]} {
			    puts $fp "\"$metadataName\", $value"
			} else {
			    # Put non-numeric value out as a quoted string.
			    puts $fp "\"$metadataName\", \"$value\""
			}
		    } else {
			# may have other formats later
		    }
		}
            }; # end foreach
        }; # end if
        for { set i $i1 } { $i <= $i2 } { incr i } {
	    if {[string equal $requestedFormat gnuplot]} {
		puts $fp "[set ::tmpX($i)] [set ::tmpY($i)]"
	    } elseif {[string equal $requestedFormat csv]} {
		puts $fp "[set ::tmpX($i)], [set ::tmpY($i)]"
	    } else {
		# may have other formats later
	    }
        }; # end for
        close $fp
        puts "Elements $i1 through $i2 of trace $current written to $fileName."
        set td(workDir) [file dirname $fileName]; # save for next time
    }; # end if    
}; # end proc


proc resetOffsetAndScale { whichSet whichTrace } {
    # Set the offsets and scales for a particular trace.
    # whichSet can take values incremental
    #                          cumulative
    global incOS cumOS

    set toffset 0.0
    set tscale  1.0
    set yoffset 0.0
    set yscale  1.0

    if { [string equal $whichSet incremental] } {
        set incOS($whichTrace) [list $toffset $tscale $yoffset $yscale]
    } else {
        set cumOS($whichTrace) [list $toffset $tscale $yoffset $yscale]
    }; # end if
}; # end proc


proc applyOffsetAndScale { whichSet whichTrace } {
    # Apply offsets and scaling to a specified trace.
    # whichSet may take values incremental
    #                          cumulative
    global incOS cumOS current

    if { $current == 0 } { return }
    puts "Apply OS set $whichSet to trace $whichTrace"

    if { [string equal $whichSet incremental] } {
        # When applying an incremental change, keep track of
        # the cumulative change.
        set toffset [lindex [set incOS($whichTrace)] 0]
        set tscale  [lindex [set incOS($whichTrace)] 1]
        set yoffset [lindex [set incOS($whichTrace)] 2]
        set yscale  [lindex [set incOS($whichTrace)] 3]
        set cum_toffset [lindex [set cumOS($whichTrace)] 0]
        set cum_tscale  [lindex [set cumOS($whichTrace)] 1]
        set cum_yoffset [lindex [set cumOS($whichTrace)] 2]
        set cum_yscale  [lindex [set cumOS($whichTrace)] 3]
        set cum_toffset \
            [expr $cum_toffset + $toffset * $cum_tscale]
        set cum_tscale [expr $cum_tscale * $tscale] 
        set cum_yoffset \
            [expr $cum_yoffset + $yoffset * $cum_yscale]
        set cum_yscale [expr $cum_yscale * $yscale]
        set cumOS($whichTrace) \
            [list $cum_toffset $cum_tscale $cum_yoffset $cum_yscale]
    } else {
        # Apply cumulative offsets and scales for this trace.
        set toffset [lindex [set cumOS($whichTrace)] 0]
        set tscale  [lindex [set cumOS($whichTrace)] 1]
        set yoffset [lindex [set cumOS($whichTrace)] 2]
        set yscale  [lindex [set cumOS($whichTrace)] 3]
    }; # end if
    ::plotXVec$whichTrace expr \
        [list (::plotXVec$whichTrace - $toffset) * $tscale ]
    ::plotYVec$whichTrace expr \
        [list (::plotYVec$whichTrace - $yoffset) * $yscale ]
}; # end proc


proc restoreOriginalTrace {} {
    global td current
    if { $current == 0 } { return }
    resetOffsetAndScale incremental $current
    resetOffsetAndScale cumulative $current
    set td(nhalf) 0
    ::tVec$current dup ::plotXVec$current
    ::vVec$current dup ::plotYVec$current
}; # end proc


proc snapToReferenceTime { t } {
    global td incOS current
    set toffset [expr $t - $td(referenceTime) ]
    set incOS($current) [list $toffset 1.0 0.0 1.0]
    applyOffsetAndScale incremental $current
}; # end proc


proc snapToReferenceValue { y } {
    global td incOS current
    set yoffset [expr $y - $td(referenceValue) ]
    set incOS($current) [list 0.0 1.0 $yoffset 1.0]
    applyOffsetAndScale incremental $current
}; # end proc

puts "End of td_br_data.tcl"
