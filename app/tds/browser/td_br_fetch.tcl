# td_br_fetch.tcl
# Procedures that get data from the td_server script.
#

proc fetchSomething { queryString } {
    global td
    # Given the query string, call upon the CGI script to
    # send the data.
    if { [string equal $td(dataSource) httpServer] } {
        set serverURL "http://$td(serverIP)"
        append serverURL "/cgi-bin/tds/td_server.tcl?"
        append serverURL $queryString

        set headerList [list "Authorization" \
            [concat "Basic " [base64::encode $td(username):$td(password)]]]

        if { $td(useProxy) == 1 } {
            http::config -proxyhost $td(proxyhost) -proxyport $td(proxyport)
        } else {
            http::config -proxyhost "" -proxyport ""
        }; # end if

        # Fetch the data
        set r [http::geturl $serverURL -headers $headerList]
        set rcode [lindex [http::code $r] 1]

        if {$rcode == 200} {
            # we received the data OK
            set rdata [http::data $r]
            # puts $rdata
        } else {
            puts "Data transfer failed with http code $rcode"
            set rdata {}
        }; # end if

        http::cleanup $r
    } else {
        # Use the server as a local script and get the data
        # via stdout.
        set server [file join $td(serverScriptHome) td_server.tcl]
        puts "Use local script: $server"
        if { [catch {exec tclsh $server $queryString } rdata ] } {
            puts "Server script failed: $rdata"
            set rdata {}
        }; # end if
    }; # end if

    if { [regexp td_server_error $rdata] > 0 } {
        puts $rdata
	set rdata {}
    }; # end if
    return $rdata
}; # end proc


proc fetchChannelNames {} {
    global td
    # Set up the request
    set queryString "facility=$td(facility)"
    append queryString "+shot=$td(shot)"
    append queryString "+channel=list"
    # get the data and return it
    return [fetchSomething $queryString]
}; # end proc


proc fetchShotList {} {
    global td
    # Set up the request
    set queryString "facility=$td(facility)"
    append queryString "+shot=list"
    # get the data and return it
    return [fetchSomething $queryString]
}; # end proc


proc fetchRunDescription {} {
    global td
    # Set up the request
    set queryString "facility=$td(facility)"
    append queryString "+shot=$td(shot)"
    append queryString "+channel=rundescription"
    # get the data and return it
    return [fetchSomething $queryString]
}; # end proc


proc fetchChannelHeader { channel_id } {
    global td
    # Set up the request
    set queryString "facility=$td(facility)"
    append queryString "+shot=$td(shot)"
    append queryString "+channel=$channel_id"
    append queryString "+part=info"

    set td(headerText) [fetchSomething $queryString]
    # puts td(header)
    foreach line [split $td(headerText) "\n"] {
        foreach {name value} [split $line] {
            set td(header.$name) $value
        }; # end foreach
    }; # end foreach
    # foreach name [array names td] { puts "$name = [set td($name)]" }
}; # end proc


proc fetchChannelData { channel_id newIndex } {
    # Fetch the channel data from the server script 
    # (either local or via http request).
    # If everything is successful, a value of 1 is returned and
    # the channel's data is contained in two new global vectors
    # ::tVec$newIndex and ::vVec$newIndex.

    global td

    # Set up the request
    set queryString "facility=$td(facility)"
    append queryString "+shot=$td(shot)"
    append queryString "+channel=$channel_id"
    append queryString "+part=data"
    set time_series_data [fetchSomething $queryString]

    if { $time_series_data == {} } {
        tk_messageBox -type ok -icon error -title td_browser \
            -message "No data was found."
        return 0
    } else {
        # Put the new data into fresh vectors
        if { [info exists ::tVec$newIndex] } {
            rbc::vector destroy ::tVec$newIndex
        }; # end if
        if { [info exists ::vVec$newIndex] } {
            rbc::vector destroy ::vVec$newIndex
        }; # end if
        rbc::vector create ::tVec$newIndex
        rbc::vector create ::vVec$newIndex
	if {[string length $td(currentTimeUnits)] == 0} {
	    set td(currentTimeUnits) $td(header.timeUnits)
	}
	set multiplyFactor [getTimeFactor $td(header.timeUnits) $td(currentTimeUnits)]
        set tlist {}
        set vlist {}
        foreach {t v} $time_series_data {
            lappend tlist [expr $t * $multiplyFactor]
            lappend vlist $v
        }; # end foreach
        ::tVec$newIndex set $tlist
        ::vVec$newIndex set $vlist
        return 1
    }; # end if

}; # end proc

