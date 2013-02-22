#!/bin/sh
#\
exec tclsh "$0" ${1+"$@"}

# File:
# td_web_client.tcl
#
# Purpose:
# This script (when called from octave, say) will go and
# ask for some data from the td_server.tcl CGI script.
# 
# Usage:
# td_web_client.tcl <facility> <shot> <channel> <part>
#
# Author: 
# Peter Jacobs
# Department of Mechanical Engineering
# The University of Queensland
#
# Revisions:
# 26-Dec-01 original code
# 29-Dec-01 extra command-line options
# 30-Dec-01 add documentation

package require http

# the facility name and shot number should have been passed to
# this script as command line arguments
if {$argc == 0} {
    # no arguments passed, set defaults
    set facility T4
    set shot 7319
    set channel 622
    set part data
} else {
    # octave seems to pass only one argument on the command line.
    # It may have spaces and so can be pulled apart as a list.
    set facility [lindex $argv 0]
    set shot [lindex $argv 1]
    set channel [lindex $argv 2]
    set part [lindex $argv 3]
}; # end if

set queryString "facility=$facility"
append queryString "+shot=$shot"
append queryString "+channel=$channel"
append queryString "+part=$part"

set serverIP "130.102.240.102"
set serverURL "http://$serverIP"
append serverURL "/cgi-bin/td_server.tcl?"
append serverURL $queryString

set matFile temporary.mat

set r [http::geturl $serverURL]
set rcode [lindex [http::code $r] 1]

if {$rcode == 200} {
    # we received the data OK
    set rdata [http::data $r]
    set fp [open $matFile w]
    puts $fp $rdata
    close $fp
} else {
    puts "Data transfer failed with http code $rcode"
}; # end if
