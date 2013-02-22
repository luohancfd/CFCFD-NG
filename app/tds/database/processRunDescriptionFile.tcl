# processRunDescriptionFile.tcl
#
# Theses procedures will add parts of the run description file 
# to the database
#
# Adapted from Aaron Brandis' code by P. Jacobs, Apr 2004.
#-----------------------------------------------------------------------------

if { $::debugLevel } {
    puts "sourcing processRunDescriptionFile.tcl"
}

# -----------------------------------------------------------------------------
# Presently, all that is left of Aaron's code is this date-conversion procedure.

proc dateConversion {dates} {
    # this proc takes dates in the format of
    # dd/mm/yy, dd/mm/yyyy, dd month yy or
    # dd month yyyy and puts in the mysql
    # date format of yyyy-mm-dd
     
    set date_string $dates
    set month_list {no_month_for_pos_0 January February March April May June \
			July August September October November December}
    # add abrev month_list
    set date_letters [split $dates]
    set k [llength $date_letters]
    set h 0
    for {set j 0} {$j < $k} {incr j 1} {
        if {[llength [lindex [split $dates] $j]] > 0} {
            set date_long($h) [lindex [split $dates] $j]
            incr h
        };# end if
    }; #end for
    for {set e 0} {$e < $h} {incr e 1} {
        for {set i 1} {$i <= 12} {incr i 1} {
            if {$date_long($e) == [lindex $month_list $i] } {
                set month $i
                set day $date_long(0)
                set year $date_long(2)
                if {[string length $year] < 3} {
                    if {$year < 69} {
                        set year [join 20$year]
                    } else {
                        set year [join 19$year]
                    }; #end if
                }; #end if
                set date_string [join $year-$month-$day]
           }; #end if
        }; #end for
    }; #end for

    set abrev_month_list {no_month_for_pos_0 Jan Feb Mar Apr May Jun \
			      Jul Aug Sep Oct Nov Dec}
    set date_letters [split $dates]
    set k [llength $date_letters]
    set h 0
    for {set j 0} {$j < $k} {incr j 1} {
        if {[llength [lindex [split $dates] $j]] > 0} {
            set date_long($h) [lindex [split $dates] $j]
            incr h
        };# end if
    }; #end for
    for {set e 0} {$e < $h} {incr e 1} {  
        for {set i 1} {$i <= 12} {incr i 1} {
            if {$date_long($e) == [lindex $abrev_month_list $i] } {
                set month $i
                set day $date_long(0)
                set year $date_long(2)
                if {[string length $year] < 3} {
                    if {$year < 69} {
                        set year [join 20$year]
                    } else {
                        set year [join 19$year]
                    }; #end if
                }; #end if 
                set date_string [join $year-$month-$day]
           }; #end if
        }; #end for
    }; #end for

    set j [llength [split $dates]]

    for {set u 0} {$u <= $j} {incr u 1} {
        set date_list($u) [lindex [split $dates] $u]
        set date_short $date_list($u)
        set date_short [split $date_short /]
        if {[llength $date_short] > 1} {
            set day [lindex $date_short 0]
            set month [lindex $date_short 1]
            set year [lindex $date_short 2]
            if {[string length $year] < 3} {
               if {$year < 69} {
                   set year [join 20$year]
               } else {
                   set year [join 19$year]
               }; #end if
            }; #end if
            set date_string [join $year-$month-$day]
        };# end if
    }; #end for
    return $date_string
}; #end proc

# --------------------------------------------------------------------------

proc removePrePostBlankLines { listOfLines } {
    # Strip blank lines from start and end of the provided list of lines.
    while { [llength $listOfLines] > 0 && \
	    [string length [string trim [lindex $listOfLines 0]]] == 0 } {
	set listOfLines [lreplace $listOfLines 0 0]
    }
    while { [llength $listOfLines] > 0 && \
	    [string length [string trim [lindex $listOfLines end]]] == 0 } {
	set listOfLines [lreplace $listOfLines end end]
    }
    return $listOfLines
}

proc fieldIsPresent { fieldName listOfPatterns line } {
    # Look for the field name in the first 15 characters of the line.
    # If found, remove it from the line and return a flag to indicate that
    # it was found together with the content of the rest of the line 
    # (i.e. the dataValue).
    # The form of a field name in MONC is
    # Some text......
    # where the number of periods could be a small as 1.
    #
    set startOfLine [string range $line 0 15]
    set found 0
    set dataValue ""
    foreach pattern $listOfPatterns {
	set fullPattern [join [list {^} $pattern {[0-9a-zA-Z _]*\.+}] ""] 
	if [regexp -nocase $fullPattern $startOfLine] {
	    # found pattern, remove it from line
	    regsub -nocase $fullPattern $line {} line
	    set restOfLine [string trim $line]
	    set sixSpaces "      "
	    set startOfSpace [string first $sixSpaces $restOfLine]
	    if { $startOfSpace > 0 } {
		# We seem to have more than the dataValue on the line
		set dataValue [string trim [string range $restOfLine 0 $startOfSpace]]
	    } else {
		set dataValue $restOfLine
	    }
	    if [string equal $fieldName "shot_number"] {
		# We only want the integer part, not any extra characters that
		# the experimenter may have appended to this value.
		regexp {[[:digit:]]+} $dataValue dataValue
		if { [string length $dataValue] == 0 } {
		    # A shot number of zero can be used as a null value
		    set dataValue 0
		}
	    }
	    set found 1
	    break
	}
    }
    return [list $found $dataValue]
}

proc processRunDescriptionFile { rundescFileName } {
    # Reads the specified run-description file and assembles the data for
    # the possible fields.
    # Whatever is left over gets entered as notes.
    # Returns 0 if no fields were found, 1 is at least one field was filled.
    #
    global pattern
    global dataValue
    global nullString

    set fieldNames [array names pattern]
    # Initially, fill field data with SQL Nulls
    foreach name $fieldNames {
	set dataValue($name) $nullString
    }
    set someFullFields 0

    if [catch {open $rundescFileName "r"} result] {
	puts "processRunDescriptionFile: $result"
	return $someFullFields
    }
    set fp $result
    set fileContent [read $fp]
    close $fp
    if { $::debugLevel >= 2 } {
	puts "--------------------"
	puts "fileContent: $fileContent"
	puts "--------------------"
    }
    set fileLines {}
    foreach line [split $fileContent "\n"] {
	lappend fileLines [string trim [textutil::tabify::untabify2 $line]]
    }
    set fileLines [removePrePostBlankLines $fileLines]
    if { $::debugLevel >= 2 } {
	puts "--------------------"
	puts "$fileLines"
	puts "--------------------"
    }

    set deleteFlags {}
    foreach line $fileLines {
	set deleteThisLine 0
	foreach name $fieldNames {
	    set result [fieldIsPresent $name [set pattern($name)] $line]
	    if [lindex $result 0] {
		set deleteThisLine 1
		set someFullFields 1
		set dataValue($name) [lindex $result 1]
		if [string equal $name date] {
		    set dataValue($name) [dateConversion [set dataValue($name)]]
		}
		if { $::debugLevel >= 2 } {
		    puts "<$name>[set dataValue($name)]</$name>"
		}
	    }
	}
	lappend deleteFlags $deleteThisLine
    }
    if { $::debugLevel >= 2 } {
	puts "deleteFlags $deleteFlags"
    }
    set notesLines {}
    for {set i 0} {$i < [llength $deleteFlags]} {incr i} {
	if [lindex $deleteFlags $i] { continue }
	lappend notesLines [lindex $fileLines $i]
    }
    set notesLines [removePrePostBlankLines $notesLines]
    if { $::debugLevel >= 2 } {
	puts "<notes>$notesLines</notes>"
    }
    set dataValue(notes) $notesLines

    return $someFullFields
}

# ------------------------------------------------------------------------------

proc writeXMLFileForShotDescription { fileName } {
    global dataValue
    global fieldNames
    global nullString
    if [catch {open $fileName "w"} result] {
	puts "writeXMLFile: $result"
	return
    } else {
	set fp $result
    }
    puts $fp "<shot_description>"
    foreach name $fieldNames {
	if [string equal $name notes] {
	    # Notes have been stored as a list of strings.
	    puts $fp "<notes>"
	    foreach line [set dataValue(notes)] {
		puts $fp "$line"
	    }
	    puts $fp "</notes>"
	} else {
	    set content [set dataValue($name)]
	    if [string equal $content $nullString] {
		# Do nothing with null value.
	    } else {
		# Write an entry for non-null a value.
		puts $fp "<$name>$content</$name>"
	    }
	}
    }
    puts $fp "</shot_description>"
    close $fp
}

# ------------------------------------------------------------------------------

proc makeFieldString { {all 0} {withType 0} } {
    # Returns a string containing the names of the data fields 
    # separated by commas.
    # If all == 1, all fields are included.  
    # By default, only fields with non-null dataValues are included.
    # If withType == 1, the SQL data type is included with the field name.
    # By default, only the field names are included.
    global sqlType dataValue fieldNames nullString
    set str ""
    foreach name $fieldNames {
	if { $all || ![string equal [set dataValue($name)] $nullString] } {
	    append str "$name"
	    if { $withType } { append str " [set sqlType($name)]" }
	    append str ","
	}
    }
    set str [string trimright $str ","]
    return $str
}

proc quoteTextValue { value } {
    # Remove all single quotes before adding a pair of single quotes
    # around whole string.
    regsub -all {'} $value {} value
    set result "'"
    append result $value "'"
    return $result
}

proc makeValueString { {all 0} } {
    # Returns a string containing the data values, comma separated.
    # If all == 1, even null values are included.
    global dataValue sqlType fieldNames nullString
    set str ""
    foreach name $fieldNames {
	set value [set dataValue($name)]
	if { $all || ![string equal $value $nullString] } {
	    if [string equal [set sqlType($name)] "text"] {
		set value [quoteTextValue $value]
	    }
	    append str "$value,"
	}
    }
    set str [string trimright $str ","]
    return $str
}

proc makeSQLCreateCmd {} {
    set cmd "CREATE TABLE shot_descriptions ([makeFieldString 1 1]);"
    return $cmd
}

proc makeSQLInsertCmd {} {
    set cmd "INSERT into shot_descriptions ([makeFieldString]) VALUES ([makeValueString]);"
    return $cmd
}

# ----------------------------- end --------------------------------------------
