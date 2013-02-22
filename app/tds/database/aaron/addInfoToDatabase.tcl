# Theses procedures will add the required shots to the database

#--------------------------------------------------------------

proc readInData {data_field line} {
    global f
    # f counts the number of data_fields
    
    if {[string first $data_field $line] != -1} {
        set a [expr 2 + [string last .. $line]]
        incr f
        set data_names($f) $data_field
	set data_field [string range $line $a end]
        set data_info($f) $data_field
    }; # end if
    set temp($f) "{$data_names($f)} $data_info($f)"
    return $temp($f)
}; # end proc

proc extractInfo {line run_no} {
    # might add [lindex $argv 1 to passed args
    # to see if X2
    global f ;# f counts the number of data fields
    # Looks for data fields and stores the required 
    # information into their corresponding variables
    # with the read_in procedure    
 
    set b [expr [string first .. $line] - 1 ] ;
    # b is the position of the end of the data field
    # it finds the first location of ".." because
    # the preceeding characters are the data fields

    set data_field [string range $line 0 $b]
    set data_field_length [string length $data_field]
    
    # run_no specific section of code
    # this section of code is required because 
    # the run number data field often
    # only has 1 "." after it
    if {[string first Run $line] != -1 && [string first numb $line] != -1 || [string first Numb $line] != -1 } {    	 
	set a [expr 1 + [string last . $line]]
	set b [expr [string first . $line] - 1 ]
	set data_field [string range $line 0 $b] 
        incr f
	set data_names_run_no($f) $data_field
        set run_no_f $f
	set run_no [string range $line $a end]
	set data_info_run_no($f) $run_no
	set temp($run_no_f) "{$data_names_run_no($run_no_f)} $data_info_run_no($run_no_f)"
	return $temp($run_no_f)
	# end of run_no specific section of code
    
    # shocktube specific section of code
    # this section of code is required because 
    # the run number data field often
    # only has 1 "." after it in X2
    } elseif {[string first Shock $line] != -1 || [string first shock $line] != -1} {; #used to have lindex argv = "X2" as well
       	set a [expr 1 + [string last . $line]]
	set b [expr [string first . $line] - 1 ]
	set data_field [string range $line 0 $b] 
        incr f
	set data_names_shock($f) $data_field
        set shock_f $f
	set shocktube [string range $line $a end]
	set data_info_shock($f) $shocktube
	set temp($shock_f) "{$data_names_shock($shock_f)} $data_info_shock($shock_f)"
	return $temp($shock_f)
	# end of shocktube specific section of code 
    } elseif {$data_field_length > 0 } {
            set result [readInData $data_field $line]
	    set data_names($f) [lindex $result 0]
	    set data_info($f) [lrange $result 1 end]
	    set temp "{$data_names($f)} $data_info($f)"
	    return $temp
    }; #end if

}; # end proc

#------------------------------------------------------------------------

proc extractInfoFromRundescAscii {arg1 spec_f} {
    global f
    # f count the number of data fields
    set f 0
    
    global line_counter
    
    global extra ; 
    # used to capture the information following
    # the data fields in the rundesc text file

    set run_no [file tail $arg1]
    set dir [file dirname $arg1]
    gets $spec_f
    set i 1; # counter
    set line_counter 1
    foreach line [split [read $spec_f] \n] {
	set extra($line_counter) $line
	set result [extractInfo $line $run_no]   
        if {$i <= $f} {
            set data_names($i) [lindex $result 0]
            set data_info($i) [lrange $result 1 end]
            set temp($i) "{$data_names($i)} $data_info($i)"
            incr i 	    
	}; #end if
        incr line_counter
    }; #end foreach
    set ascii_data ""
    for {set j 1} {$j <= $f} {incr j 1} {
        lappend ascii_data [lindex $temp($j) 0]
	# ascii data stores the rundesc information
        # in the following format: data_names(1), data_info(1), 	
        # data_names(2), data_info(2), data_names(3)
	# data_info(3) etc
        lappend ascii_data [lrange $temp($j) 1 end]
    }; #end for
    return $ascii_data 
}; # end proc


proc setTransducerCodes {no_of_channels spec_f desc_f transCodeTable i} {
    # From Descript, extract transducer Info

    gets $desc_f line
    set trans_no($i) $line
    gets $desc_f line
    set position($i) $line
    gets $desc_f line
    set sensitivity($i) $line
    gets $desc_f line
    set channel_no($i) $line
    gets $desc_f line
    set gain($i) $line

    # For each transducer number, there is a corresponding
    # code used the following code replaces the current
    # transducer number with the corresponding code which
    # is optained from the trans_code_table 
  
    set trans_no($i) [lindex $transCodeTable $trans_no($i)]

    # End of translanting transducer number to transducer code
    set temp "$trans_no($i) $position($i) $sensitivity($i) $channel_no($i) \
    $gain($i)"
    return $temp
}; # end proc


proc extractScramjetParameters {desc_f i} {
    gets $desc_f line
    set scram($i) $line
    return $scram($i)
}; # end proc


proc extractTransducerNames {no_of_channels desc_f i} {
    if {[gets $desc_f line] == 0} {
        set trans_name($i) "Unknown"
    } else { 
        set trans_name($i) $line 
    }; # end if

    return $trans_name($i)
}; # end proc


proc extractTransducerSerialNo {no_of_channels i desc_f} {
    if {[gets $desc_f line] == 0} {
        set serial_no($i) "Unknown"
    } else {
        set serial_no($i) $line
    }; # end if

    return $serial_no($i)
}; # end proc

proc dateConversion {dates} {
    # this proc takes dates in the format of
    # dd/mm/yy, dd/mm/yyyy, dd month yy or
    # dd month yyyy and puts in the mysql
    # date format of yyyy-mm-dd
     
    set date_string $dates
    set month_list {no_month_for_pos_0 January February March April May June July August September October November \
    December}
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

    set abrev_month_list {no_month_for_pos_0 Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec}
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


#--------------------------- start of main script -------------------------
#puts "Start of processing file"

global f; # f is used to count the number of data fields

global line_counter

set dir [file dirname $arg1]

# Extract Information from descript directory

if {$desc_f != "No"} {
    gets $desc_f line
    set no_of_channels $line

    set transCodeTable { no_code_4_number_0 p1 p2 q1 q2 sf1 sf2 sp ss m1 m2 d }

    set i 0
    while {$i < $no_of_channels} {
        set result [setTransducerCodes $no_of_channels $spec_f $desc_f $transCodeTable $i]
        set trans_no($i) [lindex $result 0]
        set position($i) [lindex $result 1]
        set sensitivity($i) [lindex $result 2]
        set channel_no($i) [lindex $result 3]
        set gain($i) [lindex $result 4]
        incr i
    }; #end while

    set i 0
    while {$i < 13} {
        set scram($i) [extractScramjetParameters $desc_f $i]
        incr i
    }; #end while

    set i 0
    while {$i < $no_of_channels} {
        set trans_name($i) [extractTransducerNames $no_of_channels $desc_f $i]
        incr i
    }; #end while
  
    set i 0
    while {$i < $no_of_channels} {
        set serial_no($i) [extractTransducerSerialNo $no_of_channels $i $desc_f]
        incr i
    }; #end while
}; #end if

if {$spec_f != "No"} {
    set result [extractInfoFromRundescAscii $arg1 $spec_f]
    for {set i 1} {$i <= $f} {incr i 1} {
        set data_names($i) [lindex $result [expr 2*($i-1)]] ;#2*($i-1) to get all
                                                             #the even positions
		    					 #ie 0,2,4
        set data_info($i) [lindex $result [expr 2*($i-1)+1]];#2*($i-1)+1 to get all
                                                         #the even positions
							 #ie 1,3,5						 
    };# end for

    set i 1
    set notes ""

    while {$i < $line_counter} {
        set notes [concat $notes $extra($i)];
        # stores all the extra data (which is one line) 
        # into one variable so it is easier to search
        incr i
    }; #end while

    # Begin writing txt file

    # set all previous values in sql_data_info to mysql's null character \N
    set i 1
    set max_no_data_fields 19 
    while {$i <= $max_no_data_fields} {; # this number has to be equal to the max 
        # number of data fields of all tubes.  if you add a column
        # make sure this number is still equal to the max number
        # of data fields
        set sql_data_info($i) {\N}
        incr i
    }; # end while

    # If another shock tunnel or expansion tube needs to placed into the 
    # script another data names block also needs to be added for it here

    if {[lindex $argv 1] == "T4"} {
        set i 1 
        while {$i <= $f} {
            if {[regexp -nocase {(Proj)} $data_names($i)]} {
                set sql_data_info(1) $data_info($i)
            } elseif {[regexp -nocase {(Run)} $data_names($i)] || [regexp -nocase {(Test)} $data_names($i)]} {
                set sql_data_info(2) $data_info($i)
            } elseif {[regexp -nocase {(Dat)} $data_names($i)]} {
	        set sql_data_info(3) [dateConversion $data_info($i)]
            } elseif {[regexp -nocase {(Bla)} $data_names($i)]} {
                set sql_data_info(4) $data_info($i)
            } elseif {[regexp -nocase {(Reser)} $data_names($i)]} {
                set sql_data_info(5) $data_info($i)
            } elseif {[regexp -nocase {(Driv)} $data_names($i)]} {
                set sql_data_info(6) $data_info($i)
            } elseif {[regexp -nocase {(Diap)} $data_names($i)]} {
                set sql_data_info(7) $data_info($i)
            } elseif {[regexp -nocase {(Shoc)} $data_names($i)]} {
                set sql_data_info(8) $data_info($i)
            } elseif {[regexp -nocase {(Nozz)} $data_names($i)]} {
                set sql_data_info(9) $data_info($i)
            } elseif {[regexp -nocase {(Acce)} $data_names($i)]} {
                set sql_data_info(10) $data_info($i)
            } elseif {[regexp -nocase {(Jack)} $data_names($i)]} {
                set sql_data_info(11) $data_info($i)
            } elseif {[regexp -nocase {(Brid)} $data_names($i)]} {
                set sql_data_info(12) $data_info($i)
            } elseif {[regexp -nocase {(Evac)} $data_names($i)]} {
                set sql_data_info(13) $data_info($i)
            } elseif {[regexp -nocase {(Pist)} $data_names($i)]} {
                set sql_data_info(14) $data_info($i)
            } elseif {[regexp -nocase {(Pressure)} $data_names($i)]} {
                set sql_data_info(15) $data_info($i)
            } elseif {[regexp -nocase {(Electron)} $data_names($i)]} {
                set sql_data_info(16) $data_info($i)
            } elseif {[regexp -nocase {(Skimme)} $data_names($i)]} {
                set sql_data_info(17) $data_info($i)
            } elseif {[regexp -nocase {(Condition)} $data_names($i)]} {
                set sql_data_info(18) $data_info($i)
            } else {
                set maintenance_problems [open [file join $dir/maintenance_problems_T4.txt] a]
                puts $maintenance_problems "add new column to table: $data_names($i), data from file $files and run number $run_no"
                close $maintenance_problems
            }; #end if    
            incr i
        }; #end while
        
        # If you want to add a column, just add another elseif to the list 
        # above.  Make sure you change all the sql_data_info(x) numbers -
        # the number x - to be in the right order.  Don't forget to increase
        # the number for notes as well.
    
        set max_no_T4_fields 19

        set sql_data_info($max_no_T4_fields) $notes

        set i 1
        while {$i <= $max_no_T4_fields} {; # set to number of data fields for T4
            puts -nonewline $txt_file_for_mysql "$sql_data_info($i)	"
            incr i
        };# end while
        puts $txt_file_for_mysql "\n"
    }; #end if

    if {[lindex $argv 1] == "X2"} {
        set i 1 
        while {$i <= $f} {
            if {[regexp -nocase {(Proj)} $data_names($i)]} {
                set sql_data_info(1) $data_info($i)
            } elseif {[regexp -nocase {(Run)} $data_names($i)]} {
                set sql_data_info(2) $data_info($i)
            } elseif {[regexp -nocase {(Dat)} $data_names($i)]} {
	        set sql_data_info(3) [dateConversion $data_info($i)]
            } elseif {[regexp -nocase {(Bla)} $data_names($i)]} {
                set sql_data_info(4) $data_info($i)
            } elseif {[regexp -nocase {(Reser)} $data_names($i)]} {
                set sql_data_info(5) $data_info($i)
            } elseif {[regexp -nocase {(Driv)} $data_names($i)]} {
                set sql_data_info(6) $data_info($i)
            } elseif {[regexp -nocase {(Diap)} $data_names($i)]} {
                set sql_data_info(7) $data_info($i)
            } elseif {[regexp -nocase {(Shoc)} $data_names($i)]} {
                set sql_data_info(8) $data_info($i)
            } elseif {[regexp -nocase {(Nozz)} $data_names($i)]} {
                set sql_data_info(9) $data_info($i)
            } elseif {[regexp -nocase {(Acce)} $data_names($i)] || [regexp -nocase {(Acl)} $data_names($i)]} {
                set sql_data_info(10) $data_info($i)
            } elseif {[regexp -nocase {(Opti)} $data_names($i)]} {
                set sql_data_info(11) $data_info($i)
            } else {
                set maintenance_problems [open [file join $dir/maintenance_problems_X2.txt] a]
                puts $maintenance_problems "add new column to table: $data_names($i), data from file $files and run number $run_no"
                close $maintenance_problems
            }; #end if    
            incr i
         }; #end while

        set max_no_X2_fields 12

        set sql_data_info($max_no_X2_fields) $notes

        set i 1

        while {$i <= $max_no_X2_fields} {; #set to number of data fields of X2
            puts -nonewline $txt_file_for_mysql "$sql_data_info($i)	"
            incr i
        };# end while
        puts $txt_file_for_mysql "\n"
    }; #end if
}; #end if


if {[lindex $argv 1] == "X3" && $spec_f != "No"} {
    set i 1
    while {$i <= $f} {
        if {[regexp -nocase {(Proj)} $data_names($i)]} {
            set sql_data_info(1) $data_info($i)
        } elseif {[regexp -nocase {(Run)} $data_names($i)]} {
            set sql_data_info(2) $data_info($i)
        } elseif {[regexp -nocase {(Dat)} $data_names($i)]} {
            set sql_data_info(3) [dateConversion $data_info($i)]
        } elseif {[regexp -nocase {(Bla)} $data_names($i)]} {
            set sql_data_info(4) $data_info($i)
        } elseif {[regexp -nocase {(Reser)} $data_names($i)]} {
            set sql_data_info(5) $data_info($i)
        } elseif {[regexp -nocase {(2nd D)} $data_names($i)]} {
            set sql_data_info(6) $data_info($i)
        } elseif {[regexp -nocase {(Driver)} $data_names($i)]} {
            set sql_data_info(7) $data_info($i)
        } elseif {[regexp -nocase {(Diaphragm2)} $data_names($i)] || [regexp -nocase {(Dia2)} $data_names($i)]} {
            set sql_data_info(8) $data_info($i)
        } elseif {[regexp -nocase {(Diaphragm3)} $data_names($i)] || [regexp -nocase {(Dia3)} $data_names($i)]} {
            set sql_data_info(9) $data_info($i)
        } elseif {[regexp -nocase {(Diaphragm)} $data_names($i)] || [regexp -nocase {(Diaphragm1)} $data_names($i)] || [regexp -nocase {(Dia1)} $data_names($i)]} {
            set sql_data_info(10) $data_info($i)
        } elseif {[regexp -nocase {(Shoc)} $data_names($i)]} {
            set sql_data_info(11) $data_info($i)
        } elseif {[regexp -nocase {(Acce)} $data_names($i)] || [regexp -nocase {(Atub)} $data_names($i)]} {
            set sql_data_info(12) $data_info($i)
        } elseif {[regexp -nocase {(Noz)} $data_names($i)] } {
            set sql_data_info(13) $data_info($i)
        } elseif {[regexp -nocase {(Pist)} $data_names($i)] } {
            set sql_data_info(14) $data_info($i)
        } elseif {[regexp -nocase {(1st D)} $data_names($i)] } {
            set sql_data_info(15) $data_info($i)
        } else {
            set maintenance_problems [open [file join $dir/maintenance_problems_X3.txt] a]
            puts $maintenance_problems "add new column to table: $data_names($i), data from file $files and run number $run_no"
            close $maintenance_problems
        }; #end if
        incr i
    }; #end while

    set max_no_X3_fields 16

    set sql_data_info($max_no_X3_fields) $notes
    
    set i 1
    while {$i <= $max_no_X3_fields} {; #set to number of data fields of X3
        puts -nonewline $txt_file_for_mysql "$sql_data_info($i)	"
        incr i
    }; #end while
    puts $txt_file_for_mysql "\n"    

};# end if

puts "Run = $run_no"
puts "File = $files"

if {$desc_f != "No"} {
    close $desc_f
}; #end if

if {$spec_f != "No"} {
    close $spec_f
}; #end if

#puts "End Porcessing shot"


