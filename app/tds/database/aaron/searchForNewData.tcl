#!/bin/sh
# comment to tcl \
exec tclsh "$0" ${1+"$@"}

# These procedures will find the shots the aren't already in the database
# and then call another script to add them into the database.

proc checkToSeeIfRunNoContainsLetters {run_no} {

# this procedure checks to see if the run_no contains letters
# if it does it will remove the letters from the left and right
# of the run_no.  if the run_no contains additional numbers to 
# the right of the actual run_no this will fail
# eg the following will fail if the actual run_no is 674: 
# Run2B674Execute 

set n [string index $run_no 0]

if {[regexp {^[0-9]+$} $n] == 0} {; # if the first char is a letter
                                    # perform the while loop
    while {[regexp {^[0-9]+$} $n] == 0} {; # while the first char 
        # is a letter trim the string until it becomes a number
        set run_no [string range $run_no 1 end]
        set n [string index $run_no 0]; # n will be the first
                                        # character

    }; #end while
}; #end if

set m [string index $run_no end]
if {[regexp {^[0-9]+$} $m] == 0} {;# if the last char is a letter
                                   # perform the commands and the 
                                   # while loop
    set temp $run_no
    set n [string index $run_no 0]
    set i 0
    while {[regexp {^[0-9]+$} $n] == 1} {;# while the first char 
        # is a number trim the string until it becomes a letter
        # and store the remaining string in temp
        set temp [string range $temp 1 end]
        set n [string index $temp 0]
        incr i
    }; #end while
    set n [expr [string first $temp $run_no] - 1]
    # returns the part of run_no up until the first
    # occurence of temp
    set run_no [string range $run_no 0 $n]
}; #end if
return $run_no
}; # end proc

proc x2RunAndModNames {files} {
   # this procedure is used as X2 has a diffent naming convention
    # to the other tubes at this point in time.  Eg if the shot
    # directory was s674, the mod file would be 674.mod and the
    # text file would be Run674.txt
if {[regexp {^[sS]} $files]} {
            set names_mod [string range $files 1 end]
            set names_run RUN$names_mod
            set temp "$names_mod $names_run"
            return $temp
}; #end if
}; #end proc


proc checkToSeeIfFileExistsAndOpenSpec {dir run_no spec_f files} {
    set no_file_found 0 ; 
    set result [x2RunAndModNames $files]
    set file_run [lindex $result 1]
    if ![file exists $dir/rundesc/$files.txt ] {
        if ![file exists $dir/rundesc/$file_run.txt ] {
            if ![file exists $dir/rundesc/$run_no.txt ] {
                set spec_f "No"
                if {[string length $run_no] > 2} { 
                    # so it doesn't find all files if only                                          
                    # has two characters in run_no
                    cd $dir/rundesc/
	
                    foreach files [glob -nocomplain *$run_no*] {
	                set t [expr [string first $run_no $files] - 1]
                        # finds the position where the run_no starts
                        # in the filename.  These lines are required
                        # so if the run number is say 674, files 
                        # doesn't return 674, 1674, 2674 etc... 

	                set start_of_file [string range $files 0 $t]; 
                        # obtains any characters that preceed the run_no
	                if {[regexp {^[0-9]+$} $start_of_file] == 0 && \
	                 [regexp {^[0-9]+$} [string index $start_of_file $t]] == 0} {
	                    # this checks that that preceeding chars
                            # contain letters and last char is not 
                            # a number
			     
	                    if ![file exists $dir/rundesc/$files ] {
		                set run_desc_found "No"
                                set spec_f "No"
		            } else {
         	                set run_desc_found "Yes"
		                set spec_f [open $dir/rundesc/$files r] 
                                incr no_file_found
		            }; #end if
                            cd /home/moncdata/databaseScripts
	                }; #end if
	                incr no_file_found
                    };# end foreach
                    if [file exists $dir/$run_no/$run_no.txt ] {
                        set spec_f [open $dir/$run_no/$run_no.txt r]
                        set run_desc_found "Yes"
                        incr no_file_found
                    }; #end if
                } else {
                    set spec_f "No"
                }; #end if
            } else {
                set spec_f [open $dir/rundesc/$run_no.txt r]
                set run_desc_found "Yes"
	        incr no_file_found
            }; # end if
            if {$no_file_found == 0} { ;# added
                set run_desc_found "No"    
	        set spec_f "No";
            }; #end if
	} else {
            set spec_f [open $dir/rundesc/$file_run.txt]
            incr no_file_found
            set run_desc_found "Yes"
        }; #end if
    } else {
       set spec_f [open $dir/rundesc/$files.txt]
       incr no_file_found
       set run_desc_found "Yes"
    }; #end if
    return $spec_f	
}; # end proc

proc checkToSeeIfFileExistsAndOpenDesc {dir run_no desc_f files} {
    set no_file_found 0
    set result [x2RunAndModNames $files]
    set file_mod [lindex $result 0]
    if ![file exists $dir/descript/$files.mod] {
        if ![file exists $dir/descript/$file_mod.mod ] {
            if ![file exists $dir/descript/$run_no.mod ] {
                set desc_f "No"
                if {[string length $run_no] > 2} { 
                    # so it doesn't find all files if only
                    # has two characters in run_no
                    cd $dir/descript/
                    foreach files [glob -nocomplain *$run_no*] {
	                set t [expr [string first $run_no $files] - 1]
                        # finds the position where the run_no starts
                        # in the filename. These lines are required 
                        # so if the run number is say 674, files
                        # doesn't return 674, 1674, 2674 etc ... 

	                set start_of_file [string range $files 0 $t]; 
                        # obtains any characters 
	                # that preceed the run_no
	                if {[regexp {^[0-9]+$} $start_of_file] == 0 && \
	                 [regexp {^[0-9]+$} [string index $start_of_file $t]] == 0} {				     
	                    # this checks that the preceeding chars
                            # contain letters and the last char
                            # is not a number		     
                            if ![file exists $dir/descript/$files ] {
		                set descript_found "No"
                                set desc_f "No"
		            } else {
                                set desc_f [open $dir/descript/$files r]
		                set descript_found "Yes"
                                cd /home/moncdata/databaseScripts
                                incr no_file_found
		            }; #end if
	                }; #end if
	                incr no_file_found
	            }; # end foreach
                    if {[file exists $dir/$run_no/$run_no.mod] } {
                        set desc_f [open $dir/$run_no/$run_no.mod]
                        set descript_found "Yes"
                        incr no_file_found
                    };#end if
                }; #end if     
            } else {
                set desc_f [open $dir/descript/$run_no.mod r]
	        set descript_found "Yes"
	        incr no_file_found
            }; # end if
            if {$no_file_found == 0} { ;# added
                set descript_found "No"
	        set desc_f "No"
	     }; #end if
	} else {
            set desc_f [open $dir/descript/$file_mod.mod]
            set descript_found "Yes"
            incr no_file_found
        }; #end if   
    } else {
        set desc_f [open $dir/descript/$files.mod]
        set descript_found "Yes"
        incr no_file_found
    }; #end if
    return $desc_f	
}; # end proc

proc search_sql_data {run_no} {
    set check_run_no_data [open /home/moncdata/databaseScripts/currentData.txt r]
    set check 0
    foreach line [split [read $check_run_no_data] \n] {
        if {[string first " $run_no " $line] > -1} {
	    set check 1
        }; #end if
    }; #end foreach   	   
    if {$check == 0} {
        set need_to_be_added yes
    } else {
	set need_to_be_added no
    }; #end if
    close $check_run_no_data
    return $need_to_be_added
}; #end proc

proc lowercase {file_dir} {
    cd
    cd $file_dir
    foreach files [glob *] {
        set file_end [file extension $files]
        if {$file_end == ".MOD"} {
            set files [file root $files]
	    file rename -force $files.MOD $files.mod
        }; #end if
        if {$file_end == ".TXT"} {
            set files [file root $files]
	    file rename -force $files.TXT $files.txt
        }; #end if

    }; #end foreach
}



#-----------------------------------start processing----------------------------

if {$argc == 0} {
    puts ""
    puts ""
    puts "Usage: ./searchForNewData.tcl searchSize tubeName dir ?shotName?"
    puts ""
    puts " searchSize can be set to all (add in all new shots) or set to one" 
    puts " (to add a specific shot, if this option is set, the shotName must" 
    puts " also be given) "
    puts ""
    puts " tubeName is the name of the shock tunnel / expansion tube. IE T4, X3, X2"
    puts ""
    puts " dir is the directory where the shot information is stored"
    puts ""
    puts " shotName only needs to be supplied when searchSize is set to one.  shotName is" 
    puts " just the name of the shot as it appears in its specific directory"
    puts ""
    puts "Example 1: Do all new shots for X3"
    puts "./searchForNewData.tcl all X3 /home/moncdata/X3/"
    puts "Example 2: Do one specific shot (in this case 7319) for T4"
    puts "./searchForNewData.tcl one T4 7319 /home/moncdata/T4/"
    puts ""
    puts ""
}



if {$argc > 0 && [lindex $argv 0] == "all" || [lindex $argv 0] == "one"} {
    set dir [lindex $argv 2]
    cd $dir
    lowercase $dir/descript
    lowercase $dir/rundesc
    cd
    cd $dir
    # If another shock tunnel or expansion tube is added,
    # another block of code needs to be added here to assign
    # its data text file and maintenance problems text file
    if {[lindex $argv 1] == "T4"} {
        set txt_file_for_mysql [open [file join $dir/test_data_T4.txt] w]
        # the text file must be the same name as the table name
        close $txt_file_for_mysql
        set maintenance_problems [open [file join $dir/maintenance_problems_T4.txt] w]
        close $maintenance_problems
    };#end if

    if {[lindex $argv 1] == "X2"} {
        set txt_file_for_mysql [open [file join $dir/test_data_X2.txt] w]
        # the text file must be the same name as the table name
        close $txt_file_for_mysql
        set maintenance_problems [open [file join $dir/maintenance_problems_X2.txt] w]
        close $maintenance_problems
    };#end if

    if {[lindex $argv 1] == "X3"} {
        set txt_file_for_mysql [open [file join $dir/test_data_X3.txt] w]
        # the text file must be the same name as the table name
        close $txt_file_for_mysql
        set maintenance_problems [open [file join $dir/maintenance_problems_X3.txt] w]
        close $maintenance_problems
    };#end if

    if {[lindex $argv 0] == "all"} {
        set files_to_be_searched [glob -nocomplain *]
    }; #end if

    if {[lindex $argv 0] == "one"} {
        set files_to_be_searched [glob -nocomplain *[lindex $argv 3]*]
    }; #end if

    foreach files $files_to_be_searched {
        # If another shock tunnel or expansion tube is added
        # another bit of code needs to be added to open up
        # its data text file
        if {[lindex $argv 1] == "T4"} {
            set txt_file_for_mysql [open [file join $dir/test_data_T4.txt] a] 
        }; #end if
        if {[lindex $argv 1] == "X2"} {
            set txt_file_for_mysql [open [file join $dir/test_data_X2.txt] a]
        }; #end if
        if {[lindex $argv 1] == "X3"} {
            set txt_file_for_mysql [open [file join $dir/test_data_X3.txt] a]
        }; #end if
        #puts "open sql text 2 $txt_file_for_mysql"
        cd $dir
        if {[file isdirectory $files] == 1 && [regexp {^[A-Za-z_]+$} $files] == 0 } {
            set run_no_files [checkToSeeIfRunNoContainsLetters $files]
	    set does_file_need_to_be_added [search_sql_data $run_no_files]
	    if {$does_file_need_to_be_added == "yes"} {
                set arg1 $dir/$files/
                set run_no [file tail $arg1]
	        set run_no [checkToSeeIfRunNoContainsLetters $run_no]
	        set desc_f $dir/descript/$run_no.mod
	        set desc_f [checkToSeeIfFileExistsAndOpenDesc $dir $run_no $desc_f $files]
	        set spec_f $dir/rundesc/$run_no.txt
	        set spec_f [checkToSeeIfFileExistsAndOpenSpec $dir $run_no $spec_f $files]
                source /home/moncdata/databaseScripts/addInfoToDatabase.tcl
                cd $dir/$files/
                set ascii_files [glob -nocomplain $files{A}.???]
                set has_ascii_files_check [llength $ascii_files]
	        if {$has_ascii_files_check > 0} {
	            #puts "run no $run_no already has ascii files"
	        } else {
	            puts "going to run defrock"
                    cd /home/moncdata/cfd/tds/defrock/
                    catch {exec ./defrock $dir $files} defrock_output
                    if {[regexp "End of defrock run" $defrock_output]} {
                        exec ./defrock $dir/ $files
                    }; #end if
	        }; #end if
                cd
                cd $dir
	    } else {
	     
            }; #end if
        }; #end if 
        close $txt_file_for_mysql
    }; #end foreach

    cd
    cd test
    # If another shock tunnel or expansion tube is added
    # another block needs to be added here
    if {[lindex $argv 1] == "T4"} {
        set maintenance_problems [open [file join $dir/maintenance_problems_T4.txt] r]
        set i 0
        foreach line [split [read $maintenance_problems] \n] {
            incr i
        }; #end if
        if {$i > 0} {
            puts "There appears to be some new columns that aren't in the database"
            puts "Please refer to maintenance_problems_T4.txt for details"
        }; #end if
        close $maintenance_problems
    }; #end if
    if {[lindex $argv 1] == "X2"} {
        set maintenance_problems [open [file join $dir/maintenance_problems_X2.txt] r]
        set i 0
        foreach line [split [read $maintenance_problems] \n] {
            incr i
        }; #end if
        if {$i > 0} {
            puts "There appears to be some new columns that aren't in the database"
            puts "Please refer to maintenance_problems_X2.txt for details"
        }; #end if
        close $maintenance_problems
    }; #end if
    if {[lindex $argv 1] == "X3"} {
        set maintenance_problems [open [file join $dir/maintenance_problems_X3.txt] r]
        set i 0
        foreach line [split [read $maintenance_problems] \n] {
            incr i
        }; #end if
        if {$i > 0} {
            puts "There appears to be some new columns that aren't in the database"
            puts "Please refer to maintenance_problems_X3.txt for details"
        }; #end if
        close $maintenance_problems
    }; #end if
    puts end

}; #end if argc > 0













