#!/bin/sh
# comment to tcl \
exec tclsh "$0" ${1+"$@"}

proc lowercase {dir} {
    cd
    cd $dir
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

set dir [lindex $argv 0]
lowercase $dir
