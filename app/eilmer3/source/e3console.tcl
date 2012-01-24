#!/usr/bin/wish
#
# \file e3console.tcl
# \ingroup eilmer3
# \brief GUI front-end for Eilmer3.
#
# \author P. A. Jacobs
#         School of Engineering
#         The University of Queensland
#
# \version 04-Oct-2005 : First hack adapted from plot/sourse/mb_view.tcl
# \version 20-Jun-2006 : brought over to mbcns2
# \version 03-Apr-2008 : ported to Eilmer3
# \version 15-Aug-2008 : cleaned up to leave a minimal set of options.
# \version 05-Oct-2009 : Updated for the *current* Eilmer3 code.
#
# -----------------------------------------------------------------------

proc showAboutBox {} {
    tk_messageBox -type ok -title "About e3console.tcl" \
	-message "e3console V0.4\nP.J. 05-Oct-2009"
}; # end proc showAboutBox

proc showGeneralHelp {} {
    tk_messageBox -type ok -title "General Help" -message "None available yet."
    # Should change this procedure to insert text straight into the text widget.
}; # end proc showGeneralHelp

# -----------------------------------------------------------------------

package require BWidget

# When looking at the OS name, we are expecting values such as
# "GNU/Linux" or "Cygwin"
if { [catch {exec uname -o} ::my_platform] } {
    puts "Seem to be not in a Unix-like environment."
    set ::my_platform "not_unix"
}
if {[string equal $tcl_platform(platform) "windows"] && \
	![string equal $::my_platform "Cygwin"] } {
    # we are on a plain Win32 system (without Cygwin)
    console show
}
set ::debugFlag 1
# Get the directory of the script and assume that the
# executable programs are in the same place.
if { [string length [info script]] } {
    set ::binDir [file dirname [info script]]
} else {
    set ::binDir [file dirname [info nameofexecutable]]
}
set ::workDir [pwd]
cd $::workDir
if {$::debugFlag == 1} {
    puts "bin directory         : $::binDir"
    puts "current work directory: $::workDir"
}
set ::jobName ""
set ::scriptExt ".py"
set ::subprocessId {}
set ::prepProgramName [file join $::binDir e3prep.py]
set ::runProgramName [file join $::binDir e3shared.exe]
set ::postProgramName [file join $::binDir e3post.py]
set ::viewProgramName "paraview"
set ::cygwinPython [file join / usr bin python2.4.exe]

#-----------------------------------------------------------------------------
# Default options that may be reset by the user.

set optionsPrep(svg) 0
set optionsPrep(zippedDataFiles) 1

set optionsMain(maxWallClock) ""
set optionsMain(zippedDataFiles) 1
set optionsMain(verbose) 0
set optionsMain(noComplain) 0
set optionsMain(extras) ""

set optionsPost(tindx) "all"
set optionsPost(outputFormat) "vtk-xml"
set optionsPost(zippedDataFiles) 1
set optionsPost(addMach) 0
set optionsPost(addPitotP) 0
set optionePost(addTotalP) 0
set optionsPost(extras) ""

set logFont {courier 12}

# -----------------------------------------------------------------------

proc startSubProcess { cmdString statusMessage } {
    if {$::debugFlag} { puts "cmdString is $cmdString" }
    appendToLogText "================ begin ================\n"
    set pipe [open "$cmdString"]
    fileevent $pipe readable [list Reader $pipe]
    set ::subprocessId [pid $pipe]
    updateStatusMessage $statusMessage
}; # end proc startSubProcess


proc Reader { pipe } {
    # Plumbing to enable this script to drive the external programs
    # which really do all of the work.
    if [eof $pipe] {
        catch {close $pipe}
	set ::subprocessId {}
	puts "End of pipe."
	flush stdout
	appendToLogText "================ end ==================\n"
	updateStatusMessage "Ready."
	return
    }
    gets $pipe line
    # puts $line; flush stdout
    appendToLogText "$line\n"; # reinstate the new-line character
}; # end proc Reader


proc killSubProcess {} {
    if { [string length $::subprocessId] > 0 } {
	appendToLogText "kill subprocess $::subprocessId\n"
	catch { exec kill $::subprocessId }
    }
}

# -----------------------------------------------------------------------

proc runPrep {} {
    # Construct a command for the Python interpreter.
    set inputFile [file join $::workDir $::jobName.py]
    if {[file exists $inputFile] == 0} {
	tk_messageBox -icon error -type ok \
	    -message "Could not find input script file: $inputFile"
	return
    }; # end if
    if { [string equal $::my_platform "Cygwin"] } {
	set cmd "|$::cygwinPython $::prepProgramName --job=$::jobName"
    } else {
	# Assume a unix environment.
	set cmd "|$::prepProgramName --job=$::jobName"
    }
    if { $::optionsPrep(svg) } {
	append cmd " --do-svg"
    }
    if { $::optionsPrep(zippedDataFiles) } {
	append cmd " --zip-files"
    } else {
	append cmd " --no-zip-files"
    }
    startSubProcess $cmd "Running e3prep, Please Wait..."
}; # end proc runPrep


proc helpPrep {} {
    # Construct a command for the Python interpreter.
    if { [string equal $::my_platform "Cygwin"] } {
	set cmd "|$::cygwinPython $::prepProgramName --help"
    } else {
	# Assume a unix environment.
	set cmd "|$::prepProgramName --help"
    }
    startSubProcess $cmd "Usage message for e3prep."
}; # end proc helpPrep


proc runMain {} {
    set cmd "|$::runProgramName --job=$::jobName --run"
    if { [string length [string trim $::optionsMain(maxWallClock)]] > 0 } {
	append cmd " --max-wall-clock=$::optionsMain(maxWallClock)"
    }
    if { $::optionsMain(zippedDataFiles) } {
	append cmd " --zip-files"
    } else {
	append cmd " --no-zip-files"
    }
    set extras [string trim $::optionsMain(extras)]
    if { [string length $extras] > 0 } {
	append cmd " "
	append cmd $extras
    }
    startSubProcess $cmd "Running e3shared, Please Wait..."
}; # end proc runMain


proc helpMain {} {
    set cmd "|$::runProgramName --help"
    startSubProcess $cmd "Usage message for e3shared."
}; # end proc helpMain


proc runPost {} {
    if { [string equal $::my_platform "Cygwin"] } {
	set cmd "|$::cygwinPython $::postProgramName --job=$::jobName"
    } else {
	# Assume a unix environment.
	set cmd "|$::postProgramName --job=$::jobName"
    }
    if { $::optionsPost(zippedDataFiles) } {
	append cmd " --zip-files"
    } else {
	append cmd " --no-zip-files"
    }
    if { $::optionsPost(addMach) } {
	append cmd " --add-mach"
    }
    if { $::optionsPost(addPitotP) } {
	append cmd " --add-pitot-p"
    }
    if { $::optionsPost(addTotalP) } {
	append cmd " --add-total-p"
    }
    set tindx [string trim $::optionsPost(tindx)]
    if { ([string length $tindx] > 0) } {
	# We are expecting either "all" or an integer-float pair.
	set indexValue [lindex [split $tindx] 0]
	append cmd [join [list " --tindx=" $indexValue] ""]
    }
    append cmd [join [list " --" $::optionsPost(outputFormat)] ""]
    set extras [string trim $::optionsPost(extras)]
    if { [string length $extras] > 0 } {
	append cmd " "
	append cmd $extras
    }
    startSubProcess $cmd "Running e3post, Please Wait..."
}; # end proc runPost


proc helpPost {} {
    if { [string equal $::my_platform "Cygwin"] } {
	set cmd "|$::cygwinPython $::postProgramName --help"
    } else {
	# Assume a unix environment.
	set cmd "|$::postProgramName --help"
    }
    startSubProcess $cmd "Usage message for e3post."
}; # end proc helpPost


proc runViewer {} {
    set cmd "|$::viewProgramName"
    startSubProcess $cmd "Running Viewer, close that GUI window when finished..."
}; # end proc runViewer

# --------------------------------------------------------------

proc loadScriptFile { useExtEditor } {
    set fileName [tk_getOpenFile -title "Load script File..." \
        -initialdir $::workDir ]
    puts "loadScriptFile: fileName=$fileName"
    if {[string length $fileName] > 0} {
        set shortFileName [file tail $fileName]
	set ::scriptExt [file extension $shortFileName]
	set ::jobName [file root $shortFileName]
        set ::workDir [file dirname $fileName]; # save for next time
	if { $useExtEditor == 0 } {
	    readScriptFile $fileName
	} else {
	    set cmd "|/usr/bin/emacs $fileName"
	    set pipe [open "$cmd"]
	    fileevent $pipe readable [list Reader $pipe]
	    set ::editSubprocessId [pid $pipe]
	}
    }
}; # end proc loadScriptFile


proc readScriptFile { {fileName {}} } {
    if {[string length $fileName] == 0} {
	# Try to assemble a sensible filename from the available pieces.
	set fileName $::jobName
	append fileName $::scriptExt
	set fileName [file join $::workDir $fileName]
	puts "Assembled script file name: $fileName"
    }
    set fp [open $fileName "r"]
    $::scriptTextWidget delete 1.0 end
    $::scriptTextWidget insert end [read $fp]
    close $fp
}; # end proc readScriptFile


proc saveScriptFile {} {
    set fileName [tk_getSaveFile -title "Save script File..." \
        -initialdir $::workDir ]
    puts "saveScriptFile: fileName=$fileName"
    if {[string length $fileName] > 0} {
        set shortFileName [file tail $fileName]
	set ::scriptExt [file extension $shortFileName]
	set ::jobName [file root $shortFileName]
        set ::workDir [file dirname $fileName]; # save for next time
	set fp [open $fileName "w"]
	puts $fp [$::scriptTextWidget get 1.0 end-1c]
	close $fp
    }
}; # end proc loadScriptPyFile

# --------------------------------------------------------------

proc postMainMenu {} {
   # Set up menu bar for the main window
   # m0 is the top-level menubar
   # m1 is a first-level menu
   # m2 is a second-level menu
   set m0 [menu .menubar]
   . config -menu $m0
   wm title . "e3console: Eilmer3 supervisor program."

   set m1 [menu $m0.file -tearoff 0]
   $m0 add cascade -label File -menu $m1
   $m1 add command -label "Open script file (external editor)" \
       -command { loadScriptFile 1 }
   $m1 add command -label "Open script file (internal editor)" \
       -command { loadScriptFile 0; .menubar.file entryconfigure 2 -state normal }
   $m1 add command -label "Save script file (internal editor)" \
       -command { saveScriptFile } -state disabled
   $m1 add command -label "Clear log messages" \
       -command { clearLogText }
   $m1 add separator
   $m1 add command -label Quit -command "exit"

   set m1 [menu $m0.action -tearoff 1 -title "Action"]
   $m0 add cascade -label Action -menu $m1
   $m1 add command -label "Prepare grid and initial flow field" \
       -command {runPrep}
   $m1 add command -label "Run simulation" \
       -command {runMain}
   $m1 add command -label "Extract flow field data at specified time" \
       -command {runPost}
   $m1 add command -label "Run data-viewing program" \
       -command {runViewer}
   $m1 add separator
   $m1 add command -label "Kill subprocess" -command { killSubProcess }

   $m0 add command -label "Options..." \
       -command {displayOptionsDialog}

   set m1 [menu $m0.help -tearoff 0]
   $m0 add cascade -label Help -menu $m1
   $m1 add command -label "About..." \
       -command "showAboutBox"
   $m1 add command -label "e3prep usage" \
       -command { helpPrep }
   $m1 add command -label "e3shared usage" \
       -command { helpMain }
   $m1 add command -label "e3post usage" \
       -command { helpPost }
   # $m1 add command -label "General Help" \
   #     -command "showGeneralHelp"
}; # end proc postMainMenu


proc initializeDisplay {} {
    postMainMenu

    #
    # A line indicating the current jobName and extension.
    #
    set jobNameFrame [LabelFrame .fn -text "Job Name:" -side left -borderwidth 3 -relief flat]
    set jobNameEntry [entry [$jobNameFrame getframe].jne -textvariable ::jobName -width 50]
    bind $jobNameEntry <Return> { readScriptFile }
    pack $jobNameEntry -side left -expand 1 -fill x
    set extLabel [label [$jobNameFrame getframe].el -text "Script extension:"]
    pack $extLabel -side left
    set extEntry [entry [$jobNameFrame getframe].ee -textvariable ::scriptExt -width 20 ]
    pack $extEntry -side left -expand 1 -fill x
    pack $jobNameFrame -side top -expand 1 -fill both
    #
    # A scrolling text window for the content of the user's script.
    #
    set scriptFrame [LabelFrame .tf0 -text "Input script:" -side top]
    set ::scriptTextWidget [text [$scriptFrame getframe].t -height 20 -width 70 \
				-font $::logFont -wrap char \
		       -yscrollcommand [list [$scriptFrame getframe].vsb set] ]
    set scriptScrollBar [scrollbar [$scriptFrame getframe].vsb -orient vertical \
			     -command {$::scriptTextWidget yview} ]
    pack $::scriptTextWidget -side left -expand 1 -fill both
    pack $scriptScrollBar -side left -fill y
    pack $scriptFrame -side top -expand 1 -fill both
    #
    # A scrolling text window for the text output from running the programs.
    #
    set textFrame [LabelFrame .tf -text "Log of output:" -side top]
    set ::logTextWidget [text [$textFrame getframe].t -height 10 -width 70 \
			     -font $::logFont -wrap char \
		       -yscrollcommand [list [$textFrame getframe].vsb set] ]
    set textScrollBar [scrollbar [$textFrame getframe].vsb -orient vertical \
			   -command {$::logTextWidget yview} ]
    pack $::logTextWidget -side left -expand 1 -fill both
    pack $textScrollBar -side left -fill y
    pack $textFrame -side top -expand 1 -fill both
    #
    # A status line giving some hint as to the current state.
    #
    set statusFrame [LabelFrame .fs -text "Status:" -side left -borderwidth 3 -relief flat]
    set statusEntry [entry [$statusFrame getframe].state -textvariable ::statusText ]
    pack $statusEntry -side right -expand 1 -fill x
    pack $statusFrame -side bottom -expand 1 -fill both
}; # end proc initializeDisplay

proc clearLogText {} {
    $::logTextWidget delete 1.0 end
    update idletasks
}

proc appendToLogText { {newText ""} } {
    $::logTextWidget insert end $newText
    $::logTextWidget see end
    update idletasks
}

proc updateStatusMessage { {newText ""} } {
    set ::statusText $newText
    update idletasks
}

#-----------------------------------------------------------------------------

proc checkTimesFile {} {
    set fileName $::jobName
    append fileName ".times"
    set fp [open $fileName "r"]
    set tindxList {}
    foreach line [split [read -nonewline $fp] \n] {
	set tokens [split $line]
	if { [lindex $tokens 0] != "#" } {
	    lappend tindxList [join [list [lindex $tokens 0] [lindex $tokens 1]] " "]
	}
    }
    close $fp
    lappend tindxList "all"
    $::tindxListBox configure -values $tindxList
    return
}

proc displayOptionsDialog {} {
    if { [winfo exists .optionsDialog] } {
        # Dialog widget already exists; just show it
	# after 0.5 seconds to avoid other events bringing
	# the main window to the fore.
        after 500 {raise .optionsDialog}
    } else {
        # Create the dialog widget and show it.
        toplevel .optionsDialog
        set nb [NoteBook .optionsDialog.notebook -height 200 -width 400]
        set bF [frame .optionsDialog.buttonFrame]
	#--------------------------------------------------
        set page [$nb insert end 1 -text "Prepare grids"]
	set lab0 [label $page.lab0 -text "Use e3prep.py to generate grid and initial flow-field data."]
	set frame1 [frame $page.f1 -relief groove -borderwidth 2]
	set lab1 [label $frame1.lab1 -text "Standard options:"]
        set cb2 [checkbutton $frame1.cb2 -text "Generate SVG Output (relevant for 2D only)" \
            -variable ::optionsPrep(svg)]
        set cb3 [checkbutton $frame1.cb3 \
            -text "Write compressed (zipped) grid and solution files" \
            -variable ::optionsPrep(zippedDataFiles)]
	pack $lab1 $cb2 $cb3 -anchor w
	pack $lab0 -anchor w -pady 10
	pack $frame1 -anchor w -fill x
	#--------------------------------------------------
        set page [$nb insert end 2 -text "Run simulation"]
	set lab0 [label $page.lab0 -text "Use e3shared.exe to generate subsequent flow-field data."]
        set cb0 [checkbutton $page.cb0 \
            -text "Use compressed (zipped) grid and solution files" \
            -variable ::optionsMain(zippedDataFiles)]
	set fWallClock [frame $page.fWC]
	set labWallClock [label $fWallClock.labwc -text "Maximum Wall-Clock (sec):"]
	set entryWallClock [entry $fWallClock.entrywc -width 20 \
			     -textvariable ::optionsMain(maxWallClock)]
	pack $labWallClock -side left
	pack $entryWallClock -side left
	set labExtras [label $page.labe -text "Extra options:"]
	set entryExtras [entry $page.entrye -width 60 \
			     -textvariable ::optionsMain(extras)]
	pack $lab0 -anchor w -pady 10
        pack $cb0 $fWallClock $labExtras $entryExtras -anchor w
	#--------------------------------------------------
        set page [$nb insert end 3 -text "Extract field data"]
	set lab0 [label $page.lab0 -text "Use e3post.py to extract flow-field data."]
	set frame1 [frame $page.f1]
	set labTime [label $frame1.labt -text "TimeIndex (integer value or 'all'):"]
	set ::tindxListBox [ComboBox $frame1.entryt -width 20 -entrybg white \
				-textvariable ::optionsPost(tindx) \
				-values [list "all" "9999 0.0"] ]
	set checkTime [button $frame1.butt -text "Check times file" \
			   -command checkTimesFile]
	pack $labTime -side left
	pack $::tindxListBox -side left
	pack $checkTime -side left
        set cb0 [checkbutton $page.cb0 \
            -text "Read compressed (zipped) grid and solution files" \
            -variable ::optionsPost(zippedDataFiles)]

	set buttonFrame [frame $page.bf]
	set formatFrame [LabelFrame $buttonFrame.ff -text "Output Format:" \
			     -side top -relief groove -borderwidth 2]
	set rb2 [radiobutton [$formatFrame getframe].rb2 \
		     -text "TECPLOT (block-structured mesh)" \
		     -variable ::optionsPost(outputFormat) \
		     -value "tecplot"]
	set rb3 [radiobutton [$formatFrame getframe].rb3 \
		     -text "VTK-XML (unstructured mesh)" \
		     -variable ::optionsPost(outputFormat) \
		     -value "vtk-xml"]
	set rb4 [radiobutton [$formatFrame getframe].rb4 \
		     -text "Plot3D (block-structured mesh)" \
		     -variable ::optionsPost(outputFormat) \
		     -value "plot3d"]
	grid $rb2 -row 0 -column 0 -sticky w
	grid $rb3 -row 1 -column 0 -sticky w
	grid $rb4 -row 2 -column 0 -sticky w

	set deFrame [LabelFrame $buttonFrame.lf -text "Add data elements:" \
			 -side top -relief groove -borderwidth 2]
	set cb1 [checkbutton [$deFrame getframe].cb1 \
		     -text "Mach number" \
		     -variable ::optionsPost(addMach)]
	set cb2 [checkbutton [$deFrame getframe].cb2 \
		     -text "Pitot pressure" \
		     -variable ::optionsPost(addPitotP)]
	set cb3 [checkbutton [$deFrame getframe].cb3 \
		     -text "Total pressure" \
		     -variable ::optionsPost(addTotalP)]
	grid $cb1 -row 0 -column 0 -sticky w
	grid $cb2 -row 1 -column 0 -sticky w
	grid $cb3 -row 2 -column 0 -sticky w

	pack $formatFrame -side left -anchor w
	pack $deFrame -side right -anchor e -padx 5

	set labExtras [label $page.labe -text "Extra options:"]
	set entryExtras [entry $page.entrye -width 60 \
			     -textvariable ::optionsPost(extras)]

	pack $lab0 -anchor w -pady 10 
        pack $frame1 $cb0 $buttonFrame $labExtras $entryExtras -anchor w
	#--------------------------------------------------

        NoteBook::compute_size $nb
        $nb raise 1

        set bDismiss [button $bF.dismiss -text "Lower Window" \
            -command {lower .optionsDialog .} ]
        pack $bDismiss -side left -padx 4 -pady 4

        pack $nb $bF
    }; # end if
    wm title .optionsDialog "e3console options for programs."
}; # end proc

# --------------------------------------------------------

# Everything is ready; paint the main display...
initializeDisplay
displayOptionsDialog
lower .optionsDialog .
updateStatusMessage "Ready"

# --------------- end of e3console.tcl ---------------------
