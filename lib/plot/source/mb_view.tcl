#!/usr/bin/wish
# \file mb_view.tcl
# \brief GUI front-end for mb_cont contouring program.
#
# \author P. A. Jacobs
#         Department of Mechanical Engineering
#         The University of Queensland
#
# Version...
# 29-Sep-98 : First hack
# 01-Oct-98 : Added canvas display, most bits functional
# 03-Oct-98 : Generalized variable menu
# 19-Jan-99 : small fixes
# 20-Apr-99 : extra couple of options for suppressing the colour table
#             and plot axes
# 17-Jun-01 : tidy up a little so that things work on Win32
#
# ----------------------------------------------------

# Set default values for the data stored in array mbv.
set mbv(debug) 1

# When looking at the OS name, we are expecting values such as
# "GNU/Linux" or "Cygwin"
if { [catch {exec uname -o} ::my_platform] } {
    puts "Seem to be not in a Unix-like environment."
    set ::my_platform "not_unix"
}
if { [string equal $tcl_platform(platform) "windows"] && 
     ![string equal $::my_platform "Cygwin"] } {
   # we are on a plain Win32 system
   console show
}; # end if

# Get the directory of the script and assume that the
# executable file mb_cont.exe is in the same place.
if { [string length [info script]] } {
    set mbv(binDir) [file dirname [info script]]
} else {
    set mbv(binDir) [file dirname [info nameofexecutable]]
}
set mbv(programName) [file join $mbv(binDir) mb_cont.exe]

set mbv(workDir) [pwd]
cd $mbv(workDir)

if {$mbv(debug) == 1} {
   puts "bin directory         : $mbv(binDir)"
   puts "current work directory: $mbv(workDir)"
   puts "plotting-program name : $mbv(programName)"
}; # end if

set mbv(inputFile) ""
set mbv(plotFile) ""
set mbv(varIndex) 3
set mbv(plotType) gif
set mbv(pixels_per_mm) 3.0
set mbv(gifFile) "mb_view.gif"
set mbv(psFile) "mb_view.ps"
set mbv(inputFile) ""
set mbv(colour) 1
set mbv(fill) 0
set mbv(mirror) 0
set mbv(trueShape) 1
set mbv(mesh) 0
set mbv(edge) 0
set mbv(label) 1
set mbv(table) 1
set mbv(axes) 1
set mbv(xrange) ""
set mbv(yrange) ""
set mbv(lrange) ""

# ------------------------------------------------

proc displayWaitPlacard { msgText } {
   toplevel .waitPlacard -borderwidth 4 -relief raised
   wm overrideredirect .waitPlacard 1
   after idle {
      update idletasks
      set xmax [winfo screenwidth .waitPlacard]
      set ymax [winfo screenheight .waitPlacard]
      set x0 [expr ($xmax - [winfo reqwidth .waitPlacard])/2]
      set y0 [expr ($ymax - [winfo reqheight .waitPlacard])/2]
      wm geometry .waitPlacard "+$x0+$y0"
   }
   label .waitPlacard.info -text "$msgText"
   pack .waitPlacard.info -padx 15 -pady 15
   update
}; # end proc displayWaitPlacard

proc removeWaitPlacard {} {
   update
   destroy .waitPlacard
}; # end proc

# Plumbing to enable this script to drive mb_cont
# which really does all of the work.

proc Reader { pipe } {
   if [eof $pipe] {
      catch {close $pipe}
      puts "End of pipe."
      flush stdout
      set ::EndOfPipeFlag 1
      removeWaitPlacard
      return
   }
   gets $pipe line
   puts $line
   flush stdout
}; # end proc Reader

proc runProgram {} {
   # Set up the command line that runs mb_cont
   # as a separate process.

   global mbv

   if {[file exists $mbv(inputFile)] == 0} {
      tk_messageBox -icon error -type ok \
         -message "Could not find data file: $mbv(inputFile)."
      return
   }; # end if

   if {[string compare $mbv(plotType) "gif"] == 0} { 
      set outFile $mbv(gifFile)
      set plotType "-gif"
      set pixelsOption "-pixelspermm $mbv(pixels_per_mm)"
   } else {
      set outFile $mbv(psFile)
      set plotType "-ps"
      set pixelsOption ""
   }; # end if

   if {$mbv(colour)} {
      set colourOption "-colour"
   } else {
      set colourOption ""
   }; # end if

   if {$mbv(fill)} {
      set fillOption -fill
   } else {
      set fillOption ""
   }; # end if

   if {$mbv(mirror)} {
      set mirrorOption -mirror
   } else {
      set mirrorOption ""
   }; # end if

   if {$mbv(mesh)} {
      set meshOption -mesh
   } else {
      set meshOption ""
   }; # end if

   if {$mbv(edge)} {
      set edgeOption -edge
   } else {
      set edgeOption ""
   }; # end if

   if {$mbv(trueShape)} {
      set shapeOption ""
   } else {
      set shapeOption "-notrueshape"
   }; # end if

   if {$mbv(label)} {
      set labelOption ""
   } else {
      set labelOption "-nolabel"
   }; # end if

   if {$mbv(table)} {
      set tableOption ""
   } else {
      set tableOption "-notable"
   }; # end if

   if {$mbv(axes)} {
      set axesOption ""
   } else {
      set axesOption "-noaxes"
   }; # end if

   if {[string compare $mbv(xrange) ""] == 0} {
      set xrange ""
   } else {
      set xrange "-xrange $mbv(xrange)"
   }; # end if 

   if {[string compare $mbv(yrange) ""] == 0} {
      set yrange ""
   } else {
      set yrange "-yrange $mbv(yrange)"
   }; # end if 

   if {[string compare $mbv(lrange) ""] == 0} {
      set levels ""
   } else {
      set levels "-levels $mbv(lrange)"
   }; # end if 

   set cmd "|$mbv(programName) \
      -fi $mbv(inputFile) -fo $outFile -var $mbv(varIndex) \
      $plotType $pixelsOption $colourOption $fillOption $mirrorOption \
      $meshOption $edgeOption $labelOption $shapeOption \
      $tableOption $axesOption \
      $xrange $yrange $levels"
   if {$mbv(debug)} { puts "cmd is $cmd" }

   set pipe [open "$cmd"]
   fileevent $pipe readable [list Reader $pipe]
   displayWaitPlacard "Generating Plot, Please Wait..."
   vwait ::EndOfPipeFlag
   displayGIFImage .
}; # end proc runProgram

# --------------------------------------------

proc initializeDisplay { parentFrame } {
   # The main window contains a display (canvas) for the
   # GIF image and a couple of scroll bars.

   global mbv

   set display $parentFrame; append display display
   set xsbar   $parentFrame; append xsbar xsbar
   set ysbar   $parentFrame; append ysbar ysbar

   canvas $display -width 15c -height 10c \
      -xscrollcommand "$xsbar set" \
      -yscrollcommand "$ysbar set"
   scrollbar $xsbar -orient horizontal \
      -command "$display xview"
   scrollbar $ysbar -orient vertical \
      -command "$display yview"
   grid $display -row 0 -column 0 -sticky nsew
   grid $ysbar -row 0 -column 1 -sticky ns
   grid $xsbar -row 1 -column 0 -sticky ew
   grid columnconfigure $parentFrame 0 -weight 1
   grid rowconfigure $parentFrame 0 -weight 1

   if [file exists $mbv(gifFile)] {
      set im [image create photo -file $mbv(gifFile)]
      set mbv(imageId) [$display create image 0 0 \
         -image $im -anchor nw]
      $display configure -scrollregion [$display bbox all]
      $display yview moveto 0.5
   }; # end if

}; # end proc initializeDisplay


proc displayGIFImage { parentFrame } {
   # Put the GIF image onto the display canvas.
   # If one is already displayed, delete the old and
   # put up the new.

   global mbv
   set display $parentFrame; append display display
   set xsbar   $parentFrame; append xsbar xsbar
   set ysbar   $parentFrame; append ysbar ysbar

   if [info exist mbv(imageId)] {
      # Wipe the old image before displaying the new.
      puts "Deleting old image from canvas."
      $display delete $mbv(imageId)
   }; # end if

   puts "Looking for plotFile $mbv(plotFile)."
   if {[file exists $mbv(gifFile)] == 1} {
      puts "Displaying image from file: $mbv(gifFile)"
      set im [image create photo -file $mbv(gifFile)]
      set mbv(imageId) [$display create image 0 0 \
         -image $im -anchor nw]
      $display configure -scrollregion [$display bbox all]
      $display yview moveto 0.5
   }; # end if

}; # end proc displayGIFImage

# --------------------------------------------------------

proc setInputFile {} {
   # Set the input file name and scan the beginning
   # of the file to get the variable names.

   global mbv

   set inputFileName \
      [tk_getOpenFile -title "Open a Generic Data File"]
   set mbv(inputFile) [file tail $inputFileName]
   set mbv(workDir) [file dirname $inputFileName]
   eval "cd $mbv(workDir)"
   puts "Selected data file: $mbv(inputFile)"

   if {[file exists $mbv(inputFile)] != 1} {
      tk_messageBox -title "Problem with Input File" \
         -icon warning -message "Specified Input File does not exist"
      return
   }; # end if

   set fp [open $mbv(inputFile) "r"]
   set buffer [gets $fp];  # first line is a title
   set buffer [gets $fp];  # second line is number of variables
   scan $buffer "%d" mbv(nvar)
   for {set i 0} {$i < $mbv(nvar)} {incr i} {
      set buffer [gets $fp]
      scan $buffer "%s" mbv(variable$i)
      puts "Variable $i $mbv(variable$i)"
   }; # end for

   postVariableMenu
}; # end proc setInputFile


proc setGIFFile {} {
   global mbv
   set fileName \
      [tk_getSaveFile -initialfile "mb_view.gif"]
   set mbv(gifFile) [file tail $fileName]
   puts "selected file for GIF output: $mbv(gifFile)"
}; # end proc setGIFFile

proc setPSFile {} {
   global mbv
   set fileName \
      [tk_getSaveFile -initialfile "mb_view.ps"]
   set mbv(psFile) [file tail $fileName]
   puts "selected file for PS output: $mbv(psFile)"
}; # end proc setPSFile

# --------------------------------------------------------

proc postMainMenu {} {
   # Set up menu bar for the main window
   # m0 is the top-level menubar
   # m1 is a first-level menu
   # m2 is a second-level menu

   global mbv

   set m0 [menu .menubar]
   . config -menu $m0
   wm title . "mb_view: CFD postprocessing on the cheap."

   set m1 [menu $m0.file -tearoff 0]
   $m0 add cascade -label File -menu $m1
   $m1 add command -label "Open Data File" \
      -command { setInputFile }
   $m1 add separator
   $m1 add command -label "Set PS Output File" \
      -command { setPSFile }
   $m1 add command -label "Set GIF Output File" \
      -command { setGIFFile }
   $m1 add separator
   $m1 add command -label Quit -command "exit"

   set m1 [menu $m0.plot -tearoff 1 -title "Plot Menu"]
   $m0 add cascade -label Plot -menu $m1
   $m1 add command -label "Refresh_Display" \
      -command {displayGIFImage .}
   $m1 add separator
   $m1 add command -label "Generate GIF Image" \
      -command {set mbv(plotType) gif; runProgram}
   $m1 add command -label "Generate PS File" \
      -command {set mbv(plotType) ps; runProgram}
   $m1 add command -label "Display Range Dialog" \
      -command {rangeDialog}

   set m1 [menu $m0.options -tearoff 0]
   $m0 add cascade -label Options -menu $m1
   $m1 add check -label "Colour" -variable mbv(colour)
   $m1 add check -label "Mirror Image" -variable mbv(mirror)
   $m1 add check -label "Fill" -variable mbv(fill)
   $m1 add check -label "Add Mesh" -variable mbv(mesh)
   $m1 add check -label "Add Mesh Edges" -variable mbv(edge)
   $m1 add check -label "True Shape" -variable mbv(trueShape)
   $m1 add check -label "Include Text" -variable mbv(label)
   $m1 add check -label "Include Table" -variable mbv(table)
   $m1 add check -label "Include Axes" -variable mbv(axes)

   set m1 [menu $m0.help -tearoff 0]
   $m0 add cascade -label Help -menu $m1
   $m1 add command -label "About..." \
      -command "showAboutBox"
   $m1 add command -label "General Help" \
      -command "showGeneralHelp"
}; # end proc postMainMenu


proc postVariableMenu {} {
   # Post a new menu for the variable names

   global mbv

   if [info exists mbv(variableMenuPosted)] {
      .menubar delete Variable
      destroy .menubar.variable
      unset mbv(variableMenuPosted)
   }; # end if

   set m1 [menu .menubar.variable -tearoff 0]
   .menubar insert Options cascade -label Variable -menu $m1
   for {set i 0} {$i < $mbv(nvar)} {incr i} {
      $m1 add radio -label $mbv(variable$i) \
         -variable mbv(selectedVariable) \
         -command "set mbv(varIndex) $i"
   }; # end for
   set mbv(variableMenuPosted) 1

   if {$mbv(nvar) < 3} {
      tk_messageBox -type ok -title "Problems with Data File" \
         -icon warning \
         -message "There are not enough variables in Data file."
   } else {
      # Set the default variable to contour
      set i 2
      set mbv(varIndex) $i
      set mbv(selectedVariable) $mbv(variable$i)
   }; # end if
}; # end proc postVariableMenu

# --------------------------------------------------------

# Some dialog windows...

proc showAboutBox {} {
   tk_messageBox -type ok -title "About mb_view" \
      -message "mb_view V1.01\nP.J. 20-Apr-1999"
}; # end proc showAboutBox

proc showGeneralHelp {} {
   tk_messageBox -type ok -title "General Help" \
      -message "A real CFDer would not have asked for help."
}; # end proc showGeneralHelp

proc rangeDialog {} {
   # Set up the separate window for entering the range values
   global mbv
   global global_xrange
   global global_yrange
   global global_lrange

   # Get current ranges for display
   set global_xrange $mbv(xrange)
   set global_yrange $mbv(yrange)
   set global_lrange $mbv(lrange)

   toplevel .rangeDialog
   set rF [frame .rangeDialog.rangeFrame]
   set bF [frame .rangeDialog.buttonFrame]

   set lX [label $rF.labelX -text "X; min max dX"]
   set lY [label $rF.labelY -text "Y: min max dY"]
   set lL [label $rF.labelL -text "Levels: min max dL"]
   set eX [entry $rF.entryX -textvariable global_xrange \
      -relief sunken -bg white]
   set eY [entry $rF.entryY -textvariable global_yrange \
      -relief sunken -bg white]
   set eL [entry $rF.entryL -textvariable global_lrange \
      -relief sunken -bg white]
   grid $lX $eX
   grid $lY $eY
   grid $lL $eL

   set bAccept [button $bF.accept -text Accept \
      -command {set mbv(xrange) $global_xrange; \
                set mbv(yrange) $global_yrange; \
                set mbv(lrange) $global_lrange} ]
   set bClear [button $bF.clear -text "Clear Edits" \
      -command {set global_xrange $mbv(xrange); \
                set global_yrange $mbv(yrange); \
                set global_lrange $mbv(lrange)} ]
   set bDismiss [button $bF.dismiss -text Dismiss \
      -command {destroy .rangeDialog} ]
   pack $bAccept $bClear $bDismiss -side left -padx 4 -pady 4

   pack $rF $bF
}; # end proc rangeDialog

# --------------------------------------------------------

# Everything is ready; paint the main display...
postMainMenu
initializeDisplay .

# --------------- end of mb_view.tcl ---------------------
