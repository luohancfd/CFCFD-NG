# ssp_demo.tcl
# Demonstration script for ssp graphics functions.
#
# You have two alternatives:
# (1) start wish and then source this file.
#          or
# (2) from the command line, type
#     wish ssp_demo.tcl
#
#----------------------------------------------------------

proc initializeDisplay { parentFrame } {
    puts "Entering initializeDisplay..."
    # Set up a canvas and a couple of scroll bars.

    set display $parentFrame; append display c
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

    $display configure -scrollregion [$display bbox all]
    $display configure -background white
    $display yview moveto 0.5

    return $display
}; # end proc initializeDisplay


proc postMainMenu { canvas_name } {
    # Set up menu bar for the main window
    # m0 is the top-level menubar
    # m1 is a first-level menu

    set m0 [menu .menubar]
    . config -menu $m0
    wm title . "ssp_demo: Demonstration of SSP Plotting Functions."

    set m1 [menu $m0.file -tearoff 0]
    $m0 add cascade -label File -menu $m1
    $m1 add command -label "Run" \
        -command "ssp_tk_demo $canvas_name"
    $m1 add separator
    $m1 add command -label Quit -command "exit"

    set m1 [menu $m0.help -tearoff 0]
    $m0 add cascade -label Help -menu $m1
    $m1 add command -label "About..." \
        -command "showAboutBox"
    $m1 add command -label "General Help" \
        -command "showGeneralHelp"
}; # end proc postMainMenu

# --------------------------------------------------------

# Some dialog windows...

proc showAboutBox {} {
    tk_messageBox -type ok -title "About ssp_demo" \
        -message "ssp_demo V0.1\nP.J. 28-Sep-1999"
}; # end proc showAboutBox

proc showGeneralHelp {} {
    tk_messageBox -type ok -title "General Help" \
        -message "No help text; look at the source code."
}; # end proc showGeneralHelp

#---------------------------------------------------------

puts "Begin main script for demonstration."

puts "Post menu and set up the canvas for display."
set canvas_name [initializeDisplay .]
postMainMenu $canvas_name

puts "Load the C-stuff."
load ./ssp_tk.so ssp_tk

puts "Now, let the menu system take over..."
