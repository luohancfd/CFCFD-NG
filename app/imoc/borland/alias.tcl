#
# File alias.tcl: Create a Borland compatible .LIB file using
#                 the 'alias' technique
#
# For technical background information go to
#     http://www.bcbdev.com/articles/vcdll.htm
# This script implements the procedure outlined there.
#
# Copyright (c) 2000 Helmut Giese (hgiese@ratiosoft.com)
#
# DISCLAIMER: This script is provided AS IS without warranties of
# any kind. The author assumes no risk or responsibility if you use
# it. However I would like to hear of any corrections or improvements
# you know of.
#
#
# Step 1: Create a 'definition' file from TCLxx.DLL
# -------------------------------------------------
# - Go to the .LIB directory of your Tcl installation, where you want 
#   the resulting .LIB file(s) to be placed.
# - From there we can access ..\bin\tclxx.dll (throughout the rest of 
#   this description we will talk about version 83 instead of 'xx').
# - Issue (assuming the Borland .BIN directory is on your path)
#     IMPDEF tcl83.def0 ..\bin\tcl83.dll
#   The file created by IMPDEF (tcl83.def0) lists all the functions
#   exported by tcl83.dll.
# - The resulting file looks like the following extract - except this was
#   taken from the TK DLL because it offers a broader variety of cases you
#   can expect (the lines containing LIBRARY and EXPORTS start at column 1):
#
#       LIBRARY     TK83.DLL
#
#       EXPORTS
#            [snip]
#         XmbLookupString                @721 ; XmbLookupString
#         TkWinChildProc                = _TkWinChildProc@16
#         _XInitImageFuncPtrs            @722 ; _XInitImageFuncPtrs
#
#   ATTENTION: The format for the second function (TkWinChildProc) differs
#              from the format described in
#                   http://www.bcbdev.com/articles/vcdll.htm
#              However it is returned from the current version of IMPDEF
#              (as of C++ Builder 5) working on tk83.dll, so we have to
#              work with what we've got.
#
#
# Step 2: Create a suitable .DEF file from tcl83.def0
# ------------------------------------------------------
# - Use this script to create tcl83.def from tcl83.def0
# - The syntax is
#     alias.tcl destination source
#   (to stay in tune with the calls to IMPDEF and IMPLIB) so you call
#   it like
#     alias.tcl tcl83.def tcl83.def0
#
#
# Step 3: From the .DEF file create a .LIB file to add to your project
# --------------------------------------------------------------------
# - Issue
#     IMPLIB tcl83bc55.lib tcl83.def
#   to build a .LIB file from the .DEF file created in step 2.
#   This creates a .LIB file tcl83bc55.lib. The suffix 'bc55' has been
#   added to reflect that it is a library for the Borland Compiler 5.5
#   (the same that comes with C++ Builder 5). Adapt as needed.
#


#
# How the script works:
# ---------------------
#
# - We read the file and then process line by line.
# - If the content of the line starts at column 1 or if the line is
#   empty, we just copy it as is.
# - If the line begins with a space we use regular expressions to
#   distinguish 2 cases. What we work on looks like
#       XmbLookupString                @721 ; XmbLookupString
#       TkWinChildProc                = _TkWinChildProc@16
#       _XInitImageFuncPtrs            @722 ; _XInitImageFuncPtrs
#   Line 1 and line 3 are actually of the same type - the developers
#   just chose to let the name of the function begin with an '_'.
#
#   -- The format is (white space reduced)
#         name @nnn ; name
#      Then we construct an 'alias' statement like
#         _name = name
#
#   -- The format is
#         name = _name@nn
#      Then we construct an 'alias' statement like
#         name = _name@nn
#
#   -- The name is neither of the above. Then we send a warning to
#      stderr - no information on what to do in this case.
#

set leadIn "    "               ;# pad left function definitions
set lineCnt 0
set name "\[a-zA-Z0-9_$]+"      ;# for use in regexp: function name
set ws   "\[ \\t]+"             ;# white space

proc showHelp {} {
  puts stdout ""
  puts stdout "Syntax"
  puts stdout "alias outfile infile"
  puts stdout ""
  exit
}

# check if needed arguments are given
if { $argc != 2 } showHelp
if [catch "open [lindex $argv 0] w" outF] {
  puts stdout $outF
  showHelp
}
if [catch "open [lindex $argv 1] r" inF] {
  puts stdout $inF
  showHelp
}

# read input file
while { [gets $inF line] >= 0 } { lappend input $line }
close $inF

# process input
foreach line $input {
  incr lineCnt
  if { [string index $line 0] != " " || [string trim $line] == "" } {
    # echo line as is
    puts $outF $line
  } else {
    if { $lineCnt % 50 == 0 } {
      # show progress
#      puts stderr $line
    }
    # construct alias
    set line [string trim $line]
    set lineOut $leadIn
    if [regexp -expanded \
          "($name) $ws @\[0-9]+ $ws ; $ws $name" \
          $line match fctName] {
      # case 1: name @nnn ; name
      append lineOut "_$fctName\t=\t$fctName"
    } elseif [regexp -expanded \
              "($name) $ws = $ws ($name@\[0-9]+)" \
              $line match name1 name2] {
      # case 2: name = _name@nn
      append lineOut "$name1\t=\t$name2"
    } else {
      # case 3: produce error
      append $lineOut @@@@@     ;# IMPLIB should choke on this
      puts stderr "Encountered unknown case on line $lineCnt: \"$line\""
    }
    puts $outF $lineOut
  }
}
