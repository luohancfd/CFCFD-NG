#!/usr/bin/tclsh
# poshax3-test.tcl
#
# Smoke test the Poshax3 code.
#
# DFP 14/01/2014
#
# Run some short tests of the Poshax3 program. Please add
# test cases relevant to your work.
#
# Usage:
# 1. ./poshax3-text.tcl
# 2. tclsh poshax3-test.tcl

set test_scripts [list "FireII/1634s/Panesi_comparison/FireII.test"]
# lappend test_scripts "some/other/test.test"

set original_dir [pwd]
foreach test_script $test_scripts {
    cd [file dir $test_script]
    puts "[exec date] [pwd]"
    source [file tail $test_script]
    cd $original_dir
}
puts "[exec date] Finished tests."
