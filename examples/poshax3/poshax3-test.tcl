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
# 3. tclsh poshax3-test.tcl --long-tests

set long_tests 0
for {set i 0} {$i < $argc} {incr i} {
    set arg [lindex $argv $i]
    if {[string first "long" $arg] >= 0} {
        set long_tests 1
    }
}

set test_scripts [list "FireII/1634s/Panesi_comparison/FireII.test"]
lappend test_scripts "Argon/Glass-Liu/M16.5/without_radiation/Glass-Liu-RUC.test"
if {$long_tests} {
    puts "Do long tests as well as short tests..."
    lappend test_scripts "Argon/Glass-Liu/M16.5/without_radiation/Glass-Liu-RUC.test"
} else {
    puts "Do short tests only..."
}

set original_dir [pwd]
foreach test_script $test_scripts {
    cd [file dir $test_script]
    puts "[exec date] [pwd]"
    source [file tail $test_script]
    cd $original_dir
}
puts "[exec date] Finished tests."
