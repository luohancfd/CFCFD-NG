#!/usr/bin/tclsh
# eilmer3-test.tcl
#
# Smoke test the Eilmer3 code.
#
# PJ, 11-Jan-2011, 12-Jul-2011
#
# Presently, we work through the specified directories and explicitly invoke 
# test scripts.  Of course, there must be a better way to do this using the 
# tcltest module.
# Short tests are those which may take up to 20 minutes on my workstation.
# I don't want to wait more than an hour, or so, for the set of short tests
# to run.
#
# Usage:
# 1. ./eilmer3-text.tcl
# 2. tclsh eilmer3-test.tcl
# 3. tclsh eilmer3-test.tcl --long-tests
# 4. ./eilmer3-test.tcl --dummy-run

set long_tests 0
set dummy_run 0
for {set i 0} {$i < $argc} {incr i} {
    set arg [lindex $argv $i]
    if {[string first "long" $arg] >= 0} {
        set long_tests 1
    }
    if {[string first "dummy" $arg] >= 0} {
        set dummy_run 1
    }
}
set test_scripts [list "2D/cone20-simple/cone20.test"]
lappend test_scripts "2D/sod/N2-O2/sod.test"
lappend test_scripts "2D/back-nozzle/back.test"
lappend test_scripts "2D/methane-reactor/psr.test"
lappend test_scripts "2D/estcj/estcj-pitot.test"
lappend test_scripts "2D/channel-with-bump/bump.test"
lappend test_scripts "2D/channel-with-bump/bump-marching.test"
lappend test_scripts "2D/mms_euler/mms_euler.test"
lappend test_scripts "2D/odw/odw.test"
lappend test_scripts "2D/sawada_sphere/ss3.test"
lappend test_scripts "2D/bittker-hydrogen-combustion/hydrogen.test"
lappend test_scripts "2D/radiating-cylinder/gray-gas/MC/cyl.test"
lappend test_scripts "2D/giordano/inf_cyl.test"
lappend test_scripts "3D/simple_ramp/simple_ramp.test"
lappend test_scripts "2D/radiating-cylinder/Argon/MC/cyl.test"
lappend test_scripts "2D/mms/mms.test"
lappend test_scripts "2D/nenzfr-Mach4-nozzle-eq/nozzle-eq-marching.test"
lappend test_scripts "2D/nenzfr-Mach4-nozzle-eq/nozzle-eq.test"
# lappend test_scripts "2D/shock-fitting-sphere/shock-fitting.test"
if {$long_tests} {
    puts "Do long tests as well as short tests..."
    lappend test_scripts "2D/turb-flat-plate/turb_flat_plate.test"
    lappend test_scripts "2D/nenzfr-Mach4-nozzle-noneq/nozzle-noneq.test"
    lappend test_scripts "2D/Rutowski-hemisphere/Ms_12.70/Rutowski-short.test"
} else {
    puts "Do short tests only..."
}
set original_dir [pwd]
foreach test_script $test_scripts {
    cd [file dir $test_script]
    puts "[exec date] [pwd]"
    if { !$dummy_run } {
        source [file tail $test_script]
    }
    cd $original_dir
}
puts "[exec date] Finished tests."
