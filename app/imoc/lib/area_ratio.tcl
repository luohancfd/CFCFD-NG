# area_ratio.tcl
# Area ratio - Mach number relationship 
# for quasi-one-dimensional flow of an ideal gas.
#
# PJ 04-Oct-00
#--------------------------------------------------------------

proc area_ratio { M g } {
    # Input: M : Mach number
    #        g : ratio of specific heats
    # Output: returns the area ratio A/A_star
    #
    set tr [expr 1.0 + 0.5 * ($g - 1.0) * $M * $M]
    set ex [expr 0.5 * ($g + 1.0) / ($g - 1.0)]
    set t1 [expr 2.0 / ($g + 1.0)]
    set aratio [expr (1.0 / $M) * pow( $t1 * $tr, $ex ) ]
    return $aratio
}; # end proc area_ratio
