# moc_scales.tcl

#@doc
# Routines to transform between canvas coordinates and
# world coordinates. 
#@end

package provide imoc 0.1

#-------------------------------------------------------------

proc setXYRanges { {xmin 0.0} {ymin 0.0} {xmax 1.2} {ymax 1.0} } { 
    #@proc
    #@doc
    # Set the ranges (in world units) for the x- and y-coordinates.
    #@end
    global gd

    # Re-order values if necessary
    if {$xmin > $xmax} {
        set temp $xmin
        set xmin $xmax
        set xmax $temp
    }; # end if
    if {$ymin > $ymax} {
        set temp $ymin
        set ymin $ymax
        set ymax $temp
    }; # end if

    # Make sure that we have finite ranges and set the ranges
    set gd(xmin) $xmin
    if { $xmax == $xmin } {
        set gd(xmax) [expr $xmin + 1.0]
    } else {
        set gd(xmax) $xmax
    }; # end if
    set gd(ymin) $ymin
    if { $ymax == $ymin } {
        set gd(ymax) [expr $ymin + 1.0]
    } else {
        set gd(ymax) $ymax
    }; # end if

    # Update the variables as displayed in the Range dialog.
    set gd(display_xmin) $gd(xmin)
    set gd(display_xmax) $gd(xmax)
    set gd(display_ymin) $gd(ymin)
    set gd(display_ymax) $gd(ymax)

    setXYScales
}; # end proc

proc setXYScales { } { 
    #@proc
    #@doc
    # Sets the world->canvas scales so that the full range will
    # fit on the canvas with a few pixels to spare.
    #@end
    global gd

    global gd
    set xscale [expr ($gd(canvasXsize) - $gd(canvasXoffset) - 10) / \
        ($gd(xmax) - $gd(xmin))]
    set yscale [expr ($gd(canvasYsize) - $gd(canvasYoffset) - 10) / \
        ($gd(ymax) - $gd(ymin))]

    if { $gd(sameScales) == 1 } {
        # Select the smaller scale such that the whole range fits.
        if { $xscale < $yscale } {
            set gd(xscale) $xscale
            set gd(yscale) $xscale
        } else {
            set gd(xscale) $yscale
            set gd(yscale) $yscale
        }; # end if
    } else {
        # We will allow different scales.
        set gd(xscale) $xscale
        set gd(yscale) $yscale
    }; # end if

    return { $gd(xscale) $gd(yscale) }
}; # end proc

proc setXYTics { {dx 0.2} {dy 0.2} } { 
    #@proc
    #@doc
    # Set the tic-mark spacing for the x- and y-coordinates.
    #@end
    global gd
    set gd(dx) $dx
    set gd(dy) $dy
    set gd(display_dx) $gd(dx)
    set gd(display_dy) $gd(dy)
}; # end proc

#-------------------------------------------------------------

proc canvasX { worldX } {
    #@proc
    global gd
    set xscale $gd(xscale)
    set xmin    $gd(xmin)
    return [expr (($worldX - $xmin) * $xscale) + $gd(canvasXoffset)]
}; # end proc canvasX

proc canvasY { worldY } {
    #@proc
    global gd
    set yscale $gd(yscale)
    set ymin    $gd(ymin)
    set ysize   $gd(canvasYsize)
    return [expr $ysize - $gd(canvasYoffset) - (($worldY - $ymin) * $yscale)]
}; # end proc canvasY

proc worldX { canvasX } {
    #@proc
    global gd
    set xscale $gd(xscale)
    set xmin    $gd(xmin)
    return [expr (($canvasX - $gd(canvasXoffset)) / $xscale) + $xmin]
}; # end proc canvasY

proc worldY { canvasY } {
    #@proc
    global gd
    set yscale $gd(yscale)
    set ymin    $gd(ymin)
    set ysize   $gd(canvasYsize)
    return [expr (($ysize - $canvasY - $gd(canvasYoffset)) / $yscale) + $ymin]
}; # end proc canvasY

#-------------------------------------------------------------
