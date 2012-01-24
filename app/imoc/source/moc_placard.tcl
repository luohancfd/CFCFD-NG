# moc_placard.tcl

#@doc
# Tcl procedures to post a message to the middle of the screen
# and to later remove that window.
# The user should not need to call these procedures directly.
#@end

package provide imoc 0.1

proc displayWaitPlacard { msgText } {
    #@proc
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
    #@proc
    update
    destroy .waitPlacard
}; # end proc
