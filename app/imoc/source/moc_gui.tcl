# moc_gui.tcl
#@doc Graphical User-Interface part of the MOC program.
#@end

# package provide imoc 0.1

#-----------------------------------------------------------

proc InitGUI {} {
    #@proc
    #@doc
    # Initialize the Graphical User Interface elements.
    # These elements include the scrolled canvas upon which the 
    # mesh will be drawn and the menus which provide point-and-click
    # access to the underlying procedures.
    #@end
    global gd
    global env
    global selectedNodeList

    # We need to have the IMOC binary kernel loaded 
    # prior to the GUI components.
    LoadIMOCKernel

    # We only want to initialize once...
    if [info exists gd(GUI_Initialized)] return

    # 1.1 Splash screen
    package require Tk
    displayWaitPlacard "Preparing the GUI; Please wait."

    # 1.2 Global data is stored in the array gd.
    # Set up default and initial values.

    set gd(nodeFile) ""
    set gd(WallFile) ""
    set gd(plotFile) ""

    set gd(dotsPerCm) 30; # a guess for standard screens
    # we will give plotting commands to Tk in pixel units
    set gd(vportXsize) [expr 15.0 * $gd(dotsPerCm)]
    set gd(vportYsize) [expr 10.0 * $gd(dotsPerCm)]
    # set gd(canvasXsize) [expr 2.0 * $gd(vportXsize)]
    # set gd(canvasYsize) [expr 2.0 * $gd(vportYsize)]
    set gd(canvasXsize) [expr 0.8 * [winfo screenwidth .]]
    set gd(canvasYsize) [expr 0.8 * [winfo screenheight .]]
    set gd(canvasXoffset) 40
    set gd(canvasYoffset) 40
    set gd(sameScales) 1
    setXYRanges; # default values implied
    setXYTics; # default values implied
    set gd(showNodeNumbers) 0
    set gd(showCharMesh) 1
    set gd(showStreamlines) 1
    set gd(displayDialogForCoincidentNodes) 1

    set gd(pickAction) pickNode
    set selectedNodeList {}

    # 1.3 Set up the menu and main display
    initMenu
    initDisplay
    # Put the main window near the top left of the screen
    wm geometry . +20+20
    nodelist_Clear
    refreshDisplay

    # 1.4 Everything should be loaded; 
    # we are ready for the user to take control.
    # For the cases when everything loads quickly, 
    # leave the starting placard up a while.
    showStatusMsg "Program successfully initialized."
    after 1000 [list set done 1]
    vwait done
    removeWaitPlacard
    showStatusMsg "Waiting for input."
    set gd(GUI_Initialized) "Initialized"
}; # end proc InitGUI

# --------------------------------------------------------

proc closeGUI {} {
    #@proc
    #@doc
    # Shutdown the Graphical User Interface.
    #@end
    global gd

    if [info exists gd(GUI_Initialized)] {
        nodelist_Clear
        closeDisplay
        closeMenu
        # may need to do some other things as we add features...

        unset gd(GUI_Initialized)
    }; # end if
   
}; # end proc closeGUI

# --------------------------------------------------------
