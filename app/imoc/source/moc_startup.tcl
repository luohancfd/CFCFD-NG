#!/usr/bin/wish
# moc_startup.tcl
# Startup script for the IMOC program
# PJ, 25-Jan-2000
#     19-Apr-2000 use the Windows-registry to get environment
#

if { [string compare $tcl_platform(platform) windows] == 0 } {
    package require registry 1.0

    proc winBrowserCmd {} {
        # Get the Windows Browser Command
        # PJ, 17-Apr-2000
        #
        # Adapted from the showHtml procedure
        # in the Tcl/Tk Programmer's Reference (C. Nelson)
        #
        set root HKEY_CLASSES_ROOT
        # Look up the application key for HTML files
        set appKey [registry get $root\\.html ""]
        # Get the command for opening a HTML file
        set appCmd [registry get $root\\$appKey\\shell\\open\\command ""]
        # Strip out some unwanted elements, if present
        regsub {%1} $appCmd {} appCmd
        regsub -all {\"} $appCmd {} appCmd
        regsub -- {-nohome} $appCmd {} appCmd
        # Double-up the backslashes for later use with with eval
        regsub -all {\\} $appCmd {\\\\} appCmd
        return $appCmd
    }; # end proc winBrowserCmd
}; # end if

#----------------------------------------------------------------
# Let's Begin...
#----------------------------------------------------------------

if {[ string compare $tcl_platform(platform) "windows" ] == 0 } {
    console show
    puts "Windows Platform."
} else {
    puts "Assume UNIX platform."
}; # end if

# Start with the determination of our environment.
# Global data is stored in the array gd.
#
# For a Unix platform, assume that
# the following environment variables have been set.
#     IMOC_HOME
#     BROWSER
# If they haven't been set, say so and stop.
#
if [ catch { set env(IMOC_HOME) } gd(IMOC_HOME) ] {
    if { [string compare $tcl_platform(platform) windows] == 0 } {
	# On Windows, first option is to use the registry 
	# to find the install directory.
	# Failing that, use the script location and work back.
	if [catch {
                set keyname "HKEY_LOCAL_MACHINE\\SOFTWARE\\cfcfd\\imoc\\0.1"
                puts "registry keyname = $keyname"
                set gd(IMOC_HOME) [registry get $keyname HOME]
                puts "From windows registry: IMOC_HOME = $gd(IMOC_HOME)"
	} result] {
	    puts "Registry entry not found."
	    set gd(IMOC_HOME) ""
	}
    } else {
	# On *nix, try the user's home directory for the imoc_bin subdir.
	if [catch { set env(HOME) } result] {
	}
	set test_name [file join $result "imoc_bin"]
	if [file isdirectory $test_name] {
	    set gd(IMOC_HOME) $test_name
	} else {
	    puts "imoc_bin directory not found."
	    set gd(IMOC_HOME) ""
	}
    }; # end if
}; # end if

if [string equal $gd(IMOC_HOME) ""] {
    # One last place to look for IMOC_HOME...
    puts "Try to use the script info to locate the imoc package."
    set startingScriptName [info script]
    set startingScriptDir [file dir $startingScriptName]
    set currentDir [pwd]
    cd $startingScriptDir
    set startingScriptDir [pwd]
    set gd(IMOC_HOME) \
	[eval file join [lreplace [file split $startingScriptDir] end end]]
    cd $currentDir
}

if [ catch { set env(BROWSER) } gd(htmlViewer) ] {
    if { [string compare $tcl_platform(platform) windows] == 0 } {
        set gd(htmlViewer) [winBrowserCmd]
        puts "From windows registry: HTML Browser is $gd(htmlViewer)"
    } else {
	# On *nix, look for a potential browser.
	if [file exists /usr/bin/firefox] {
	    set gd(htmlViewer) /usr/bin/firefox
	} else {
	    puts stderr "The BROWSER environment variable has not been set"
	    puts stderr "and I couldn't find the firefox browser."
	    puts stderr "Set BROWSER to contain the command for a HTML Browser."
	    exit
	}
    }; # end if
}; # end if

set gd(workDir)    [pwd]
set scriptName     [info script]
puts "IMOC_HOME directory   : $gd(IMOC_HOME)"
puts "Current work directory: $gd(workDir)"
puts "Script Name           : $scriptName"
puts "HTML Viewer           : $gd(htmlViewer)"

# lappend auto_path [file join $gd(IMOC_HOME) source]
# package require imoc
# Manually source the package files.
set packageFiles [list moc_kernel.tcl moc_gui.tcl moc_placard.tcl \
		  moc_plot.tcl moc_nodelist.tcl moc_scales.tcl \
		  moc_unitproc.tcl moc_syn_cmds.tcl moc_menu.tcl]
foreach tclFile $packageFiles {
    source [file join $gd(IMOC_HOME) source $tclFile]
}

if {$tcl_interactive == 1} {
    puts "The console will accept commands."
    puts "Do you want to start the GUI (yes/no)? "
    set answer [gets stdin]
    # puts "answer = $answer"
    set answer [string trim $answer]
    if {[regexp -nocase {^y} $answer] == 1} {
        puts "Start the GUI."
        InitGUI
    } else {
        puts "Load the computational kernel only."
        LoadIMOCKernel
    }; # end if
} else {
    puts "The console will not accept commands; you must work from the GUI."
    InitGUI
}; # end if
