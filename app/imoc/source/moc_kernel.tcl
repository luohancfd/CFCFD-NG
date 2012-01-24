# moc_kernel.tcl
#@doc
# Tcl part of the IMOC kernel.
#@end
# PJ
# 1998,1999,2000
#

# package provide imoc 0.1

proc LoadIMOCKernel {} {
    #@proc
    #@doc
    # Loads the MOC computational kernel as an extension to
    # the Tcl interpreter.
    #@end
    global gd env tcl_platform

    # We should olny do this procedure once.
    if { [info exists gd(kernelAlreadyLoaded)] } return

    puts "Load IMOC module..."
    set suffix [info sharedlibextension]
    puts "Tcl platform is $tcl_platform(platform)"
    if { [string compare $tcl_platform(platform) windows] == 0 } {
        # For windows, we have to have a specific DLL for 
        # each major version of Tcl/Tk (unless we start using stubs)
        # These DLLs have names "imocXX.dll" where XX are the digits
        # of the major and minor Tcl version numbers.
        set tclver [info tclversion]
        regsub {\.} $tclver {} tclver
        set sharedLibName1 [file join $gd(IMOC_HOME) source imoc$tclver$suffix]
        set sharedLibName2 [file join $gd(IMOC_HOME) borland imoc$tclver$suffix]
    } else {
        # For Unix, one shared object fits all Tcl versions 8.x ...	
        set sharedLibName1 [file join $gd(IMOC_HOME) source imoc$suffix]
        set sharedLibName2 [file join $gd(IMOC_HOME) unix imoc$suffix]
    }; # end if
    if [file exists $sharedLibName1] {
	set sharedLibName $sharedLibName1
    } elseif [file exists $sharedLibName2] {
	set sharedLibName $sharedLibName2
    } else {
	puts "Cannot find module: $sharedLibName1, $sharedLibName2"
	set sharedLibName ""
    }
    puts "Loadable module (shared library) is $sharedLibName"
    if [ catch {load $sharedLibName Imoc} load_result ] {
        bell
        puts "Failed to load imoc shared library."
        puts "Result: $load_result"
        puts "Will exit after 3 seconds..."
        after 3000 [list set done 1]
        vwait done
        exit
    } else {
        puts "Shared object/DLL loaded OK, $load_result"
    }; # end if

    set mocStatus [InitMOC]
    if {$mocStatus == 0} {
        puts "MOC kernel module initialized OK."
    } else {
        bell
        puts "MOC kernel module failed to initialize."
        puts "Will exit after 3 seconds..."
        after 3000 [list set done 1]
        vwait done
        exit
    }; # end if

    set mocStatus [InitWall]
    if {$mocStatus == 0} {
        puts "MOC wall module initialized OK."
    } else {
        bell
        puts "MOC wall module failed to initialize."
        puts "Will exit after 3 seconds..."
        after 3000 [list set done 1]
        vwait done
        exit
    }; # end if


    SetDebugLevel 0
    SetGamma 1.4
    SetAxiFlag 0

    set gd(StreamStepSize) 0.05
    set gd(echoCommands) 0

    # Last thing... 
    # set a flag to indicate that we have already passed this way.
    set gd(kernelAlreadyLoaded) 1
}; # end proc LoadIMOCKernel

# 1.3 Define some Useful functions

proc SetGamma { { gamma 1.4 } } {
    #@proc
    #@doc
    # Is a wrapper for SetGamma_C()
    #@end
    global gd
    SetGamma_C $gamma
    set gd(gamma) $gamma
    set gd(dispaly_gam) $gamma
}; # end function SetGamma

proc SetAxiFlag { { axi_flag 0 } } {
    #@proc
    #@doc
    # Is a wrapper for SetAxiFlag_C()
    #@end
    global gd
    SetAxiFlag_C $axi_flag
    set gd(axiFlag) $axi_flag
}; # end proc SetAxiFlag

proc SetDebugLevel { { dbg_level 0 } } {
    #@proc
    #@doc
    # Is a wrapper for SetDebugLevel_C()
    #@end
    global gd
    SetDebugLevel_C $dbg_level
    set gd(debug) $dbg_level
}; # end proc SetDebugLevel

proc GetNodeData { id args } {
    #@proc
    #@doc
    # Gets all of the data values for a specified node. <BR>
    # <BR>
    # Input: <BR>
    # id   : Node index <BR>
    # args : names of the requested fields <BR>
    #  <BR>
    # Returns: <BR>
    # (1) an empty string if the node doesn't exist OR <BR>
    # (2) All field names and values if args is empty OR <BR>
    # (3) The specific fields requested in name-value pairs. <BR>
    # <BR>
    # The returned string is suitable for sending to SetNodeData
    # so long as a suitable argument list is built from it and
    # eval is used to concatenate the arguments into a single list.
    #@end
    #
    set validFields [list X Y Mach Nu Theta P0 T0 \
        CPlusUp CMinusUp CZeroUp CPlusDown CMinusDown CZeroDown]
    set resultString ""
    set separator ""
    if { [ValidNode $id] } {
        if {[llength $args] == 0} {
            set selectedFields $validFields
        } else {
            set selectedFields $args
        }; # end if
        foreach field $selectedFields {
            set value [GetNodeData_C $id $field]
            append resultString $separator
            set separator " "
            append resultString $field $separator $value
        }; # end foreach
    }; # end if

    return $resultString
}; # end proc SetNodeData


proc SetNodeData { id args } {
    #@proc
    #@doc
    # Gets the data values for a specified node. <BR>
    # <BR>
    # Input: <BR>
    # id   : Node index <BR>
    # args : a list of names of the fields and values to be set <BR>
    #        These data are assumed to be in field-value pairs. <BR>
    # <BR>
    # Returns: <BR>
    #  0 if all was OK <BR>
    # -1 on failure <BR>
    #@end
    #
    if { [ValidNode $id] } {
        if {[GetDebugLevel] >= 1} {
            puts "args = $args"
            puts "length of list = [llength $args]"
        }; # end if
        foreach {field value} $args {
            if {[GetDebugLevel] >= 1} {
                puts "Node $id field $field value $value"
            }; # end if
            SetNodeData_C $id $field $value
        }; # end foreach
        return 0
    } else {
        return -1
    }; # end if

}; # end proc SetNodeData


proc GetNodeDataValue { id field } {
    #@proc
    #@doc
    # Gets the data for one particular field.
    #@end
    set value [lindex [split [GetNodeData $id $field]] 1]
    return $value
}; # end proc GetNodeDataValue


proc DeleteAll {} {
    #@proc
    #@doc
    # Deletes all nodes and walls. (Effectively starting anew.)
    #@end
    for {set iw 0} {$iw <= 1} {incr iw} {
        if { [WallIsPresent $iw] } {
            WallDeletePoints $iw
        }; # end if
    }; # end for

    set nodeid [GetNextNodeId -1]
    while { $nodeid >= 0 } {
        DeleteNode $nodeid
        set nodeid [GetNextNodeId -1]
    }; # end while
}; # end proc
