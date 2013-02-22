# readTFile.tcl
# Procedures to read the set of reduced HEG data files (aka T-files)
# and write a set of ASCII files that are compatible with the 
# td_browser.
#
# PJ @HEG.DLR.Goettingen
#    August 2002
# -------------------------------------------------------------------
#
# Sample use from within tclsh...
# % source readTFile.tcl
# % readTFile ../moncdata/HEG ../moncdata/HEG/T0049500 T0049500
#
# -------------------------------------------------------------------

proc readTFile { srcDir destDir shotName } {
    # Reads the data-description (IXC) file to get the channel
    # parameters and then reads the binary data file to get the
    # measured data, one channel at a time.
    # Input...
    # srcDir   : directory in which T Files are located
    # destDir  : directory into which new format files are written
    # shotName : basic name of the data files 
    #            (usually made from the shot number)

    set f [open [file join $srcDir $shotName.IXC] r]
    set fileContent [read -nonewline $f]
    close $f
    set linesFromIXCFile [split $fileContent "\n"]
    set lineCount1 [llength $linesFromIXCFile]

    set dataFile [open [file join $srcDir $shotName.DAT] r]
    fconfigure $dataFile -translation binary

    if { [file isdirectory $destDir] == 0 } {
        file mkdir $destDir
    }
    set listFileName [join [list $shotName A.LST] ""]
    set channelListFile [open [file join $destDir $listFileName] "w"]

    # Put the lines for a single channel together and then
    # extract the data describing parameters.
    set nChannels 0
    set compoundLine ""
    set lineIsComplete 0
    set nSegments 0
    set lineCount2 0
    foreach thisLine $linesFromIXCFile {
        # puts $thisLine
        incr lineCount2

        # We know that we have a complete channel description 
        # when we see that the first character is not a backslash
        # and we have already been accumulating lines.
        if { [string first "\\" $thisLine] != 0 &&
             [string length $compoundLine] > 0 } {
            set lineIsComplete 1
        }

        if { $lineIsComplete == 0 } {
            append compoundLine $thisLine
            # Every line that begins with a backslash indicates
            # a segment of samples for the current channel. 
            if { [string first "\\" $thisLine] == 0 } {
                incr nSegments
            }
            # Look for the next line only if we know that there
            # are more to process.
            if { $lineCount2 < $lineCount1 } {
                continue
            }
        } 

        # At this point, we should have a complete channel description.
        incr nChannels
        set listOfParameters [split $compoundLine "\\"]

        set ch(id) [string trim [lindex $listOfParameters 0]]
        regsub {:} $ch(id) {_} ch(id)
        regsub {#} $ch(id) {_} ch(id)
        set ch(name) [string trim [lindex $listOfParameters 1]]
        set ch(nSegments) $nSegments
        set ch(yFactor) [string trim [lindex $listOfParameters 2]]
        set ch(yOffset) [string trim [lindex $listOfParameters 3]]
        set ch(yUnits) [string trim [lindex $listOfParameters 4]]
        puts -nonewline "$ch(id) $ch(name) time-segments $ch(nSegments): "
        # puts "$ch(yFactor) $ch(yOffset) $ch(yUnits)"

        set p 4; # starting offset for indexing into parameter list
        for {set i 0} {$i < $ch(nSegments)} {incr i} {
            set ch(N,$i) [string trim [lindex $listOfParameters [incr p]]] 
            set ch(dt,$i) [string trim [lindex $listOfParameters [incr p]]] 
            set ch(t0,$i) [string trim [lindex $listOfParameters [incr p]]] 
            set ch(tUnits,$i) [string trim [lindex $listOfParameters [incr p]]] 
            # puts -nonewline "segment $i: N=[set ch(N,$i)] "
            # puts -nonewline "dt=[set ch(dt,$i)] t0=[set ch(t0,$i)] "
            # puts "units=[set ch(tUnits,$i)]"

            # Now, attempt to read the actual binary data.
            # The original data was written by an Intel-type PC
            # so we guess that the format is 16-bit little-endian.
            set ch(listOfValues,$i) {}
            for {set k 0} {$k < [set ch(N,$i)]} {incr k} {
                set bytes [read $dataFile 2]
                binary scan $bytes "s" signedValue
                lappend ch(listOfValues,$i) $signedValue
            }
            puts -nonewline " read [llength [set ch(listOfValues,$i)]] values,"
        }
        puts ""
        writeChannelData $destDir $shotName ch
        puts $channelListFile "$ch(id) $ch(name)"

        # Reset for the next channel.
        set compoundLine $thisLine
        set lineIsComplete 0
        set nSegments 0
    }; # end foreach

    close $dataFile
    close $channelListFile

    # Some more house-keeping.
    set saveDir [pwd]
    cd $destDir
    puts "gzip ASCII files (if any)"
    set targetPattern [join [list $shotName A *] ""]
    foreach targetFile [glob -nocomplain $targetPattern] {
        catch { exec gzip -f $targetFile } gzipOutput
        # puts "; file $targetFile, gzipOutput is $gzipOutput"
    }; # end for
    # Tidy up, just in case we have dumped core here.
    file delete -force core
    cd $saveDir

    return $nChannels
}; # end proc readIXCFile


proc writeChannelData {destDir shotName channelArray} {
    upvar $channelArray ch

    # Use the unique channel id as the file extension.
    set fileName [join [list $shotName A. [set ch(id)]] ""]
    set f [open [file join $destDir $fileName] "w"]

    # For the moment, 
    # try to make the header compatible with the T4 data.
    puts $f "\# dataSource HEG via T-file and readTFile.tcl"
    puts $f "\# shotName $shotName"
    puts $f "\# channelId [set ch(id)]"
    puts $f "\# transducerName [set ch(name)]"
    puts $f "\# withTimeColumn yes"
    set totalPoints 0
    for {set i 0} {$i < $ch(nSegments)} {incr i} {
        incr totalPoints [set ch(N,$i)]
    }
    puts $f "\# dataPoints $totalPoints"
    puts $f "\# dataType scaled"
    puts $f "\# dataUnits [set ch(yUnits)]"
    puts $f "\# timeStart [set ch(t0,0)]"
    puts $f "\# timeInterval [set ch(dt,0)]"
    puts $f "\# timeUnits [set ch(tUnits,0)]"
    puts $f "\# transducerSensitivity [set ch(yFactor)]"
    puts $f "\# transducerSensitivityUnits none"
    puts $f "\# transducerLocation 0.0"
    puts $f "\# transducerSerialNumber 0"
    puts $f "\# gain 1.0"
    puts $f "\# qfluxgain 1.0"
    puts $f ""

    set y0 [set ch(yOffset)]
    set yFactor [set ch(yFactor)]
    for {set i 0} {$i < $ch(nSegments)} {incr i} {
        # Pick out the current segment
        set t0 [set ch(t0,$i)]
        set dt [set ch(dt,$i)]
        set vList [set ch(listOfValues,$i)]
        for {set k 0} {$k < [set ch(N,$i)]} {incr k} {
            # Reconstruct the time and floating-point data value.
            set t [expr $t0 + $k * $dt]
            set intValue [lindex $vList $k]
            set y [expr $y0 + $intValue * $yFactor]
            puts $f "$t $y"
        }; # end for k
    }; # end for i

    close $f  
}; # end proc writeChannelData
    