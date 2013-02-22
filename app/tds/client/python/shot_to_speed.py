#!/usr/bin/env python
#----------------------------------------------------------------------------
# Header Information
#----------------------------------------------------------------------------

# file: shot_to_speed.py
#
# author: Carolyn Jacobs
# version: 01-Feb-07

# Updated 13-Apr-07 by Carolyn Jacobs
# -- includes the line to set the env to python

#----------------------------------------------------------------------------
# Preamble
#----------------------------------------------------------------------------

"""
This code combines modified versions of two previous codes
- read_tunnel_data.py and shock_time.py.

The code can be used to do the following:
1. Extract the data from the gzipped experimental data files
2. Plot the various channels using Gnuplot (the code to use pylab is still shown)
3. Write the data into a Gnuplot file of time and pressure columns
4. Extract the shock speed as it moves down the tunnel

The shock capturing uses a trigger pressure level to tell if the shock has
arrived at the transducer. As such, if the plots show a strong spike in the
transducer signal ahead of the shock arrival, the time for that transducer
will be incorrect and affect the average shock speed.

The code requires the following Python packages: Numeric, Gnuplot (or pylab if
preferable - make sure the code is altered), os and gzip.
"""

from Numeric import *
import os
import gzip

#------------------------------------------------------------------------------
# Define Functions
#------------------------------------------------------------------------------

# This function finds the shock arrival time for each transducer by using a set
# pressure level as a trigger. At the present time, the trigger is set to 20kPa.
def findshock(t, p):
    found = 0
    p_trigger = 20.0    # 20 kPa is the threshold for the trigger
    p_old = 0.0
    t_old = 0.0
    for linestep in xrange(len(v)):
        if found == 0:
             if ( p[linestep] > p_trigger ):
                 frac = (p_trigger - p_old) / (p[linestep] - p_old)
                 shockt = t_old + frac * (t[linestep] - t_old)
                 found = 1
             t_old = t[linestep]
             p_old = p[linestep]
    return shockt    

# This function extracts the information contained in the gzipped data file
# header section.
def split_metadata(line):
    line = line[1:]  # throw away the sharp character
    line = line.strip()
    words = line.split()
    name = words[0]
    try:
        stringvalue = words[1]
    except:
        stringvalue = "unknown"
    return name, stringvalue

# This function reads the channel data from a directory on the local machine.
def read_channel_from_file(shotDir, fileName):
    """
    Get the channel data from a file on the local disk.

    Returns the metadata as a dictionary and the signal data
    as a pair of arrays containing times and corresponding values.
    """
    f = gzip.open(os.path.join(shotDir, fileName), "rb")
    # Process the header...
    metadata = {}
    line = f.readline().strip()
    while line[0] == "#":
        # process a line of metadata
        name, stringvalue = split_metadata(line)
        metadata[name] = stringvalue
        # Read the following line for the next iteration, if any.
        line = f.readline().strip()
        if len(line) == 0: break
    # We should have come across the blank line separating the
    # file header from the file data.
    nsample = int(metadata["dataPoints"])
    dt = float(metadata["timeInterval"])
    t0 = float(metadata["timeStart"])
    time_data = zeros((nsample,), Float)
    scaled_data = zeros((nsample,), Float)
    for i in range(nsample):
        line = f.readline().strip()
        try:
            if metadata["withTimeColumn"] == "no":
                # Only one number on the line.
                sample = float(line)
                t = t0 + i * dt
            else:
                listOfPieces = line.split()
                t = float(listOfPieces[0])
                sample = float(listOfPieces[1])
            time_data[i] = t
            scaled_data[i] = sample
        except:
            break
    f.close()
    # print "Number of samples read", len(scaled_data)
    return (metadata, time_data, scaled_data)

# This function writes the plotted data into a file in Gnuplot
#format for later use
def write_plot_data(fileName, ylabelText, facilityName, shotName, signalName, t, v):
    # Open the file for writing
    newfileName = shotName + "A." + signalName + ".out"
    f = open(newfileName, "w")

    # Write the headers for the new file
    header_string1 = "# Gnuplot output file for tunnel data. \n"
    header_string2 = "# Data extracted from facility: " + facilityName + "\n"
    header_string3 = "# Data extracted from shot: " + shotName + "\n"
    header_string4 = "# Data extracted from file: " + fileName + "\n"
    header_string5 = "# " + xlabelText + ";" + ylabelText + "\n\n"
    f.write(header_string1)
    f.write(header_string2)
    f.write(header_string3)
    f.write(header_string4)
    f.write(header_string5)

    # Write the data for the file
    for stepData in xrange(len(v)):
        newLine = str(t[stepData]) + " " +  str(v[stepData]) + "\n"
        f.write(newLine)
    
    # Close the file
    f.close()
    return 

#------------------------------------------------------------------------------
# Main Code
#------------------------------------------------------------------------------

if __name__ == "__main__":

    # Print to the screen the details of what is happening.
    print '-----------------------------------------------------------------------'
    print ''
    print 'Python code to calculate the average shock speed down the tunnel '
    print 'using the experimental data files.'
    print ''
    print "author: Carolyn Jacobs"
    print "version: 01-Feb-07"
    print ''
    print '-----------------------------------------------------------------------'

    # Insert the details of the facility and shot of interest
    print 'Starting calculation...'
    facilityName = raw_input("Facility name:")
    shotName = raw_input("Shot name:")
    print ''

    # Select the directory in which the data is stored. This is the level at which
    # the folder containing the shot data files is kept, not the data files
    # themselves. 
    print "Default directory containing facility data: ../"
    testDir = raw_input("Is this correct? (y/n)")
    if testDir == "n":
        rootDir = raw_input("Directory relative to this containing facility data:")
    else:
        rootDir = "../"

    # Decide here whether or not average values of the pressures over the test time
    # are desired
    recAvg = raw_input("Average the transducer pressures over the test time? (y/n):")

    # Decide here whether or not to save the plots of the channels
    recPlot = raw_input("Save the signal plots to a file? (y/n):")

    # Decide here whether or not to write the data to Gnuplot format files
    recFiles = raw_input("Write the data to Gnuplot format files? (y/n):")
    
    # The data files typically have a name constructed as so...
    shotDir = rootDir + facilityName + "/" + shotName

    # Define the channels to be extracted
    # These and the locations of the transducers need to be defined differently for
    # each facility. This is not an input option at the moment as most people only
    # require one facility setup.
    # For X2:
    #  - can't use at6 because location is the same as at5 (divides by zero)
    at1 = "310"
    at2 = "330"
    at3 = "320"
    at4 = "410"
    at5 = "420"
    numSignals = 5
    signals = [ at1, at2, at3, at4, at5 ]
    locations = [ 8.765, 9.015, 10.829, 12.675, 12.855 ]
    flowtime = ravel(zeros((1,numSignals), Float))
    shock_speed = ravel(zeros((1,numSignals-1), Float))
    
    # Loop through the channels to extract the data for each automatically
    for channels in xrange(numSignals):
        
        # The data files typically have a name constructed as so...
        signalName = signals[channels]
        fileName = shotName + "A." + signalName + ".gz"
        metadata, t, v = read_channel_from_file(shotDir, fileName)
        ylabelText = metadata["transducerName"] + " in " + metadata["dataUnits"]
        xlabelText = "time in " + metadata["timeUnits"]

        try:
            import Gnuplot
            # remove the persist=1 below if you do not wish to view the plots
            g = Gnuplot.Gnuplot(persist=1)  
            g.xlabel(xlabelText)
            g.ylabel(ylabelText)
            g.title("Data for facility " + facilityName + ", shot " + metadata['shotName'])
            g.plot(Gnuplot.Data(t, v, with='lines'))
            if recPlot == "y":
                plotname = metadata['shotName'] + metadata["transducerName"] + ".eps"
                g.hardcopy(filename=plotname, terminal='postscript', enhanced=1, mode='eps', color=1, fontname='Times-Roman', fontsize=14)
   # The commented section contains the options for viewing pylab plots of the data if
   # that is preferable. It does not save the files though and the recPlot option will
   # not have any effect.
   #    try:
   #        from pylab import *
   #        plot(t, v)
   #        title("Data for facility " + facilityName + ", shot " + metadata['shotName'])
   #        xlabel(xlabelText)
   #        ylabel(ylabelText)
   #        show()
   #        # junk = raw_input("Press enter...")
        except:
            print ""
            print "Plotting doesn't work but we did get the following stuff..."
            print ylabelText
            print xlabelText
            print "Length of time-series data=", len(v)

        # This gives an average value of pressure for a time interval of choice
        if recAvg == "y":
            print "Averaging channel data through test time for channel", metadata["transducerName"]
            tmin = int(raw_input("Lower limit of test time (microseconds):"))
            tmax = int(raw_input("Upper limit of test time (microseconds):"))
            tt = logical_and(t > tmin, t < tmax)
            avge = average(v, weights=tt)
            print "Average=", avge, metadata["dataUnits"], "over", sum(tt), "samples"
            print ""

        # This writes the data extracted from the gzipped files into Gnuplot
        # format text files
        if recFiles == "y":
            outfile = write_plot_data(fileName, ylabelText, facilityName, shotName, signalName, t, v)

        # Call the function to determine when the shock arrives at the transducer
        flowtime[channels] = findshock(t, v)
        
    # Use the stored data to calculate the shock speed between each point.
    for step in xrange(numSignals-1):
        xstep = locations[step+1] - locations[step]
        tstep = flowtime[step+1] - flowtime[step]
        shock_speed[step] = xstep / tstep * 1.0e6
        print 'Shock speed between', locations[step], 'm and', locations[step+1], 'm is', shock_speed[step], 'm/s'
            
    # Calculate the average shock speed.
    avg_speed = average(shock_speed)
            
    # Print pretty output to the screen.
    print ''
    print 'The average shock speed down tunnel is', avg_speed, 'm/s'
    print '                                                ...Calculation finished'
    print ''
    print '-----------------------------------------------------------------------'

#------------------------------------------------------------------------------
