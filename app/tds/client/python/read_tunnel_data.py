# file: read_tunnel_data.py
# 
# author: Peter J.
# version: 16-Dec-05
#
"""
This module provides functions to read the tunnel data from files on your computer.
"""
import Numeric
import os
import gzip

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
    time_data = Numeric.zeros((nsample,), Numeric.Float)
    scaled_data = Numeric.zeros((nsample,), Numeric.Float)
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

#----------------------------------------------------------------------------

if __name__ == "__main__":
    print "Begin demonstration of the data-file reading function..."
    facilityName = "T4"; shotName = "7319"; signalName = "110"
    # The data files typically have a name constructed as so...
    shotDir = "../../moncdata/" + facilityName + "/" + shotName
    fileName = shotName + "A." + signalName + ".gz"
    metadata, t, v = read_channel_from_file(shotDir, fileName)
    ylabelText = metadata["transducerName"] + " in " + metadata["dataUnits"]
    xlabelText = "time in " + metadata["timeUnits"]
    try:
        from pylab import *
        plot(t, v)
        title("Data for facility "+facilityName+", shot "+metadata['shotName'])
        xlabel(xlabelText); ylabel(ylabelText)
        show()
        # junk = raw_input("Press enter...")
    except:
        print "Plotting doesn't work but we did get the following stuff..."
        print ylabelText
        print xlabelText
        print "length of time-series data=", len(v)
    print "Average through test time:"
    tt = Numeric.logical_and(t > 7500, t < 8000)
    avge = Numeric.average(v, weights=tt)
    print "average=", avge, "over", sum(tt), "samples"
    print "Done."

        
