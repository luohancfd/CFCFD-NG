#! /usr/bin/env python
## \file tek2som.py
## \brief Tektronics CSV file to Son-of-Monc formatted files.
## \author Peter Jacobs
## \version 1.0, 30-Sep-2005 adapted from pieces of dbox_data.py.
##

import sys
import os
import gzip
try:
    from numarray import *
except:
    try:
        from Numeric import *
    except:
        print "Failed to import either numarray or Numeric"
import mx.DateTime

def print_usage():
    print "Tektronics-CSV to Son-of-Monc Converter..."
    print "Usage:"
    print "$ tek2som.py <shot-directory>"

#-----------------------------------------------------------------

class SignalInfo:
    "Information for a particular data signal."
    def __init__(self, name="noname",
                 card_id=0,
                 channel_id=0,
                 subchannel_id=0,
                 external_gain=1.0,
                 sensitivity=1.0,
                 units="volts",
                 position=0.0,
                 serial_number = "unknown",
                 transducer_type = "unknown",
                 print_it=0):
        self.name          = name
        self.card_id       = card_id
        self.channel_id    = channel_id
        self.subchannel_id = subchannel_id
        self.external_gain = external_gain
        self.sensitivity   = sensitivity
        self.units         = units
        self.position      = position
        self.serial_number = serial_number
        self.transducer_type = transducer_type
        self.time_data     = []
        self.raw_data      = []
        self.scaled_data   = []
        self.FS            = 0.0
        self.offset        = 0.0
        self.dt            = 0.0
        self.t0            = 0.0
        if print_it:
            print "Signal name=%s card_id=%d channel_id=%d subchannel=%d" % \
                  (self.name, self.card_id, self.channel_id, self.subchannel_id)
            print "    external_gain=%f sensitivity=%e units=%s" % \
                  (self.external_gain, self.sensitivity, self.units)
            print "    position=%f serial_number=%s transducer_type=%s" % \
                  (self.position, self.serial_number, self.transducer_type)
        return

signalList = []

def scan_signal_config_file(configFileName):
    """
    Scans the information for each configured signal from text form.

    There should be one line for each signal to be configured.
    Lines starting with a hash (or sharp) character are comments
    (and are ignored).
    """
    global signalDict
    print "Scan Signal Configuration..."
    fp = open(configFileName, "r")
    configText = fp.read()
    fp.close()
    for line in configText.split("\n"):
        line = line.strip()
        if len(line) == 0: continue
        if line[0] == "#": continue
        # process the noncomment line
        # The definition of what is expected on a line is given
        # by the following lines.
        # Always check here for the latest definition.
        elements = line.split()
        signal = SignalInfo( name          = elements[0],
                             card_id       = int(elements[1]),
                             channel_id    = int(elements[2]),
                             subchannel_id = int(elements[3]),
                             external_gain = float(elements[4]),
                             sensitivity   = float(elements[5]),
                             units         = elements[6],
                             position      = float(elements[7]),
                             serial_number = elements[8],
                             transducer_type = elements[9],
                             print_it      = 1)
        signalList.append(signal)
    return


def scale_signal_data(signal):
    "Scale the raw data using the sensitivity and the initial level as zero."
    N = 20
    start_level = sum(signal.raw_data[5:5+N])/float(N)
    signal.offset = start_level
    signal.scaled_data = (signal.raw_data - start_level) / signal.sensitivity \
                         / signal.external_gain
    return


def write_tek_data_to_file(signal, f, extn, shotId):
    "Send the formatted data to already-opened file f."
    f.write("# dataSource tektronics DSO processed by tek2som.py\n")
    f.write("# dateTime %s\n" % mx.DateTime.now())
    f.write("# shotName %s\n" % shotId)
    f.write("# channelId %s\n" % extn)
    f.write("# withTimeColumn no\n")
    nsample = len(signal.scaled_data)
    f.write("# dataPoints %d\n" % nsample)
    f.write("# dataType scaled\n")
    f.write("# dataUnits %s\n" % signal.units)
    f.write("# timeStart %e\n" % signal.t0)
    f.write("# timeAverageWindow 0.0\n")
    f.write("# timeInterval %e\n" % signal.dt)
    f.write("# timeUnits microseconds\n")
    f.write("# transducerSensitivity %e\n" % signal.sensitivity)
    f.write("# transducerSensitivityUnits %s\n" % signal.units)
    f.write("# transducerName %s\n" % signal.name)
    f.write("# transducerLocation %e\n" % signal.position)
    f.write("# transducerSerialNumber %s\n" % signal.serial_number)
    f.write("# transducerType %s\n" % signal.transducer_type)
    f.write("# gain %e\n" % signal.external_gain)
    f.write("# qfluxgain 1.0\n")
    f.write("# fullScaleVolts %e\n" % signal.FS)
    f.write("# offsetVolts %e\n" % signal.offset)
    f.write("\n")
    for i in range(nsample):
        f.write("%e\n" % signal.scaled_data[i])
    return

#--------------------------------------------------------------
# Start of processing...

if len(sys.argv) < 2:
    print_usage()
    sys.exit()
    
# assume that we have been given the shot-directory name
shotDir = sys.argv[1]
print "Shot directory=", shotDir
if os.path.exists(shotDir) and os.path.isdir(shotDir):
    print "Shot directory exists."
    if shotDir[-1] == "/":
        # Remove the trailing slash so that split works its magic.
        shotDir = shotDir[:-1]
    facilityDir, shotId = os.path.split(shotDir)
    print "facilityDir=", facilityDir, "shotId=", shotId
else:
    print "Shot directory doesn't exist; quitting."
    sys.exit()

descriptFileName = os.path.join(shotDir, shotId+".txt")
if not os.path.isfile(descriptFileName):
    print "Text run-description is missing (but we can proceed)."

csvFileName = os.path.join(shotDir, shotId+".CSV")
if not os.path.isfile(csvFileName):
    print "CSV file is missing, quitting."
    sys.exit()

configFileName = os.path.join(shotDir, shotId+".config")
if not os.path.isfile(configFileName):
    print "Config file is missing; quitting."
    sys.exit()

scan_signal_config_file(configFileName)
print "Found", len(signalList), "signal definitions."

print "Begin scanning CSV file..."
fp = open(csvFileName, "r")
csvText = fp.read()
fp.close()
for line in csvText.split("\n"):
    line = line.strip()
    if len(line) == 0: continue
    if line[0] == "#": continue
    if line.find("Volts") >= 0: continue
    # process the noncomment line
    items = line.split(',')
    for j in range(len(signalList)):
        signalList[j].time_data.append(float(items[j*2]))
        signalList[j].raw_data.append(float(items[j*2+1]))
# Although assembled in lists, we really need the data in Numeric arrays
for signal in signalList:
    signal.time_data = array(signal.time_data)
    signal.raw_data = array(signal.raw_data)
print "Read", len(signalList[0].time_data), "samples."

print "Scale raw data..."
for signal in signalList:
    scale_signal_data(signal)
    # Convert times from (presumed) seconds to microseconds.
    signal.t0 = signal.time_data[0] * 1.0e6
    signal.dt = (signal.time_data[1] - signal.time_data[0]) * 1.0e6
    
print "Begin writing Son-of-Monc files..."
# Write the data files, one signal at a time.
# For each data file, add an entry to the list file.
file_for_list = shotId + "A.LST.gz"
flist = gzip.open(os.path.join(shotDir, file_for_list), "wb")
for signal in signalList:
    extn = str(signal.card_id) + str(signal.channel_id) + \
           str(signal.subchannel_id)
    flist.write("%s %s\n" % (extn, signal.name) )
    fileName = shotId + "A." + extn + ".gz"
    print "Writing data file for signal", extn, "fileName", fileName
    fdata = gzip.open(os.path.join(shotDir, fileName), "wb")
    write_tek_data_to_file(signal, fdata, extn, shotId)
    fdata.close
flist.close()

print "Done."
