#! /usr/bin/env python
## \file lvm2som.py
## \brief LabView Measurement file to Son-of-Monc formatted files.
## \author Peter Jacobs
## \version 1.0, 13-Mar-2007 adapted from tek2som.py.
##

import sys
import os
import gzip
try:
    from numpy import *
except:
    try:
        from Numeric import *
    except:
        print "Failed to import either numpy or Numeric"
import mx.DateTime

def print_usage():
    print "LabView Measurement to Son-of-Monc Converter..."
    print "Usage:"
    print "python lvm2som.py <shot-directory>"
    print ""
    print "Note that all of the files are expected to be in the shot directory."
    print "These files include:"
    print "   the LabView Measurement file (.lvm)"
    print "   the run description (.txt)"
    print "   the signal configuration file (.config)"

DEBUG_FLAG = 0

#-----------------------------------------------------------------

class SignalInfo:
    """
    Information for a particular data signal.
    """
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
                             print_it      = DEBUG_FLAG)
        signalList.append(signal)
    return


def scale_signal_data(signal):
    """
    Scale the raw data using the sensitivity and the initial level as zero.
    """
    N = 20
    start_level = sum(signal.raw_data[5:5+N])/float(N)
    signal.offset = start_level
    signal.scaled_data = (signal.raw_data - start_level) / signal.sensitivity \
                         / signal.external_gain
    return


def write_lvm_data_to_file(signal, f, extn, shotId, record_date, record_time):
    """
    Send the formatted data to already-opened file f.
    """
    f.write("# dataSource LabVIEW Measurement file processed by lvm2som.py\n")
    if record_date == "":
        f.write("# dateTime %s\n" % mx.DateTime.now())
    else:
        record_date = record_date.replace('/', '-')
        f.write("# dateTime %s %s\n" % (record_date, record_time))
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
# Configuration data is stored in dictionaries
main_header = {
    "Multi_Headings": "No",
    "Operator": "",
    "Project": "",
    "Separator": "Tab",
    "Time": "",
    "Date": "",
    "X_Columns": "One",
    "Time_Pref": "Relative"
    }
def update_main_header_dictionary(line_of_text):
    items = line_of_text.split(separator)
    for entry in main_header.keys():
        if entry == items[0]:
            main_header[entry] = items[1]
            if DEBUG_FLAG: print entry, main_header[entry]
    return

section = {
    "Channels": 0,
    "Samples": [],
    "Date": [],
    "Time": [],
    "Y_Dimension": [],
    "Y_Unit_Label": [],
    "X_Dimension": [],
    "X_Unit_Label": [],
    "X0": [],
    "Delta_X": []
    }
def update_section_dictionary(line_of_text):
    items = line_of_text.split(separator)
    if items[0] == "Channels":
        section["Channels"] = int(items[1])
        if DEBUG_FLAG: print "Channels", section["Channels"]
    if items[0] == "Samples":
        for chan in range(section["Channels"]):
            section["Samples"].append(int(items[chan+1]))
        if DEBUG_FLAG: print "Samples", section["Samples"]
    if items[0] == "Date":
        for chan in range(section["Channels"]):
            section["Date"].append(items[chan+1])
        if DEBUG_FLAG: print "Date", section["Date"]
    if items[0] == "Time":
        for chan in range(section["Channels"]):
            section["Time"].append(items[chan+1])
        if DEBUG_FLAG: print "Time", section["Time"]
    if items[0] == "Y_Unit_Label":
        for chan in range(section["Channels"]):
            section["Y_Unit_Label"].append(items[chan+1])
        if DEBUG_FLAG: print "Date", section["Y_Unit_Label"]
    if items[0] == "X_Dimension":
        for chan in range(section["Channels"]):
            section["X_Dimension"].append(items[chan+1])
        if DEBUG_FLAG: print "X_Dimension", section["X_Dimension"]
    if items[0] == "X0":
        for chan in range(section["Channels"]):
            section["X0"].append(float(items[chan+1]))
        if DEBUG_FLAG: print "X0", section["X0"]
    if items[0] == "Delta_X":
        for chan in range(section["Channels"]):
            section["Delta_X"].append(float(items[chan+1]))
        if DEBUG_FLAG: print "Delta_X", section["Delta_X"]
    return
def fill_unspecified_section_elements():
    if len(section["X_Unit_Label"]) == 0 and section["X_Dimension"] == "Time":
        for chan in range(section["Channels"]):
            section["X_Unit_Label"].append("Seconds")
    else:
        print "Don't know what to do with X_Dimension."
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
    if DEBUG_FLAG: print "Shot directory exists."
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

lvmFileName = os.path.join(shotDir, shotId+".lvm")
if not os.path.isfile(lvmFileName):
    print "LabView Measurement file is missing, quitting."
    sys.exit()

configFileName = os.path.join(shotDir, shotId+".config")
if not os.path.isfile(configFileName):
    print "Config file is missing; quitting."
    sys.exit()

scan_signal_config_file(configFileName)
print "Found", len(signalList), "signal definitions."

print "Begin scanning LVM file..."
fp = open(lvmFileName, "r")
lvmText = fp.read()
fp.close()
lvmLines = lvmText.split("\n")

# First, determine the field separator.
for line in lvmLines:
    if line.find("Separator") >= 0:
        if line.find("Tab"):
            separator = "\t"
            print "Separator is Tab"
        elif line.find("Comma"):
            separator = ","
            print "Separator is Comma"
        else:
            print "Separator is not comma or tab; quitting."
            sys.exit()
        break

print "Start scanning main header."
current_section = "main_header"
sample_count = 0
for line in lvmLines:
    line = line.strip()
    if len(line) == 0: continue
    if line.find("LabVIEW Measurement") >= 0: continue
    if current_section == "main_header":
        if line.find("***End_of_Header***") >= 0:
            print "End of main header."
            print "Start scanning section header."
            current_section = "section_header"
            continue
        update_main_header_dictionary(line)
    if current_section == "section_header":
        if line.find("***End_of_Header***") >= 0:
            print "End of section header."
            print "Start scanning data."
            current_section = "section_data"
            continue
        update_section_dictionary(line)
    if current_section == "section_data":
        # Process the a line from the data section.
        # Only pull out data for signals that appear in the configuration file.
        # Reconstruct the time from the header data because the columns seem to
        # have only 6 decimal digits and cannot record sub-microsecond increments
        # reliably.
        items = line.split(separator)
        if line.find("X_Value") >= 0: continue
        if main_header["X_Columns"] == "One":
            for j in range(len(signalList)):
                x0 = section["X0"][j]
                dt = section["Delta_X"][j]
                t = x0 + sample_count * dt
                signalList[j].time_data.append(t)
                signalList[j].raw_data.append(float(items[j+1]))
        elif main_header["X_Columns"] == "No":
            for j in range(len(signalList)):
                x0 = section["X0"][j]
                dt = section["Delta_X"][j]
                t = x0 + sample_count * dt
                signalList[j].time_data.append(t)
                signalList[j].raw_data.append(float(items[j]))
        elif main_header["X_Columns"] == "Multi":
            for j in range(len(signalList)):
                x0 = section["X0"][j]
                dt = section["Delta_X"][j]
                t = x0 + sample_count * dt
                signalList[j].time_data.append(t)
                signalList[j].raw_data.append(float(items[j*2+1]))
        else:
            print "Don't know what format each line has."
            print "X_Columns=", main_header["X_Columns"]
            sys.exit()
        sample_count += 1

# Although assembled in lists, we really want the data in Numeric arrays
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
    write_lvm_data_to_file(signal, fdata, extn, shotId,
                           main_header["Date"], main_header["Time"])
    fdata.close
flist.close()

print "Done."
