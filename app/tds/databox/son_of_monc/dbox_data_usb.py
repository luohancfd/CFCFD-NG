## \file dbox_data.py
## \brief Data definitions and services for the BCD Databox Viewer program
## \author Peter Jacobs
## \version 1.0, 22-Oct-2004, separated out of dbox_view.py
##          1.01, 4-Apr-2005, Rainer's more complete CSV output is included
##          1.02,10-Jan-2007, Small changes to accommodate the new USB module
##

import sys
import Tkinter
import tkFileDialog
import tkSimpleDialog
import Pmw
import mx.DateTime
try:
    from Numeric import *
except:
    try:
        from numpy import *
    except:
        print "Failed to import either numpy or Numeric"
from dbox_services_usb import *
from my_progress import myProgressDialog

#------------------------------------------------------------------------

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
        self.raw_data      = None
        self.scaled_data   = None
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

class CardInfo:
    "Information for a particular A/D converter card."
    def __init__(self, card_id):
        self.id_number = card_id
        self.used = 0
        self.channels_used = {1:0, 2:0, 3:0}
        # Each A/D channel may have a number of subchannels
        # as specified in the following dictionary.
        self.subchannels = {1:0, 2:0, 3:0}
        self.dt_sample = 0.0
        return

class DataBoxInfo:
    "Management information for the whole databox"
    def __init__(self, dbox_service):
        self.shotId = "shot"
        self.dataDir = "data"
        self.backupDataDir = "data-backup"
        fname = "run_description.template"
        if os.path.exists(fname):
            f = open(fname, "r")
            self.runDescription = f.read()
            f.close()
        else:
            self.runDescription = \
"""Project........
Run number.....
Date...........%s
Blame..........
""" % mx.DateTime.now()
        fname = "signal.config"
        if os.path.exists(fname):
            f = open(fname, "r")
            self.signalConfigText = f.read()
            f.close()
        else:
            self.signalConfigText = ""
        self.signalDict = {}
        self.cardDict = {}
        self.lock = 0
        try:
            if dbox_service.databox_is_present():
                self.databoxIsPresent = 1
            else:
                self.databoxIsPresent = 0
        except Exception:
            print "Something unusual is wrong with the Databox or its device driver."
            junk = raw_input("Press RETURN to continue...")
            # Continue on to allow non-databox related activities such as
            # data-file management, etc
            self.databoxIsPresent = 0
        scan_signal_config(self, dbox_service)
        return

    # We are going to use the lock to allow the check_status functions
    # to be automatically scheduled and, if the system is busy doing something
    # else, we will be able see that and skip on the check_status.
    def wait_to_acquire_lock(self):
        "Block until the lock is available."
        while self.lock > 0:
            pass
        self.lock += 1
        return self.lock

    def release_lock(self):
        "Give back the lock."
        if self.lock > 0:
            self.lock -= 1
        return self.lock

    def lock_is_free(self):
        return self.lock == 0

#------------------------------------------------------------------------

def do_save_all_data(dbox_info):
    """
    Write the data to a gzip-files in a subdirectory named after the shotId.

    This also writes the current run description text.
    Should check for the existence of old data files.
    """
    import gzip
    for dataDir in [dbox_info.dataDir, dbox_info.backupDataDir]:
        shotId = dbox_info.shotId
        shotDir = os.path.join(dataDir, shotId)
        if not os.path.isdir(dataDir):
            os.mkdir(dataDir)
        if not os.path.exists(shotDir):
            os.mkdir(shotDir)

        # Run description text is added to the same directory as
        # the data files.
        file_for_runDesc = shotId + ".txt"
        frd = open(os.path.join(shotDir, file_for_runDesc), "wt")
        frd.write(dbox_info.runDescription)
        frd.close()

        # Signal configuration text is added to the same directory as
        # the data files.
        file_for_signalConfig = shotId + ".config"
        fsc = open(os.path.join(shotDir, file_for_signalConfig), "wt")
        fsc.write(dbox_info.signalConfigText)
        fsc.close()

        # Write the data files, one signal at a time.
        # For each data file, add an entry to the list file.
        file_for_list = shotId + "A.LST.gz"
        flist = gzip.open(os.path.join(shotDir, file_for_list), "wb")
        pd = myProgressDialog(maxValue=len(dbox_info.signalDict.keys()),
                              title="dbox_view",
                              message="Writing data to files.")
        for signalId in dbox_info.signalDict.keys():
            signal = dbox_info.signalDict[signalId]
            extn = str(signalId[0]) + str(signalId[1]) + str(signalId[2])
            flist.write("%s %s\n" % (extn, signal.name) )
            fileName = shotId + "A." + extn + ".gz"
            print "Writing data file for signal", signalId, "fileName", fileName
            fdata = gzip.open(os.path.join(shotDir, fileName), "wb")
            write_scaled_data_to_file(dbox_info, signal, fdata, extn)
            fdata.close
            pd.progress += 1
        flist.close()
        pd.close()
        print "Finished writing data and list files to directory=", dataDir
    #
    return

def do_export_csv_data(dbox_info):
    """
    Write the data to a CSV file in a subdirectory named after the shotId.

    Rainer Kirchhartz' arrangement for writing everything to the CSV file.
    """
    print "Begin exporting data to CSV file."
    dataDir = dbox_info.dataDir
    shotId = dbox_info.shotId
    shotDir = os.path.join(dataDir, shotId)
    file_for_CSV = shotId + ".csv"
    f = open(os.path.join(shotDir, file_for_CSV), "wt")
    signalIdList = dbox_info.signalDict.keys()
    signalIdList.sort()
    f.write("dateTime,")
    for signalId in signalIdList:
        f.write("%s," % mx.DateTime.now())
    f.write("\nshot_id,")
    for signalId in signalIdList:
        f.write("%s," % (dbox_info.shotId,))
    f.write("\nsignal_id,")
    for signalId in signalIdList:
        extn = str(signalId[0]) + str(signalId[1]) + str(signalId[2])
        f.write("%s," % (extn,))
    f.write("\nwithtimecolumn,")
    for signalId in signalIdList:
        f.write("no,")
    f.write("\ndatapoints,")
    max_nsample = 0
    for signalId in signalIdList:
        signal = dbox_info.signalDict[signalId]
        nsample = len(signal.scaled_data)
        f.write("%d," % (nsample,) )
        max_nsample = max(max_nsample, nsample)
    f.write("\ndataType,")
    for signalId in signalIdList:
        f.write("scaled," )
    f.write("\ndataunits,")
    for signalId in signalIdList:
        signal = dbox_info.signalDict[signalId]
        f.write("%s," % (signal.units,) )
    f.write("\ntimestart,")
    for signalId in signalIdList:
        signal = dbox_info.signalDict[signalId]
        f.write("%f," % (signal.t0,) )
    f.write("\ntimeAverageWindow,")
    for signalId in signalIdList:
        f.write("0.0,")
    f.write("\ntimeInterval,")
    for signalId in signalIdList:
        signal = dbox_info.signalDict[signalId]
        f.write("%f," % (signal.dt,) )
    f.write("\ntimeUnits,")
    for signalId in signalIdList:
        f.write("microseconds,")
    f.write("\ntransducerSensitivity,")
    for signalId in signalIdList:
        signal = dbox_info.signalDict[signalId]
        f.write("%f," % (signal.sensitivity,) )
    f.write("\ntransducerSensitivityUnits,")
    for signalId in signalIdList:
        signal = dbox_info.signalDict[signalId]
        f.write("%s," % (signal.units,) )
    f.write("\ntransducer_name,")
    for signalId in signalIdList:
        signal = dbox_info.signalDict[signalId]
        f.write("%s," % (signal.name,) )
    f.write("\nposition,")
    for signalId in signalIdList:
        signal = dbox_info.signalDict[signalId]
        f.write("%f," % (signal.position,) )
    f.write("\ntransducerSerialNumber,")
    for signalId in signalIdList:
        signal = dbox_info.signalDict[signalId]
        f.write("%s," % (signal.serial_number,) )
    f.write("\ntransducerType,")
    for signalId in signalIdList:
        signal = dbox_info.signalDict[signalId]
        f.write("%s," % (signal.transducer_type,) )
    f.write("\ngain,")
    for signalId in signalIdList:
        signal = dbox_info.signalDict[signalId]
        f.write("%s," % (signal.external_gain,) )
    f.write("\nqfluxgain,")
    for signalId in signalIdList:
        signal = dbox_info.signalDict[signalId]
        f.write("1.0,")
    f.write("\nfullScaleVolts,")
    for signalId in signalIdList:
        signal = dbox_info.signalDict[signalId]
        f.write("%s," % (signal.FS,) )
    f.write("\noffsetVolts,")
    for signalId in signalIdList:
        signal = dbox_info.signalDict[signalId]
        f.write("%s," % (signal.offset,) )
    f.write("\n")
    f.write("\n")
    for i in range(max_nsample):
        f.write("%d" % i)
        f.write(",")
        for signalId in signalIdList:
            signal = dbox_info.signalDict[signalId]
            if i < len(signal.scaled_data):
                f.write("%f," % (signal.scaled_data[i],) )
            else:
                f.write(",")
        f.write("\n")
    f.close()
    print "Finished exporting data to CSV file."
    #
    return

def do_export_csv_data_old_version(dbox_info):
    """
    Write the data to a CSV file in a subdirectory named after the shotId.
    """
    print "Begin exporting data to CSV file."
    dataDir = dbox_info.dataDir
    shotId = dbox_info.shotId
    shotDir = os.path.join(dataDir, shotId)
    file_for_CSV = shotId + ".csv"
    f = open(os.path.join(shotDir, file_for_CSV), "wt")
    signalIdList = dbox_info.signalDict.keys()
    signalIdList.sort()
    for signalId in signalIdList:
        extn = str(signalId[0]) + str(signalId[1]) + str(signalId[2])
        f.write("%s," % (extn,))
    f.write("signal_id\n")
    for signalId in signalIdList:
        signal = dbox_info.signalDict[signalId]
        f.write("%s," % (signal.name,) )
    f.write("signal_name\n")
    for signalId in signalIdList:
        signal = dbox_info.signalDict[signalId]
        f.write("%s," % (signal.units,) )
    f.write("units\n")
    for signalId in signalIdList:
        signal = dbox_info.signalDict[signalId]
        f.write("%f," % (signal.position,) )
    f.write("position\n")
    max_nsample = 0
    for signalId in signalIdList:
        signal = dbox_info.signalDict[signalId]
        nsample = len(signal.scaled_data)
        f.write("%d," % (nsample,) )
        max_nsample = max(max_nsample, nsample)
    f.write("nsample\n")
    for signalId in signalIdList:
        signal = dbox_info.signalDict[signalId]
        f.write("%d," % (signal.t0,) )
    f.write("t0\n")
    for signalId in signalIdList:
        signal = dbox_info.signalDict[signalId]
        f.write("%d," % (signal.dt,) )
    f.write("dt\n")
    for i in range(max_nsample):
        for signalId in signalIdList:
            signal = dbox_info.signalDict[signalId]
            if i < len(signal.scaled_data):
                f.write("%f," % (signal.scaled_data[i],) )
            else:
                f.write(",")
        f.write("%d\n" % i)
    f.close()
    print "Finished exporting data to CSV file."
    #
    return

#------------------------------------------------------------------------

def scale_signal_data(signal):
    "Scale the raw data using the sensitivity and the initial level as zero."
    N = 20
    start_level = sum(signal.raw_data[5:5+N])/float(N)
    signal.offset = start_level
    signal.scaled_data = (signal.raw_data - start_level) / signal.sensitivity \
                         / signal.external_gain
    return

def reconstruct_raw_voltages(signal):
    "Presumably we have loaded (old) scaled data from disc."
    signal.raw_data = signal.scaled_data * signal.sensitivity * signal.external_gain \
                      + signal.offset
    return

#------------------------------------------------------------------------

def write_scaled_data_to_file(dbox_info, signal, f, extn):
    "Send the formatted data to already-opened file f."
    f.write("# dataSource dbox_view_usb v1.2 2007\n")
    f.write("# dateTime %s\n" % mx.DateTime.now())
    f.write("# shotName %s\n" % dbox_info.shotId)
    f.write("# channelId %s\n" % extn)
    f.write("# withTimeColumn no\n")
    try:
        nsample = len(signal.scaled_data)
    except:
        nsample = 0
        print "WARNING: channel", extn, "seems to have no data."
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

def read_data_from_file(dbox_info, signal, f, extn):
    "Get the formatted data from the already-opened file f."
    # Default data, in case the file does not tell us all...
    withTimeColumn = 0
    nsample = 0
    signal.units = "Volts"
    signal.t0 = 0.0
    signal.dt = 1.0
    signal.name = "None"
    signal.position = 0.0
    signal.serial_number = "1234"
    signal.transducer_type = "unknown"
    signal.external_gain = 1.0
    signal.FS = 5.0
    signal.offset = 0.0
    # Process the header
    line = f.readline().strip()
    while line[0] == "#":
        # process a line of metadata
        name, stringvalue = split_metadata(line)
        print "metadata: name=", name, " value=", stringvalue
        # if name == "shotName": dbox_info.shotId = stringvalue
        if name == "withTimeColumn": withTimeColumn = stringvalue
        if name == "dataPoints": nsample = int(stringvalue)
        if name == "dataUnits": signal.units = stringvalue
        if name == "timeStart": signal.t0 = float(stringvalue)
        if name == "timeInterval": signal.dt = float(stringvalue)
        if name == "transducerSensitivity": signal.sensitivity = float(stringvalue)
        if name == "transducerSensitivityUnits": signal.units = stringvalue
        if name == "transducerName": signal.name = stringvalue
        if name == "transducerLocation": signal.position= float(stringvalue)
        if name == "transducerSerialNumber": signal.serial_number = stringvalue
        if name == "transducerType": signal.transducer_type = stringvalue
        if name == "gain": signal.external_gain = float(stringvalue)
        if name == "fullScaleVolts": signal.FS = float(stringvalue)
        if name == "offsetVolts": signal.offset = float(stringvalue)
        # Read the following line for the next iteration, if any.
        line = f.readline().strip()
        if len(line) == 0: break
    # We should have come across the blank line separating the
    # file header from the file data.
    print "Read sampled data"
    try:
        # For the older Numeric
        signal.scaled_data = zeros((nsample,), Float)
    except:
        # For the newer numpy
        signal.scaled_data = zeros((nsample,), float)
    for i in range(nsample):
        line = f.readline().strip()
        if withTimeColumn == "no":
            # Only one number on the line.
            samplestring = line
        else:
            # Ignore the timestamp.
            samplestring = line.split()[1]
        try:
            signal.scaled_data[i] = float(samplestring)
        except:
            break
    print "Number of samples read", len(signal.scaled_data)
    return

#------------------------------------------------------------------------

def scan_signal_config(dbox_info, dbox_service):
    """
    Scans the information for each configured signal from text form.

    There should be one line for each signal to be configured.
    Lines starting with a hash (or sharp) character are comments
    (and are ignored).
    """
    # Clear out old config data
    dbox_info.cardIdList = []
    dbox_info.cardDict = {}
    dbox_info.signalDict = {}
    
    print "----------------------------"
    print "Scan Signal Configuration..."
    for line in dbox_info.signalConfigText.split("\n"):
        line = line.strip()
        if len(line) == 0: continue
        if line[0] == "#": continue
        # process the noncomment line
        # The definition of what is expected on a line is given by the
        # following lines. Always check here for the latest definition.
        elements = line.split()
        try:
            signal = SignalInfo( name            = elements[0],
                                 card_id         = int(elements[1]),
                                 channel_id      = int(elements[2]),
                                 subchannel_id   = int(elements[3]),
                                 external_gain   = float(elements[4]),
                                 sensitivity     = float(elements[5]),
                                 units           = elements[6],
                                 position        = float(elements[7]),
                                 serial_number   = elements[8],
                                 transducer_type = elements[9],
                                 print_it        = 1)
        except:
            print "Something is wrong with this signal config line:"
            print line
            print "Will ignore this signal definition and continue with the next."
            continue
        # Each signal is identified by the tuple.
        signalId = (signal.card_id, signal.channel_id, signal.subchannel_id)
        # Check that we don't configure the same combination of
        # card, channel and subchannel twice.
        if dbox_info.signalDict.has_key(signalId):
            print "Configuration error: already done", signalId
            print "Will ignore this redundant definition and continue with the next."
            continue
        else:
            dbox_info.signalDict[signalId] = signal
        if not dbox_info.cardDict.has_key(signal.card_id):
            # New card, add it to the dictionary of cards.
            card = CardInfo(signal.card_id)
            card.used = 1
            dbox_info.cardDict[signal.card_id] = card
        else:
            # Card already exists, get a reference to it.
            card = dbox_info.cardDict[signal.card_id]
        card.channels_used[signal.channel_id] = 1
        if signal.subchannel_id > 0:
            # This signal is part of a multiplexed set.
            card.subchannels[signal.channel_id] += 1
        # end of processing noncomment line
    # end of loop body for each line
    print "Configured cards and channels:"
    for id in dbox_info.cardDict.keys():
        card = dbox_info.cardDict[id]
        print "Card:", card.id_number, \
              "Channels used:", card.channels_used, \
              "SubChannel counts:", card.subchannels
        if not dbox_service.card_is_present(id):
            print "WARNING: this card is not present in the databox."
    print "----------------------------------"
    return

def rescan_signal_config_for_new_scales(dbox_info):
    """
    Keep the old signals but update the signal specifications such as
    transducer sensitivity, etc.
    """
    print "----------------------------"
    print "Rescan Signal Configuration for new sensitivities, etc."
    for line in dbox_info.signalConfigText.split("\n"):
        line = line.strip()
        if len(line) == 0: continue
        if line[0] == "#": continue
        # process the noncomment line
        elements = line.split()
        # Each signal is identified by the tuple.
        card_id       = int(elements[1])
        channel_id    = int(elements[2])
        subchannel_id = int(elements[3])
        signalId = (card_id, channel_id, subchannel_id)
        signal = dbox_info.signalDict[signalId]
        signal.external_gain = float(elements[4])
        signal.sensitivity   = float(elements[5])
        signal.units         = elements[6]
        signal.position      = float(elements[7])
        signal.serial_number = elements[8]
        signal.transducer_type = elements[9]
    print "----------------------------"
    return

#------------------------------------------------------------------------
