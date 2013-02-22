#! /usr/bin/env python
## \file dbox_services_usb.py
## \brief Service functions for the BCD databox connected via USB.
## \author Peter Jacobs
## \version 1.0 08-Jan-2007
## \version 1.1 19-Jan-2007 version 3.1 firmware, checksum, try several ports
##
## This module provides a collection of service functions that make
## interaction with the BCD databox via Luke Hilliard's USB
## supervisory card.  This is far more convenient than directly reading
## from and writing to the registers via the ISA bus extender.
##

import os, sys, time, serial
from copy import copy
try:
    from Numeric import *
except:
    try:
        from numpy import *
    except:
        print "Failed to import either numpy or Numeric"


#-------------------------------------------------------------------

def cleanUpString(str):
    "Clean out new-line, carriage-return, extra spaces and null characters."
    str = str.replace('\n', '')
    str = str.replace('\r', '')
    str = str.replace('\x00', '')
    while str.find('  ') >= 0:
        str = str.replace('  ', ' ')
    return str

def rotate_right(c):
    "Rotate bits in a 16-bit word as per Alex Martelli's 2001 message."
    if c & 1:
        return (c >> 1) | 0x8000
    else:
        return c >> 1

#-------------------------------------------------------------------
    
def open_connection(deviceName, baudrate):
    """
    Opens a USB serial connection and attempts to talk to the databox.
    """
    try:
        print "Try serial port", deviceName, "at", baudrate, "baud"
        # We don't want our program to hang if there is no response.
        sp = serial.Serial(port=deviceName, baudrate=baudrate,
                           bytesize=7, parity='O',
                           stopbits=2, timeout=2.5 )
        print "    Successfully opened serial port", sp.portstr
    except:
        sp = None
        print "    Failed to open serial port", deviceName
        return sp
    if databox_is_present(sp):
        print "    Found databox."
        sp.write('B1'); response = sp.readline()
        sp.flushInput()
        # The verbose response to version command may be several lines.
        sp.write('v'); response = sp.readline()
        while response:
            print cleanUpString(response)
            response = sp.readline()
        sp.write('B0'); response = sp.read(1)
        sp.flushInput()
    else:
        print "    Did not find databox."
        sp.close()
        sp = None
    return sp

def close_connection(sp):
    if sp: sp.close()
    return

def databox_is_present(sp):
    """
    Check for presence of databox by sending a simple command
    and checking the response.
    """
    if sp == None: return 0
    sp.flushInput()
    sp.write('B0'); response = sp.read(1) # should be 'F'
    sp.write('y'); response = sp.read(1)
    if len(response) == 0 or response == '0':
        print "There was no response from the databox."
        print "Check that it is turned on and the cable is connected."
        return 0
    else:
        if response.find('1') >= 0:
            print "Databox responds OK."
            sp.flushInput()
            return 1
        else:
            print "Databox responds strangely."
            print "    response was", response, " but should have been 1"
            sp.flushInput()
            return 0

#-------------------------------------------------------------------

class BCDDataBox:
    "Provides the basic services for dealing with the BCD databox."
    
    def __init__(self):
        """
        Attempts to open a connection and see if the databox is present.
        """
        # The following name will hold the connection to the serial port.
        self.sp = None
        # Assume that USB ports start at COM5.
        if sys.platform == 'win32':
            deviceList = [('COM1:', 115200),
                          ('COM2:', 115200),
                          ('COM3:', 115200),
                          ('COM4:', 115200),
                          ('COM5:', 230400),
                          ('COM6:', 230400),
                          ('COM7:', 230400),
                          ('COM8:', 230400)]
        else:
            deviceList = [('/dev/ttyS0', 115200),
                          ('/dev/ttyS1', 115200),
                          ('/dev/ttyUSB0', 230400),
                          ('/dev/ttyUSB1', 230400)]
        for deviceName,baudrate in deviceList:
            self.sp = open_connection(deviceName, baudrate)
            if self.sp: break
        return

    def finalize(self):
        close_connection(self.sp)
        return

    def databox_is_present(self):
        return databox_is_present(self.sp)
    
    def reset_cards_with_list_and_arm(self, cardList):
        """
        Synchronise specified cards and arm the box.
        
        Card list is the list of card numbers that are
        to be synchronised prior to arming the box.
        This is probably the better service to use because the
        Python functions (synchronise_card and arm) are likely
        to be rather slow.
        """
        for card in cardList:
            if card > 7 or card < 1: continue
            self.sp.write('N'+str(card))
        self.sp.write('A1')
        self.sp.flushInput()
        return
    
    def trigger(self):
        """
        Send a trigger signal to all timebases.
        """
        self.sp.write('T1')
        self.sp.flushInput()
        return

    def enable_threshold_reporting(self, i):
        """
        Write 'z1' to enable reporting of trigger-level threshold.
              'z0' to suppress reporting of triger-level threshold.
              
        When reporting is suppressed, dummy values are returned.
        It seems that the T4 databox does not tolerate any fiddling with
        the k-register when the levels are such that a trigger will occur
        if the box is armed.
        """
        self.sp.write('z'+str(i)); response = self.sp.read(1)
        self.sp.flushInput()
        return
    
    def get_trigger_unit_settings(self, unit):
        """
        Return a dictionary indicating the status of the trigger unit.
        """
        if unit > 3: unit = 3
        if unit < 1: unit = 1
        coupling = 'Unknown'
        slope = 'Unknown'
        threshold = -1
        self.sp.flushInput()
        self.sp.write('c'+str(unit)); response = self.sp.read(1)
        if response == 'D': coupling = 'DC'
        if response == 'A': coupling = 'AC'
        self.sp.write('s'+str(unit)); response = self.sp.read(1)
        if response == 'R': slope = '+'
        if response == 'F': slope = '-'
        self.sp.write('k'+str(unit)); response = self.sp.read(3)
        # The following calculation converts from the range -99..99
        # to something closer to the 5..995 range on the dial.
        threshold = (int(cleanUpString(response)) + 99)*5
        self.sp.flushInput()
        return {"threshold": threshold, "slope": slope, "coupling": coupling}

    def card_is_present(self, ncard):
        """
        Returns 1 is a card is present.

        The test is to assume that for an existing card,
        the latest 15-bit value will be something other
        than 77777(octal).  If a card is not present, the
        resistor networks will pull the bus high.
        """
        if self.sp == None: return 0
        self.sp.flushInput()
        self.sp.write('x'+str(ncard)); response = self.sp.read(1)
        self.sp.flushInput()
        return response == '1'

    def card_is_still_sampling(self, ncard):
        """
        Returns 1 if the card is still sampling, 0 otherwise.
        """
        self.sp.flushInput()
        self.sp.write('a'+str(ncard)); response = self.sp.read(1)
        self.sp.flushInput()
        if response == 'T' or response == '1':
            return 1
        else:
            return 0

    def still_sampling(self, cardList):
        """
        Returns True if at least one card in the supplied list
        is still sampling.
        """
        count = 0
        for id in cardList:
            if self.card_is_still_sampling(id): count += 1
        return count > 0
        
    def get_card_status(self, ncard):
        """
        Return a dictionary containing the card status.
        """
        sampling = self.card_is_still_sampling(ncard)
        self.sp.flushInput()
        self.sp.write('N'+str(ncard))
        self.sp.write('r'); response = self.sp.read(4)
        first_word = int(cleanUpString(response), 16)
        self.sp.write('t'); response = self.sp.read(1)
        time_base = int(response)
        if time_base == 0:
            # this is an assumption built into the hardware
            time_base = 1
            time_scale = 4
        else:
            # time_base is already correct
            time_scale = 1
        self.sp.flushInput()
        return {'sampling': sampling, 'time_base': time_base,
                'time_scale': time_scale, 'first_word': first_word}

    def get_timebase_settings(self, unit, print_it=0):
        """
        Returns a dictionary containing the timebase settings.
        """
        self.sp.flushInput()
        self.sp.write('d'+str(unit)); response = self.sp.read(4)
        pretrigger_samples = int(response, 10)
        self.sp.write('b'+str(unit)); response = self.sp.read(1)
        if response.isdigit():
            buffer_size = int(response) * 1024
        else:
            buffer_size = 0
        self.sp.write('p'+str(unit)); response = self.sp.read(2)
        try:
            sample_period = int(response)
        except:
            sample_period = 0
        self.sp.write('m'+str(unit)); response = self.sp.read(1)
        if response == '1':
            multiplier = 100
        else:
            multiplier = 1
        self.sp.write('u'+str(unit)); response = self.sp.read(1)
        trigger_unit = int(response)
        #
        info = {'buffer_size': buffer_size,
                'pretrigger_samples': pretrigger_samples,
                'sample_period': sample_period*multiplier,
                'trigger_unit': trigger_unit}
        if print_it:
            print "Timebase", unit, ":", info
        self.sp.flushInput()
        return info

    def get_timebase_status(self, print_it=0):
        """
        Returns the details for all timebases in a dictionary.
        """
        resultDict = {}
        for time_base in [1,2,3]:
            resultDict[time_base] = self.get_timebase_settings(time_base)
        return resultDict
    
    def read_data_buffer(self, ncard, channel):
        """
        Read the data buffer for a particular card and channel.
        """
        self.sp.flushInput()
        self.sp.write('N'+str(ncard))
        self.sp.write('D'+str(channel))
        expected_nchar = 20 + 16*1024 + 8
        text_data = self.sp.read(expected_nchar)
        self.sp.flushInput()
        received_checksum = text_data[-8:-4]
        reserved_trailer = text_data[-4:]
        if len(text_data) != expected_nchar or reserved_trailer != 'zzzz':
            print "Response contains", len(text_data), \
                  "characters, expected", expected_nchar
            print "checksum:", received_checksum, "trailer:", reserved_trailer
            return None
        else:
            return text_data

    def decode_data_buffer(self, text_data, first_word=0,
                           buffer_size=8192, mux=0, print_it=0):
        """
        Returns a list of signals decoded from the data_buffer bytes.
        Each signal consists of a list of voltage samples.

        text_data   : string of bytes from the USB controller
        first_word  : pointer to the oldest word in the ring buffer
                      Note that the USB controller sends the data through
                      in correct order so a value of first_word=0 is default.
        buffer_size : number of sampled data points
        mux         : number of multiplexed channels
                      0 == no multiplexing
                      1, 2, 3, 4 == this many signals through the multiplexer
        The convention for decoding multiplexed signals is that
        the first subchannel will have a higher starting level
        than the others.
        This will have been manually set on the multiplexer box and
        is not communicated trhough the electronics.
        """
        # Split the text into sections.
        header = text_data[0:20]
        data_string = text_data[20:20+(16*1024)]
        received_checksum = int(text_data[-8:-4], 16)
        print "Card:", header[0], "Channel:", header[1], "header:", header
        #
        # Pull a couple items out of the header
        FS = float(header[17:20])
        ground_input = int(header[12]) # raw slope and coupling
        #
        # First, get all of the sampled data into one array.
        try:
            # For the newer numpy
            voltage = zeros((buffer_size,), float)
        except:
            # For the older Numeric
            voltage = zeros((buffer_size,), Float)
        chksum = 0
        for i in range(buffer_size):
            word_pointer = (first_word + i) % 8192
            low_byte = ord(data_string[word_pointer*2])
            high_byte = ord(data_string[word_pointer*2 + 1])
            int_value = (high_byte-0x30)*64 + (low_byte-0x30)
            chksum = rotate_right(chksum)
            chksum += int_value
            chksum &= 0xffff
            sample_volts = (float(int_value)/2048.0 - 1.0) * FS
            voltage[i] = sample_volts
            if i < 3 and print_it > 0:
                print "i=", i, "bytes=", low_byte, high_byte, \
                      "volts=", sample_volts
        if chksum == received_checksum:
            print "Check sums matched. OK."
        else:
            print "Check sums did not match:", received_checksum, chksum
        if mux <= 1:
            # One signal only, even if it goes through the multiplexer.
            v_demux_sorted = [voltage,]
        else:
            # For two or more multiplexed signals,
            # separate and then decide which is first.
            # We will end up with a list of signals.
            v_demux = []
            for m in range(mux):
                v_demux.append(voltage[m::mux])
            # It is expected that the first subchannel has the
            # biggest initial voltage level.
            biggest = 0  # Starting guess only.
            biggest_sumv = sum(v_demux[0][5:25])
            for m in range(1,mux):
                sumv = sum(v_demux[m][5:25])
                if sumv > biggest_sumv:
                    biggest = m
                    biggest_sumv = sumv
            # Put in order, starting with the biggest signal.
            v_demux_sorted = []
            for m in range(mux):
                v_demux_sorted.append(v_demux[(biggest+m) % mux])
        return (v_demux_sorted, ground_input, FS)
    

#-----------------------------------------------------------------------

if __name__ == '__main__':
    print "Test dbox_services_usb.py..."
    def delay_seconds(n=1):
        import time
        t0 = time.time()
        while time.time() < t0+n: pass
        return
    
    dbox = BCDDataBox()
    if dbox.sp == None:
        print "Did not find a databox; exiting."
        sys.exit()

    for unit in [1,2,3]:
        print "trigger unit", unit, ":", \
              dbox.get_trigger_unit_settings(unit)

    for unit in [1,2,3]:
        print "timebase unit", unit, ":", \
              dbox.get_timebase_settings(unit)

    card_list = []
    for n in range(1,8):
        if dbox.card_is_present(n):
            print "Card", n, "is present"
            card_list.append(n)
        else:
            print "Card", n, "is missing"

    for card in card_list:
        print "Card", card, "status:", dbox.get_card_status(card)

    print "still_sampling:", dbox.still_sampling(card_list)
    print "Synchronize cards, arm and then trigger box"
    dbox.reset_cards_with_list_and_arm([1,7])
    print "still_sampling:", dbox.still_sampling(card_list)
    delay_seconds(3)
    print "send trigger signal"
    dbox.trigger()
    delay_seconds(3)
    print "still_sampling:", dbox.still_sampling(card_list)

    text_data = dbox.read_data_buffer(ncard=1, channel=1)
    v_sorted, gnd_input, FS = dbox.decode_data_buffer(text_data, print_it=1)
    print "gnd_input:", gnd_input, "FS:", FS
    print "v_sorted:"
    
    dbox.finalize()
