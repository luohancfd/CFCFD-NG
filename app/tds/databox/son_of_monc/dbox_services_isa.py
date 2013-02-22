#! /usr/bin/env python
## \file dbox_services.py
## \brief Service functions for the BCD databox
## \author Peter Jacobs
## \version 1.0 19-Oct-2004
##
## This module provides a collection of service functions that make
## interaction with the BCD databox more convenient than directly reading from
## and writing to the registers directly.
## The databox registers are accessed via a separate user-space device driver.
##

import os, time
from copy import copy
try:
    from numarray import *
except:
    try:
        from Numeric import *
    except:
        print "Failed to import either numarray or Numeric"


class BCDDataBox:
    "Provides the basic services for dealing with the BCD databox."
    
    def __init__(self):
        "Does a little checking of the user-space device driver."
        if self.device_driver_is_running():
            print "Databox device driver appears to be running properly."
            if self.databox_is_present():
                print "Databox appears to be present and operating."
            else:
                print "Databox does not appear to be present."
                print "Check that it is turned on and the cable is connected."
        else:
            print "Cannot find device files for databox."
            print "Check that the device driver is running."
        return

    def device_driver_is_running(self):
        "Look for one of the device files."
        return os.path.exists('/dev/databox/trigger_control')

    def databox_is_present(self):
        """
        Returns 1 if the databox is present and turned on.
        
        We assume that ISA ports that do not have actual hardware
        on the far end will have all bits held high.
        Our actual test is to look at all 6 timebase status registers
        because they should always be present and at least one should
        not have the value 0xFF.
        """
        f=open("/dev/databox/timebase_status", "r")
        text_data = f.read()
        f.close()
        stringList = text_data.split(" ")
        all_FF = 1
        for i in range(6):
            register_value = int(stringList[i], 16)
            if register_value != 0xFF: all_FF = 0
        return not all_FF
    
    def arm(self):
        """
        Arms the trigger units by writing 1 to bit 5.

        BCD's procedure for synchronising the sampling should be done
        before arming the box.
        """
        f=open("/dev/databox/trigger_control","w")
        f.write("0x20")
        f.close()
        return

    def synchronise_card(self, card):
        """
        BCD's procedure for synchronising sampling.
        
        Should only be done if there is no active sampling.
        It may be too slow to coordinate this from python;
        maybe build it into the device driver.
        """
        f = open("/dev/databox/card_control", "w")
        f.write("%x 0x80" % card)
        f.close()
        f = open("/dev/databox/card_control", "w")
        f.write("%x 0x00" % card)
        f.close()
        return

    def reset_cards_with_list_and_arm(self, cardList):
        """
        Synchronise specified cards and arm the box.
        
        Card list is the list of card numbers that are
        to be synchronised prior to arming the box.
        This is probably the better service to use because the
        Python functions (synchronise_card and arm) are likely
        to be rather slow.
        """
        flagList = [0, 0, 0, 0, 0, 0, 0]
        for card in cardList:
            if card > 7 or card < 1: continue
            flagList[card-1] = 1
        f = open("/dev/databox/reset_cards_and_arm", "w")
        f.write("%x %x %x %x %x %x %x" % tuple(flagList) )
        f.close()
        return
    
    def trigger(self):
        "Send a trigger signal to all timebases by writing 1 to bit 4."
        f=open("/dev/databox/trigger_control","w")
        f.write("0x10")
        f.close()
        return

    def read_trigger_status(self, unit):
        "Read the status register for the specified trigger unit"
        if unit > 3: unit = 3
        if unit < 1: unit = 1
        control_value = 0x04 | unit
        f=open("/dev/databox/trigger_control", "w")
        f.write("0x%x" % control_value)
        f.close()
        steady = 0
        while not steady: 
            f=open("/dev/databox/trigger_status", "r")
            text_data = f.read()
            f.close()
            low_byte = int(text_data.split()[0], 16)
            steady = low_byte < 0x80   # bit-7 is zero when steady
        return text_data

    def decode_trigger_status(self, text_data):
        "Decode the two hexadecimal values for the trigger unit registers."
        byteList = text_data.split()
        low_byte = int(byteList[0], 16)
        high_byte = int(byteList[1], 16)
        frac = float(high_byte) / 255.0
        # Display the threshold in range 0-1000 to match dial on databox
        K = (1000.0 * frac) + (0.0 * (1.0 - frac))
        if (low_byte > 5) & 0x01 :
            s3="+"
        else:
            s3="-"
        if (low_byte >> 4) & 0x01 :
            c3="DC"
        else:
            c3="AC" 
        if (low_byte >> 3) & 0x01 :
            s2="+"
        else:
            s2="-"
        if (low_byte >> 2) & 0x01 :
            c2="DC"
        else:
            c2="AC"  
        if (low_byte >> 1) & 0x01 :
            s1="+"
        else:
            s1="-"
        if low_byte & 0x01 :
            c1="DC"
        else:
            c1="AC"
        info = {"threshold": K,
                "slope": {1: s1, 2: s2, 3: s3},
                "coupling": {1: c1, 2: c2, 3: c3}}
        return info
        
    def read_card_status(self, ncard, channel=1):
        "Read the card status register."
        f=open("/dev/databox/card_status","w")
        reset = 0
        message = "%x %x %x" % (ncard, channel, reset)
        f.write(message)
        f.close()
        f=open("/dev/databox/card_status","r")
        text_data = f.read()
        f.close()
        return text_data

    def card_is_present(self, ncard):
        """
        Returns 1 if the specified card is present in the databox.
        
        We assume that a card which is present will not have all bits high
        in the status registers.  In the event that the device-driver is not
        running or there is something unusual wrong with the system, return
        that the card is not present.
        """
        try:
            text_data = self.read_card_status(ncard)
            stringList = text_data.split(" ")
            regA = int(stringList[0], 16)
            regB = int(stringList[1], 16)
            card_flag = not(regA == 0xFF and regB == 0xFF)
        except Exception:
            card_flag = 0
        return card_flag
    
    def decode_card_status(self, text_data, print_it=0):
        "Decode the card status to return a dictionary with the information."
        stringList = text_data.split(" ")
        regA = int(stringList[0], 16)
        regB = int(stringList[1], 16)
        sampling = (regB & 0x80) >> 7   # just bit 7
        time_base = (regB & 0x60) >> 5  # extract bits 5 and 6
        if time_base == 0:
            time_base = 1
            time_scale = 4
        else:
            time_scale = 1
        B = regB & 0x1F      # bits 0-4, MSB of buffer pointer
        b = regA    # least-significant byte of buffer pointer
        first_word = B * 256 + b
        info = {"sampling": sampling,
                "time_base": time_base,
                "time_scale": time_scale,
                "first_word": first_word}
        if print_it > 0:
            print "sampling=%d, time_base=%d, first_word=%d" % \
                  (sampling, time_base, first_word)
        return info

    def card_is_still_sampling(self, card):
        """
        Return the sampling status of the card.

        1 == sampling
        0 == hold
        """
        card_info = self.decode_card_status(self.read_card_status(card))
        return card_info["sampling"]

    def still_sampling(self, cardList):
        """
        Returns True if at least one card in the supplied list
        is still sampling.
        """
        count = 0
        for id in cardList:
            if self.card_is_still_sampling(id): count += 1
        return count > 0
    
    def read_timebase_status(self):
        "Read the timebase status registers."
        f=open("/dev/databox/timebase_status","r")
        text_data = f.read()
        f.close()
        return text_data

    def decode_timebase_status(self, text_data, print_it=0):
        "Decode the timebase status to return a nested dictionary with the data."
        stringList = text_data.split(" ")
        resultDict = {}
        for time_base in [1, 2, 3]:
            # Each timebase has two registers
            # and the data is in hexadecimal strings.
            regA = int(stringList[2*(time_base-1)], 16)
            regB = int(stringList[2*(time_base-1)+1], 16)
            # pretrigger_samples
            ccc = (regB & 0x70) >> 4  # bits 4 thru 6
            dddd = regB & 0x0F  # bits 0-3, are least significant BCD digit
            pretrigger_samples = 1000 * ccc + dddd * 100 + 92
            # buffer_size
            bb = (regA & 0xc0) >> 6  # bits 6,7 are buffer-size selector
            sizeList = [0, 2048, 8192, 4096]
            buffer_size = sizeList[bb]
            # trigger_unit
            unitList = [1, 2, 0, 3]
            tt = (regA & 0x30) >> 4  # bits 4,5
            trigger_unit = unitList[tt]
            # sample_period
            m = (regA & 0x08) >> 3  # bit 3
            multiplier = 1
            if m == 1: multiplier *= 100
            rrr = (regA & 0x07)  # bits 0,1,2
            periodList = [0, 50, 20, 10, 5, 2, 1, 0]
            sample_period = periodList[rrr] * multiplier
            # Collect the information for this timebase.
            info = {"buffer_size": buffer_size,
                    "pretrigger_samples": pretrigger_samples,
                    "sample_period": sample_period,
                    "trigger_unit": trigger_unit}
            if print_it > 0:
                print "Timebase", time_base
                print "  buf_size=", buffer_size, "words " \
                      "pretrig=", pretrigger_samples, \
                      "sample_period=", sample_period, "us " \
                      "trigger_unit=", trigger_unit
            resultDict[time_base] = info
        return resultDict

    def read_data_buffer(self, ncard, channel):
        "Read the data buffer for a particular card and channel."
        f=open("/dev/databox/data_buffer","w")
        message = "%x %x" % (ncard, channel)
        f.write(message)
        f.close()
        f=open("/dev/databox/data_buffer","r")
        binary_data = f.read()
        f.close()
        return binary_data

    def decode_data_word(self, low_byte, high_byte):
        """
        Returns a tuple of decoded data from one sample word.
        """
        # The 1000.0 values are obviously not available on the databox dial
        # but we want a non-zero full scale value to check that the ground
        # input really does give close to zero values.
        FSList = [1000.0, 5.0, 2.0, 1.0, 0.5, 0.2, 0.1, 1000.0]
        L = (low_byte & 0xF0) >> 4   # bits 4,5,6 and 7
        H = high_byte
        S = (low_byte & 0x0e) >> 1   # bits 1,2 and 3
        FS = FSList[S]
        sample_volts = ((16.0 * H + L)/2048.0 - 1.0) * FS
        C = (low_byte & 0x01)
        return (sample_volts, S, FS, C)

    def decode_data_buffer(self, text_data, first_word,
                           buffer_size=8192, mux=0, print_it=0):
        """
        Returns a list of signals decoded from the data_buffer bytes.
        Each signal consists of a list of voltage samples.

        text_data   : string of 16384 bytes from the ring buffer
        first_word  : pointer to the oldest word in the ring buffer
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
        # First, get all of the sampled data into one array.
        voltage = zeros((buffer_size,), Float)
        for i in range(buffer_size):
            word_pointer = (first_word + i) % 8192
            low_byte = ord(text_data[word_pointer*2])
            high_byte = ord(text_data[word_pointer*2 + 1])
            sample_volts, S, FS, C = self.decode_data_word(low_byte, high_byte)
            voltage[i] = sample_volts
            if i < 3 and print_it > 0:
                print "i=", i, "bytes=", low_byte, high_byte, \
                      "S=", S, "FS=", FS, "C=", C, "volts=", sample_volts
        ground_input = (S == 7 and C == 0)
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
                # The following code *should* work but it fails
                # -- don't know why.
                # v_demux.append(copy(voltage[m::mux]))
            # Following code is for a list of lists.
            # for m in range(mux):
            #     v_demux.append([])
            # i = 0
            # while i < buffer_size:
            #     for m in range(mux):
            #         v_demux[m].append(voltage[i+m])
            #     i += mux
            
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
    print "Test dbox_services.py..."
    dbox = BCDDataBox()

    print "Look at timebases:"
    text_data = dbox.read_timebase_status()
    print "Timebase text_data:", text_data
    timebase_status = dbox.decode_timebase_status(text_data, 1)

    for i in range(1,4):
        text_data = dbox.read_trigger_status(i)
        print "trigger", i, "text_data=", text_data
        print "    info=", dbox.decode_trigger_status(text_data)
    
    # In the test box, we have card 4 only.
    n = 4
    cardList = [n,]

    print "Prepare to sample:"
    while dbox.still_sampling(cardList):
        time.sleep(0.2)
    # dbox.synchronise_card(n)
    # dbox.arm()
    dbox.reset_cards_with_list_and_arm(cardList)
    print "Have armed databox; waiting 2 seconds before trigger..."
    time.sleep(2)
    dbox.trigger()
    print "Have sent trigger signal."

    print "Look at card status but don't wait for sampling to stop."
    text_data = dbox.read_card_status(n)
    print "Card", n, "text_data=", text_data
    card_status = dbox.decode_card_status(text_data, 1)

    print "Look at card status again. This time, do wait for sampling to stop."
    while dbox.still_sampling(cardList):
        time.sleep(0.2)
    text_data = dbox.read_card_status(n)
    print "Card", n, "text_data=", text_data
    card_status = dbox.decode_card_status(text_data, 1)
    print "Look again. This should return exactly the same card status."
    text_data = dbox.read_card_status(n)
    print "Card", n, "text_data=", text_data
    card_status = dbox.decode_card_status(text_data, 1)

    print "Collect the sampled data."
    c = 1
    text_data = dbox.read_data_buffer(n, c)
    print "data buffer: type=", type(text_data), " size=", len(text_data)
    voltage, ground_input, FS = \
        dbox.decode_data_buffer(text_data, card_status["first_word"], print_it=1)
    print "The first few voltage values are:"
    for i in range(10):
        print voltage[i],
    print ""

    try:
        import Gnuplot
        # Construct the time stamps for the samples
        sample_period = timebase_status[card_status["time_base"]]["sample_period"]
        sample_times = []
        for i in range(len(voltage)):
            sample_times.append(i * sample_period)
        
        g1 = Gnuplot.Gnuplot(debug=1)
        d1 = Gnuplot.Data(sample_times, voltage, with="lines")
        g1.xlabel("time, microseconds")
        g1.ylabel("voltage")
        g1("set yrange [%f:%f]" % (-FS,FS) )
        g1.title("Card %d, Channel %d" % (n, c) )
        g1.plot(d1)
        junk = raw_input("Press return to continue...")
        g1.hardcopy("sample_plot.eps", terminal="postscript", mode="eps",
                    fontsize=20)
    except:
        print "Gnuplot appears to be unavailable."
        
    print "Done."
    
#------------------------------------------------------------------------
