# file: fetch_tunnel_data.py
# 
# author: Peter J.
# version: 15-Dec-05
#
"""
This module provides functions to fetch the tunnel data from the web server.
"""
import Numeric
import urllib2

# The following recipe for interacting with the web server
# is lifted from Martelli's "Python in a Nutshell".
# It allows us to get past the Basic Authentication challenge
# from the server.
server = "http://www.mech.uq.edu.au/"
x = urllib2.HTTPPasswordMgrWithDefaultRealm()
x.add_password(None, server, "tdguest", "tdpasswd")
auth = urllib2.HTTPBasicAuthHandler(x)
opener = urllib2.build_opener(auth)
urllib2.install_opener(opener)

def fetch_channel_data(facilityName="T4", shotName="7319", channelName="110"):
    """
    Get the channel data as a pair of arrays containing times and
    corresponding values.
    """
    text = fetch_text_from_server(facilityName, shotName, channelName, "data")
    times = []; values = []
    for line in text:
        words = line.split()
        times.append(float(words[0]))
        values.append(float(words[1]))
    return (Numeric.array(times), Numeric.array(values))

def fetch_channel_header(facilityName="T4", shotName="7319", channelName="110"):
    """
    Fetch the metadata from the head of the data file and return as a dictionary.
    """
    text = fetch_text_from_server(facilityName, shotName, channelName, "info")
    headerDict = {}
    for line in text:
        key, value = line.split()
        headerDict[key] = value
    return headerDict

def fetch_text_from_server(facilityName, shotName, channelName, partName):
    """
    Build the HTTP request string, request the specific data.
    Returns the text as a list of lines.
    """
    cgiScript = "cgi-bin/tds/td_server.tcl"
    requestString = "?facility=" + facilityName + "+shot=" + shotName
    if len(channelName) > 0: requestString += "+channel=" + channelName
    if len(partName) > 0: requestString += "+part=" + partName
    fullRequest = server+cgiScript+requestString
    f = urllib2.urlopen(fullRequest)
    text_received = f.readlines()
    return text_received

#----------------------------------------------------------------------------

if __name__ == "__main__":
    print "Begin demonstration of the data fetching functions..."
    facility = "T4"
    shot = "7319"
    chan = "110"
    # print fetch_text_from_server(facility, shot, chan, "info")
    t, v = fetch_channel_data(facility, shot, chan)
    metadata = fetch_channel_header(facility, shot, chan)
    ylabelText = metadata["transducerName"] + " in " + metadata["dataUnits"]
    xlabelText = "time in " + metadata["timeUnits"]
    try:
        from pylab import *
        plot(t, v)
        title("Data retrieved for "+facility+" shot "+shot)
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

        
