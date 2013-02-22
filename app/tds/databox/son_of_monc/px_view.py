#! /usr/bin/env python
## \file px_view.py
## \brief Example program to read shot data and display the spatial distribution.
## \author Peter Jacobs
## \version 1.0, 20-Jan-2005
##

versionString = "1.0"

import sys
import os
import gzip
import Tkinter
import tkFileDialog
import tkSimpleDialog
import tkMessageBox
import Pmw
import math
try:
    from numarray import *
except:
    try:
        from Numeric import *
    except:
        print "Failed to import either numarray or Numeric"
from my_progress import myProgressDialog

#------------------------------------------------------------------------

class Signal(object):
    "Somewhere convenient for storage of the signal data items"
    def __init__(self, name="unknown"):
        self.N = 0
        self.dt = 0.0
        self.t0 = 0.0
        self.tmax = 0.0
        self.position = 0.0
        self.name = name
        self.data_units = "kPa"
        self.onplot = Tkinter.IntVar()
        self.onplot.set(1)  # default to showing signal
        return

# Some global data for the application.
signalDict = {}
nominal_time = 0.0
time_window = 0.0
time_min = 0.0
time_max = 1.0
time_inc = 0.01
shotId = "unknown"

#------------------------------------------------------------------------

def split_metadata(line):
    "Split a header line into a name and value pair."
    line = line[1:]  # throw away the sharp character
    line = line.strip()
    words = line.split()
    name = words[0]
    try:
        stringvalue = words[1]
    except:
        stringvalue = "unknown"
    return name, stringvalue

def read_data_from_file(signal, f):
    "Get the formatted signal data from the already-opened file f."
    # Process the header
    metadata = {}
    line = f.readline().strip()
    while line[0] == "#":
        # process a line of metadata
        name, stringvalue = split_metadata(line)
        metadata[name] = stringvalue
        # print "metadata: name=", name, " value=", stringvalue
        # Read the following line for the next iteration, if any.
        line = f.readline().strip()
        if len(line) == 0: break
    # We should have come across the blank line separating the
    # file header from the file data.
    signal.N = int(metadata["dataPoints"])
    signal.dt = float(metadata["timeInterval"])
    signal.t0 = float(metadata["timeStart"])
    signal.position = float(metadata["transducerLocation"])
    signal.data_units = metadata["dataUnits"]
    # signal.name = metadata["transducerName"]
    
    print "Read sampled data: expect %d samples" % (signal.N,)
    signal.sample_time = zeros((signal.N,), Float)
    signal.scaled_data = zeros((signal.N,), Float)
    if metadata["withTimeColumn"] == "no":
        # Only one number on the line.
        for i in range(signal.N):
            line = f.readline().strip()
            samplestring = line
            signal.sample_time[i] = (signal.dt * i) + signal.t0
            signal.scaled_data[i] = float(samplestring)
    else:
        for i in range(signal.N):
            line = f.readline().strip()
            timestamp, samplestring = line.split()
            signal.sample_time[i] = float(timestamp)
            signal.scaled_data[i] = float(samplestring)
    print "Number of samples read", len(signal.scaled_data)
    signal.tmax = signal.sample_time[-1]
    print "time from %f to %f" % (signal.t0, signal.tmax)
    return

def signal_statistics(signal, nominal_time, time_window):
    "Returns the mean and standard deviation in the specified time window."
    t1 = nominal_time - 0.5 * time_window
    t2 = nominal_time + 0.5 * time_window
    mask = logical_and(greater(signal.sample_time, t1),
                       less(signal.sample_time, t2))
    selection = compress(mask,signal.scaled_data)
    N = len(selection)
    if N > 0:
        mean = sum(selection) / N
    else:
        mean = 0.0
    devsqr = power(selection - mean, 2)
    if N > 1:
        stddev = math.sqrt(sum(devsqr) / (N - 1))
    else:
        stddev = 0.0
    print "signal: %s, x=%f, mean=%f, stddev=%f, over %d points" % \
          (signal.name, signal.position, mean, stddev, N)
    return mean, stddev

#------------------------------------------------------------------------

def setup_GUI_menus(root):
    "Establish the menu options and their function calls"
    menuBar = Tkinter.Menu(root)
    
    fileMenu = Tkinter.Menu(menuBar)
    fileMenu.add_command(label="Read shot data", command=import_shot_dialog)
    fileMenu.add_command(label="Save plot (as postscript)",
                         command=save_postscript_plot_dialog)
    fileMenu.add_command(label="Quit", command=quit_dialog)
    menuBar.add_cascade(label="File", menu=fileMenu)

    optionMenu = Tkinter.Menu(menuBar)
    optionMenu.add_checkbutton(label="Hold y-axis at zero", variable=force_y_zero)
    optionMenu.add_checkbutton(label="Include std-deviation bars",
                               variable=include_stddev_bars)
    optionMenu.add_checkbutton(label="Connect the dots", variable=connect_the_dots)
    menuBar.add_cascade(label="Options", menu=optionMenu)

    root.config(menu=menuBar)
    return

def quit_dialog():
    sys.exit()
    return
    
def import_shot_dialog():
    "Load previously recorded data for all signals in a shot."
    global root, signalDict, shotId
    shotDir = tkFileDialog.askdirectory(title="Choose a shot directory")
    if shotDir:
        busy()
        print "Import shot description and data from %s" % shotDir
        (dataDir, shotId) = os.path.split(shotDir)
        print "shotId=", shotId, "dirName=", shotDir
        fileName = os.path.join(shotDir, shotId + "A.LST.gz")
        signalDict = {}
        print "Import signal list from", fileName
        f = gzip.open(os.path.join(shotDir, fileName), "rb")
        # read signal list
        line = f.readline().strip()
        while line:
            signalId, signalName = line.split()
            signalDict[signalId] = Signal(signalName)
            print "signalId %s --> signalName %s" % (signalId, signalName)
            line = f.readline().strip()
        f.close

        pd = myProgressDialog(maxValue=len(signalDict.keys()),
                              title="px_view",
                              message="Importing data from files.")
        for signalId in signalDict.keys():
            root.update()
            signal = signalDict[signalId]
            extn = signalId
            fileName = shotId + "A." + signalId + ".gz"
            print "Reading data file for signal", signalId, "fileName", fileName
            fdata = gzip.open(os.path.join(shotDir, fileName), "rb")
            read_data_from_file(signal, fdata)
            fdata.close
            pd.progress += 1
        pd.close()
        not_busy()
    print "Finished importing data."
    dialog = Pmw.MessageDialog(root, title="Importing data.",
                               defaultbutton=0, buttons=("OK",),
                               message_text="Finished importing data.")
    result = dialog.activate()
    root.update()
    wait_a_while(200) # seem to need a bit of time to clear dialog
    refresh_signal_selection_page()
    refresh_time_range()
    return

def save_postscript_plot_dialog():
    "Offer the option of saving the plot to a file in postscript format."
    global dbox_info, canvasWidget
    fileName = tkFileDialog.asksaveasfilename()
    if fileName:
        c = canvasWidget
        c.postscript(file=fileName, pagewidth='12c')
    return

#------------------------------------------------------------------------

def wait_a_while(ms):
    """Wait for the specified number of milliseconds.

    This function is called in various places to allow time
    for the GUI elements to be updates.  If there are now
    pauses in some of the processor-intensive procedures,
    the GUI becomes disconcertingly unresponsive.
    """
    global root, wait_flag
    print "Waiting...",
    wait_flag = Tkinter.StringVar()
    root.after(ms, lambda : wait_flag.set("done"))
    root.wait_variable(wait_flag)
    print "done."
    return

def busy():
    "Display the watch cursor to indicate that the program is busy."
    global root
    root.config(cursor="watch")
    root.update()
    return

def not_busy():
    "Go back to the default cursor to indicate that the program is not busy."
    global root
    root.config(cursor="")
    root.update()
    return

#------------------------------------------------------------------------

def layout_signal_selection_page(parent):
    "Packs a frame for later filling with check-buttons"
    global signalSelectionFrame
    f = Tkinter.Frame(parent, borderwidth=1)
    f.pack(expand=1, fill="both")
    signalSelectionFrame = f
    return

def refresh_signal_selection_page():
    "Put the check-buttons, one for each signal, into the frame."
    global signalSelectionFrame
    f = signalSelectionFrame
    signals_per_row = 5
    count = 0
    signalList = signalDict.keys()
    signalList.sort()
    for signalId in signalList:
        signal = signalDict[signalId]
        row = count / signals_per_row
        col = count % signals_per_row
        b = Tkinter.Checkbutton(f, text=signal.name,
                                variable=signal.onplot,
                                anchor='w')
        b.grid(row=row, column=col, sticky='w')
        count += 1
    print "Put up %d check boxes" % count
    b = Tkinter.Button(f, text="Refresh time range",
                       command=refresh_time_range)
    b.grid(row=row+1, column=0, columnspan=5, sticky="w")
    return

def refresh_time_range(event=None):
    "The time range should be relevant to the selected signals."
    global signalDict, time_min, time_max, time_inc
    global timeSliderWidget
    count = 0
    for signalId in signalDict.keys():
        signal = signalDict[signalId]
        if signal.onplot.get() == 1:
            count += 1
            print "check signal", count
            if count == 1:
                time_min = signal.t0
                time_max = signal.tmax
            else:
                time_min = min(signal.t0, time_min)
                time_max = max(signal.tmax, time_max)
    time_max = max(time_max, (time_min + 1.0))
    time_inc = max((time_max - time_min)/100.0, 1.0)
    print "new time range from %e to %e" % (time_min, time_max)
    s = timeSliderWidget
    s.configure(from_=time_min, to=time_max, resolution=time_inc)
    return

def refresh_plot(event=None):
    "Completely redraw the plot page with the latest selections."
    global signalDict, force_y_zero, shotId
    global canvasWidget, timeWindowEntryWidget, nominalTimeEntryWidget
    busy()
    # Get signal means in selected period
    nominal_time = float(nominalTimeEntryWidget.get())
    time_window = float(timeWindowEntryWidget.get())
    xyylist = []
    data_units = None
    for signalId in signalDict.keys():
        signal = signalDict[signalId]
        if signal.onplot.get() == 1:
            (mean, stddev) = signal_statistics(signal, nominal_time, time_window)
            xyylist.append((signal.position, mean, stddev))
            if not data_units: data_units = signal.data_units
    xyylist.sort()
    if len(xyylist) == 0:
        not_busy()
        print "No selected signals."
        return
    xlist = [t[0] for t in xyylist]
    ylist = [t[1] for t in xyylist]
    y2list = [t[2] for t in xyylist]
    # Sort through data to get an idea of magnitudes
    if force_y_zero.get() == 1:
        ymin = 0.0
    else:
        ymin = min(ylist)
    ymax = max(max(ylist), ymin+1.0)
    xmin = min(xlist)
    xmax = max(max(xlist), xmin+1.0)
    # Erase all previous items on the canvas   
    c = canvasWidget
    c.delete("boxes")
    c.delete("text")
    c.delete("axes")
    c.delete("datalines")
    c.delete("datapoints")
    # Plot size and position
    xsize = 600
    ysize = 380
    xleft = 10
    ytop  = 10
    # outline box with identifying title
    c.create_rectangle(xleft, ytop, xleft+xsize, ytop+ysize, tags="boxes")
    plotLabel = "Shot %s, time=%.2f microsec" % (shotId, nominal_time)
    c.create_text(xleft+50, ytop+5, text=plotLabel, anchor="nw", tags="boxes")
    # axes
    x1 = xleft + int(0.10 * xsize)
    x2 = xleft + int(0.90 * xsize)
    xaxis_length = x2 - x1
    y0 = ytop + int(0.90 * ysize)
    yaxis_length = int(0.80 * ysize)
    c.create_line(x1, y0-yaxis_length, x1, y0, tags="axes")
    c.create_line(x1, y0, x2, y0, tags="axes")
    c.create_text(x1-5, y0-yaxis_length, text=("%.2f" % ymax), anchor='e', tags="axes")
    c.create_text(x1-5, y0, text=("%.2f" % ymin), anchor='e', tags="axes")
    c.create_text(x1, y0+5, text=("%.2f" % xmin), anchor='n', tags="axes")
    c.create_text(x2, y0+5, text=("%.2f" % xmax), anchor='n', tags="axes")
    c.create_text(x1-5, y0-yaxis_length/2, text=data_units, anchor='e', tags="axes")
    c.create_text(x1+xaxis_length/2, y0+5, text="position", anchor="n", tags="axes")

    if len(signalDict.keys()) > 0:
        r = 4
        coordList = []
        for i in range(len(xlist)):
            xp = x1 + int(xaxis_length * (xlist[i] - xmin) / (xmax - xmin))
            yp = y0 - int(yaxis_length * (ylist[i] - ymin) / (ymax - ymin))
            # Copy the data point for adding the line through the points.
            coordList.append(xp)
            coordList.append(yp)
            if include_stddev_bars.get() == 1:
                # Put the variation bars in place.
                yd = int(yaxis_length * y2list[i] / (ymax - ymin))
                c.create_line(xp-r, yp-yd, xp+r, yp-yd, tags="datalines", fill="red")
                c.create_line(xp-r, yp+yd, xp+r, yp+yd, tags="datalines", fill="red")
                c.create_line(xp, yp-yd, xp, yp+yd, tags="datalines", fill="red")
            # A big dot for the mean value.
            c.create_arc(xp-r, yp-r, xp+r, yp+r, tags="datapoints",
                         extent=359, style="chord", fill="red")
        if len(coordList) > 2 and connect_the_dots.get() == 1:
            # only when there is more than one data point
            c.create_line(coordList, tags="datalines", fill="red")
            c.tag_raise("datapoints", "datalines")
    not_busy()
    return

def set_nominal_time(event=None):
    "Get the nominal time from the slider widget."
    global timeSliderWidget, nominalTimeEntryWidget
    nominal_time = float(timeSliderWidget.get())
    e = nominalTimeEntryWidget
    e.delete(0, Tkinter.END)
    e.insert(0, str(nominal_time))
    print "nominal_time now is", nominal_time
    refresh_plot()
    return

def layout_plot_page(parent):
    "Sets up a canvas for later drawing of the individual plots."
    global canvasWidget, timeWindowEntryWidget, nominalTimeEntryWidget
    global timeSliderWidget, time_min, time_max, time_inc
    f = Tkinter.Frame(parent, borderwidth=1)
    f.pack(expand=1, fill="both")
    c = Tkinter.Canvas(f, bg="gray85", width=620, height=400)
    c.pack()
    canvasWidget = c
    #
    # Entry widgets get a frame of their own.
    fe = Tkinter.Frame(parent, borderwidth=1)
    #
    row = 0
    l = Tkinter.Label(fe, text="nominal_time:")
    l.grid(row=row, column=0, sticky='w')
    f2 = Tkinter.Frame(fe, borderwidth=1)
    e = Tkinter.Entry(f2, bg="white", relief='sunken', width=20)
    e.insert(0, str(nominal_time))
    e.bind('<Return>', refresh_plot)
    e.pack(side="left")
    nominalTimeEntryWidget = e
    s = Tkinter.Scale(f2, orient="horizontal", from_=time_min, to=time_max,
                      resolution=time_inc, command=set_nominal_time,
                      showvalue=0)
    s.pack(side="left")
    timeSliderWidget = s
    f2.grid(row=row, column=1, sticky='w')
    row = 2
    l = Tkinter.Label(fe, text="time_window:")
    l.grid(row=row, column=0, sticky='w')
    e = Tkinter.Entry(fe, bg="white", relief='sunken', width=20)
    e.insert(0, str(time_window))
    e.grid(row=row, column=1, sticky='w')
    e.bind('<Return>', refresh_plot)
    timeWindowEntryWidget = e

    fe.pack(expand=1, fill='both')
    return

#------------------------------------------------------------------------

if __name__ == '__main__':
    global root, wait_flag, force_y_zero
    root = Tkinter.Tk()
    Pmw.initialise()
    
    root.title("px_view Data Viewer")
    busy()
    Pmw.aboutversion(versionString)
    Pmw.aboutcopyright("Centre for Hypersonics 2005")
    Pmw.aboutcontact("Peter Jacobs\nemail: peterj@mech.uq.edu.au")
    if len(sys.argv) > 1:
        first_cmd_arg = sys.argv[1]
    else:
        first_cmd_arg = ""
    if first_cmd_arg.rfind("v") >= 0:
        # Only advertise if the first command line argument is something
        # like --version -version -v version v ...
        about = Pmw.AboutDialog(root, applicationname="px_view Example Data Viewer")
        root.lower(about)  # for advertising...
        root.after(3000, lambda : about.lower())

    force_y_zero = Tkinter.IntVar()
    force_y_zero.set(1)
    include_stddev_bars = Tkinter.IntVar()
    include_stddev_bars.set(1)
    connect_the_dots = Tkinter.IntVar()
    connect_the_dots.set(1)

    setup_GUI_menus(root)
    nb = Pmw.NoteBook(root)
    pSignal = nb.add("Signal Selection")
    pPlot = nb.add("Plotted Data")
    layout_signal_selection_page(pSignal)
    layout_plot_page(pPlot)
    nb.pack(fill=Tkinter.BOTH, expand=1)
    nb.setnaturalsize()
    not_busy()

    # Give control to the GUI elements.
    root.mainloop()
        
#-------------------------------------------------------------------------
