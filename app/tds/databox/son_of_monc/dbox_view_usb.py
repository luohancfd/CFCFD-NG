#! /usr/bin/env python
## \file dbox_view_usb.py
## \brief Manager/viewer for the BCD Databox
## \author Peter Jacobs
## \version 1.0, 20-Oct-2004
## \version 1.1, 29-Jan-2007, configuration file, small UI tweaks
## \version 1.2, 05-Feb-2007, select threshold-reporting behaviour
##

versionString = "1.2"

import sys
import os
import gzip
import Tkinter
import tkFileDialog
import tkSimpleDialog
import tkMessageBox
import Pmw
try:
    from Numeric import *
except:
    try:
        from numpy import *
    except:
        print "Failed to import either numpy or Numeric"
from dbox_services_usb import *
from dbox_data_usb import *
from my_progress import myProgressDialog

# Global data
global dbox_info, dbox_service, dbox_lock

print "Initialize..."
dbox_service = BCDDataBox()
dbox_info = DataBoxInfo(dbox_service)

#------------------------------------------------------------------------

def setup_GUI_menus(root):
    "Establish the menu options and their function calls"
    menuBar = Tkinter.Menu(root)
    
    fileMenu = Tkinter.Menu(menuBar)
    fileMenu.add_command(label="Select data directory",
                         command=select_data_dir_dialog)
    fileMenu.add_command(label="Select backup directory",
                         command=select_backup_data_dir_dialog)
    fileMenu.add_separator()
    fileMenu.add_command(label="Read signal config",
                         command=read_config_file_dialog)
    fileMenu.add_command(label="Save signal config",
                         command=save_config_file_dialog)
    fileMenu.add_separator()
    fileMenu.add_command(label="Read run description",
                         command=read_run_description_dialog)
    fileMenu.add_command(label="Save run description",
                         command=save_run_description_dialog)
    fileMenu.add_separator()
    fileMenu.add_command(label="Read shot data from disc", command=import_shot_dialog)
    fileMenu.add_command(label="Save data", command=save_data_dialog)
    fileMenu.add_command(label="Export data to CSV file", command=export_csv_dialog)
    fileMenu.add_separator()
    fileMenu.add_command(label="Save all plots (as postscript)",
                         command=save_postscript_allplot_dialog_askName)
    fileMenu.add_command(label="Save plots on screen (as postscript)",
                         command=save_postscript_plot_dialog_askName)
    fileMenu.add_command(label="Save plots on screen as shot-N.ps",
                         command=save_postscript_plot_dialog_automatic)
    fileMenu.add_separator()
    fileMenu.add_command(label="Quit", command=quit_dialog)
    menuBar.add_cascade(label="File", menu=fileMenu)

    root.config(menu=menuBar)
    return

def quit_dialog():
    global root, dataLoadedButNotSaved
    if dataLoadedButNotSaved:
        dialog = Pmw.MessageDialog(root, title="Quit Databox Manager",
                     defaultbutton=0, buttons=("Quit anyway", "Oops"),
                     message_text="You have loaded data which has not been saved.")
        result = dialog.activate()
    else:
        result = "Quit anyway"
    if result == "Quit anyway":
        sys.exit()
    return
    
def read_config_file_dialog(fileName = None):
    """
    Offer the option of selecting a configuration file.

    """
    global dbox_info, dbox_service, signalConfigTextWidget
    if fileName == None:
        fileName = tkFileDialog.askopenfilename()
    if fileName:
        f = open(fileName, "r")
        dbox_info.signalConfigText = f.read()
        f.close
        T = signalConfigTextWidget
        T.delete('1.0', Tkinter.END)
        T.insert(Tkinter.END, dbox_info.signalConfigText)
        scan_signal_config(dbox_info, dbox_service)
        print "Finished reading config file."
    return

def save_config_file_dialog():
    """
    Offer the option of saving a configuration file.

    """
    global dbox_info, dbox_service, signalConfigTextWidget
    fileName = tkFileDialog.asksaveasfilename()
    if fileName:
        T = signalConfigTextWidget
        dbox_info.signalConfigText = T.get('1.0',Tkinter.END)
        scan_signal_config(dbox_info, dbox_service)
        f = open(fileName, "wt")
        f.write(dbox_info.signalConfigText)
        f.close
        print "Finished writing config file."
    return

def read_run_description_dialog(fileName = None):
    """
    Look for an existing run description to import.
    """
    global dbox_info, runDescTextWidget
    if fileName == None:
        fileName = tkFileDialog.askopenfilename()
    if fileName:
        f = open(fileName, "r")
        dbox_info.runDescription = f.read()
        f.close()
        T = runDescTextWidget
        T.delete('1.0', Tkinter.END)
        T.insert(Tkinter.END, dbox_info.runDescription)
    return

def save_run_description_dialog():
    """
    Offer the option of saving the run description text to a file.
    """
    global dbox_info, runDescTextWidget
    fileName = tkFileDialog.asksaveasfilename()
    if fileName:
        T = runDescTextWidget
        dbox_info.runDescText = T.get('1.0',Tkinter.END)
        f = open(fileName, "wt")
        f.write(dbox_info.runDescText)
        f.close
    return

def save_postscript_allplot_dialog_askName():
    """
    Offer the option of saving all plots to a file in postscript format.
    """
    global dbox_info
    fileName = tkFileDialog.asksaveasfilename()
    if fileName:
        save_postscript_plot(fileName, 1)
    return

def save_postscript_plot_dialog_askName():
    """
    Offer the option of saving the screen view to a file in postscript format.
    """
    global dbox_info
    fileName = tkFileDialog.asksaveasfilename()
    if fileName:
        save_postscript_plot(fileName, 0)
    return

def save_postscript_plot_dialog_automatic():
    """
    Save the plot to an automatically-named file in postscript format.
    """
    global dbox_info
    # Try to build a unique file name for the plot.
    shotId = dbox_info.shotId
    N = 1
    fileName = shotId + "-"+ str(N) + ".ps"
    while os.path.exists(fileName):
        N = N + 1
        fileName = shotId + "-"+ str(N) + ".ps"
    if fileName:
        save_postscript_plot(fileName, 0)
    return

def save_postscript_plot(fileName, allPlots=0):
    """
    Do the actual work of saving the plot in postscript format.
    """
    global canvasWidget, canvasWidth, canvasHeight
    c = canvasWidget
    if allPlots:
        # Render all plots to postscript.
        try:
            x1, y1, x2, y2 = c.bbox("boxes")
        except:
            # in case there are no items found on the canvas
            x1 = 0; x2=canvasWidth; y1=0; y2=canvasHeight
        print "bbox: x1=", x1, "y1=", y1, "x2=", x2, "y2=", y2
        width = x2 - x1
        height = y2 - y1
        if height > 2 * width:
            c.postscript(file=fileName, y=y1, height=height, pageheight='24c')
        else:
            c.postscript(file=fileName, y=y1, height=height, pagewidth='12c')
    else:
        # Render just the plots that are presently seen on screen.
        c.postscript(file=fileName, pagewidth='12c')
    return

def arm_databox_dialog():
    """
    After asking for confirmation
    (because it will destroy any data in the databox buffers),
    do synchronisation on all used cards and arm (i.e. start sampling)
    """
    global root
    dialog = Pmw.MessageDialog(root, title="Preparing to arm databox",
             defaultbutton=0, buttons=("Arm Databox", "Cancel"),
             message_text="Beware that arming will wipe data in databox buffers.")
    result = dialog.activate()
    if result == "Arm Databox":
        dbox_info.wait_to_acquire_lock()
        dbox_service.reset_cards_with_list_and_arm(dbox_info.cardDict.keys())
        dbox_info.release_lock()
    return

def trigger_databox_dialog():
    "Send trigger signal to databox without asking for confirmation."
    dbox_info.wait_to_acquire_lock()
    dbox_service.trigger()
    dbox_info.release_lock()
    return

def select_data_dir_dialog():
    "Bring up a GUI dialog for directory selection."
    global dbox_info, dataDirEntryWidget
    dataDir = tkFileDialog.askdirectory()
    if dataDir:
        dbox_info.dataDir = dataDir
        e = dataDirEntryWidget
        e.delete(0, Tkinter.END)
        e.insert(0, dataDir)
    return

def select_backup_data_dir_dialog():
    "Bring up a GUI dialog for directory selection."
    global dbox_info, backupDataDirEntryWidget
    dataDir = tkFileDialog.askdirectory()
    if dataDir:
        dbox_info.backupDataDir = dataDir
        e = backupDataDirEntryWidget
        e.delete(0, Tkinter.END)
        e.insert(0, dataDir)
    return

#------------------------------------------------------------------------

def do_collect_data():
    """
    Fetch the data from the databox buffers for the actively-used cards.

    If a card is still sampling when we get to it, wait a little and
    check again every so often.
    """
    global root, dbox_info, dbox_service, dataLoadedButNotSaved
    dbox_info.wait_to_acquire_lock()
    # Do a check of all of the configured cards
    # to see if they have stopped sampling.
    cardList = dbox_info.cardDict.keys()
    ncss = 0
    for card_id in cardList:
        card = dbox_info.cardDict[card_id]
        print "Collect data for card:", card.id_number
        card_status = dbox_service.get_card_status(card_id)
        ncss = card_status["sampling"]
    if ncss > 0:
        print "Cards still sampling; collect failed."
        dialog = Pmw.MessageDialog(root, title="Collect data failed",
                     defaultbutton=0, buttons=("OK", "This is Crap"),
                     message_text="You have tried to collect the sampled data\n" +
                                   "before all of the cards have finished sampling.")
        result = dialog.activate()
        root.update()
        wait_a_while(200) # seem to need a bit of time to clear dialog
        print "Ignored request to collect data."
        dbox_info.release_lock()
        return
    else:
        # Do the collection.
        busy()
        timebase_status = dbox_service.get_timebase_status(1)
        pd = myProgressDialog(maxValue=len(cardList), title="dbox_view",
                              message="Collecting data from the A/D cards.")
        for card_id in cardList:
            root.update()
            card = dbox_info.cardDict[card_id]
            print "Collect data for card:", card.id_number
            card_status = dbox_service.get_card_status(card_id)
            card_dt = float(timebase_status[card_status["time_base"]]["sample_period"]) / \
                      float(card_status["time_scale"])
            buffer_size = timebase_status[card_status["time_base"]]["buffer_size"]
            for channel_id in card.channels_used.keys():
                mux = card.subchannels[channel_id]
                print "Channel:", channel_id, "SubChannel count:", mux
                attempt_count = 0
                text_data = None
                while text_data == None and attempt_count < 3:
                    text_data = dbox_service.read_data_buffer(card_id, channel_id)
                    attempt_count += 1
                if text_data == None:
                    print "Failed to collect data in three attempts."
                    print "Proceeding to the next channel."
                    continue
                voltage_demux, ground_input, FS = \
                    dbox_service.decode_data_buffer(text_data, 0, buffer_size, mux, 1)
                if mux < 1:
                    # No multiplexing
                    subchannel_id = 0
                    signalId = (card_id, channel_id, subchannel_id)
                    if not dbox_info.signalDict.has_key(signalId): continue
                    try:
                        signal = dbox_info.signalDict[signalId]
                        signal.dt = card_dt
                        signal.raw_data = voltage_demux[0]
                        signal.FS = FS
                        print "Signal", signalId, "dt_sample=", signal.dt, \
                              "number-of-samples=", len(signal.raw_data)
                        scale_signal_data(signal)
                    except Exception:
                        print "Problems trying to collect signal", signalId
                else:
                    # Multiplexed signals are labelled 1..mux
                    subchannel_list = range(1,mux+1)
                    for subchannel_id in subchannel_list:
                        signalId = (card_id, channel_id, subchannel_id)
                        try:
                            signal = dbox_info.signalDict[signalId]
                            signal.dt = card_dt * len(subchannel_list)
                            signal.raw_data = voltage_demux[subchannel_id - 1]
                            signal.FS = FS
                            print "Signal", signalId, "dt_sample=", signal.dt, \
                                  "number-of-samples=", len(signal.raw_data)
                            scale_signal_data(signal)
                        except Exception:
                            print "Problems trying to collect signal", signalId
                    # end of loop body for multiplexed signals for one card
                # end of loop body for collecting data for one channel
            pd.progress += 1
            # end of loop body for collecting data for one card
        pd.close()
        not_busy()
        dataLoadedButNotSaved = 1
    # We are done.
    dbox_info.release_lock()
    print "Begin refreshing plot page."
    refresh_plot_page()
    print "Done refreshing plot page."
    # Rainer wants the plot page to be brought up straight after collecting the data.
    global nb
    nb.selectpage("Plotted Data")
    return
    
def import_shot_dialog():
    "Load previously recorded data from file."
    # load run description
    # load signals
    # load data
    # put this data into the GUI widgets
    global root, dbox_info
    global signalConfigTextWidget
    global shotIdEntryWidget, dataDirEntryWidget
    shotDir = tkFileDialog.askdirectory(title="Choose a shot directory")
    if shotDir:
        busy()
        print "Import shot description and data from %s" % shotDir
        (dataDir, shotId) = os.path.split(shotDir)
        print "shotId=", shotId, "dirName=", shotDir
        e = dataDirEntryWidget
        e.delete(0, Tkinter.END)
        e.insert(0, dataDir)
        dbox_info.dataDir = dataDir
        e = shotIdEntryWidget
        e.delete(0, Tkinter.END)
        e.insert(0, shotId)
        dbox_info.shotId = shotId
        fileName = os.path.join(shotDir, shotId + ".txt")
        if os.path.exists(fileName):
            print "Import run description from", fileName
            read_run_description_dialog(fileName)
        fileName = os.path.join(shotDir, shotId + ".config")
        if os.path.exists(fileName):
            print "Import signal configuration from", fileName
            read_config_file_dialog(fileName)
        else:
            noConfigDialog = Pmw.MessageDialog(root, title="Importing data.",
                               defaultbutton=0, buttons=("OK",),
                               message_text="There is no signal config file.")
            result = noConfigDialog.show()
            not_busy()
            root.update()
            wait_a_while(200) # seem to need a bit of time to clear dialog
            return
        pd = myProgressDialog(maxValue=len(dbox_info.signalDict.keys()),
                              title="dbox_view",
                              message="Importing data from files.")
        for signalId in dbox_info.signalDict.keys():
            root.update()
            signal = dbox_info.signalDict[signalId]
            extn = str(signalId[0]) + str(signalId[1]) + str(signalId[2])
            fileName = shotId + "A." + extn + ".gz"
            print "Reading data file for signal", signalId, "fileName", fileName
            try:
                fdata = gzip.open(os.path.join(shotDir, fileName), "rb")
                read_data_from_file(dbox_info, signal, fdata, extn)
                fdata.close
                reconstruct_raw_voltages(signal)
            except:
                print "Failed to read data file:", fileName
            pd.progress += 1
        pd.close()
        not_busy()
        print "Begin refreshing plot page."
        refresh_plot_page()
        print "Done refreshing plot page."
    print "Finished importing data."
    dataLoadedButNotSaved = 0  # memory copy is same as disc copy
    global nb
    nb.selectpage("Plotted Data")
    return

def apply_current_scales_to_data():
    """
    Assuming that the on-screen signal configuration has been updated,
    rescale the in-memory raw data to get new scaled data.
    """
    global dbox_info, dataLoadedButNotSaved, signalConfigTextWidget
    dbox_info.signalConfigText = signalConfigTextWidget.get('1.0',Tkinter.END)
    rescan_signal_config_for_new_scales(dbox_info)
    for signalId in dbox_info.signalDict.keys():
        signal = dbox_info.signalDict[signalId]
        scale_signal_data(signal)
    dataLoadedButNotSaved = 1
    return

def save_data_control_s(event):
    """
    User pressed Control-s.
    """
    print "Activating save_data_dialog procedure."
    save_data_dialog()
    return

def save_data_dialog():
    """
    Write out all of the in-memory data,
    making sure that we have the latest information from
    the GUI display widgets.
    """
    global root, dbox_info, dataLoadedButNotSaved
    global runDescTextWidget, signalConfigTextWidget
    setShotIdFromEntryWidget(None)
    setDataDirFromEntryWidget(None)
    setBackupDataDirFromEntryWidget(None)
    dbox_info.runDescription = runDescTextWidget.get('1.0',Tkinter.END)
    dbox_info.signalConfigText = signalConfigTextWidget.get('1.0',Tkinter.END)
    targetDir = os.path.join(dbox_info.dataDir, dbox_info.shotId)
    if os.path.exists(targetDir):
        dialog = Pmw.MessageDialog(root, title="Save Data",
                     defaultbutton=0, buttons=("Proceed", "Cancel"),
                     message_text="You are about to overwrite existing data files.")
        result = dialog.activate()
        root.update()
        wait_a_while(200) # seem to need a bit of time to clear dialog
    else:
        result = "Proceed"
    if result == "Proceed":
        busy()
        try:
            do_save_all_data(dbox_info)
            dataLoadedButNotSaved = 0
            dialog_message_text = "Finished saving data."
        except:
            dialog_message_text = "Failed to save data."
        not_busy()
        print dialog_message_text
        dialog = Pmw.MessageDialog(root, title="Saving data.",
                     defaultbutton=0, buttons=("OK",),
                     message_text=dialog_message_text)
        result = dialog.show()
        root.update()
        wait_a_while(200) # seem to need a bit of time to clear dialog
    return


def export_csv_dialog():
    """
    Write out all of the in-memory data,
    making sure that we have the latest information from
    the GUI display widgets.
    """
    global root, dbox_info
    global runDescTextWidget, signalConfigTextWidget
    setShotIdFromEntryWidget(None)
    setDataDirFromEntryWidget(None)
    setBackupDataDirFromEntryWidget(None)
    dbox_info.runDescription = runDescTextWidget.get('1.0',Tkinter.END)
    dbox_info.signalConfigText = signalConfigTextWidget.get('1.0',Tkinter.END)
    targetDir = os.path.join(dbox_info.dataDir, dbox_info.shotId)
    busy()
    do_export_csv_data(dbox_info)
    not_busy()
    dialog = Pmw.MessageDialog(root, title="Export data to CSV.",
                               defaultbutton=0, buttons=("OK",),
                               message_text="Finished exporting data.")
    result = dialog.show()
    root.update()
    wait_a_while(200) # seem to need a bit of time to clear dialog
    return

#------------------------------------------------------------------------

def setShotIdFromEntryWidget(event):
    global dbox_info, shotIdEntryWidget
    dbox_info.shotId = shotIdEntryWidget.get()
    print "Shot identity now set to ", dbox_info.shotId
    return

def setDataDirFromEntryWidget(event):
    global dbox_info, dataDirEntryWidget
    dbox_info.dataDir = dataDirEntryWidget.get()
    print "Data directory now set to ", dbox_info.dataDir
    return

def setBackupDataDirFromEntryWidget(event):
    global dbox_info, backupDataDirEntryWidget
    dbox_info.backupDataDir = backupDataDirEntryWidget.get()
    print "Backup data directory now set to ", dbox_info.backupDataDir
    return

def layout_run_description_page(parent):
    "Layout the display of the data location and shot identity."
    global dbox_info, shotIdEntryWidget, runDescTextWidget
    global backupDataDirEntryWidget, dataDirEntryWidget
    f = Tkinter.Frame(parent, borderwidth=1)
    #
    row = 0
    l = Tkinter.Label(f, text="Shot Identity:")
    l.grid(row=row, column=0, sticky='w')
    e = Tkinter.Entry(f, bg="white", relief='sunken', width=70)
    e.insert(0, dbox_info.shotId)
    e.grid(row=row, column=1, sticky='w')
    e.bind('<Return>', setShotIdFromEntryWidget)
    shotIdEntryWidget = e
    #
    row = 1
    l = Tkinter.Label(f, text="Data Directory:")
    l.grid(row=row, column=0, sticky='w')
    e = Tkinter.Entry(f, bg="white", relief='sunken', width=70)
    e.insert(0, dbox_info.dataDir)
    e.grid(row=row, column=1, sticky='w')
    e.bind('<Return>', setDataDirFromEntryWidget)
    dataDirEntryWidget = e
    #
    row = 2
    l = Tkinter.Label(f, text="Backup Directory:")
    l.grid(row=row, column=0, sticky='w')
    e = Tkinter.Entry(f, bg="white", relief='sunken', width=70)
    e.insert(0, dbox_info.backupDataDir)
    e.grid(row=row, column=1, sticky='w')
    e.bind('<Return>', setBackupDataDirFromEntryWidget)
    backupDataDirEntryWidget = e
    #
    # Display the run description in a scrollable text widget
    row = 3
    l = Tkinter.Label(f, text="Run Description:")
    l.grid(row=row, column=0, sticky='nw')
    twf = Tkinter.Frame(f)
    s = Tkinter.Scrollbar(twf)
    T = Tkinter.Text(twf, bg="white")
    s.pack(side='right', fill='both')
    T.pack(side='left', fill='both')
    s.config(command=T.yview)
    T.config(yscrollcommand=s.set)
    T.insert(Tkinter.END, dbox_info.runDescription)
    twf.grid(row=row, column=1)
    runDescTextWidget = T
    #
    f.pack(expand=1, fill='both')
    return

def layout_signal_config_page(parent):
    "Layout the display of the signal configuration file."
    global dbox_info, signalConfigTextWidget 
    f = Tkinter.Frame(parent, borderwidth=1)
    #
    # Display the run description in a scrollable text widget
    row = 0
    l = Tkinter.Label(f, text="Signals:")
    l.grid(row=row, column=0, sticky='nw')
    twf = Tkinter.Frame(f)
    s = Tkinter.Scrollbar(twf)
    T = Tkinter.Text(twf, bg="white")
    s.pack(side='right', fill='both')
    T.pack(side='left', fill='both')
    s.config(command=T.yview)
    T.config(yscrollcommand=s.set)
    T.insert(Tkinter.END, dbox_info.signalConfigText)
    twf.grid(row=row, column=1)
    signalConfigTextWidget = T
    #
    rescaleB = Tkinter.Button(f, text="Rescale in-memory data",
                              command=apply_current_scales_to_data)
    rescaleB.grid(row=1,column=0)
    f.pack(expand=1, fill='both')
    return

#------------------------------------------------------------------------

def layout_databox_status_page(parent):
    "Layout the display of the databox status information."
    global dbox_info
    global cardStatusEntryWidget
    global timebaseStatusEntryWidget
    global triggerStatusEntryWidget
    global stateOfPlayLabel
    #
    # Put the whole collection into a frame that doesn't expand
    # so that we don't end up with all the entries spread out.
    gStatus = Pmw.Group(parent, tag_text="Status")
    gStatus.pack()
    
    # Right side of status group will contain trigger units and timebase units
    fRight = Tkinter.Frame(gStatus.interior(), borderwidth=1)
    fRight.pack(side='right', expand=1, fill='both')
    #
    g1 = Pmw.Group(fRight, tag_text="Trigger Units")
    g1.pack(expand=1, fill='both')
    triggerStatusEntryWidget = {}
    for i in range(3):
        l = Tkinter.Label(g1.interior(), text=("%d:" % (i+1)) )
        l.grid(row=i, column=0, sticky='w')
        e = Tkinter.Entry(g1.interior(), relief='sunken', width=50)
        e.grid(row=i, column=1, sticky='w')
        triggerStatusEntryWidget[i+1] = e
    #
    g2 = Pmw.Group(fRight, tag_text="Timebase Units")
    g2.pack(expand=1, fill='both')
    timebaseStatusEntryWidget = {}
    for i in range(3):
        l = Tkinter.Label(g2.interior(), text=("%d:" % (i+1)) )
        l.grid(row=i, column=0, sticky='w')
        e = Tkinter.Entry(g2.interior(), relief='sunken', width=50)
        e.grid(row=i, column=1, sticky='w')
        timebaseStatusEntryWidget[i+1] = e
    #
    # Left side of status group will contain status of cards
    fLeft = Tkinter.Frame(gStatus.interior(), borderwidth=1)
    fLeft.pack(side='left', expand=1, fill='both')
    g3 = Pmw.Group(fLeft, tag_text="Cards")
    g3.pack(expand=1, fill='both')
    cardStatusEntryWidget = {}
    for i in range(7):
        l = Tkinter.Label(g3.interior(), text=("%d:" % (i+1)) )
        l.grid(row=i, column=0, sticky='w')
        e = Tkinter.Entry(g3.interior(), relief='sunken', width=50)
        e.grid(row=i, column=1, sticky='w')
        if dbox_service.card_is_present(i+1):
            e.configure(bg='gray')
        cardStatusEntryWidget[i+1] = e
    #
    gAction = Pmw.Group(parent, tag_text="Action")
    gAction.pack()
    refreshB = Tkinter.Button(gAction.interior(), text="Refresh status",
                              command=refresh_databox_status_page)
    armB = Tkinter.Button(gAction.interior(), text="Arm databox",
                          command=arm_databox_dialog)
    triggerB = Tkinter.Button(gAction.interior(), text="Send trigger",
                           command=trigger_databox_dialog)
    collectB = Tkinter.Button(gAction.interior(), text="Collect data",
                           command=do_collect_data)
    refreshB.pack(side='left')
    armB.pack(side='left')
    triggerB.pack(side='left')
    collectB.pack(side='left')
    if not dbox_info.databoxIsPresent:
        refreshB.configure(state="disabled")
        armB.configure(state="disabled")
        triggerB.configure(state="disabled")
        collectB.configure(state="disabled")
    #
    return

def refresh_databox_status_page(rescheduleCheckWhenDone=0):
    """
    Do status checking without upsetting any other activities
    of the databox.
    """
    global dbox_info
    global cardStatusEntryWidget
    global timebaseStatusEntryWidget
    global triggerStatusEntryWidget
    global stateOfPlayLabel
    #
    if not dbox_info.databoxIsPresent:
        print "Cannot refresh databox status page if the databox is not present."
        return
    
    if dbox_info.lock_is_free():
        dbox_info.wait_to_acquire_lock()
        for i in range(1,4):
            triggerStatusEntryWidget[i].delete(0, Tkinter.END)
        for i in range(1,4):
            timebaseStatusEntryWidget[i].delete(0, Tkinter.END)
        for i in range(1,8):
            cardStatusEntryWidget[i].delete(0, Tkinter.END)
        #
        for i in range(1,4):
            trs = dbox_service.get_trigger_unit_settings(i)
            display_data = "slope=%s  coupling=%s  threshold=%6.0f" % \
                           (trs["slope"], trs["coupling"], trs["threshold"])
            triggerStatusEntryWidget[i].insert(Tkinter.END, display_data)
        #
        tbs = dbox_service.get_timebase_status()
        for i in range(1,4):
            display_data = "dt=%6dus  pretrigger=%4d  trigger_unit=%d  buffer_size=%d" \
                           % (tbs[i]['sample_period'], tbs[i]['pretrigger_samples'],
                              tbs[i]['trigger_unit'], tbs[i]['buffer_size'])
            timebaseStatusEntryWidget[i].insert(Tkinter.END, display_data)

        for card_id in dbox_info.cardDict.keys():
            card = dbox_info.cardDict[card_id]
            cs = dbox_service.get_card_status(card_id)
            display_text = "sampling=%d timebase=%d time_scale=%d first_word=%d" % \
                           (cs["sampling"], cs["time_base"],
                            cs["time_scale"], cs["first_word"])
            ew = cardStatusEntryWidget[card_id]
            ew.insert(Tkinter.END, display_text)
            if cs["sampling"] > 0:
                ew.configure(bg='yellow')
            else:
                ew.configure(bg='gray')
        # We are done updating the status page
        dbox_info.release_lock()
    else:
        # We have been locked out by some other function.
        pass
    #
    if rescheduleCheckWhenDone > 0:
        root.after(1000, refresh_databox_status_page,1)
    return

#------------------------------------------------------------------------

def layout_plot_page(parent):
    "Sets up a scrolled canvas for later drawing of the individual plots."
    global dbox_info, canvasWidget, plotsPerRowComboBox
    global canvasWidth, canvasHeight, configData
    # Put a drop-down list of options at the top of the page.
    f0 = Tkinter.Frame(parent, borderwidth=1)
    f0.pack(expand=1, fill="both")
    Tkinter.Label(f0, text="Plots per row:").pack(side="left")
    plotsPerRowComboBox = Pmw.ComboBox(f0, scrolledlist_items=("1","2","3","4"),
                                       entry_width=5,
                                       selectioncommand=update_plot_page)
    plotsPerRowComboBox.selectitem(configData['plots-per-row'])
    plotsPerRowComboBox.pack(side="left")
    f = Tkinter.Frame(parent, borderwidth=1)
    f.pack(expand=1, fill="both")
    # We want a large canvas that can be scrolled.
    # Since the Python Megawidgets version seems to be broken,
    # build our own with standard Tkinter widgets.
    scf = Tkinter.Frame(f)
    scf.pack(expand=1, fill="both")
    s = Tkinter.Scrollbar(scf)
    canvasWidth = int(configData['plot-window-width'])
    canvasHeight = int(configData['plot-window-height']) * 2
    # The screen view of the canvas with be shorter than its full height.
    c = Tkinter.Canvas(scf, bg="gray85", width=canvasWidth,
                       height=int(configData['plot-window-height']),
                       scrollregion=(0,0,canvasWidth,canvasHeight) )
    s.pack(side='right', fill='y')
    c.pack(side='left', fill='both')
    s.config(command=c.yview)
    c.config(yscrollcommand=s.set)
    canvasWidget = c
    return

def plot_a_signal(signal, xleft, ytop, xsize, ysize):
    "Plot the raw voltages."
    global canvasWidget
    label = "%s (%d %d %d)" % (signal.name, signal.card_id, signal.channel_id,
                               signal.subchannel_id)
    c = canvasWidget
    # outline box
    c.create_rectangle(xleft, ytop, xleft+xsize, ytop+ysize, tags="boxes")
    # signal label at top middle
    ymid = ytop + ysize
    xmid = xleft + xsize / 2
    c.create_text(xmid, ytop+10, text=label, tags="text")
    # axes
    x1 = xleft + int(0.10 * xsize)
    x2 = xleft + int(0.90 * xsize)
    xaxis_length = x2 - x1
    y0 = ytop + int(0.50 * ysize)
    yaxis_length = int(0.40 * ysize)
    c.create_line(x1, y0-yaxis_length, x1, y0+yaxis_length, tags="axes")
    c.create_line(x1, y0, x2, y0, tags="axes")
    c.create_text(x1, y0-yaxis_length, text=("%.1f" % signal.FS),
                  anchor='e', tags="axes")
    c.create_text(x1, y0+yaxis_length, text=("%.1f" % -signal.FS),
                  anchor='e', tags="axes")
    try:
        N = len(signal.raw_data)
    except:
        N = 0
    tmax = int(signal.dt * N)
    c.create_text(x2, y0, text=("%d" % tmax), anchor='n', tags="axes")
    if N > 0:
        # now, plot every nth data point
        nth = 1
        x = x1 + (xaxis_length * arange(0, N, nth)) / N
        y = y0 - (signal.raw_data[0::nth] / signal.FS * yaxis_length)
        coordList = transpose(array([x,y])).tolist()
        c.create_line(coordList, tags="datalines", fill="red")
    return

def update_plot_page(arg):
    "Callback function for the ComboBox."
    global dbox_info
    if len(dbox_info.signalDict.keys()) > 0:
        refresh_plot_page()
    return

def refresh_plot_page():
    "Draw all of the signals onto the plot page."
    global canvasWidget, plotsPerRowComboBox
    busy()
    c = canvasWidget
    c.delete("boxes")
    c.delete("text")
    c.delete("axes")
    c.delete("datalines")
    plotsPerRow = int(plotsPerRowComboBox.get())
    xsize = int(600 / plotsPerRow)
    ysize = 100
    x0 = 10
    y0 = 0
    plotsToLeft = 0
    ytop = y0
    signalIdList = dbox_info.signalDict.keys()
    signalIdList.sort()
    for i in signalIdList:
        signal = dbox_info.signalDict[i]
        xleft = x0 + xsize * plotsToLeft
        plot_a_signal(signal, xleft, ytop, xsize, ysize)
        plotsToLeft += 1
        if plotsToLeft >= plotsPerRow:
            ytop += ysize
            plotsToLeft = 0
    xmax = x0 + plotsPerRow * xsize
    ymax = ytop + ysize
    c.configure(scrollregion=(0,0,xmax,ymax))
    not_busy()
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

def initializeConfigData():
    # We will keep the program's configuration data as a dictionary of strings.
    global configData
    configData = {
        'plots-per-row': '1',
        'plot-window-width': '620',
        'plot-window-height': '400',
        'enable-threshold-reporting': '0',
        'data-dir': 'data',
        'backup-data-dir': 'data-backup',
        }
    iniFileName = "dbox_view.ini"
    if os.path.isfile(iniFileName):
        print "Read user configData from", iniFileName
        iniFile = open(iniFileName, "r")
        for line in iniFile:
            # Remove new-line, carriage-return and redundant white-space.
            str = line.replace('\n', '')
            str = line.replace('\r', '')
            while line.find('  ') >= 0:
                line = line.replace('  ', ' ')
            line = line.strip()
            if len(line) == 0: continue # nothing left of the line
            if line[0] == '#': continue # ignore comment line
            itemList = line.split()
            try:
                print "    setting \'%s\' to value \'%s\'" % (itemList[0],itemList[1])
                configData[itemList[0]] = itemList[1]
            except:
                pass
    return

#------------------------------------------------------------------------

if __name__ == '__main__':
    global root, wait_flag
    root = Tkinter.Tk()
    Pmw.initialise()
    initializeConfigData()
    dbox_info.dataDir = configData['data-dir']
    dbox_info.backupDataDir = configData['backup-data-dir']
    root.title("BCD DataBox Manager")
    busy()
    Pmw.aboutversion(versionString)
    Pmw.aboutcopyright("Centre for Hypersonics 2004,2007")
    Pmw.aboutcontact("Peter Jacobs\nemail: peterj@mech.uq.edu.au")
    if len(sys.argv) > 1:
        first_cmd_arg = sys.argv[1]
    else:
        first_cmd_arg = ""
    if first_cmd_arg.rfind("v") >= 0:
        # Only advertise if the first command line argument is something
        # like --version -version -v version v ...
        about = Pmw.AboutDialog(root, applicationname="BCD Databox Manager")
        root.lower(about)  # for advertising...
        root.after(3000, lambda : about.lower())
   
    setup_GUI_menus(root)
    nb = Pmw.NoteBook(root)
    pRun = nb.add("Run Information")
    pSignal = nb.add("Signal Config")
    pStatus = nb.add("Databox")
    pPlot = nb.add("Plotted Data")
    layout_run_description_page(pRun)
    layout_signal_config_page(pSignal)
    layout_databox_status_page(pStatus)
    layout_plot_page(pPlot)
    nb.pack(fill=Tkinter.BOTH, expand=1)
    nb.setnaturalsize()
    root.bind("<Control-s>", save_data_control_s) # single key-stroke to save data

    dataLoadedButNotSaved = 0
    if dbox_info.databoxIsPresent:
        dbox_service.enable_threshold_reporting(int(configData['enable-threshold-reporting']))
        refresh_databox_status_page(1)
    else:
        dialog = Pmw.MessageDialog(root, title="Databox Status Note.",
            defaultbutton=0, buttons=("OK",),
            message_text="Cannot see the databox.\nIs it plugged in and turned on?")
        result = dialog.show()
        root.update()
    not_busy()

    # Give control to the GUI elements.
    root.mainloop()
        
#-------------------------------------------------------------------------
