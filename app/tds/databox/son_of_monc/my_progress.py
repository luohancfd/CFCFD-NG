## \file my_progress.py
## \brief A class to make the use of the PmwContribD.ProgressMeter
##        a bit more convenient.
## \author P.A. Jacobs
## \version 1.0, 15-Jan-2004
##

import Tkinter
import Pmw
from PmwContribD import ProgressMeter as PM

#-----------------------------------------------------------------------------

class myProgressDialog(object):
    """
    Put up a new Toplevel window with a progress meter and a message.
    Allow the progress meter to be updated by assigning to the progress property.
    """
    def __init__(self, maxValue=100,
                 title="Progress Meter",
                 message="Progress of task"):
        "Create a dialog window with progress meter in the middle of the screen."
        self.dialog = Pmw.Dialog(title="Progress Meter",
                                 buttons=("Dismiss",))
        self.dialog.withdraw()
        self.message = Tkinter.Label(self.dialog.interior(),
                                     text=message)
        self.message.pack()
        self.meter = PM.ProgressMeter(self.dialog.interior(),
                                      finishvalue=maxValue,
                                      progresscolor="lightblue")
        self._maxValue = maxValue
        self._progress = 0
        self.meter.updateProgress(self._progress)
        self.meter.pack()
        # Shift the dialog to be right in front of your view
        xleft = (self.dialog.winfo_screenwidth() - self.dialog.winfo_width())/2
        ytop = (self.dialog.winfo_screenheight() - self.dialog.winfo_height())/3
        self.dialog.geometry("+%d+%d" % (xleft,ytop))
        self.dialog.show()
        return

    def getProgressValue(self):
        return self._progress
    def setProgressValue(self, value):
        if value < self._maxValue:
            self._progress = value
        else:
            self._progress = self._maxValue
        try:
            self.meter.updateProgress(self._progress)
            self.dialog.update()
        except:
            print "Failed to increment progress meter."
        return self._progress
    progress = property(getProgressValue, setProgressValue)

    def finished(self):
        return self._progress == self._maxValue
    
    def close(self):
        "Clear away the toplevel window."
        try:
            self.dialog.withdraw()
            del self.dialog
        except:
            pass
        return

    def __del__(self):
        "Cleanup before garbage collection."
        self.close()
        return
    
#-----------------------------------------------------------------------------
if __name__ == "__main__":
    root = Tkinter.Tk()
    Pmw.initialise(root)

    def startProgressMeter():
        global pd
        pd = myProgressDialog(maxValue=3, title="Test", message="Save Data Files")
        return
    
    def incrProgressMeter():
        global pd
        pd.progress += 1
        return

    def disposeProgressMeter():
        global pd
        pd.close()
        return

    def periodicIncrement():
        global pd, root
        pd.progress += 1
        if not pd.finished():
            root.after(1000, periodicIncrement)
        return
    
    b2 = Tkinter.Button(root, text="Start non-modal Progress Meter",
                        command=startProgressMeter)
    b2.pack()
    b3 = Tkinter.Button(root, text="Increment Progress",
                        command=incrProgressMeter)
    b3.pack()
    b4 = Tkinter.Button(root, text="Dispose of Progress",
                        command=disposeProgressMeter)
    b4.pack()

    # Try out periodic updates with progress meter
    pd = myProgressDialog(maxValue=20, title="Test1", message="Count Up")
    periodicIncrement()
    
    root.mainloop()

