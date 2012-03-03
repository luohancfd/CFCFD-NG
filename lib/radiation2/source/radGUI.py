#! /usr/bin/env python
#----------------------------------------------------------------------
# A GUI for the spectral radiation models library.
# Based on the simple.py wxPython example. 
#----------------------------------------------------------------------

import wx

from librad2 import *
from gaspy import *

class MyFrame(wx.Frame):
    """
    This is MyFrame.  It just shows a few controls on a wxPanel,
    and has a simple menu.
    """
    def __init__(self, parent, title):
        wx.Frame.__init__(self, parent, -1, title,
                          pos=(150, 150), size=(350, 200))

        # Create the menubar
        menuBar = wx.MenuBar()

        # and a menu 
        menu = wx.Menu()

        # add an item to the menu, using \tKeyName automatically
        # creates an accelerator, the third param is some help text
        # that will show up in the statusbar
        menu.Append(wx.ID_EXIT, "E&xit\tAlt-X", "Exit this simple sample")

        # bind the menu event to an event handler
        self.Bind(wx.EVT_MENU, self.OnTimeToClose, id=wx.ID_EXIT)

        # and put the menu on the menubar
        menuBar.Append(menu, "&File")
        self.SetMenuBar(menuBar)

        self.CreateStatusBar()
        

        # Now create the Panel to put the other controls on.
        panel = wx.Panel(self)

        # and a few controls
        text = wx.StaticText(panel, -1, "Photaura GUI version 0.1")
        text.SetFont(wx.Font(14, wx.SWISS, wx.NORMAL, wx.BOLD))
        text.SetSize(text.GetBestSize())
        btn = wx.Button(panel, -1, "Close")
        create_btn = wx.Button(panel, -1, "Create spectral model")
        calc_btn = wx.Button(panel, -1, "Compute integrated emission")

        # bind the button events to handlers
        self.Bind(wx.EVT_BUTTON, self.OnTimeToClose, btn)
        self.Bind(wx.EVT_BUTTON, self.CreateModels, create_btn)
        self.Bind(wx.EVT_BUTTON, self.CalculateIntegratedEmission, calc_btn)

        # Use a sizer to layout the controls, stacked vertically and with
        # a 10 pixel border around each
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(text, 0, wx.ALL, 10)
        sizer.Add(btn, 0, wx.ALL, 10)
        sizer.Add(create_btn, 0, wx.ALL, 10)
        sizer.Add(calc_btn, 0, wx.ALL, 10)
        panel.SetSizer(sizer)
        panel.Layout()


    def OnTimeToClose(self, evt):
        """Event handler for the button click."""
        print "See ya later!"
        self.Close()

    def CreateModels(self, evt):
        """Event handler for the button click."""
        self.psm = Photaura("rad-model.lua")
        self.gm = create_gas_model("gas-model.lua")
        
    def CalculateIntegratedEmission(self, evt):
        """Event handler for the button click."""
        nsp = self.gm.get_number_of_species()
        ntm = self.gm.get_number_of_modes()
        
        Q = Gas_data(self.gm)
        
        for iT in range(ntm):
        	Q.T[iT] = 1.0e4
        Q.p = 1.0e5
        Q.massf[0] = 1.0
        Q.rho = Q.p / 286. / Q.T[0]
        j_total = self.psm.radiative_integrated_emission_for_gas_state( Q, True )
        print "j_total = %e W/m3-sr" % j_total


class MyApp(wx.App):
    def OnInit(self):
        frame = MyFrame(None, "Simple wxPython App")
        self.SetTopWindow(frame)

        print "Print statements go to this stdout window by default."

        frame.Show(True)
        return True
        
def main():
	app = MyApp(redirect=True)
	app.MainLoop()
	
if __name__=="__main__":
	main()
