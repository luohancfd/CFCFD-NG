#! /usr/bin/env python
#----------------------------------------------------------------------
# A GUI for the spectral radiation models library.
# Based on the simple.py wxPython example. 
#----------------------------------------------------------------------

from librad import *
from gaspy import *
from cfpylib.gasdyn.cea2_gas import *

import os
import pprint
import random
import wx

# The recommended way to use wx with mpl is with the WXAgg
# backend. 
#
import matplotlib
matplotlib.use('WXAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar
from matplotlib import rc
    
def make_cea2_dictionary( species_list ):
    nsp = len(species_list)
    reactants = dict()
    for sp in species_list:
        # replace names containing '_plus' with '+' 
        if ( sp.find("_plus")>=0 ): sp = sp[0:sp.find("_plus")] + "+"
        # replace names containing '_minus' with '-' 
        if ( sp.find("_minus")>=0 ): sp = sp[0:sp.find("_minus")] + "-"
        reactants.setdefault(sp,0.0)
    return reactants
    
def get_cea2_composition( sp, species_data ):
    # replace names containing '_plus' with '+' 
    if ( sp.find("_plus")>=0 ): sp = sp[0:sp.find("_plus")] + "+"
    # replace names containing '_minus' with '-' 
    if ( sp.find("_minus")>=0 ): sp = sp[0:sp.find("_minus")] + "-"
    if sp in species_data.keys():
        return species_data[sp]
    else:
        return 0.0


class BarsFrame(wx.Frame):
    """ The main frame of the application
    """
    title = 'Photaura GUI version 0.2'
    
    def __init__(self):
        wx.Frame.__init__(self, None, -1, self.title)
        
        # to enable use of symbols
        rc('text', usetex=True)
        
        self.load_models()
        
        self.create_menu()
        self.create_status_bar()
        self.create_main_panel()
        
        self.U_box.SetValue(str(10000.))
        self.p_box.SetValue(str(10.))
        self.T_box.SetValue(str(300.))
        self.draw_figure()
        
    def load_models(self):
        # radiation and gas-model
        self.psm = Photaura("rad-model.lua")
        self.gm = create_gas_model("gas-model.lua")
        
        # cea2 gas interface
        nsp = self.gm.get_number_of_species()
        self.species = []
        for isp in range(nsp):
            self.species.append( self.gm.species_name(isp) )
        reactants = make_cea2_dictionary( self.species )
        self.cea = Gas( reactants, with_ions=True, trace=1.0e-15 )
        
        # also initialise some other useful objects
        self.Q = Gas_data(self.gm)
        self.X = CoeffSpectra(self.psm)
        self.I = SpectralIntensity(self.psm)
        
    def calculate_spectra(self, ds=0.1, hwhm=5):        
        s0 = 0.0; T0 = 0.0
        s1 = ds; T1 = 0.0
        
        LOS = LOS_data(self.psm,1,T0,T1)
        divq = new_doublep()
        LOS.set_rad_point(0,self.Q,divq,s0+0.5*ds,ds)
        for inu in range(self.psm.get_spectral_points()):
            self.I.I_nu[inu] = 0.0
        I_total = LOS.integrate_LOS(self.I)
        print "I_total = %e W/m3-sr" % I_total
        
        self.X = LOS.get_rpoint_pointer(0).X_
        
    def prepare_spectra(self, hwhm):
        if hwhm>0.0:
            self.I.apply_apparatus_function( hwhm )
        
        self.lambda_nm = []
        self.I_lambda = []
        
        for inu in range(self.psm.get_spectral_points()):
            nu = self.I.nu[inu]
            self.lambda_nm.append( nu2lambda(nu) )
            self.I_lambda.append( self.I.I_nu[inu] * nu**2 / RC_c_SI *1.0e-10 )

    def create_menu(self):
        self.menubar = wx.MenuBar()
        
        menu_file = wx.Menu()
        m_expt = menu_file.Append(-1, "&Save plot\tCtrl-S", "Save plot to file")
        self.Bind(wx.EVT_MENU, self.on_save_plot, m_expt)
        menu_file.AppendSeparator()
        m_exit = menu_file.Append(-1, "E&xit\tCtrl-X", "Exit")
        self.Bind(wx.EVT_MENU, self.on_exit, m_exit)
        
        menu_help = wx.Menu()
        m_about = menu_help.Append(-1, "&About\tF1", "About the Photaura GUI")
        self.Bind(wx.EVT_MENU, self.on_about, m_about)
        
        self.menubar.Append(menu_file, "&File")
        self.menubar.Append(menu_help, "&Help")
        self.SetMenuBar(self.menubar)

    def create_main_panel(self):
        """ Creates the main panel with all the controls on it:
             * mpl canvas 
             * mpl navigation toolbar
             * Control panel for interaction
        """
        self.panel = wx.Panel(self)
        
        # Create the mpl Figure and FigCanvas objects. 
        # 5x4 inches, 100 dots-per-inch
        #
        self.dpi = 100
        self.fig = Figure((5.0, 4.0), dpi=self.dpi)
        self.canvas = FigCanvas(self.panel, -1, self.fig)
        
        # Since we have only one plot, we can use add_axes 
        # instead of add_subplot, but then the subplot
        # configuration tool in the navigation toolbar wouldn't
        # work.
        #
        
        self.axes = self.fig.add_subplot(111)
        
        # Bind the 'pick' event for clicking on one of the bars
        #
        self.canvas.mpl_connect('pick_event', self.on_pick)
        
        self.U_box = wx.TextCtrl(
            self.panel, 
            size=(200,-1),
            style=wx.TE_PROCESS_ENTER)
        
        self.p_box = wx.TextCtrl(
            self.panel, 
            size=(200,-1),
            style=wx.TE_PROCESS_ENTER)
        
        self.T_box = wx.TextCtrl(
            self.panel, 
            size=(200,-1),
            style=wx.TE_PROCESS_ENTER)
        
        self.computebutton = wx.Button(self.panel, -1, "Compute")
        self.Bind(wx.EVT_BUTTON, self.on_compute_button, self.computebutton)

        self.cb_grid = wx.CheckBox(self.panel, -1, 
            "Show Grid",
            style=wx.ALIGN_RIGHT)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, self.cb_grid)

        # Create the navigation toolbar, tied to the canvas
        #
        self.toolbar = NavigationToolbar(self.canvas)
        
        #
        # Layout with box sizers
        #
        
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        self.vbox.AddSpacer(10)
        
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        self.hbox.Add(self.p_box, 0, border=3, flag=flags)
        self.hbox.Add(self.T_box, 0, border=3, flag=flags)
        self.hbox.Add(self.U_box, 0, border=3, flag=flags)
        self.hbox.Add(self.computebutton, 0, border=3, flag=flags)
        self.hbox.Add(self.cb_grid, 0, border=3, flag=flags)
        self.hbox.AddSpacer(30)
        
        self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        
        self.panel.SetSizer(self.vbox)
        self.vbox.Fit(self)
    
    def create_status_bar(self):
        self.statusbar = self.CreateStatusBar()
        
    def calculate_gas_state(self):
        nsp = self.gm.get_number_of_species()
        ntm = self.gm.get_number_of_modes()
        
        # velocity
        str = self.U_box.GetValue()
        tks = str.split()
        if len(tks)!=1:
            print "Temperature: %s not understood." % str
            return
        U_inf = float(tks[0])
        
        # pressure
        str = self.p_box.GetValue()
        tks = str.split()
        if len(tks)!=1:
           print "Pressure: %s not understood." % str
           return
        p_inf = float(tks[0])
        
        # temperature
        str = self.T_box.GetValue()
        tks = str.split()
        if len(tks)!=1:
            print "Pressure: %s not understood." % str
            return
        T_inf = float(tks[0])
        
        # composition
        self.cea.reactants["N2"] = 0.767
        self.cea.reactants["O2"] = 0.233
        
        # compute equilibrium post-shock gas state
        self.cea.set_pT(p_inf,T_inf)
        self.cea.shock_process( U_inf )
        
        # fill out the full gas-state with EOS call
        for iT in range(ntm):
            self.Q.T[iT] = self.cea.T
        self.Q.rho = self.cea.rho
        for isp,sp in enumerate(self.species):
            self.Q.massf[isp] = get_cea2_composition(sp,self.cea.species)
        self.gm.eval_thermo_state_rhoT(self.Q)
        self.Q.print_values(False)

    def draw_figure(self):
        """ Draws the figure with new data
        """
        self.calculate_gas_state()
        self.calculate_spectra(ds=0.1)
        self.prepare_spectra(hwhm=5)
        self.redraw_figure()
        
    def redraw_figure(self):
        """ Redraws the figure using old data
        """
        # clear the axes and redraw the plot anew
        #
        self.axes.clear()
        
        self.axes.grid(self.cb_grid.IsChecked())
        self.axes.set_xlabel("Wavelength, $\lambda$ (nm)")
        self.axes.set_ylabel("Spectral radiance, I$_\lambda$ { W/cm$^2$-$\mu$m-sr")
        self.axes.set_yscale('log')
        self.plot_data = self.axes.plot(
            self.lambda_nm,
            self.I_lambda,
            linewidth=1,
            color=(0, 0, 0),
            )[0]
        
        self.canvas.draw()
    
    def on_cb_grid(self, event):
        self.redraw_figure()
    
    def on_compute_button(self, event):
        self.draw_figure()
    
    def on_pick(self, event):
        # The event received here is of the type
        # matplotlib.backend_bases.PickEvent
        #
        # It carries lots of information, of which we're using
        # only a small amount here.
        # 
        box_points = event.artist.get_bbox().get_points()
        msg = "You've clicked on a bar with coords:\n %s" % box_points
        
        dlg = wx.MessageDialog(
            self, 
            msg, 
            "Click!",
            wx.OK | wx.ICON_INFORMATION)

        dlg.ShowModal() 
        dlg.Destroy()       

    def on_save_plot(self, event):
        file_choices = "PNG (*.png)|*.png"
        
        dlg = wx.FileDialog(
            self, 
            message="Save plot as...",
            defaultDir=os.getcwd(),
            defaultFile="plot.png",
            wildcard=file_choices,
            style=wx.SAVE)
        
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.canvas.print_figure(path, dpi=self.dpi)
            self.flash_status_message("Saved to %s" % path)
        
    def on_exit(self, event):
        self.Destroy()
        
    def on_about(self, event):
        msg = """ A GUI for the photaura spectral radiation model
        
         * Use the matplotlib navigation bar
         * Enter a pressure, temperature and velocity into the text
         * box and click "Draw!"
         * Show or hide the grid
         * Save the plot to a file using the File menu
         * Click on a bar to receive an informative message
        """
        dlg = wx.MessageDialog(self, msg, "About", wx.OK)
        dlg.ShowModal()
        dlg.Destroy()
    
    def flash_status_message(self, msg, flash_len_ms=1500):
        self.statusbar.SetStatusText(msg)
        self.timeroff = wx.Timer(self)
        self.Bind(
            wx.EVT_TIMER, 
            self.on_flash_status_off, 
            self.timeroff)
        self.timeroff.Start(flash_len_ms, oneShot=True)
    
    def on_flash_status_off(self, event):
        self.statusbar.SetStatusText('')


if __name__ == '__main__':
    app = wx.PySimpleApp()
    app.frame = BarsFrame()
    app.frame.Show()
    app.MainLoop()
