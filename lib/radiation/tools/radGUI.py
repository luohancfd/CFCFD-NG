#! /usr/bin/env python
#----------------------------------------------------------------------
# A GUI for the spectral radiation models library.
# Based on the simple.py wxPython example. 
#----------------------------------------------------------------------

import os
import pprint
import random
try:
    import wx
except:
    print "radGUI uses the wxPython library to create the graphical user-interface."
    print "Install wxPython and try again."
    sys.exit()

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

from radpy import *
from radmodel import *
from radiator_library import *

from gaspy import *
from cfpylib.gasdyn.cea2_gas import *

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
        
        # Set some default parameters
        self.p_box.SetValue(str(10.))
        self.T_box.SetValue(str(300.))
        self.U_box.SetValue(str(10000.))
        self.slab_width_box.SetValue(str(0.1))
        self.lambda_min_box.SetValue(str(self.psm.get_lambda_min()))
        self.lambda_max_box.SetValue(str(self.psm.get_lambda_max()))
        self.spectral_points_box.SetValue(str(self.psm.get_spectral_points()))
        self.apparatus_function_type = "Voigt"
        self.gaussian_box.SetValue(str(5.))
        self.lorentzian_box.SetValue(str(5.))
        self.problem_type = "shock"
        self.gas_mixture = "air"
        self.draw_figure()
        
    def create_default_radiation_input_file(self):
        gdata = GlobalRadData()
        
        # first the spectral data
        gdata.spectral_model = "photaura"
        gdata.lambda_min = 50.0
        gdata.lambda_max = 2000.0
        gdata.spectral_points = 195000

        # now the radiators
        radiators = available_radiators.keys()
        for rad_name in radiators:
            if rad_name.find("+")>=0 or rad_name.find("-")>=0: continue
	    rad = gdata.request_radiator(rad_name)
	    rad.isp = self.species.index(rad_name)
	    rad.iT = 0
	    rad.iTe = 3
	    rad.default_data()
	    if rad.type == "atomic_radiator":
                rad.line_set = rad.available_line_sets["all_lines"]
            if rad.type == "diatomic_radiator":
                rad.iTr = 1
                rad.iTv = 2
                
        # create the lua file
        gdata.write_LUA_file( "rad-model.lua", "radGUI.py" )
        
    def load_models(self):
        # species list
        self.species = [ 'Ar', 'Ar_plus', 'C', 'C_plus', 'C2', 'CN', 'CN_plus', 'CO2', 'CO', 'CO_plus',
                         'H', 'H_plus', 'H2', 'N', 'N_plus', 'N2', 'N2_plus', 'NCO', 'NO', 'NO_plus',
                         'O', 'O_plus', 'O2', 'O2_plus', 'Xe', 'Xe_minus', 'e_minus' ]

        # radiation model
        self.create_default_radiation_input_file()
        self.psm = Photaura("rad-model.lua")
        
        # cea2 gas interface
        reactants = make_reactants_dictionary( self.species )
        self.cea = Gas( reactants, with_ions=True, trace=1.0e-15 )
        
        # also initialise some other useful objects
        self.nsp = len(self.species)
        self.ntm = 4
        self.Q = Gas_data(self.nsp,self.ntm)
        self.X = CoeffSpectra(self.psm)
        self.I = SpectralIntensity(self.psm)
        
    def calculate_spectra(self):
        s0 = 0.0; T0 = 0.0
        s1 = self.slab_width; T1 = 0.0
        
        LOS = LOS_data(self.psm,1,T0,T1)
        divq = new_doublep()
        j_total = LOS.set_rad_point(0,self.Q,divq,s0+0.5*self.slab_width,self.slab_width)
        print "j_total = %e W/m3-sr" % j_total
        self.I.compute_spectral_distribution( self.psm )
        self.I.reset_intensity_vectors()
        I_total = LOS.integrate_LOS(self.I)
        print "I_total = %e W/m2-sr" % I_total
        
        self.X = LOS.get_rpoint_pointer(0).X_.clone()
        
    def prepare_spectra(self):
        if self.A:
            self.I.apply_apparatus_function( self.A ) 
        
        self.lambda_nm = []
        self.I_lambda = []
        
        for inu in range(self.I.nu.size()):
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
            size=(100,-1),
            style=wx.TE_PROCESS_ENTER)
        
        self.p_box = wx.TextCtrl(
            self.panel, 
            size=(100,-1),
            style=wx.TE_PROCESS_ENTER)
        
        self.T_box = wx.TextCtrl(
            self.panel, 
            size=(100,-1),
            style=wx.TE_PROCESS_ENTER)
            
        self.slab_width_box = wx.TextCtrl(
            self.panel, 
            size=(100,-1),
            style=wx.TE_PROCESS_ENTER)
            
        self.problem_choice = wx.Choice(
            self.panel, 
            size=(100, -1),
            choices = [ "shock", "pT" ])
        self.Bind(wx.EVT_CHOICE, self.select_problem_type, self.problem_choice)
        
        self.gas_choice = wx.Choice(
            self.panel, 
            size=(100, -1),
            choices = [ "air", "Mars", "nitrogen", "argon" ])
        self.Bind(wx.EVT_CHOICE, self.select_gas_mixture, self.gas_choice)
        
        self.lambda_min_box = wx.TextCtrl(
            self.panel, 
            size=(100,-1),
            style=wx.TE_PROCESS_ENTER)
            
        self.lambda_max_box = wx.TextCtrl(
            self.panel, 
            size=(100,-1),
            style=wx.TE_PROCESS_ENTER)
            
        self.spectral_points_box = wx.TextCtrl(
            self.panel, 
            size=(100,-1),
            style=wx.TE_PROCESS_ENTER)
        
        self.apparatus_choice = wx.Choice(
            self.panel, 
            size=(100, -1),
            choices = [ "Voigt", "SQRT_Voigt", "none" ])
        self.Bind(wx.EVT_CHOICE, self.select_apparatus_function, self.apparatus_choice)
        
        self.gaussian_box = wx.TextCtrl(
            self.panel, 
            size=(100,-1),
            style=wx.TE_PROCESS_ENTER)
            
        self.lorentzian_box = wx.TextCtrl(
            self.panel, 
            size=(100,-1),
            style=wx.TE_PROCESS_ENTER)
                       
        self.computebutton = wx.Button(self.panel, -1, "Compute")
        self.Bind(wx.EVT_BUTTON, self.on_compute_button, self.computebutton)
        
        self.writebutton = wx.Button(self.panel, -1, "Write to file")
        self.Bind(wx.EVT_BUTTON, self.on_write_button, self.writebutton)

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

        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL        
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        # self.ptype = wx.StaticText(self, -1, "Problem type: ")
        # self.hbox.Add(self.ptype, 0, border=3, flag=flags)
        self.hbox.Add(self.problem_choice, 0, border=3, flag=flags)
        self.hbox.Add(self.gas_choice, 0, border=3, flag=flags)
        
        self.hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox2.Add(self.p_box, 0, border=3, flag=flags)
        self.hbox2.Add(self.T_box, 0, border=3, flag=flags)
        self.hbox2.Add(self.U_box, 0, border=3, flag=flags)
        self.hbox2.AddSpacer(10)
        
        self.hbox3 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox3.Add(self.slab_width_box, 0, border=3, flag=flags)
        self.hbox3.AddSpacer(10)
        
        self.hbox4 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox4.Add(self.lambda_min_box, 0, border=3, flag=flags)
        self.hbox4.Add(self.lambda_max_box, 0, border=3, flag=flags)
        self.hbox4.Add(self.spectral_points_box, 0, border=3, flag=flags)
        self.hbox4.AddSpacer(10)
        
        self.hbox5 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox5.Add(self.apparatus_choice, 0, border=3, flag=flags)
        self.hbox5.Add(self.gaussian_box, 0, border=3, flag=flags)
        self.hbox5.Add(self.lorentzian_box, 0, border=3, flag=flags)
        self.hbox5.AddSpacer(10)
        
        self.hbox6 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox6.Add(self.computebutton, 0, border=3, flag=flags)
        self.hbox6.Add(self.writebutton, 0, border=3, flag=flags)
        self.hbox6.Add(self.cb_grid, 0, border=3, flag=flags)
        self.hbox6.AddSpacer(10)
        
        self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        self.vbox.Add(self.hbox2, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        self.vbox.Add(self.hbox3, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        self.vbox.Add(self.hbox4, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        self.vbox.Add(self.hbox5, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        self.vbox.Add(self.hbox6, 0, flag = wx.ALIGN_LEFT | wx.TOP)
                
        self.panel.SetSizer(self.vbox)
        self.vbox.Fit(self)
        
    def select_problem_type(self,event):
        print 'select_problem_type: %s\n' % event.GetString()
        
        if event.GetString() == 'pT':
            print 'pressure-temperature problem selected.'
            self.problem = "pT"
        elif event.GetString() == 'shock':
            print 'shock problem selected.'
            self.problem = "shock"
        else:
            print 'problem type: %s not recognised.' % event.GetString()
            sys.exit()
            
    def select_apparatus_function(self,event):
        print 'select_apparatus_function: %s\n' % event.GetString()
        
        if event.GetString() == 'Voigt':
            print 'Voigt apparatus function selected.'
            self.apparatus_function_type = "Voigt"
        elif event.GetString() == 'SQRT_Voigt':
            print 'Square-root Voigt apparatus function selected.'
            self.apparatus_function_type = "SQRT_Voigt"
        elif event.GetString() == 'none':
            print 'No apparatus function selected.'
            self.apparatus_function_type = "none"
        else:
            print 'apparatus function type: %s not recognised.' % event.GetString()
            sys.exit()
            
    def select_gas_mixture(self,event):
        print 'select_problem_type: %s' % event.GetString()
        
        if event.GetString() == 'air':
            print 'air selected. mass-fractions: 23.3% O2, 76.7% N2.'
            self.gas_mixture = "air"
        elif event.GetString() == 'Mars':
            print 'Mars gas selected. mass-fractions: 97% CO2, 3% N2.'
            self.gas_mixture = "Mars"
        elif event.GetString() == 'nitrogen':
            print 'nitrogen gas selected. mass-fractions: 100% N2.'
            self.gas_mixture = "nitrogen"
        elif event.GetString() == 'argon':
            print 'argon gas selected. mass-fractions: 100% Ar.'
            self.gas_mixture = "argon"
        else:
            print 'gas mixture: %s not recognised.' % event.GetString()
            sys.exit()

    def create_status_bar(self):
        self.statusbar = self.CreateStatusBar()
        
    def calculate_gas_state(self):      
        # compute equilibrium post-shock gas state
        self.cea.set_pT(self.p,self.T)
        if self.problem_type=="shock":
            self.cea.shock_process( self.U_inf )
        
        # fill out the full gas-state with EOS call
        for iT in range(self.ntm):
            self.Q.T[iT] = self.cea.T
        self.Q.rho = self.cea.rho
        for isp,sp in enumerate(self.species):
            self.Q.massf[isp] = get_species_composition(sp,self.cea.species)
        self.Q.p = self.cea.p
        if "e_minus" in self.species:
            ie = self.species.index("e_minus")
            self.Q.p_e = self.Q.massf[ie] * self.Q.rho / RC_m_SI * RC_k_SI * self.Q.T[self.ntm-1]
        self.Q.print_values(False)
        
    def read_parameter_boxes(self):
        """ Reads in the user-defined parameters from the GUI boxes
        """     
        # pressure
        str = self.p_box.GetValue()
        tks = str.split()
        if len(tks)!=1:
           print "Pressure: %s not understood." % str
           return
        self.p = float(tks[0])
        
        # temperature
        str = self.T_box.GetValue()
        tks = str.split()
        if len(tks)!=1:
            print "Temperature: %s not understood." % str
            return
        self.T = float(tks[0])
        
        # velocity
        if self.problem_type=="shock":
            str = self.U_box.GetValue()
            tks = str.split()
            if len(tks)!=1:
                print "Velocity: %s not understood." % str
                return
            self.U_inf = float(tks[0])
        
        # composition
        for isp in range(len(self.cea.reactants)):
            self.cea.reactants[isp] = 0.0
        if self.gas_mixture=="air":
            self.cea.reactants["N2"] = 0.767
            self.cea.reactants["O2"] = 0.233
        elif self.gas_mixture=="Mars":
            self.cea.reactants["CO2"] = 0.97
            self.cea.reactants["O2"] = 0.03
        elif self.gas_mixture=="nitrogen":
            self.cea.reactants["N2"] = 1.0
        elif self.gas_mixture=="argon":
            self.cea.reactants["Ar"] = 1.0
        else:
            print "Gas mixture: %s not understood." % self.gas_mixture
            return
        
        # read and set the spectral parameters
        # lambda_min
        str = self.lambda_min_box.GetValue()
        tks = str.split()
        if len(tks)!=1:
            print "Lambda min: %s not understood." % str
            return
        lambda_min = float(tks[0])
        # lambda_max
        str = self.lambda_max_box.GetValue()
        tks = str.split()
        if len(tks)!=1:
            print "Lambda max: %s not understood." % str
            return
        lambda_max = float(tks[0])
        # spectral_points
        str = self.spectral_points_box.GetValue()
        tks = str.split()
        if len(tks)!=1:
            print "spectral points: %s not understood." % str
            return
        spectral_points = int(tks[0])
        # set the new parameters in the spectral model (spectral blocks is assumed to be one)
        self.psm.new_spectral_params(lambda_min,lambda_max,spectral_points,1)
        print "lambda_min = ", self.psm.get_lambda_min()
        print "lambda_max = ", self.psm.get_lambda_max()
        print "spectral_points = ", self.psm.get_spectral_points()
    
        # read the slab width
        str = self.slab_width_box.GetValue()
        tks = str.split()
        if len(tks)!=1:
            print "Slab width: %s not understood." % str
            return
        self.slab_width = float(tks[0])
            
        # lorentzian width
        str = self.lorentzian_box.GetValue()
        tks = str.split()
        if len(tks)!=1:
            print "Lorentzian HWHM: %s not understood." % str
            return
        lorentzian_hwhm = float(tks[0])
        
        # gaussian width
        str = self.gaussian_box.GetValue()
        tks = str.split()
        if len(tks)!=1:
            print "Gaussian HWHM: %s not understood." % str
            return
        gaussian_hwhm = float(tks[0])
        
        # FIXME: nu_sample
        nu_sample = 10
    
        if lorentzian_hwhm + gaussian_hwhm > 0:
            if self.apparatus_function_type=="Voigt":
                self.A = Voigt( lorentzian_hwhm, gaussian_hwhm, nu_sample )
            elif self.apparatus_function_type=="SQRT_Voigt":
                self.A = SQRT_Voigt( lorentzian_hwhm, gaussian_hwhm, nu_sample )
            elif self.apparatus_function_type=="none":
                self.A = None
        else:
            self.A = None

    def draw_figure(self):
        """ Draws the figure with new data
        """
        self.read_parameter_boxes()
        self.calculate_gas_state()
        self.calculate_spectra()
        self.prepare_spectra()
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
        
    def on_write_button(self, event):
        print "Writing intensity spectra to file: intensity-spectra.txt"
        self.I.write_to_file("intensity-spectra.txt")
        print "Writing coefficient spectra to file: intensity-spectra.txt"
        self.X.write_to_file("coefficient-spectra.txt")
    
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
