#
# fit.py (c) Stuart B. Wilkins 2008
#
# $Id: fit.py 220 2011-07-31 14:33:31Z tombeale $
# $HeadURL: https://pyspec.svn.sourceforge.net/svnroot/pyspec/trunk/pyspec/fit.py $
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Part of the "pyspec" package
#
"""FIT - Part of the pyspec routines for fitting of data

These routines / classes provide a method for fitting of data using mostly least
squares methods. There are two main methods here. The "fit" class provides a
class for fitting of data. The "fitdata" subroutine serves as a wrapper around
the "fit" class. Here the data is taken from the current selected figure.
Basically this means that you can fit any set of data which is on the current figure.

Examples:

f = fit(x = [xdata], y = [ydata], funcs = [func1, func2])
f.run()

result = f.result
errors = f.stdev

"""
from scipy import *
import pylab
import cmd
import time
import warnings

try:
    from scipy.odr import Model, Data, ODR
except ImportError:
    warnings.warn("NO scipy.odr package ... unable to use ODR regression module")
try:
    from scipy.optimize import leastsq
except ImportError:
    warnings.warn("ERROR: NO scipy.leastsq package ... unable to use leastsq regression module")
try:
    import pyspec.mpfit as mpfit
except ImportError:
    warnings.warn("ERROR: NO mpfit package ... unable to use MPFIT regression module")
try:
    import pyspec.pylevmar as levmar
except ImportError:
    warnings.warn("NO levmar package ... unable to use LEVMAR regression module")

__version__   = "$Revision: 220 $"
__author__    = "Stuart B. Wilkins <stuwilkins@mac.com>"
__date__      = "$LastChangedDate: 2011-08-01 00:33:31 +1000 (Mon, 01 Aug 2011) $"
__id__        = "$Id: fit.py 220 2011-07-31 14:33:31Z tombeale $"

class FitPlotCmd(cmd.Cmd):
    """Simple command processor for interactive fitting"""

    prompt = 'fit> '
    figure = None
    figuresep = None

    def __init__(self, fitclass, userplotfunc = None):
        """Initialize command line interface"""
        cmd.Cmd.__init__(self)
        self.fit = fitclass

        if userplotfunc is not None:
            self._plot = userplotfunc

    def do_op(self, line):
        """op {optimizer}
        Display or select the optimizer to use."""

        if not line:
            print "Optimizer = '%s'" % self.fit.optimizer
        else:
            self.fit.optimizer = line
            print "Optimizer set to = '%s'" % self.fit.optimizer

    def do_ch(self, line):
        """ch [num] [val]
        Change guess parameter"""

        argv = line.split()
        if len(argv) == 2:
            num = int(line.split()[0])
            if (num >= 0) and (num < len(self.fit.guess)):
                self.fit.guess[num] = float(line.split()[1])
                self.do_plot()
            else:
                print "Invalid parameter number %d\n" % num
        else:
            print "Syntax : ch [num] [newval]\n"

    def do_fx(self, line):
        """fx [num] [1/0]
        Fix parameter number [number] 1 = fixed 0 = floating"""

        argv = line.split()
        if len(argv) == 2:
            num = int(argv[0])
            if (num >= 0) and (num < len(self.fit.guess)):
                if self.fit.ifix is None:
                    self.fit.ifix = zeros(len(self.fit.guess))
                if int(argv[1]) == 1:
                    self.fit.ifix[num] = 1
                else:
                    self.fit.ifix[num] = 0
                self.do_plot()
            else:
                print "Invalid parameter number %d\n" % num

        else:

            print "Syntax : fx [num] [1/0]\n"

    def do_lm(self, line):
        """lm [num] [-ve] [+ve]
        Set parameter limits (X = no limit).
        Hint : using lm [num] X X cancels all limits"""

        argv = line.split()
        if len(argv) == 3:
            if self.fit.ilimited is None:
                self.ilimited = zeros((len(self.guess), 2), dtype = int)
                self.ilimits = zeros((len(self.guess), 2))

            p = int(argv[0])
            for n in range(2):
                if argv[n+1] == 'X':
                    self.fit.ilimited[p][n] = 0
                    self.fit.ilimits[p][n] = 0.0
                else:
                    self.fit.ilimited[p][n] = 1
                    self.fit.ilimits[p][n] = float(argv[n+1])

        else:
            print "Syntax : lm [num] [-ve] [+ve]"

    def do_cp(self, line):
        """cp
        Copy current fit values to guess values"""
        if self.fit.result is None:
            print "Error : There are no fitted values. Please fit data first.\n"
        else:
            self.fit.guess = self.fit.result.copy()
            print "\nCopied results of fit to new initial guess.\n"

    def do_plotsep(self, line):
        """plotsep [list]
        Plot the separate components of the fit.
        [list] contains a list of the components to plot."""

        if line == '':
            toplot = range(len(self.fit.funcs))
        else:
            toplot = eval(line)

        self.figuresep = pylab.figure()

        for l in reversed(self._hasplot('Guess', self.figuresep)):
            del self.figure.axes[0].lines[l]

        x, y = self.fit.evalfitfunc(nxpts = 200, p = self.fit.guess, mode = 'sep')
        for i in toplot:
            self._plot(x,y[:,i],'b-', label = 'Guess')

        if self.fit.result is not None:
            for l in reversed(self._hasplot('Fit', self.figuresep)):
                del self.figure.axes[0].lines[l]

            x, y = self.fit.evalfitfunc(nxpts = 200, p = self.fit.result, mode = 'sep')
            for i in toplot:
                self._plot(x,y[:,i],'g-', label = 'Fit')

        pylab.legend()

        for i in toplot:
            print "Plotting %d" % (i)


    def do_replot(self, plt = None):
        """replot
        Replot the graph on a new figure window"""
        self.do_plot(True)

    def do_plot(self, plt = None):
        """plot
        Plot the data and current fit"""

        if plt or (self.figure is None):
            self.figure = pylab.figure()
        else:
            pylab.figure(self.figure.number)

        if self._hasplot('Data', self.figure) == []:
            self._plot(self.fit.datax, self.fit.datay, 'r*', label = 'Data')

        for l in reversed(self._hasplot('Guess', self.figure)):
            del self.figure.axes[0].lines[l]

        x, y = self.fit.evalfitfunc(nxpts = 200, p = self.fit.guess)
        self._plot(x,y,'b-', label = 'Guess')

        if self.fit.result is not None:
            for l in reversed(self._hasplot('Fit', self.figure)):
                del self.figure.axes[0].lines[l]

            x, y = self.fit.evalfitfunc(nxpts = 200, p = self.fit.result)
            self._plot(x,y,'g-', label = 'Fit')

        pylab.legend()

    def _plot(self, *args, **kwargs):
        """Routine to overload to plot data"""
        pylab.plot(*args, **kwargs)

    def do_fit(self, line):
        """fit
        Run the optimizer and fit data."""
        self.fit.go(interactive = True)

        x, y = self.fit.evalfitfunc(nxpts = 200)

        for l in reversed(self._hasplot('Fit', self.figure)):
            del self.figure.axes[0].lines[l]

        self._plot(x,y,'g-', label = 'Fit')
        pylab.legend()

    def _hasplot(self, name, f = figure):
        if f is None:
            return []
        if f.axes == []:
            return []

        all = []
        for num in range(len(f.axes[0].lines)):
            line = f.axes[0].lines[num]
            if line.get_label() == name:
                all.append(num)
        return all


    def do_show(self, line):
        """show
        Show the current fit parameters"""

        k = 0
        sep = "----------------------------------------------------------------------------------------"
        yesno = ['no','yes']

        for i in range(len(self.fit.funcs)):
            pnames = self.fit.funcs[i](0,0,'params')
            print sep
            print "Func : ", self.fit.funcs[i](0, 0, 'name')
            print sep
            print " No : Parameter       :  Guess      :  Fit         : Fix :   -ve Limit  :   +ve Limit  :"
            print sep
            for j in range(len(pnames)):
                print "%3d : %-15s : % 6.4e :" % (k, pnames[j], self.fit.guess[k]),
                if self.fit.result is not None:
                    print " % 6.4e :" %  self.fit.result[k],
                else:
                    print "             :",
                if self.fit.ifix is not None:
                    print "%3s :" % yesno[int(self.fit.ifix[k])],
                else:
                    print "    :",

                if self.fit.ilimited is not None:
                    for n in range(2):
                        if self.fit.ilimited[k,n]:
                            print " % 6.4e :" % self.fit.ilimits[k,n],
                        else:
                            print "    NONE     :",
                else:
                    print "        :        :",

                print ""
                k += 1
            print sep

    def do_quit(self):
        """quit
        Quit the interactive fitting program."""
        return True

    def emptyline(self):
        pass # Do nothing

    def preloop(self):
        print "Interactive Fitting"
        print "-------------------"
        print
        print "Type \"help\" for list of commands"
        print

    def postloop(self):
        print

    def do_EOF(self, line):
        return True

def fitdata(funcs, ax = None, showguess = False, xdata=None, ydata=None, *args, **kwargs):
    """Force fitting of a function to graphed data

    Parameters
    ----------
    funcs    : list
       list of the fit functions
    ax       : matplotlib axis instance
       axis to fit (if none the current axis is used)

    All the additional *args and **kwargs are passed onto
    the fit class (see class pyspec.fit.fit for details).

    Returns
    -------

    pyspec.fit.fit instance with result.
    """

    if ax is None:
        l = pylab.gca().lines
    else:
        l = ax.lines

    if len(l) == 0:
        print "No data on graph"
        return None

    # If we have 3 lines .. error bar plot

    if len(l) == 3:
        a = array([l[0].get_ydata()[0], l[1].get_ydata()[0], l[2].get_ydata()[0]])
        a = a.argsort()
        n = a[1]
    else:
        n = 0
	
	if xdata==None:             
		xdata = l[n].get_xdata()
    if ydata==None:
	    ydata = l[n].get_ydata()

    f = fit(x = xdata, y = ydata, funcs = funcs, *args, **kwargs)
    f.go()

    # Now plot the results

    pylab.hold('on')
    if showguess:
        x, y = f.evalfitfunc(nxpts = 200, p = f.guess, x = xdata)
        pylab.plot(x, y, 'r-')

    x, y = f.evalfitfunc(nxpts = 200, x = xdata)
    pylab.plot(x, y, 'g-')

    return f

class fit:
    """General class to perform non-linear least squares fitting of data.

    This class serves as a wrapper around various fit methods, requiring a standard
    function to allow for easy fitting of data.

    Parameters:
    -----------
            x : ndarray
                x data
            y : ndarray
                y data
            e : ndarray
                y error
            funcs : list
                A list of functions to be fitted.
                These functions should be formed like the standard
                fitting functions in the pyspec.fitfuncs.
            guess : array or string
                Array of the initial guess, or string for guess type:
                    'auto' : Auto guess of parameters
            quiet : bool
                If True then don't print to the console any results.
            ifix : ndarray
                An array containing a '1' for fixed parameters
            xlimits : ndarray
                An (n x 2) array of the limits in x to fit.
            xlimitstype : string
                Either 'world' or 'index'
            optimiser : string      :
                Optimizer to use for fit: 'mpfit', 'leastsq', 'ODR'
            interactive : bool
                True/False launch interactive fitting mode
            r2min : float
                If r^2 is less than this value, drop into interactive mode

    After a sucsesful fit, the class will contain members for the following.
        {fit}.result
           result of fit
        {fit}.stdev
           errors on results (sigma)
        {fit}.fit_result
           an array of both results and stdev for convinience
        {fit}.covar
           covarience matrix of result
        {fit}.corr
           correlation matrix of result
        {fit}.r2
           R^2 value for fit.

    """
    def __init__(self, x = None, y = None, e = None,
                 funcs = None, guess = None,
                 quiet = False, optimizer = 'mpfit',
                 ifix = None,
                 ilimited = None, ilimits = None,
                 xlimits = None, xlimitstype = 'world',
                 interactive = False,
                 r2min = 0.0, debug = 0):
        self.debug = debug

        self.optimizer = optimizer
        self.quiet = quiet

        self.fitInteractive = interactive
        self.r2min = r2min

        self.guess = self._checkForArray(guess)
        self.result = None
        self.fit_result = None

        self.ifix = self._checkForArray(ifix)
        self.ilimited = self._checkForArray(ilimited)
        self.ilimits = self._checkForArray(ilimits)

        if xlimits != None:
            self.xlimits = self._checkForArray(xlimits)
        else:
            self.xlimits = None

        self.xlimitstype = xlimitstype

        self._niter = 1
        self._lastRunTime = None

        self.setFuncs(funcs)
        self.setData(x, y, e)

    def _checkForArray(self, a):
        if type(a) == list:
            return array(a)
        else:
            return a

    def setData(self, x, y, e):

        self.datax = x
        self.datay = y
        if e is None:
            self.datae = self.datay * 0 + 1.0
        else:
            self.datae = e

        if self.xlimits is not None:

            # For fitting we use the self._xdata so apply the limits

            if len(self.xlimits.shape) == 1:
                self.xlimits = array([self.xlimits])
            if len(self.xlimits.shape) != 2:
                raise Exception("Invalid array of x limits.")

            if self.xlimitstype == 'world':
                self._datax = array([])
                self._datay = array([])
                self._datae = array([])
                for limit in self.xlimits:
                    mask = (self.datax > limit[0])
                    mask = mask & (self.datax < limit[1])
                    self._datax = concatenate((self._datax, self.datax[mask]))
                    self._datay = concatenate((self._datay, self.datay[mask]))
                    self._datae = concatenate((self._datae, self.datae[mask]))

            elif self.xlimitstype == 'index':
                self._datax = array([])
                self._datay = array([])
                self._datae = array([])
                for limit in self.xlimits:
                    self._datax = concatenate((self._datax, self.datax[limit[0]:limit[1]]))
                    self._datay = concatenate((self._datay, self.datay[limit[0]:limit[1]]))
                    self._datae = concatenate((self._datae, self.datae[limit[0]:limit[1]]))
            else:
                raise Exception('Invalid x limits mode %s.' % self.xlimitstype)
        else:
            self._datax = self.datax
            self._datay = self.datay
            self._datae = self.datae

    def setFuncs(self, funcs):
        self.funcs = funcs

    def result(self):
        return self.fit_result

    def resultDict(self):
        return None

##
## The functions called by the fitting routines
##

    def _residuals(self, p):
        """Residuals function for scipy leastsq optimizer"""
        f = self._evalfunc(self._toFullParams(p), x = self._datax)
        return ravel(self._datay - f)

    def _residualsMPFIT(self, p, fjac = None):
        """Residuals function for MPFIT optimizer"""
        f = self._evalfunc(self._toFullParams(p), x = self._datax)
        rsid = (self._datay - f) / self._datae
        return 0, ravel(rsid)

    def _modelODR(self, p = None, x = None):
        """Model function for ODR"""
        return(self._evalfunc(self._toFullParams(p), x = self._datax))

    def _modelLEVMAR(self, estimate = None, measurement = None, data = None):
        """Model function for LEVMAR"""
        f = ravel(self._evalfunc(self._toFullParams(estimate), x = self._datax))
        return f

    def _toFullParams(self, p):
        """Return the full parameter list

        Sets the full parameter list, substituting the guess for
        the fixed parameters
        """
        g = self.guess.copy()
        g[self.ifix == 0] = array(p)
        return g

##
## Functions to evaluate the fitting functions
##

    def evalfitfunc(self, nxpts = None, p = None, x = None):
        """Evaluate the fit functions with the fesult of a fit

        Parameters
        ----------
        nxpts : int
            Number of x data points if using the range of the input data.
            If none then the x points of the dataset are used.
        p : ndarray
            Parameters of function. If None, use current fit result.
        x : ndarray
            Evaluate fit function at each point defined by the ndarray.

        Returns
        -------
        returns f(x) : ndarray

        """
        if x is None:
            x = self._datax

        if x.ndim == 1:
            if nxpts is not None:
                step = float( x.max() - x.min() ) / nxpts
                x = arange(x.min(), x.max(), step)

        f = self._evalfunc(x = x, p = p)
        return x, f

    def _evalfunc(self, p = None, x = None, mode = 'sum'):
        """Evaluate fitting function, not to be called by user"""
        
        if mode == 'sum':
            f = 0.0
        else:
            f = array([])

        ps = 0
        if x is None:
            x = self._datax

        if p is None:
            if self.result is not None:
                p = self.result
            else:
                p = self.guess

        for i in range(len(self.funcs)):
            pe = ps + len(self.funcs[i](0,0,'params'))
            if mode == 'sum':
                f += self.funcs[i](x , p[ps:pe])
            else:
                f = concatenate((f, self.funcs[i](x, p[ps:pe])))

            ps = pe

        if mode != 'sum':
            f = f.reshape(len(self.funcs),-1)
            f = f.T

        return f

    def _fitguess(self, mode = 'guess'):
        """Evaluate the fit functions to get the initial guess"""
        out = array([])
        for i in range(len(self.funcs)):
            gp = self.funcs[i](self._datax, self._datay, mode = mode)
            out = concatenate((out, gp), 1)
        return out

##
## Routines to run the different optimizers
##

    def _run_odr(self):
        """Run an ODR regression"""
        linear = Model(self._modelODR)
        mydata = Data(ravel(self._datax), ravel(self._datay), 1)
        myodr = ODR(mydata, linear,
                    beta0 = self._guess,
                    maxit = 10000)

        myoutput = myodr.run()

        self._result = myoutput.beta
        self._stdev = myoutput.sd_beta
        self._covar = myoutput.cov_beta
        self._odr = myoutput

    def _run_mpfit(self):
        """Run an MPFIT regression"""
        # Make up the PARINFO list

        parinfo = []
        for i in range(len(self._guess)):
            pdict = {'value' : self._guess[i].copy(),
                     'fixed' : 0,
                     'limited' : self._ilimited[i],
                     'limits' : self._ilimits[i]}
            parinfo.append(pdict)

        quiet = 1
        if self.debug:
            quiet = 0

        m = mpfit.mpfit(self._residualsMPFIT, self._guess,
                        parinfo = parinfo, quiet = quiet, debug = self.debug)

        self._result = m.params
        self._niter = m.niter

        if m.perror is not None:
            self._stdev = m.perror
        else:
            self._stdev = zeros(len(self._guess))

        if m.covar is not None:
            self._covar = m.covar
        else:
            self._covar = zeros((len(self._guess), len(self._guess)))

        self._mpfit = m

    def _run_leastsq(self):
        """Run a scipy leastsq regression"""
        plsq = leastsq(self._residuals, self._guess, Dfun = None,  full_output = 1, factor = 0.1)
        self._result = plsq[0]  # Make the stored guess the last value
        self._stdev = sqrt(diag(plsq[1].T))
        self._covar = plsq[1].T
        self._leastsq = plsq
        self._niter = plsq[2]['nfev']

    def _run_levmar(self):
        """Run a pylavmar regression"""
        #opts = array([1e-3, 1e-5, 1e-5, 1e-5, 1e-3])
        result, covar, iterations, run_info = levmar.ddif(self._modelLEVMAR,
                                                          self._guess,
                                                          ravel(self._datay), 10000)

        if result is not None:
            self._result = result
            self._covar = covar
            self._stdev = sqrt(diag(covar))
        else:
            self._result = self._guess.copy()
            self._stdev = zeros(len(self._result))
            self._covar = zeros((len(self._result), len(self._result)))

        self._niter = iterations
        self._levmar_run_info = run_info
        print "Finished LEVMAR"

##
## Run the optimization
##

    def run(self, interactive = False):
        """Start the fit

        Parameters
        ----------
        interactive : bool
           If True, start the fit in interactive mode.
        """
        self.go()

    def interactive(self):
        cli = FitPlotCmd(self)
        cli.cmdloop()

    def go(self, interactive = False):
        """Start the fit

        Parameters
        ----------
        interactive : bool
           If True, start the fit in interactive mode.
        """

        # Get the initial guess by calling the functions if we
        # have not supplied a guess.

        if not interactive:

            if self.guess is None:
                self.guess = self._fitguess()
            elif self.guess == 'graph':
                self.guess = self._fitguess('graphguess')
            elif self.guess == 'interactive':
                self.interactive()

            if self.fitInteractive:
                self.interactive()

    # Now we have the guess remove the fixed parameters from the fit

        if self.ifix is None:
            self.ifix = zeros(len(self.guess))
            self._guess = self.guess.copy()
        else:
            # Limit fitting parameters to only thoes which are nessesary
            self._guess = self.guess[nonzero(self.ifix == 0)].copy()

        if self.ilimited is None:
            self._ilimited = zeros((len(self._guess), 2), dtype = int)
            self._ilimits = zeros((len(self._guess), 2))
        else:
            self._ilimited = self.ilimited[nonzero(self.ifix == 0)].copy()
            self._ilimits = self.ilimits[nonzero(self.ifix == 0)].copy()

        # Here we will time the run for informtion
        t1 = time.time()
        self._niter = 1

        if self.optimizer == 'ODR':
            self._run_odr()
        elif self.optimizer == 'mpfit':
            self._run_mpfit()
        elif self.optimizer == 'leastsq':
            self._run_leastsq()
        elif self.optimizer == 'levmar':
            self._run_levmar()
        else:
            raise Exception("Unknown fitting optimizer '%s'" % self.optimizer)

        self._lastRunTime = time.time() - t1

        # If we have fixed parameters then return the full parameters list
        # and return the full list of errors

        self.result = self.guess.copy()
        self.result[nonzero(self.ifix == 0)] = self._result

        self.ilimited = zeros((len(self.guess), 2), dtype = int)
        self.ilimited[nonzero(self.ifix == 0)] = self._ilimited

        self.ilimits = zeros((len(self.guess), 2))
        self.ilimits[nonzero(self.ifix == 0)] = self._ilimits

        self.stdev = zeros(len(self.guess))
        self.stdev[nonzero(self.ifix == 0)] = self._stdev

        # Now deal with the covarience matrix

        if self.ifix.sum():
            self.covar = zeros((len(self.result), len(self.result)))
            ii = 0
            for i in range(len(self.stdev)):
                jj = 0
                for j in range(len(self.stdev)):
                    if self.ifix[i] or self.ifix[j]:
                        # We have a fixed paremeter
                        self.covar[i,j] = 0.0
                    else:
                        # Non fixed
                        self.covar[i,j] = self._covar[ii,jj]
                        jj+=1
                if not self.ifix[i]:
                    ii+=1
        else:
            self.covar = self._covar

        ##
        ## Calculate some useful
        ##

        # Correlation matrix
        self.corr = self.covar * 0.0
        for i in range(len(self.corr)):
            for j in range(len(self.corr)):
                self.corr[i,j] = self.covar[i,j] / sqrt(self.covar[i,i] * self.covar[j,j])

        # Chi^2
        self.chiSquared()

        # R^2
        ssReg = pow(self._evalfunc() - self._datay, 2).sum()
        ymean = self._datay.sum() / len(self._datay)
        ssTot = pow(ymean - self._datay, 2).sum()
        self.r2 = 1.0 - (ssReg / ssTot)

        # Make the dictionaryies

        self.resultsDict = self._makeDict(self.result)

        # Print out the result to console

        if self.quiet == False:
            print self.textResult()

        self.fit_result = vstack((self.result, self.stdev))

        # Here we can check for the r^2 and see if it is ok,
        # if not we can set the interactive bit and go again!

        if not interactive:
            if self.r2 < self.r2min:
                self.interactive()

        return self.fit_result

    def _makeDict(self, values):
        """Make and return a dictionary of the parameters or errors"""

        alld = []
        # First make up the
        for f in self.funcs:
            d = {}
            pnames = f(0,0,'params')
            for (name, val) in zip(pnames, values):
                d[name] = val
            alld.append(d)

        return alld

    def __str__(self):
        """Returns the fit summary"""
        return self.textResult()

    def textResult(self):
        """Return a string containing the text results of the fit"""
        p = ""

        if self.optimizer == 'ODR':
            p += "Fitted with ODRPACK\n"
            p += "--------------------\n"
        elif self.optimizer == 'mpfit':
            p += "Fitted with MPFIT\n"
            p +=  "------------------\n"
        elif self.optimizer == 'leastsq':
            p += "Fitted with 'scipy' leastsq\n"
            p += "----------------------------\n"
        elif self.optimizer == 'levmar':
            p += "Fitted with 'pylevmar' LEVMAR module\n"
            p += "------------------------------------\n"

        ## Put in here the x limits

        p += "Fit results to function(s)\n"
        for i in range(len(self.funcs)):
            p += "    %s\n" % self.funcs[i](0,0,'name')
        sep = "--------------------------------------------------------------------------------------------------\n"
        p  += sep
        p  += "Parameter       :  Init. Guess :  Value      :  Error     :  Fixed     :  -ve Limit :  +ve Limit :\n"
        p  += sep

        yesno = ['NO', 'YES']

        k = 0
        for i in range(len(self.funcs)):
            pnames = self.funcs[i](0,0,'params')
            for j in range(len(pnames)):
                p += "%-15s : % 6.4e : % 6.4e : % 6.4e :  %-6s    :" % (pnames[j], self.guess[k], self.result[k], self.stdev[k], yesno[int(self.ifix[k])])
                for n in range(2):
                    if self.ilimited[k, n]:
                        p += " %6.4e :" % self.ilimits[k,n]
                    else:
                        p += "  %-6s    :" % 'None'
                p += "\n"
                k += 1
            p += sep

        # Display the normalized covarience matrix

        allpnames = []
        for f in self.funcs:
            allpnames += f(0,0,'params')

        p += "\nCorrelation :\n"
        p += sep
        rsidsdev = self.residualSDev()
        for i in range(len(self.covar)):
            for j in range(i,len(self.covar)):
                if i != j:
                    if self.covar[i,j] != 0:
                        p += "%-15s <=> %-15s = %6.4e\n" % (allpnames[i], allpnames[j], self.corr[i,j])

        # Display the fit parameters

        p += "\nGoodness of fit:\n"
        p += sep
        p += "%-10s = %6.4e\n" % ("Chi**2", self.chi2)
        p += "%-10s = %f\n" % ("R^2", self.r2)

        p += sep

        # Now display optimizer specific code

        if self.optimizer == 'ODR':
            for reason in self._odr.stopreason:
                p += 'ODR Stop = %s\n' % reason
        if self.optimizer == 'mpfit':
            p += 'MPFIT Status = %s\n' % self._mpfit.statusNiceText[self._mpfit.status]
            p += 'MPFIT Warning = %s\n' % self._mpfit.errmsg
            p += 'MPFIT computed in %d iterations\n' % self._mpfit.niter
        if self.optimizer == 'levmar':
            p += "\nLEVMAR Output :\n"
            p += sep
            for x in self._levmar_run_info.iteritems():
                p += 'LEVMAR %-20s : %g\n' % x

        if self._lastRunTime is not None:
            p += "\nRuntime\n"
            p += sep
            p += "Fitted in %d iterations\n" % self._niter
            p += "Fit took %0.3f seconds (%0.3f sec per iteration)\n" % (self._lastRunTime,
                                                                       self._lastRunTime / self._niter)
            p += sep
        return p

    def chiSquared(self, norm = True, dist = 'poisson'):
        """ Return the chi-squared value for the fit

        Calculate the chi-squared value for the fit. This is defined as
           :math:`\chi^2 = \sum_N (x_{d,n} - x_{m,n})`

        Where d is the data and m is the model. The normalized chi-squared is given by
           :math:`\chi^2_{norm} = \chi^2 / M`

        where :math:`M = N - P`, where P is the number of parameters.

        If dist is 'poisson' then the data is divided by the model answer.
        i.e.
           :math:`\chi^2_{poisson} = \sum_N {(x_{d,n} - x_{m,i})} / {x_{m,i}}`

        Parameters
        ----------

        norm : bool
           If true calculate the normalized chi-squared. 
        dist : string
           The distribution to use, currently only 'poisson'

        Returns
        -------
        chi-squared value

        """

        N = len(self._datax)
        P = len(self.result)
        self.chi2 = pow(self._evalfunc() - self._datay, 2)
        if dist == 'poisson':
            self.chi2 = self.chi2 / self._evalfunc()

        self.chi2 = self.chi2.sum()
        if norm:
            self.chi2 = self.chi2 / (N - P)

        return self.chi2

    def residualSDev(self):
        """Calculate the sandard deviation of the residuals"""
        resid = self._evalfunc() - self._datay
        mean = resid.sum() / len(resid)
        stdev = sqrt(pow(resid - mean, 2).sum() / len(resid))
        self.resid_sdev = stdev
        return stdev

if __name__ == "__main__":
    import fitfuncs

    x = pylab.arange(-5, 5, 0.5)
    y = x * 10 - 3
    y = y + rand(len(y)) * 10

    f = fit(x = x, y = y, funcs = [fitfuncs.linear], guess = 'interactive')
    f.go()
