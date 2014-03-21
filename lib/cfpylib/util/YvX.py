# filename: YvX.py
# description: A class describing some property Y against some property X
# ..Author: Dan Potter
# ..Version: 04-Mar-2013: simple edits by PJ, mainly format

import sys
from numpy import *
try:
    from scipy import interpolate
except:
    # print "Cannot import scipy.interpolate."
    # print "Interpolation will not be functional"
    with_scipy = False
else:
    with_scipy = True
from datetime import datetime
try:
    import pylab
except:
    # print "Cannot import pylab."
    # print "Plotting will not be functional"
    with_pylab = False
else:
    with_pylab = True
    
def eval_gaussian(x,gamma):
    return 1.0 / ( gamma * sqrt(2*pi) ) * exp( - x**2 / ( 2.0 * gamma**2 ) )
    
labels = []

class YvX:
    """
    Structure containing some property Y against some property X
    """
    def __init__( self, argA, argB=0, argC=1, argD=True ):
        if isinstance(argA,str):
            # initialise from input file
            self.create_from_file( argA, argB, argC, argD )
        else:
            # initialise from lists or arrays
            self.create_from_lists( argA, argB, argD )

    def create_from_lists(self, x_list, y_list, with_spline=True ):
        if len(x_list)!=len(y_list):
            print "YvX: x and y lists are not the same length"
            sys.exit()
        if x_list[0] > x_list[-1]:
            x_list.reverse()
            y_list.reverse()
        self.x_array = array(x_list)
        self.y_array = array(y_list)
        if with_scipy:
            if with_spline:
                self.spline_fit = interpolate.splrep( self.x_array, self.y_array )
            else:
                self.spline_fit = None
            self.linear_fit = interpolate.interp1d( self.x_array, self.y_array )
        print "Successfully created a YvX class from provided x and y lists of length: ", len(x_list)

    def create_from_file(self, infile_name, x_col=0, y_col=1, with_spline=True ):
        iFile = open( infile_name, "r" )
        x_list = []; y_list = []
        while 1:
             line = iFile.readline()
             tks = line.split()
             if len(tks)==0: break
             if line[0]=="#" or line[0]=='"': continue
             x_list.append(float(tks[x_col]))
             y_list.append(float(tks[y_col]))
        if x_list[0] > x_list[-1]:
            print "create_from_file: reversing data"
            x_list.reverse()
            y_list.reverse()
        self.x_array = array(x_list)
        self.y_array = array(y_list)
        if with_scipy:
            if with_spline:
                self.spline_fit = interpolate.splrep( self.x_array, self.y_array )
            else:
                self.spline_fit = None
            self.linear_fit = interpolate.interp1d( self.x_array, self.y_array )
        print "Successfully created a YvX class of size: %d from file: %s\n" % (len(x_list), infile_name)

    def y_from_x(self, x_val, use_spline=True, extrapolate_endpoints=True):
        # NOTE: x_val can be an array
        if not with_scipy:
            print "scipy is required for y_from_x to function"
            sys.exit()
        if isinstance( x_val, ndarray ):
            if use_spline:
                    y_vals = interpolate.splev( x_val, self.spline_fit, der=0 )
            else:
                    y_vals = self.linear_fit( x_val )
            for i,x in enumerate(x_val):
                if extrapolate_endpoints:
                    if x<self.x_array[0]:
                        y_vals[i] = self.y_array[0]
                    elif x_val>self.x_array[-1]:
                        y_vals[i] = self.y_array[-1]
                else:
                    if self.check_range( x )==False: y_vals[i] = 0.0
            return y_vals
        else:
            if extrapolate_endpoints:
                if x_val<self.x_array[0]:
                    return self.y_array[0]
                elif x_val>self.x_array[-1]:
                    return self.y_array[-1]
            else:
                if self.check_range( x_val )==False: return 0.0
            if use_spline:
                    return float(interpolate.splev( x_val, self.spline_fit, der=0 ))
            else:
                    return float(self.linear_fit( x_val ))

    def check_range(self, x_val):
        if x_val<self.x_array[0] or x_val>self.x_array[-1]: return False
        else: return True

    def write_to_file( self, outfile_name="YvX.txt", header="", include_integral=False ):
        print "Writing %d datapoints to file: %s" % ( len(self.x_array), outfile_name )
        oFile = open( outfile_name, 'w' )
        oFile.write("# filename: %s\n" % outfile_name )
        oFile.write("#  created: %s\n" % datetime.now() )
        oFile.write( header )
        for i in range(len(self.x_array)):
            if include_integral:
                oFile.write("%e \t %e \t %e\n" % ( self.x_array[i], self.y_array[i], self.integral[i] ) )
            else:
                oFile.write("%e \t %e\n" % ( self.x_array[i], self.y_array[i] ) )
        oFile.close()
        print "Done"

    def plot_data( self, title="", xlabel="", ylabel="", label="YvX_data", 
                   new_plot=True, show_plot=False, include_integral=False, 
                   rep='-', logscale_y=False, xrange=None, yrange=None ):
        if not with_pylab:
            print "Cannot plot data as pylab was not loaded."
            sys.exit()
        if new_plot: pylab.figure()
        if include_integral: pylab.subplot(211)
        pylab.title(title)
        pylab.xlabel(xlabel)
        pylab.ylabel(ylabel)
        if logscale_y: pylab.yscale('log')
        if xrange: pylab.xlim(xrange[0],xrange[1])
        if yrange: pylab.ylim(yrange[0],yrange[1])
        pylab.plot( self.x_array, self.y_array, rep )
        global labels
        labels.append( label )
        if include_integral:
            pylab.subplot(212)
            if xrange: pylab.xlim(xrange[0],xrange[1])
            self.integrate()
            pylab.plot( self.x_array, self.integral, rep )
            if show_plot:
                labels_ = []
                for label in labels:
                    labels_.append( label + " integral" ) 
                    pylab.legend( labels_, loc="best" )
        if show_plot:
            pylab.legend( labels, loc="best" ) 
            pylab.show()

    def plot_spline( self, title="", xlabel="", ylabel="", label="YvX_data", 
                     new_plot=True, show_plot=False, include_integral=False, 
                     rep='-', logscale_y=False, xrange=None, yrange=None  ):
        if not with_pylab:
            print "Cannot plot data as pylab was not loaded."
            sys.exit()
        if new_plot: pylab.figure()
        if include_integral: pylab.subplot(211)
        pylab.title(title)
        pylab.xlabel(xlabel)
        pylab.ylabel(ylabel)
        if xrange: pylab.xlim(xrange[0],xrange[1])
        if yrange: pylab.ylim(yrange[0],yrange[1])
        if logscale_y: pylab.yscale('log')
        x_spline = []; y_spline = []
        nx = 10*len(self.x_array)
        dx = float( self.x_array[-1] - self.x_array[0] ) / float( nx - 1 )
        for ix in range(nx):
            x_spline.append( self.x_array[0] + ix*dx )
            y_spline.append( self.y_from_x( self.x_array[0] + ix*dx ) )
        pylab.plot( x_spline, y_spline, rep )
        global labels
        labels.append( label )
        if include_integral:
            pylab.subplot(212)
            if xrange: pylab.xlim(xrange[0],xrange[1])
            self.integrate()
            pylab.plot( self.x_array, self.integral, rep )
            labels.append( label + " integral" )
        if show_plot:
            pylab.legend( labels, loc="best" ) 
            pylab.show()
            
    def show_plot(self, labels=None):
    	if labels!=None:
    		pylab.legend( labels, loc="best" ) 
    	pylab.show() 

    def limit_x_range( self, x_min, x_max ):
        # NOTE: x_array must contain ascending ordered values for this to work
        ix_min = searchsorted( self.x_array, x_min )
        ix_max = searchsorted( self.x_array, x_max )
        self.x_array = self.x_array[ix_min:ix_max]
        self.y_array = self.y_array[ix_min:ix_max]
        self.recompute_spline()
        print "Requested x_min = %f, new x_min = %f" % ( x_min, self.x_array[0] )
        print "Requested x_max = %f, new x_max = %f" % ( x_max, self.x_array[-1] )
        
    def integrate( self, x_min=None, x_max=None):
        if x_min==None: x_min = self.x_array[0]
        if x_max==None: x_max = self.x_array[-1]
        if x_min > x_max:
            tmp = x_min
            x_min = x_max
            x_max = tmp
        # simple trapezoidal integration
        self.integral = [0.0]; integral = 0.0
        for ix,x in enumerate(self.x_array):
            if x < x_min or x > x_max or ix==0: continue
            integral += 0.5 * ( self.y_array[ix] + self.y_array[ix-1] ) * abs( x - self.x_array[ix-1] )
            self.integral.append( integral )
        print "integral = ", integral
        return integral

    def recompute_spline(self, s=0):
        if not with_scipy:
            print "scipy is required for interpolation"
            sys.exit()
        if self.spline_fit:
            self.spline_fit = interpolate.splrep( self.x_array, self.y_array, s=s )
        if self.linear_fit:
            self.linear_fit = interpolate.interp1d( self.x_array, self.y_array )

    def apply_filter(self,name,params):
        if name=="gaussian":
            gamma = params["gamma"]
            x_tmp = []; y_tmp = []
            dx0 = self.x_array[1] - self.x_array[0]
            x = -gamma * 5
            while x < 0:
                x_tmp.append(x); x+= dx0
            for x in self.x_array:
                x_tmp.append(x)
            for x in x_tmp:
                y = 0.0
                for i,x_dash in enumerate(self.x_array):
                    delta_x = x_dash - x
                    if i==0: dx = 0.5 * ( self.x_array[1] - self.x_array[0] )
                    elif i==(len(self.x_array) - 1): dx = 0.5 * ( self.x_array[-1] - self.x_array[-2] )
                    else: dx = 0.5 * ( self.x_array[i+1] - self.x_array[i-1] )
                    f = self.y_array[i]
                    g = eval_gaussian(delta_x,gamma)
                    y += f*g*dx
                y_tmp.append(y)
            self.x_array = array(x_tmp)
            self.y_array = array(y_tmp)
        else:
            print "filter %d not implemented yet" % name
            sys.exit()
        
                    
            
            
        
