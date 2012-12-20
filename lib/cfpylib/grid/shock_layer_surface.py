"""
shock_layer_surface.py: A class for generating a 2D surface between two
                        curves ensuring that the cells adjacent to the 
                        east curve are normal to the surface.

Author: DFP

Version: 30-Oct-2012 : initial repository version, ported from my local version.
         19-Dec-2012 : addition of parametric surface function
"""


try:
    from libprep3 import *
except:
    print "Could not import libprep3."
    sys.exit()
from math import *
from cfpylib.nm.zero_solvers import bisection
from cfpylib.gasdyn.billig import *
from cfpylib.util.YvX import *
try:
    import matplotlib.pyplot as plt
except:
    print "Could not import matplotlib - plotting disabled."
    with_mpl = False
else:
    with_mpl = True
try:
    from scipy.optimize import *
except:
    print "Could not import scipy - shock fitting disabled."
    with_scipy = False
else:
    import scipy
    if float(scipy.version.version[:-2])<0.1: 
        print "WARNING: difficulties with the optimisation routines may be experienced with scipy versions less than 0.10.0"
    with_scipy = True
try:
    from numpy import *
except:
    print "Could not import numpy - shock fitting disabled."
    with_numpy = False
else:
    with_numpy = True

def RCF(a,b,beta):
    return RobertsClusterFunction(a, b, beta)

def HCF(dL0,dL1):
    return HypertanClusterFunction(dL0,dL1)

def BRCF(a,b,c,d,beta0,beta1,gamma):
    RCF_a = RCF(a,b,beta0)
    RCF_b = RCF(c,d,beta1)
    return DiscontinuousUnivariateFunction(gamma,RCF_a,RCF_b)

def TRCF( a, b, c, d, e, f, beta0, beta1, beta2, gamma0, gamma1 ):
    rcf0 = RCF(a,b,beta0)
    rcf1 = RCF(c,d,beta1)
    rcf2 = RCF(e,f,beta2)
    brcf01 = DiscontinuousUnivariateFunction(gamma0,rcf0,rcf1)
    return DiscontinuousUnivariateFunction(gamma1,brcf01,rcf2)

def BHCF( dL0, dL1, dL2, gamma ):
    hcf0 = HCF(dL0,dL1)
    hcf1 = HCF(dL1,dL2)
    return DiscontinuousUnivariateFunction(gamma,hcf0,hcf1)

def BHRCF( beta, dL0, dL1, gamma ):
    rcf = RCF(0,1,beta)
    hcf = HCF(dL0,dL1)
    return DiscontinuousUnivariateFunction(gamma,rcf,hcf)

def VCF(dL0, dL1, L, n):
    return ValliammaiFunction(dL0,dL1,L,n)

def BVCF(beta0,beta1,beta2,gamma,R,n):
    VCF_a = VCF(beta0*R,beta1*R,R*gamma,int(n*gamma))
    VCF_b = VCF(beta1*R,beta2*R,R*(1-gamma),int(n*(1-gamma)))
    return DiscontinuousUnivariateFunction(gamma,VCF_a,VCF_b)

class ShockLayerSurface:
    def __init__(self, east, west):
	self.east = east
	self.west = west

    def eval(self, r, s):
    	# eval the location on the east face
    	ep = self.east.eval(s)
    	# eval normal angle on the east face
    	dpdt = self.east.dpdt(s)
    	if dpdt.x == 0.0:
            theta = 0.0
        else:
      	    theta = -atan(dpdt.y / dpdt.x) + pi / 2.	
    	#  find the west point
    	if s == 0.0: t = 0
        elif s == 1.0: t = 1
    	elif theta != 0.0:
            def f(t):
    	        wp = self.west.eval(t)
    	        L = (wp.y - ep.y) / sin(theta) 
    	        lp = Vector3(ep.x - L * cos(theta), ep.y + L * sin(theta)) 
    	        return wp.x - lp.x
    	    t = bisection(f, 0.0, 1.0)
    	else:
            def f(t):
    	        wp = self.west.eval(t)
    	        L = (ep.x - wp.x) / cos(theta) 
    	        lp = Vector3(ep.x - L * cos(theta), ep.y + L * sin(theta)) 
    	        return wp.y - lp.y
    	    t = bisection(f, 0.0, 1.0)
    	wp = self.west.eval(t)
            # make the normal line
    	line = Line(ep, wp)
    	# eval the point
    	return line.eval(1 - r)

    def str():
        return "ShockLayerSurface"

def make_parametric_surface(bx_scale=1.0, by_scale=1.0, M_inf=1.0, R=1.0, axi=0, east=None, shock=None, f_s=1.1):  
    # model geometry
    if not east:
        # assume a quarter circle with radius R
        o = Vector3(0.0, 0.0)
        a = Vector3(-R, 0.0)
        b = Vector3(0.0, R)
        east = Arc(a, b, o)
    # else: the east boundary has been provided

    # inflow boundary
    if shock == None:
        x_limit = b.x
        inflow_nodes = []
        np = 100
        y_top = by_scale * y_from_x(-x_limit / bx_scale, M_inf, theta=0.0, axi=axi, R_nose=R)
        dy = y_top / (np - 1)
        for iy in range(np):
            y = dy * iy
            x = -bx_scale * x_from_y(y / by_scale, M_inf, theta=0.0, axi=axi, R_nose=R)
            inflow_nodes.append(Vector3(x, y))
        # create the inflow spline
        west = Spline(inflow_nodes)    
    else:  
        # create inflow nodes by extending the distances between the wall and the shock by some factor f_s
        inflow_nodes = []
        N = 100
        for i in range(N):
            s = float(i) / float(N - 1)
            ep = east.eval(s)
            dpdt = east.dpdt(s)
            gamma = -atan(dpdt.y / dpdt.x) + pi / 2.
            #  find the west point
            if s == 0.0: t = 0
            elif s == 1.0: t = 1
            else:
                def f(t):
                    wp = shock.eval(t)
                    L = (wp.y - ep.y) / sin(gamma) 
                    lp = Vector3(ep.x - L * cos(gamma), ep.y + L * sin(gamma)) 
                    return wp.x - lp.x
                t = bisection(f, 0.0, 1.0)
            # print "t = ", t
            wp = shock.eval(t)
            L = vabs(wp - ep)
            # print "i = %d, L = %e" % (i, L)
            sp = Vector3(ep.x - L * f_s * cos(gamma), ep.y + L * f_s * sin(gamma)) 
            inflow_nodes.append(sp)
        # create the inflow spline
        west = Spline(inflow_nodes)

    return ShockLayerSurface(east, west)

def extract_shock_coords( sol ):
    # find the shock (S=1) coordinates from an ExistingSolution
    shock_points = []
    shock_y_coords = []
    shock_x_coords = []
    for flow in sol.flow:
        for i in range(flow.ni):
            for j in range(flow.nj):
                for k in range(flow.nk):
                    if flow.data["S"][i][j][k]==1:
                        shock_points.append( Vector3( flow.data["pos.x"][i][j][k], flow.data["pos.y"][i][j][k], flow.data["pos.z"][i][j][k] ) ) 
                        shock_x_coords.append( flow.data["pos.x"][i][j][k] )
                        shock_y_coords.append( flow.data["pos.y"][i][j][k] )
    print "Found %d shock points" % len(shock_points)
    shock_x_coords, shock_y_coords = zip(*sorted(zip(shock_x_coords, shock_y_coords)))
    return shock_x_coords, shock_y_coords

def fit_billig2shock( sol, axi, M_inf, R, body=None ):
    # see if the required packages are available
    if not with_numpy or not with_scipy:
        print "Error: numpy and scipy are required for shock fitting"
        sys.exit()
    if not body:
        # assume a quarter circle with radius R
        o = Vector3(0.0, 0.0)
        a = Vector3(-R, 0.0)
        b = Vector3(0.0, R)
        body = Arc(a, b, o)
    shock_x_coords, shock_y_coords = extract_shock_coords( sol )
    # make a best fit to the shock location point cloud
    def shock_x_from_y( p, y ):
        print p
        _bx_scale, _by_scale, _M_inf, _Rn = p
        x = []
        for _y in y:
            _x = - _bx_scale * x_from_y(_y/_by_scale, _M_inf, theta=0.0, axi=axi, R_nose=_Rn)
            x.append( _x )
        result = array(x)
        return result

    def shock_y_from_x( p, x ):
        print p
        _bx_scale, _by_scale, _M_inf, _Rn = p
        y = []
        for _x in x:
            _y = _by_scale * y_from_x(-_x*_bx_scale, _M_inf, theta=0.0, axi=axi, R_nose=_Rn)
            y.append( _y )
        result = array(y)
        return result

    def residuals( p, x, y ):
        x_dash = shock_x_from_y(p,y)
        result = sqrt(sum((array(x) - x_dash)**2)/len(x))
        return result
    
    p0 = [ 1.0, 1.0, M_inf, R ]
    plsq = fmin_slsqp(residuals, p0, args=(shock_x_coords, shock_y_coords), bounds=[(1.0e-1,1e1),(1.0e-1,1e1),(1.0e0,1e2),(1.0e-4,1e0)], fprime=None)
    # p = p0
    p = plsq
    fit_x = []
    fit_y = copy(shock_y_coords)
    fit_y.sort()
    insert(fit_y,0,0.0)
    fit_x = shock_x_from_y(p,fit_y)
    
    # make a plot if matplotlib is available
    if with_mpl:
        plt.plot(fit_x,fit_y,"b-",label="fit")
        plt.plot(shock_x_coords,shock_y_coords,"g.",label="points")
        plt.grid()
        plt.legend()
        plt.show()  

    # sample along the billig curve to make a spline
    x_limit = body.eval(1).x
    np = 100
    y_top = p[1] * y_from_x(-x_limit/p[0], p[2], theta=0.0, axi=axi, R_nose=p[3])
    dy = y_top / ( np - 1 )
    shock_nodes = []
    for iy in range(np):
        y = dy * iy
        x = - p[0] * x_from_y(y/p[1], p[2], theta=0.0, axi=axi, R_nose=p[3])
        shock_nodes.append( Node(x,y) )
    shock_spline = Spline( shock_nodes )
    
    # find intersection of surface normal with the shock spline at t=1
    # eval the location on the east face
    ep = body.eval(1)
    # eval normal angle on the east face
    dpdt = body.dpdt(1)
    if dpdt.x == 0.0:
        theta = 0.0
    else:
        theta = -atan(dpdt.y / dpdt.x) + pi / 2.    
    #  find the west point
    if theta != 0.0:
        def f(t):
            wp = shock_spline.eval(t)
            L = (wp.y - ep.y) / sin(theta) 
            lp = Vector3(ep.x - L * cos(theta), ep.y + L * sin(theta)) 
            return wp.x - lp.x
        t = bisection(f, 0.0, 1.0)
    else:
        def f(t):
            wp = shock_spline.eval(t)
            L = (ep.x - wp.x) / cos(theta) 
            lp = Vector3(ep.x - L * cos(theta), ep.y + L * sin(theta)) 
            return wp.y - lp.y
        t = bisection(f, 0.0, 1.0)
    wp = shock_spline.eval(t)
    
    # sample along the billig curve to make a spline
    np = 100
    y_top = wp.y
    dy = y_top / ( np - 1 )
    shock_nodes = []
    for iy in range(np):
        y = dy * iy
        x = - p[0] * x_from_y(y/p[1], p[2], theta=0.0, axi=axi, R_nose=p[3])
        shock_nodes.append( Node(x,y) )
    shock_spline = Spline( shock_nodes )
    
    return shock_spline
