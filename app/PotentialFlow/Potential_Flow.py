## \Potential_Flow.py
#
""" 
Script to create Potential Flow Flow-Fields

Author: Ingo Jahn
Last modified: 11/08/2015 
"""

import numpy as np
import matplotlib.pyplot as plt

class PlotPotentialFlow:
    """
    class for PotentiaFlow-fields
    """
    def __init__(self):
        self.size()
    ##
    def size(self, x0=0.0, x1=1.0, y0=0.0, y1 = 1.0):
        self.x0 = x0
        self.x1 = x1
        self.y0 = y0
        self.y1 = y1
    ##
    def calc(self,A,n=100):
        self.A = A
        # create mesh
        xx = self.x0 + np.arange(n) * (self.x1-self.x0) / float(n-1)
        yy = self.y0 + np.arange(n) * (self.y1-self.y0) / float(n-1)
        self.X,self.Y = np.meshgrid(xx,yy)
        self.PSI = np.zeros([n,n])
        self.UU = np.zeros([n,n])
        self.VV = np.zeros([n,n])      
        # calculate stream functions and velocities
        for i in range(n):
            x = xx[i]
            for j in range(n):
                y = yy[j]
                psi = 0.
                U = 0.
                V = 0
                for it in range(len(A)):
                    psi = psi + A[it].evalP(x,y)
                    u,v = A[it].eval(x,y)
                    U = U + u
                    V = V + v
                self.PSI[j,i] = psi
                self.UU[j,i] = U
                self.VV[j,i] = V
    def evalP(self,x,y):
	# calculate Psi at a point
	PSI = 0.
        for it in range(len(self.A)):
	    PSI = PSI + self.A[it].evalP(x,y)
	return PSI
    ##
    def eval(self,x,y):
        # calculate U and V at a point
	U = 0.
        V = 0
	for it in range(len(self.A)):
	    u,v = self.A[it].eval(x,y)
            U = U + u
            V = V + v
	return U, V
    ##
    def evalPressure(self,x,y,rho):
        # calculate pressure reduction
	u,v = self.eval(x,y)
	Umag2 = u**2 + v**2
	dP = - 0.5 * rho * Umag2
	return dP
    ##
    def LinevalU(self,x0,y0,x1,y1,n=100,plot_flag=0):
        # calculate u-velocity at N points linearly spaced between point 0 and 1
        xx = x0 + np.arange(n) * (x1-x0) / float(n-1)
        yy = y0 + np.arange(n) * (y1-y0) / float(n-1)
        UU = np.zeros(n)
        for i in range(n):
            u,v = self.eval(xx[i],yy[i])
	    UU[i] = u
        if plot_flag == 1:
            plt.figure()
            plt.plot(UU)
	    plt.title('U-velocity')
	    plt.xlabel('Position along line')
	    plt.ylabel('U-velocity')
        return UU
    ##
    def LinevalV(self,x0,y0,x1,y1,n=100,plot_flag=0):
        # calculate u-velocity at N points linearly spaced between point 0 and 1
        xx = x0 + np.arange(n) * (x1-x0) / float(n-1)
        yy = y0 + np.arange(n) * (y1-y0) / float(n-1)
	VV = np.zeros(n)
	for i in range(n):
	    u,v = self.eval(xx[i],yy[i])
	    VV[i] = v
	if plot_flag == 1:
            plt.figure()
	    plt.plot(VV)
	    plt.title('V-velocity')
	    plt.xlabel('Position along line')
	    plt.ylabel('V-velocity')
	return VV			
    ##
    def LinevalPressure(self,x0,y0,x1,y1,rho,n=100,plot_flag=0):
        # calculate u-velocity at N points linearly spaced between point 0 and 1
        xx = x0 + np.arange(n) * (x1-x0) / float(n-1)
        yy = y0 + np.arange(n) * (y1-y0) / float(n-1)
	PP = np.zeros(n)
	for i in range(n):
	    u,v = self.eval(xx[i],yy[i])
	    PP[i] = - 0.5 * rho * (v**2 + u**2)
	if plot_flag == 1:
            plt.figure()
	    plt.plot(PP)
	    plt.title('Pressure')
	    plt.xlabel('Position along line')
	    plt.ylabel('Pressure')
	return PP		
    ##
    def plot_Streamlines(self):
        plt.figure()
        plt.streamplot(self.X,self.Y,self.UU,self.VV, density = 2, linewidth = 1, arrowsize=2, arrowstyle='->')
        plt.title('Streamlines')
        plt.xlabel('X-Position')
        plt.ylabel('Y-Position')
        plt.gca().set_aspect('equal')
        plt.gca().set_xlim([self.x0,self.x1])
        plt.gca().set_ylim([self.y0,self.y1])
    ##
    def plot_Streamlines_magU(self,magU_max = 100,levels=20):
        plt.figure()
        magU = (self.VV**2 + self.UU**2)**0.5
        magU[magU>magU_max] = magU_max 
        CS = plt.contourf(self.X, self.Y, magU,levels)
        plt.colorbar(CS)
        plt.streamplot(self.X,self.Y,self.UU,self.VV, density = 2, linewidth = 1, arrowsize=2, arrowstyle='->',color='k')
        plt.title('Streamlines and Velocity Magnitude')
        plt.xlabel('X-Position')
        plt.ylabel('Y-Position')
        plt.gca().set_aspect('equal')
        plt.gca().set_xlim([self.x0,self.x1])
        plt.gca().set_ylim([self.y0,self.y1])
    ##
    def plot_U(self,U_min = -100., U_max = 100., levels=20):
        U = self.UU
        U[U<U_min] = U_min
        U[U>U_max] = U_max  
        self.plot_cf(U,levels=levels,label="U-velocity") 
    ##
    def plot_V(self,V_min = -100., V_max = 100., levels=20):
        V = self.VV
        V[V<V_min] = V_min
        V[V>V_max] = V_max  
        self.plot_cf(V,levels=levels,label="V-velocity") 
    ##
    def plot_magU(self,magU_max = 100, levels=20):
        magU = (self.VV**2 + self.UU**2)**0.5
        magU[magU>magU_max] = magU_max 
        self.plot_cf((self.VV**2 + self.UU**2)**0.5,levels=levels,label="Velocity Magnitude") 
    ##
    def plot_cf(self,Z,levels=20,label="Label"):
        plt.figure()
        CS = plt.contourf(self.X, self.Y, Z,levels)
        # set graph details       
        plt.colorbar(CS)
        plt.title(label)
        plt.xlabel('X-Position')
        plt.ylabel('Y-Position')
        plt.legend
        plt.gca().set_aspect('equal')
        plt.gca().set_xlim([self.x0,self.x1])
        plt.gca().set_ylim([self.y0,self.y1])
    ## 
    def plot_P(self,P_inf = 0., rho=1.225, P_min=-100., P_max=100., levels=20):
        # limit pressure to 
        P = P_inf - 0.5 * rho * (self.VV**2 + self.UU**2)
        P[P<P_min] = P_min
        P[P>P_max] = P_max    
        self.plot_cf(P,levels=levels,label="Pressure")
    ##
    def plot_Cp(self, U_inf = 0., rho=1.225, Cp_min = -5., Cp_max = 5., levels=20):
        if float(U_inf) == 0.:
            print "For case with U_inf = 0., Cp becomes infinite everywhere"
        else:
            # Limit CP to account for localised high velocities
            Cp = (0.5* rho*U_inf**2 - 0.5* rho * (self.VV**2 + self.UU**2)) / (0.5 * rho * U_inf**2)
            Cp[Cp < Cp_min] = Cp_min
            Cp[Cp > Cp_max] = Cp_max
            self.plot_cf(Cp,levels=levels,label="Pressure Coefficient - Cp (Note limited to +/-5.)")
    ##
    def plot_Psi(self, levels=20):
        # Plot stream function Psi
        self.plot_cf(self.PSI,levels=levels,label="Streamfunction PSI (Care required with sources/sinks)") 


## Definition of classes used as Building Blocks 
class UniformFlow:
    """
    class that creates a uniform flow field for potential flow
    UniformFlow(u,v)
    u  -  x-component of velocity
    v  -  y-component of velocity
    """
    def __init__(self,u,v):
        self.u = u
        self.v = v
    ##
    def evalP(self,x,y):
        P = self.u*y - self.v*x
        return P        
    ##
    def eval(self,x,y):
        u = self.u
        v = self.v
        return u,v        

class Source:
    """
    class that creates a source for potential flow.
    Source(x0,y0,m)
    x0  -  x-position of Source
    y0  -  y-position of Source
    m   -  total flux generated by source (for sink set -ve)
    """
    def __init__(self,x0,y0,m):
        self.x0 = x0
        self.y0 = y0
        self.m = m
    ##
    def evalP(self,x,y):
        theta = np.arctan2(y-self.y0,x-self.x0)
        P = theta * self.m /(2*np.pi)
        return P    
    ##
    def eval(self,x,y):
        r = ( (x-self.x0)**2 + (y-self.y0)**2 )**0.5
        u = self.m / (2*np.pi) * (x - self.x0) / (r**2)
        v = self.m / (2*np.pi) * (y - self.y0) / (r**2)
        return u,v  
 
class Vortex:
    """
    class that creates an irrotational vortex for potential flow.
    Vortex(x0,y0,K)
    x0  -  x-position of Vortex core
    y0  -  y-position of Vortex core
    K   -  Strength of Vortex
    """
    def __init__(self,x0,y0,K):
        self.x0 = x0
        self.y0 = y0
        self.K = K
    ##
    def evalP(self,x,y):
        r = ( (x-self.x0)**2 + (y-self.y0)**2 )**0.5
        P = - self.K * np.log(r)
        return P  
    ##
    def eval(self,x,y):
        r = ( (x-self.x0)**2 + (y-self.y0)**2 )**0.5
        u = self.K / (2*np.pi) * (y - self.y0) / (r**2)
        v = - self.K / (2*np.pi) * (x - self.x0) / (r**2)
        return u,v  

class Doublet:
    """
    class that creates a doublet. 
    If combined with uniform flow of veloicty U_inf in the +x direction, this creates the flow around a cylidner.
    Doublet(x0,y0,a,U_inf)
    x0  -  x-position of Vortex core
    y0  -  y-position of Vortex core
    a   -  radius of cylinder generated if superimposed to Uniform Flow
    U-inf  -  Strength of uniform flow 
    """		
    def __init__(self,x0,y0,a,U_inf):
        self.x0 = x0
        self.y0 = y0
        self.a = a
        self.U_inf = U_inf
    def evalP(self,x,y):
        P = self.U_inf * (y-self.y0) * ( - self.a**2 / ((y-self.y0)**2 + (x-self.x0)**2) ) 
        # set to zero inside circle
        #if ((y-self.y0)**2 + (x-self.x0)**2) < self.a**2:
        #    P = np.nan
        return P
    ##
    def eval(self,x,y):
        u = self.U_inf * self.a**2 * - ((x-self.x0)**2 - (y-self.y0)**2) / ( ((x-self.x0)**2 + (y-self.y0)**2)**2 ) 
        v = self.U_inf * self.a**2 * -2. * (x-self.x0) * (y-self.y0) / ( ((x-self.x0)**2 + (y-self.y0)**2)**2 )
        # set to zero inside circle
        #if ((y-self.y0)**2 + (x-self.x0)**2) < self.a**2:
        #    u = np.nan
        #    v = np.nan
        return u,v
		
class User_Defined:
    """
    Special Userdefined building block (flow through 90 degree corner)
    User_Defined(x0,y0,A)
    x0  -  x-position of corner
    y0  -  y-position of corner
    A   -  Strength of Flow 
    """
    def __init__(self,x0,y0,A):
        self.x0 = x0
        self.y0 = y0
        self.A = A
    ##
    def evalP(self,x,y):
        P = self.A * (x-self.x0) * (y-self.y0)
        return P 
    ##
    def eval(self,x,y):
        u = self.A * (x-self.x0) 
        v = - self.A * (y-self.y0)
        return u,v 

class Name:
    """
    Template for user generated building blocks
    Name(x0,y0,Var1,Var2,Var3)
    x0  -  x-position
    y0  -  y-position
    Var1  -  Variable1
    Var2  -  Variable2
    Var3  -  Variable3
    """
    def __init__(self,x0,y0,Var1,Var2,Var3):
        self.x0 = x0
        self.y0 = y0
        self.Var1 = Var3
        self.Var2 = Var2
        self.Var3 = Var1 
    ##
    def evalP(self,x,y):
        ## Function for streamfunction goes here
        P = 0
        return P 
    ##
    def eval(self,x,y):
        ## Functions for u and v velocity go here. I.e. differentiate stream function with respect to x and y 
        u = 0 
        v = 0
        return u,v 
		
		

## Main section of the code, executed if running the file
if __name__ == "__main__":

    # List of Building Blocks 
    # Uniform Flows  
    A1 = UniformFlow(5.,0.)
    A2 = UniformFlow(0.,0.)
    # Vortices
    C1 = Vortex(0.0,0.0,5.)
    C2 = Vortex(0.0,-0.5,0.5)
    C3 = Vortex(0.0,0.0,10.0)
    # Sources
    D1 = Source(0.5,0., 5.)
    D2 = Source( 0.1,0.,-100.)
    D3 = Source(0.1,0.,0.24)
    # User Defined functions
    U1 = User_Defined(0.,0.,5)
    DD = Doublet(0.,0.,0.1,5.)


    # Initialise instance of Plotting Function
    T = PlotPotentialFlow()   	# create instance of the PotentialFlow-field class
    # Set dimensions of Plotting area
    T.size(-2.0, 2.0, -1.0 ,1.0)   #(x_min, x_max, y_min, y_max)

    # Evaluate PotentialFlow-field over a grid
    T.calc([A1,D1],n=100)   # ([List of elements], level of discretisation)

    # plot Data over flow-field area
    T.plot_Streamlines()    # create Streamline plot.
    #T.plot_magU(magU_max = 10., levels = 20)         # create plot of Velocity magnitude    
    #T.plot_Streamlines_magU(magU_max = 10.,levels = 20) # crete plot of Streamlines + velocity magnitude
    #T.plot_U(U_min = -10., U_max = 10., levels =20)            # create plot of U-velocity contours 
    #T.plot_V(V_min = -10., V_max = 10., levels =20)            # create plot of V-velocity contours
    #T.plot_P(P_inf = 0., rho=1.225, P_min=-100., P_max=200.)   # create plot of Pressure contours
    #T.plot_Cp(U_inf = 5., rho=1.225, Cp_min = -5., Cp_max = 5., levels=20) # create plot of Pressure coefficient contours
    # T.plot_Psi(levels = 20)                                    	# cretae plot of Psi contours 

    # extract data at points
    # print 'Psi = ', T.evalP(0.,0.)
    # print '(u, v) = ', T.eval(0.,0.)
    # print 'dP = ', T.evalPressure(0.,0.,rho = 1.225)

    # extact data along lines 
    # lines are defined as x0,y0,x1,y1
    # T.LinevalU(-0.5,-0.5,-0.5,0.5,plot_flag=1)
    # T.LinevalV(-0.5,-0.5,-0.5,0.5,plot_flag=1)
    # T.LinevalPressure(-0.5,-0.5,-0.5,0.5,rho = 1.225, plot_flag=1)

    # Make sure plots are displayed on the screen
    plt.show()
