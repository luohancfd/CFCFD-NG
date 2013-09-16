import sys
from scipy.optimize import *
from numpy import *

from radpy import *

def pow(x,n):
    return x**n

def eval_width( x,p ):
    T_e = x
    gamma_S0,n = p
    return gamma_S0 * pow( (T_e / 1.0e4), n ) 

def residuals ( p,y,x ):
    return y - eval_width(x,p)

class Transition:
    def __init__( self, data_block ):
        # element name and wavelength on line 4
        print data_block
        line = data_block[3]
        tks = line.split()
        self.element = tks[2]
        self.lambda_ang = float(tks[5])

        # multiplet number of line 5
        line = data_block[4]
        tks = line.split()
        self.multiplet = int(tks[2])


        # transition label on line 6
        print data_block[5]
        line = data_block[5]
        tks = line.split("=")
        tks = tks[1].split(" / ")
        conf = tks[0]
        state = tks[1]
        tks = conf.split("-")
        self.lower_conf = tks[0]
        self.upper_conf = tks[1]
        tks = state.split("-")
        self.lower_state = tks[0]
        self.upper_state = tks[1]

        # now the data table
        self.T_list = []
        self.width_list = []
        self.shift_list = []
        self.alpha_list = []
        self.beta_list = []
        for line in data_block[8:12]:
            tks = line.split()
            print tks
            self.T_list.append( float(tks[0]) )
            self.width_list.append( float(tks[1].replace("-0","e-0").replace("+0","e+0") ) )
            self.shift_list.append( float(tks[2].replace("-0","e-0").replace("+0","e+0") ) )
            self.alpha_list.append( float(tks[3].replace("-0","e-0").replace("+0","e+0") ) )
            self.beta_list.append( float(tks[4].replace("-0","e-0").replace("+0","e+0") ) )


        # make the curve fit
        p0 = [ self.width_list[1], 0.3333 ]
        plsq = leastsq(residuals,p0,args=(array(self.width_list),array(self.T_list)))
        self.gamma_S0 = plsq[0][0]
        self.n = plsq[0][1]
        
        print "gamma_S0 = %e, n = %e" % ( self.gamma_S0, self.n )
        
    def get_transition_label(self):
        return "Upper config and state: %s - %s\nLower config and state: %s - %s\n" % ( self.upper_conf, self.upper_state, self.lower_conf, self.lower_state )

    def calculate_stark_width(self,Te,Ne):
        gamma_S = self.gamma_S0 * pow( (Te / 1.0e4), self.n ) * ( Ne / 1.0e16 )

        return gamma_S 

def approx_stark_width(Te,Ne,I,E_u,lambda_cl,model="johnston"):
    
    if model=="johnston":
        constA = 1.69e10
        constB = 2.623
    elif model=="cowley":
        constA = 9.27e07
        constB = 2.000
    elif model=="arnold":
        constA = 4.20e07
        constB = 2.000
    elif model=="park":
        constA = 5.0e7
        constB = 2.0

    constC =  0.5 * constA / pow( fabs(I - E_u), constB )
    gamma_S0 = constC * RC_c
    gamma_S = gamma_S0 * pow( (Te / 1.0e4), 0.33 ) * ( Ne / 1.0e16 )
    
    nu_cl = lambda2nu( lambda_cl/10. )
    
    # return the Stark HWHM in units of Ang
    return gamma_S / nu_cl * lambda_cl


def read_Griem_file( filename ):
    # open the data file
    infile = open( filename, "r" )
    lines = infile.readlines()
    infile.close()

    data_blocks = []
    data_block = []
    for line in lines:
        if line.find("Element :") >= 0:
            if len(data_block) > 0:
                data_blocks.append( data_block )
            data_block = []
        else:
            data_block.append(line)

    print "found %d data blocks" % len(data_blocks)

    transitions = []
    for data_block in data_blocks:
        transitions.append( Transition(data_block) )
        print transitions[-1].get_transition_label()

    # first test the temperature scaling model 
    for transition in transitions:
        rms_error = 0
        for i,T in enumerate(transition.T_list):
            scaled_width = transition.calculate_stark_width(T,1.0e16)
            print "T = %f, width = %f, scaled width = %f" % ( T, transition.width_list[i], scaled_width ) 
            rms_error += ((scaled_width- transition.width_list[i])/transition.width_list[i])**2
        print "RMS error = ", rms_error

    return transitions

def test_approx_models( transition, I ):
    print "wavelength = %f ang" % ( transition.lambda_ang )
    E = 1.0e8 / transition.lambda_ang
    for i,T in enumerate(transition.T_list):
        gamma_j = approx_stark_width(T,1.0e16,I,E,transition.lambda_ang,"johnston")
        gamma_c = approx_stark_width(T,1.0e16,I,E,transition.lambda_ang,"cowley")
        gamma_a = approx_stark_width(T,1.0e16,I,E,transition.lambda_ang,"arnold")
        gamma_p = approx_stark_width(T,1.0e16,I,E,transition.lambda_ang,"park")
        print "gamma_griem = %e, gamma_j = %e, gamma_c = %e, gamma_a = %e, gamma_p = %e" % ( transition.width_list[i], gamma_j, gamma_c, gamma_a, gamma_p )

