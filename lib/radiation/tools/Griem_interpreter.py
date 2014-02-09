import sys
from scipy.optimize import *
from numpy import *
import itertools

from radpy import *

def pow(x,n):
    return x**n

def eval_width( x,p ):
    T_e = x
    gamma_S0,n = p
    return gamma_S0 * pow( (T_e / 1.0e4), n ) 

def residuals ( p,y,x ):
    return y - eval_width(x,p)
    
def get_NIST_ASD_confs_and_states( line ):
    conf_u = line.conf_u.split(".")[-1].upper()[:2]
    conf_l = line.conf_l.split(".")[-1].upper()[:2]
    term = line.term_u
    if "[" in term:
        state_u = term[term.index("[")+1:term.index("]")]
    else:
        state_u = term
    term = line.term_l
    if "[" in term:
        state_l = term[term.index("[")+1:term.index("]")]
    else:
        state_l = term
    return conf_u, conf_l, state_u, state_l
    
def get_TOPBase_confs_and_states( line ):
    conf_u = line.conf_u.split(".")[-1].upper()[:2]
    conf_l = line.conf_l.split(".")[-1].upper()[:2]
    term = line.term_u
    if "[" in term:
        state_u = term[term.index("[")+1:term.index("]")]
    else:
        state_u = term
    term = line.term_l
    if "[" in term:
        state_l = term[term.index("[")+1:term.index("]")]
    else:
        state_l = term
    return conf_u, conf_l, state_u, state_l
    
def process_Griem_float( tk ):
    return float(tk.replace("-0","e-0").replace("+0","e+0") )

class GriemStarkWidth:
    def __init__( self, data_block, verbose=False ):
        # element name and wavelength on line 4
        if verbose:
            print data_block
        line = data_block[3]
        tks = line.split()
        self.element = tks[2].replace(" ","")
        self.lambda_ang = float(tks[5])

        # multiplet number of line 5
        line = data_block[4]
        tks = line.split()
        self.multiplet = int(tks[2])


        # transition label on line 6
        if verbose:
            print data_block[5]
        line = data_block[5]
        tks = line.split("=")
        tks = tks[1].split(" / ")
        conf = tks[0]
        state = tks[1]
        tks = conf.split("-")
        self.lower_conf = tks[0].replace(" ","").replace("\n","")
        self.upper_conf = tks[1].replace(" ","").replace("\n","")
        tks = state.split("-")
        self.lower_state = tks[0].replace(" ","").replace("\n","")
        self.upper_state = tks[1].replace(" ","").replace("\n","")

        # now the data table
        self.T_list = []
        self.width_list = []
        self.shift_list = []
        self.alpha_list = []
        self.beta_list = []
        for line in data_block[8:12]:
            tks = line.split()
            if verbose: print tks
            self.T_list.append( process_Griem_float(tks[0]) )
            self.width_list.append( process_Griem_float(tks[1]) )
            self.shift_list.append( process_Griem_float(tks[2]) )
            self.alpha_list.append( process_Griem_float(tks[3]) )
            self.beta_list.append( process_Griem_float(tks[4]) )


        # make the curve fit
        p0 = [ self.width_list[1], 0.3333 ]
        plsq = leastsq(residuals,p0,args=(array(self.width_list),array(self.T_list)))
        self.gamma_S0 = plsq[0][0]
        self.n = plsq[0][1]
        
        if verbose:
            print "gamma_S0 = %e, n = %e" % ( self.gamma_S0, self.n )
        
    def get_transition_label(self):
        return "Upper config and state: %s - %s\nLower config and state: %s - %s\n" % ( self.upper_conf, self.upper_state, self.lower_conf, self.lower_state )

    def get_string(self):
        return "%s, lambda_ang = %e, multiplet = %d, self.lower_conf = %s, self.upper_conf = %s, self.lower_state = %s, self.upper_state = %s, " % ( self.element, self.lambda_ang, self.multiplet, self.lower_conf, self.upper_conf, self.lower_state, self.upper_state )

    def calculate_stark_width(self,Te,Ne):
        gamma_S = self.gamma_S0 * pow( (Te / 1.0e4), self.n ) * ( Ne / 1.0e16 )

        return gamma_S 

    def find_associated_line(self,lines,line_source,tol,allow_inexact_matches,verbose):
        # set the get_confs_and_states function pointer
        if line_source=="NIST_ASD":
            get_confs_and_states = get_NIST_ASD_confs_and_states
        elif line_source=="TOPBase":
            get_confs_and_states = get_TOPBase_confs_and_states
        else:
            print "GriemStarkWidth::find_associated_line()"
            print "line_source not recognised"
            sys.exit()
        associated_line = None
        # firstly try the closest line with matching confs and states
        delta_lambda_min = 9.9e9
        closest_line = None
        found = False
        if verbose: 
            print "------------------------------------------------------------"
        for line in lines:
            # check if this line has the same conf and term
            conf_u, conf_l, state_u, state_l = get_confs_and_states( line )
            if conf_u==self.upper_conf and conf_l==self.lower_conf and state_u==self.upper_state and state_l==self.lower_state:
                # see how close this line is 
                delta_lambda = abs(line.lambda_ul - self.lambda_ang/10.0)
                # print "this line is %e percent away and matches: %s" % ( delta_lambda, line.get_string() )
                if delta_lambda < delta_lambda_min:
                    delta_lambda_min = delta_lambda
                    closest_line = line
        if closest_line:
            delta_lambda = abs(closest_line.lambda_ul - self.lambda_ang/10.0)
            error = (delta_lambda/(self.lambda_ang/10.0))
            if error <= tol:
                if verbose:
                    print "a matching line within %f percent was found" % ( tol*100 )
                associated_line = closest_line
                found = True
            elif verbose:
                print "the closest matching line is more than %f percent away: %e percent" % ( tol*100, error*100 )
        elif verbose:
            print "no matching line was found"
        if not found and allow_inexact_matches:
            # use the best matching line within the tolerance
            close_lines = []
            close_errors = []
            for line in lines:
                error = abs(line.lambda_ul - self.lambda_ang/10.0) / (self.lambda_ang/10.0)
                if error < tol:
                    close_lines.append(line)
                    close_errors.append(error)
                    # print "this line is %e percent away: %s" % ( delta_lambda, line.get_string() )
            best_matches = 0
            best_line = None
            best_error = 9.9e9
            for line,error in itertools.izip(close_lines,close_errors):
                if verbose: print line.get_string()
                # check what conf's and term's match
                conf_u, conf_l, state_u, state_l = get_confs_and_states( line )
                number_of_matches = 0
                if conf_u==self.upper_conf:
                    #print "conf_u: %s != upper_conf: %s" % ( conf_u, self.upper_conf )
                    number_of_matches += 1
                if conf_l==self.lower_conf:
                    #print "conf_l: %s != lower_conf: %s" % ( conf_l, self.lower_conf )
                    number_of_matches += 1
                if state_u==self.upper_state:
                    #print "conf_u: %s != upper_conf: %s" % ( state_u, self.upper_state )
                    number_of_matches += 1
                if state_l==self.lower_state:
                    #print "conf_l: %s != lower_conf: %s" % ( state_l, self.lower_state )
                    number_of_matches += 1
                if number_of_matches > best_matches:
                    best_line = line
                    best_matches = number_of_matches
                    best_error = error
                if number_of_matches == best_matches and error < best_error:
                    best_line = line
                    best_matches = number_of_matches
                    best_error = error
            if best_line:
                if verbose: 
                    print "A line within %e percent and with %d matches will be used..." % ( best_error*100, best_matches )
                found = True
                associated_line = best_line
            elif verbose:
                print "No line out of %d close lines was adequate for this stark parameter: " % len(close_lines)
                print self.get_string()
                
        return associated_line

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


def read_Griem_file( filename, verbose=False ):
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

    transitions = []
    for data_block in data_blocks:
        transitions.append( GriemStarkWidth(data_block) )
        if verbose:
            print transitions[-1].get_transition_label()

    # first test the temperature scaling model 
    for transition in transitions:
        rms_error = 0
        for i,T in enumerate(transition.T_list):
            scaled_width = transition.calculate_stark_width(T,1.0e16)
            if verbose:
                print "T = %f, width = %f, scaled width = %f" % ( T, transition.width_list[i], scaled_width ) 
            rms_error += ((scaled_width- transition.width_list[i])/transition.width_list[i])**2
        if verbose: print "RMS error = ", rms_error

    print "Found Griem Stark width models for %d transitions" % len(transitions)

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

def add_Stark_width_parameters_to_lines( Stark_widths, lines, line_source="NIST_ASD", tol=0.02, allow_inexact_matches=False, verbose=False ):
    found_count = 0
    for Stark_width in Stark_widths:
        associated_line = Stark_width.find_associated_line(lines,line_source,tol,allow_inexact_matches,verbose)
        if associated_line:
            associated_line.n = Stark_width.n
            associated_line.gamma_S0 = Stark_width.gamma_S0
            found_count += 1
    print "Found associated lines for %d out of %d provided Stark width models" % ( found_count, len(Stark_widths) )
    