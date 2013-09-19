import sys

from rl_defs import AtomicLevel, AtomicLine, TOPBasePICSLevel
from radpy import RC_Ry

def get_SLP( SLP ):
    S = (int(SLP[0])-1)/2
    L = int(SLP[1])
    parity = int(SLP[2])
    # use NIST Parity indices
    if parity==0: parity=2
    
    return S,L,parity

def make_term( S, L, parity ):
    # make NIST-like term string
    term = "%d" % (2*S+1)
    if   L==0: term += "S"
    elif L==1: term += "P"
    elif L==2: term += "D"
    elif L==3: term += "F"
    elif L==4: term += "G"
    elif L==5: term += "H"
    else:
        print "L is greater than 5!"
        sys.exit()
    if parity==1: term += "*"
    
    return term
    
def make_conf( CONF ):
    # make NIST-like conf string
    tmp = CONF.replace("(",".(").replace(")",").").replace("O","*")
    tks = tmp.split()
    conf = tks[0]
    for tk in tks[1:]: conf += "." + tk
    conf = conf.replace("..",".")
    
    return conf
    
def get_l( conf ):
    # l should be the last letter (s,p,d or f)
    l = -1
    for c in conf:
        if   c=="s": l = 0
        elif c=="p": l = 1
        elif c=="d": l = 2
        elif c=="f": l = 3
        elif c=="g": l = 4
        elif c=="h": l = 5

    if l < 0:
        print "l is greater than 5!"
        print conf
        sys.exit()

    return l
    
def get_n( conf ):
    # n should be last number in the string after a dot
    ns = []
    for i in range(50):
        ns.append(str(i))
    for i in range(len(conf)-1):
        if conf[i]=="." and conf[i+1] in ns: dot_index = i
        else: dot_index = -1
        
    for i in range(dot_index+2,len(conf)):
        n_maybe = conf[dot_index+1:i]
        if n_maybe in ns: n = int(n_maybe) 
        
    return n

def make_level_from_TOPBase_string( raw_level_string ):
    # extract the sub strings
    tks = raw_level_string[:22].split()
    ilev = tks[0]
    NZ = tks[1]
    NE = tks[2]
    iSLP = tks[3]
    iLV = tks[4]
    iCONF = raw_level_string[23:34]
    tks = raw_level_string[35:].split()
    E_RYD = tks[0]
    TE_RYD = tks[1]
    gi = tks[2]
    QD = tks[3]
    EQN = tks[4]
    RL_NS = tks[5]
        
    # make the easy data
    E = float(TE_RYD) * RC_Ry  # convert Ry -> cm-1
    g = int(float(gi))
    
    # make the more complicated data
    S,L,parity = get_SLP( iSLP )
    term = make_term( S, L, parity )
    conf = make_conf( iCONF )
    l = get_l( conf )
    n = get_n( conf )
        
    return AtomicLevel(n,E,g,l,L,S,parity,conf,term)
    
def make_line_from_TOPBase_string( raw_level_string ):
    # extract the sub strings
    # NOTE: here i is upper level and j is lower level
    tks = raw_level_string[:31].split()
    iline = tks[0]
    NZ = tks[1]
    NE = tks[2]
    iSLP = tks[3]
    jSLP = tks[4]
    iLV = tks[5]
    jLV = tks[6]
    iCONF = raw_level_string[32:44]
    jCONF = raw_level_string[48:59]
    tks = raw_level_string[64:].split()
    gF = tks[0]
    gA = tks[1]
    WL_A = tks[2]
    gi = tks[3]
    gj = tks[4]

    # make the easy data
    g_u = int(float(gi))
    E_u = -1
    g_l = int(float(gj))
    E_l = -1
    A_ul = float(gA) / float(g_u)
    lambda_ul = float(WL_A) / 10.0
    
    # make the more complicated data
    S,L,parity = get_SLP( iSLP )
    term_u = make_term( S, L, parity )
    conf_u = make_conf( iCONF )
    S,L,parity = get_SLP( jSLP )
    term_l = make_term( S, L, parity )
    conf_l = make_conf( jCONF )

    return AtomicLine(g_u, E_u, g_l, E_l, A_ul, lambda_ul, conf_u, term_u, conf_l, term_l, acc="TB", type_str=""  )

def make_PICS_from_TOPBase_strings( header_line, data_lines ):
    # first the header line
    tks = header_line.split()
    NZ = int(tks[1])
    NE = int(tks[2])
    ISLP = tks[3]
    ILV = int(tks[4])
    E_RYD = float(tks[5])
    NP = int(tks[6])
    
    # decode some of this stuff
    S,L,parity = get_SLP( ISLP )
    term = make_term( S, L, parity )
    E_cm = E_RYD * RC_Ry       # convert Ry -> cm-1
    ilevTB = ILV 
    
    # check the data_lines size is correct
    if len(data_lines)!=NP:
        print "len(data_lines) = %d, NP = %d!" % ( len(data_lines), NP )
        print data_lines
        sys.exit()
        
    # make the list of energies and cross-sections (currently, unsure of units: probably Ry and cm2)
    E_au_list = []
    sigma_au_list = []
    for line in data_lines:
        tks = line.split()
        E_au_list.append( float(tks[0]) )
        sigma_au_list.append( float(tks[1]) )
        
    return TOPBasePICSLevel( E_cm, term, ilevTB, E_au_list, sigma_au_list )

def read_level_file( filename, include_psuedocontinuum_levels, echo_result=False ):
    # read in the file
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()
   
    # check the header
    if "==================================================" not in lines[0] or \
    "i NZ NE iSLP iLV iCONF                 E(RYD)      TE(RYD)   gi         QD     EQN    RL(NS)" not in lines[1] or \
     "==================================================" not in lines[2]:
        print "TOPBase level file doesn't have expected format"
        sys.exit()

    # read in levels, one line at a time
    raw_level_data = []
    for line in lines:
        tks = line.split()
        if len(tks)<=1: continue
        try: ilev = int(tks[0]) - 1
        except: continue
        else: raw_level_data.append( line )
        
    print "Found %d TOPBase levels in file %s" % ( len(raw_level_data), filename)
        
    # pull out the data strings and make AtomicLevel instances
    levels = []
    for raw_level_string in raw_level_data:
        levels.append( make_level_from_TOPBase_string( raw_level_string ) )
        if echo_result: print levels[-1].get_string()
        
    # sort the levels by ascending energy
    sorted_levels = []
    nlevels = len(levels)
    while len(sorted_levels)<nlevels:
        E_min = 1.0e12
        for ilev,level in enumerate(levels):
            if level.E < E_min:
                ilev_min = ilev
                E_min = level.E
        sorted_levels.append(levels[ilev_min])
        levels.pop(ilev_min)
        
    # FIXME: apply psuedocontinuum level filter

    return sorted_levels
    
def read_line_file( filename, echo_result=False ):
    # read in the file
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()
   
    # check the header
    if "==================================================" not in lines[0] or \
    "      i NZ NE iSLP jSLP iLV jLV iCONF           jCONF                  gF   gA(S-1)     WL(A)   gi   gj" not in lines[1] or \
     "==================================================" not in lines[2]:
        print "TOPBase line file doesn't have expected format"
        sys.exit()
        
    # read in lines, one line at a time
    raw_line_data = []
    for line in lines:
        tks = line.split()
        if len(tks)<=1: continue
        try: ilev = int(tks[0]) - 1
        except: continue
        else: raw_line_data.append( line )
        
    print "Found %d TOPBase lines in file %s" % ( len(raw_line_data), filename)
        
    # pull out the data strings and make AtomicLine instances
    lines = []
    for raw_line_string in raw_line_data:
        lines.append( make_line_from_TOPBase_string( raw_line_string ) )
        if echo_result: print lines[-1].get_string()
        
    # make all transitions in the positive (absorption) direction
    for line in lines:
        if line.A_ul < 0.0:
            # swap upper and lower level data
            # degeneracies
            tmp = line.g_u
            line.g_u = line.g_l
            line.g_l = tmp
            # energies
            tmp = line.E_u
            line.E_u = line.E_l
            line.E_l = tmp
            # configuration
            tmp = line.conf_u
            line.conf_u = line.conf_l
            line.conf_l = tmp
            # term
            tmp = line.term_u
            line.term_u = line.term_l
            line.term_l = tmp
            # lev
            tmp = line.ilev_u
            line.ilev_u = line.ilev_l
            line.ilev_l = tmp
            # and finally make the transition probability positive
            line.A_ul *= -1.0
        
    return lines
    
def read_PICS_file( filename, echo_result=False ):
    # read in the file
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()
   
    # check the header
    if "===========================" not in lines[0] or \
    "I  NZ  NE  ISLP  ILV        E(RYD)      NP" not in lines[1] or \
     "============================" not in lines[2]:
        print "TOPBase PICS file doesn't have expected format"
        sys.exit()
        
    # read in lines, one line at a time
    PICSs = []
    data = []
    header = ""
    for line in lines[3:]:
        tks = line.split()
        if len(tks)<=1: continue
        try: I = int(tks[0])
        except: data.append(line)
        else:
            if len(header) > 0:
                PICSs.append( make_PICS_from_TOPBase_strings( header, data ) )
                header = ""
                data = [] 
            header = line

    print "Found %d TOPBase PICS in file %s" % ( len(PICSs), filename)

    # reset the energies to be referenced from 0 at the ground state
    E0 = PICSs[0].E
    for PICS in PICSs:
        PICS.E -= E0
    
    if echo_result:
        for PICS in PICSs:
            print PICS.get_string()
        
    return PICSs
    
def check_line_levels( lines, levels ):
    # find the upper and lower levels of the lines in the list
    
    return
    
def add_level_data_to_lines( lines, levels, exit_when_not_found=False ):
    # find the upper and lower level energies of the lines in the list
    filtered_lines = []
    for line in lines:
        u_found = False
        l_found = False
        for ilev,level in enumerate(levels):
            if line.conf_u == level.conf and line.term_u == level.term:
                line.E_u = level.E
                line.ilev_u = ilev
                u_found = True
            if line.conf_l == level.conf and line.term_l == level.term:
                line.E_l = level.E
                line.ilev_l = ilev
                l_found = True
        if not u_found and exit_when_not_found:
            print "The upper state for the line below was not found in the provided list of levels!"
            print line.get_string()
            sys.exit()
        if not l_found and exit_when_not_found:
            print "The lower state for the line below was not found in the provided list of levels!"
            print line.get_string()
            sys.exit()
        if u_found and l_found:
            filtered_lines.append( line )
        # FIXME: currently setting type (optically allowed/forbidden??) to empty
        #        eventually use level data to set this correctly
        line.type_str = ""
            
    print "Found upper and lower level energies for %d out of %d lines" % ( len(filtered_lines), len(lines) )
        
    return filtered_lines       

def get_PICS_with_level_indices_and_datapoints( levels, PICSs, tol=1.0e-6, require_term_match=True, verbose=False ):
    new_PICSs = []
    for PICS in PICSs:
        # check that this PICS has datapoints
        if len(PICS.sigma_list)==0: continue
        found = False
        # find the levels with matching terms (if required)
        if not require_term_match:
            filtered_ilevs = range(len(levels))
        else:
            filtered_ilevs = []
            for ilev,level in enumerate(levels):
                if PICS.term == level.term:
                    filtered_ilevs.append( ilev )
        if verbose: print "found %d levels meeting the term requirement" % len(filtered_ilevs)
        # find which of these levels meet the energy requirement
        further_filtered_ilevs = []
        for ilev in filtered_ilevs:
            level = levels[ilev]
            error = abs(level.E - PICS.E)
            if error < tol: further_filtered_ilevs.append(ilev)
        if verbose: print "found %d levels meeting the energy requirement" % len(further_filtered_ilevs)
        # find the level with the closest energy
        closest_ilev = None
        smallest_error = 9.9e9
        for ilev in further_filtered_ilevs:
            level = levels[ilev]
            error = abs(level.E - PICS.E)
            if error < smallest_error:
                smallest_error = error
                closest_ilev = ilev
        if verbose: print "The closest level has an energy %e away" % smallest_error
        # check that this error is good enough
        if smallest_error < tol:
            PICS.ilev = closest_ilev
            new_PICSs.append( PICS) 
            found = True
        
        if not found and verbose:
            print "WARNING: Did not find a level for PICS with E = %e, term = %s" % ( PICS.E, PICS.term )

    print "Found level indices for %d out of %d PICS" % ( len(new_PICSs), len(PICSs) )

    return new_PICSs
