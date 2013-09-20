import sys

from rl_defs import AtomicLevel, AtomicLine

# Some useful functions

def remove_spaces( tk ):
    return tk.replace(" ","")

def remove_braces( tk ):
    return tk.replace("(","").replace(")","").replace("[","").replace("]","")
    
def get_quantum_numbers( conf, term ):
    # l should be the last letter (s,p,d or f) in the string
    for c in conf:
        try: int(c)
        except:
            if   c=="s": l = 0
            elif c=="p": l = 1
            elif c=="d": l = 2
            elif c=="f": l = 3
            else: l = -1
        
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
    
    # parity is odd (1) if last char in term is '*', otherwise even (2)
    if term[-1]=="*": parity =  1
    else            : parity =  2
    
    # S and L should be first and second characters in term
    # NOTE: actual formula according to LAUX is 2S+1
    if term=="*":
        S = -1
        L = -1
    else: 
        S = (int(term[0])-1)/2
        L_str = term[1]
        if   L_str=="S": L = 0
        elif L_str=="P": L = 1
        elif L_str=="D": L = 2
        elif L_str=="F": L = 3
        else: L = -1
    
    return n,l,L,S,parity

# level classes

class IonizationLimit:
    def __init__(self,label,E):
        self.label = label
        self.E = E

class GroupedLevel:
    def __init__(self, level_tks ):
        # extract config and term info
        tks = level_tks[0]
        self.config = remove_spaces(tks[0])
        self.term = remove_spaces(tks[1])
        # get the g and E values for each individual level
        self.g_list = []
        self.E_list = []
        for tks in level_tks:
            if remove_spaces(tks[3])!="":
                g_tks = tks[2]
                if "/2" in g_tks:
                    self.g_list.append( int(remove_spaces(g_tks).split("/")[0])+1 )
                else:
                    self.g_list.append( int(remove_spaces(g_tks))*2+1 )
                self.E_list.append( float(remove_braces(remove_spaces(tks[3]))) )
        self.nlevels = len(self.g_list)
        # create multiplet level data
        self.g = sum( self.g_list )
        tmp = 0.0
        for i in range(self.nlevels):
            tmp += self.g_list[i] * self.E_list[i]
        self.E = tmp / self.g
        
    def get_python_string( self ):
        n,l,L,S,parity = get_quantum_numbers( self.config, self.term )
        return "%2d  %9.2f  %3d  %2d  %2d  %2d  %2d" % ( n, self.E, self.g, l, L, S, parity )

    def get_lua_string(self ):
        n,l,L,S,parity = get_quantum_numbers( self.config, self.term )
        return "%2d,  %9.2f,  %3d,  %2d,  %2d,  %2d,  %2d" % ( n, self.E, self.g, l, L, S, parity ) 

    def convert_grouped_data_to_AtomicLevel(self):
        n,l,L,S,parity = get_quantum_numbers( self.config, self.term )
        return AtomicLevel(n=n,E=self.E,g=self.g,l=l,L=L,S=S,parity=parity,conf=self.config,term=self.term)

    def convert_individual_data_to_AtomicLevel(self,i):
        n,l,L,S,parity = get_quantum_numbers( self.config, self.term )
        return AtomicLevel(n=n,E=self.E_list[i],g=self.g_list[i],l=l,L=L,S=S,parity=parity,conf=self.config,term=self.term)

class LevelData:
    def __init__(self):
        self.levels = []
        self.ionization_limits = []

    def make_levels_from_raw_data( self, raw_level_data ):
        count = 0
        for level_tks in raw_level_data:
            self.levels.append( GroupedLevel( level_tks ) )
            count += self.levels[-1].nlevels
        print "created %d multiplet levels, totally %d individual levels" % ( len(self.levels), count )
        
    def make_python_table( self, ofile_name, species, label ):
        full_label = species + "_" + label
        ofile = open( ofile_name, "w" )
        ofile.write( "%s = AtomicLevelSet()\n" % full_label )
        ofile.write( "%s.levels = [ '' ] * %d\n" % ( full_label, len(self.levels) ) )
        ofile.write( "%s.comments = '# grouped levels from NIST ASD'\n" % full_label )
        ofile.write( "# ----------------------------------------------------------------------\n" )
        ofile.write( "#              No.                 n    E(cm-1)    g   l   L   S  parity\n" )
        ofile.write( "# ----------------------------------------------------------------------\n" )
        for ilev,level in enumerate(self.levels):
            spaces = " "*(4-len(str(ilev)))
            ofile.write( "%s.levels[%d]%s=  '%s'\n" % ( full_label, ilev, spaces, level.get_python_string() ) )
        ofile.write( "#-----------------------------------------------------------------------\n" )
        ofile.write( "%s.available_level_sets[%s] = %s\n" % ( species, label, full_label ) )
        ofile.close()       
        
    def make_lua_table( self, ofile_name, species ):
        ofile = open( ofile_name, "w" )
        ofile.write( "%s.electronic_levels = {\n" % species )
        ofile.write( "   n_levels =  %d,\n" % ( len(self.levels) ) )
        ofile.write( "   ref = 'Grouped levels from NIST ASD'\n")
        ofile.write( "   -- ===========================================================\n" )
        ofile.write( "   --  No.       n      E(cm-1)    g    l    L    S    parity\n")
        ofile.write( "   -- ===========================================================\n" )
        for ilev,level in enumerate(self.levels):
            spaces = " "*(4-len(str(ilev)))
            ofile.write( "   ilev_%d%s= { %s },\n" % ( ilev, spaces, level.get_lua_string() ) )
        ofile.write( "   -- ===========================================================\n" )
        ofile.write( "}\n" )
        ofile.close() 

# line classes

class IndividualLine:
    def __init__(self, config_i, config_k, term_i, term_k, line_tks ):
        self.config_i = config_i
        self.config_k = config_k
        self.term_i = term_i
        self.term_k = term_k
        # get the parameters for this line
        tmp = line_tks[3].split("-")
        self.Ei = float( remove_spaces(tmp[0]) )
        self.Ek = float( remove_spaces(tmp[1]) )
        tmp = line_tks[4].split("-")
        self.gi = int( remove_spaces(tmp[0]) )
        self.gk = int( remove_spaces(tmp[1]) )
        A_str = remove_spaces(line_tks[5])
        if A_str=="":
            self.A = 0.0
        else:
            self.A = float( A_str )
        self.acc = remove_spaces( line_tks[6] )
        self.type = remove_spaces( line_tks[7] )

    def convert_to_AtomicLine( self ):
        return AtomicLine(g_u=self.gk, E_u=self.Ek, g_l=self.gi, E_l=self.Ei, A_ul=self.A, lambda_ul=-1, conf_u=self.config_k, term_u=self.term_k, conf_l=self.config_i, term_l=self.term_i, acc=self.acc, type_str=self.type, ilev_u=-1, ilev_l=-1 )

class MultipletLine:
    def __init__(self, prev_line, grouped_line_tks ):
        line_tks = grouped_line_tks[0]
        # extract config and term info
        tks = line_tks[1].split("-")
        if remove_spaces(line_tks[1])=="":
            # use prev line config data
            self.config_i = prev_line.config_i
            self.config_k = prev_line.config_k
        else:
            self.config_i = remove_spaces(tks[0])
            self.config_k = remove_spaces(tks[1])
        tks = line_tks[2].split("-")
        self.term_i = remove_spaces(tks[0])
        self.term_k = remove_spaces(tks[1])
        # check for multiplet data
        tmp = line_tks[3]
        if remove_spaces(tmp)=="":
            has_multiplet_data = False
        else:
            has_multiplet_data = True
            tks = line_tks[3].split("-")
            self.Ei = float( remove_spaces( tks[0] ) )
            self.Ek = float( remove_spaces( tks[1] ) )
            tks = line_tks[4].split("-")
            self.gi = int( remove_spaces( tks[0] ) )
            self.gk = int( remove_spaces( tks[1] ) )
            tk = line_tks[5]
            A_str = remove_spaces( tk )
            if A_str=="":
                self.A = 0.0
            else:
                self.A = float( remove_spaces( tk ) )
            tk = line_tks[6]
            self.acc = remove_spaces( tk )
            tk = line_tks[7]
            self.type = remove_spaces( tk )
        # create the individual lines
        self.lines = []
        for line_tks in grouped_line_tks[2:-1]:
            self.lines.append( IndividualLine( self.config_i, self.config_k, self.term_i, self.term_k, line_tks ) )
        self.nlines = len(self.lines)
        if not has_multiplet_data:
            # calculate multiplet data ourselves
	    icopies = []; kcopies = []
	    for iml in range(len(self.lines)):
	        # find unique levels
	        for jml in range(iml+1,len(self.lines)):
		    if self.lines[iml].Ei==self.lines[jml].Ei and iml not in icopies:
		        icopies.append(jml)
		    if self.lines[iml].Ek==self.lines[jml].Ek and iml not in kcopies:
		        kcopies.append(jml)
	    # i and kcopies should hold list of repeated states
            self.Ei = 0; self.Ek = 0
            self.gi = 0; self.gk = 0
            self.A = 0
            for i,line in enumerate(self.lines):
                if i not in icopies:
                    self.gi += line.gi
                    self.Ei += line.gi*line.Ei
                if i not in kcopies:
                    self.gk += line.gk
                    self.Ek += line.gk*line.Ek
                self.A += line.gk*line.A
            self.Ei /= self.gi
            self.Ek /= self.gk
            self.A /= self.gk

class LineData:
    def __init__(self):
        self.lines = []
        self.n_multiplet_lines = 0
        self.n_individual_lines = 0 

    def make_lines_from_raw_data( self, raw_line_data ):
        count = 0
        prev_line = None
        for line_tks in raw_line_data:
            self.lines.append( MultipletLine( prev_line, line_tks ) )
            prev_line = self.lines[-1]
            count += prev_line.nlines
        self.n_multiplet_lines = len(self.lines)
        self.n_individual_lines = count
        print "created %d multiplet lines, totally %d individual lines" % ( self.n_multiplet_lines, self.n_individual_lines )


# read functions

def read_level_file( filename, omit_psuedocontinuum_levels=True, use_individual_levels=False, echo_result=False ):
    # read in the file
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()

    # check the header
    tks = remove_spaces(lines[1]).split("|")
    if tks[0]!="Configuration" or tks[1]!="Term" or tks[2]!="J" or tks[3]!="Level":
        print "The header in file %s does not appear to match the expected format:" % filename
        print "Configuration          | Term   |   J |                 Level    |"
        sys.exit()

    # read in levels, one line at a time
    level_data = LevelData()
    raw_level_data = []
    level_tks = [] 
    for line in lines[3:]:
        tks = line.split("|")
        if len(tks)<=1 or "-----------" in line: continue
        if remove_spaces(tks[1])=="Limit":
            E_limit = float(remove_braces(tks[3]))
            print "Encountered %s ionization limit at %f" % ( tks[0], E_limit )
            level_data.ionization_limits.append( IonizationLimit( tks[0], E_limit ) )
            continue
        if ( remove_spaces(tks[0])=="" and remove_spaces(tks[3])=="" ) or tks[0]=="-----------------------":
            if len(level_tks)>0:
                raw_level_data.append( level_tks )
            level_tks = []
            continue
        else:
            level_tks.append( tks )

    # convert the tks to data
    level_data.make_levels_from_raw_data( raw_level_data )

    levels = []
    omission_count = 0
    for grouped_level in level_data.levels:
        if omit_psuedocontinuum_levels and grouped_level.E > level_data.ionization_limits[0].E:
            omission_count += 1
            continue
        elif use_individual_levels:
            for i in range(len(grouped_level.E_list)):
                levels.append( grouped_level.convert_individual_data_to_AtomicLevel(i) )
        else:
            levels.append( grouped_level.convert_grouped_data_to_AtomicLevel() )
        
    print "omitted %d grouped levels above the first ionization limit" % omission_count
    
    # sort the levels by ascending energy
    levels.sort(key=lambda x: x.E, reverse=False)

    return levels

def read_line_file( filename ):
    # read in the file
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()

    # check the header
    tks = remove_spaces(lines[1]).split("|")
    if tks[0]!="Mult." or tks[1]!="Configurations" or tks[2]!="Terms":
        print "The header in file %s does not appear to match the expected format:" % filename
        print "Mult. |           Configurations                 |     Terms    |       Ei           Ek       | gi   gk |    Aki    | Acc. |Type|"
        sys.exit()
        
    # read in lines, one block at a time
    raw_line_data = []
    line_tks = []
    for line in lines[5:]:
        tks = line.split("|")
        if len(tks)<=1: continue
        if remove_spaces(tks[0])!="":
            if len(line_tks)>0:
                raw_line_data.append( line_tks )
            line_tks = [ tks ]
            continue
        else:
            line_tks.append( tks )

    # convert the tks to data
    line_data = LineData()
    line_data.make_lines_from_raw_data( raw_line_data )

    lines = []
    for multiplet_line in line_data.lines:
        for line in multiplet_line.lines:
            if line.A!=0.0:
                lines.append( line.convert_to_AtomicLine() )

    return lines

def add_level_data_to_lines( lines, levels, exit_when_not_found=False, tol=1.0e-3 ):
    # find the upper and lower level indices of the lines in the list
    filtered_lines =[]
    for line in lines:
        u_found = False
        l_found = False
        for ilev,level in enumerate(levels):
            if line.conf_u == level.conf and line.term_u == level.term and abs(line.E_u - level.E)<tol:
                line.ilev_u = ilev
                u_found = True
            if line.conf_l == level.conf and line.term_l == level.term and abs(line.E_l - level.E)<tol:
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
            
    print "Found upper and lower level energies for %d out of %d lines" % ( len(filtered_lines), len(lines) )
        
    return filtered_lines 
    

