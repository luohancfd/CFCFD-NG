#!/usr/bin/env python
# extract_NIST_data.py

# Expected NIST format:
# -------------------------------------------------------------
# Configuration      | Term     |    J |             Level    |
# -------------------|----------|------|----------------------|
#                    |          |      |                      |
# 2s2.2p3            | 4S*      |  3/2 |             0        |

import sys, os
from copy import copy
from datetime import datetime 
import matplotlib.pyplot as plt
from radpy import *
from math import exp

def energy_group( E ):
    for i,E_limit in enumerate(energy_group_limits[:-1]):
        if E>=E_limit and E<energy_group_limits[i+1]: return i
    return len(energy_group_limits)-1

class Level:
    def __init__(self,n=-1,E=-1,g=-1,l=-1,L=-1,S=-1,parity=-1,conf="",term=""):
        self.n = n
        self.E = E
        self.g = g
        self.l = l
        self.L = L
        self.S = S
        self.parity = parity
        self.conf = conf
        self.term = term
    def string(self):
        lev_string = "n = " + str(self.n)
        lev_string += ", E = " + str(self.E)
        lev_string += ", g = " + str(self.g)
        lev_string += ", l = " + str(self.l)
        lev_string += ", L = " + str(self.L)
        lev_string += ", S = " + str(self.S)
        lev_string += ", parity = " + str(self.parity) + "\n"
        return lev_string
        
def get_quantum_numbers( conf, term ):
    # l should be the last letter (s,p,d or f) in the string
    for c in conf:
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
        S = (int(term[0]))/2
        L_str = term[1]
        if   L_str=="S": L = 0
        elif L_str=="P": L = 1
        elif L_str=="D": L = 2
        elif L_str=="F": L = 3
        else: L = -1
    
    return n,l,L,S,parity
    
def grouping_test( i_base, i ):
    # test level i to see if it should be grouped with level i_base
    if test_type==0:
        # energy level proximity
        if levels[i_base].E > individual_level_limit and abs(levels[i].E - levels[i_base].E)/levels[i].E < fE: return True
        else: return False
    elif test_type==1:
        # customized: Johnston N groupings
        if levels[i_base].n == 4 and levels[i].n == 4 and levels[i_base].l==2 and levels[i].l==2:
             # group all 4d states
             return True
        elif levels[i_base].n == 4 and levels[i].n == 4 and levels[i_base].l==3 and levels[i].l==3:
             # group all 4f states
             return True
        elif levels[i_base].n == levels[i].n and levels[i_base].n>4:
             # group all shells greater than 4 
             return True
        else:
             # an independent level
             return False
    elif test_type==2:
        # customized: Johnston O groupings
        if levels[i_base].n == levels[i].n and levels[i_base].n>5 and levels[i_base].n<10:
             # group all shells greater than 5 and less than 10
             return True
        elif levels[i_base].n>=10 and levels[i].n>=10:
             # group all remaining high-lying states
             return True
        else:
             # an independent level
             return False
    elif test_type==3:
        # group by shell number
        if levels[i_base].n == levels[i].n:
             # group all shells
             return True
        else:
             # an independent level
             return False
    elif test_type==4:
        # custom grouping
        if energy_group(levels[i_base].E)==energy_group(levels[i].E):
            return True
        else:
            return False
    else:
        return False
        
             
def find_energy_level_group( conf, term, groups ):
    for igrp, group in enumerate(groups):
        for ilev in group:
            if conf==levels[ilev].conf and term==levels[ilev].term: return igrp
    
    # if we get here, the search failed
    # print "search for conf: %s, term: %s failed." % ( conf, term )
    # for igrp, group in enumerate(groups):
        # for ilev in group:
            # print "levels[ilev].conf = %s, and term==levels[ilev].term = %s" % ( levels[ilev].conf, levels[ilev].term )
    # sys.exit()
    
    # Instead -> return -1 to indicate this level is super-ionized
    return -1
    
def extract_and_combine_levels( NIST_level_file ):
    infile = open(NIST_level_file, 'r')
    
    global levels
    levels = []

    global E_ionized
    E_ionized = 1.0e99

    # skip 3 header lines
    for i in range(4):
        line = infile.readline()

    count = -1
    n_ionized_levels = 0
    
    while line:
        # get on with some work
        line = infile.readline()
        print line
        # print line
        tks = line.split()
        
        # check if this is the very last line, will be of length 1
        # print len(tks)
        if len(tks)==0:
            break
        elif len(tks)==1:
            continue
        elif tks[0]=="|" and tks[1]=="|":
            print "a gap line"
            continue
        elif tks[1]=="II" or tks[1]=="III":
            # start of super-ionization levels
            tmp = tks[9].replace("[","").replace("]","")
            if "(" in tmp:
                E_ionized = float(tmp[:tmp.index("(")])
            else:
                E_ionized = float(tmp)
            continue
            
        count += 1
        
        # this should always be a 'composite' header
        term = tks[2]
        conf = tks[0]
        n,l,L,S,parity = get_quantum_numbers( conf, term )
        g = []; E = []
        
        if "/" in tks[4]:
            i_slash = tks[4].index("/")
            g.append(int(tks[4][:i_slash]) + 1)
        else:
            g.append(int(tks[4])*2 + 1)
            
        if tks[6]=="|":
            # missing energy level, record as -1
            E.append(-1)
        else:
            E.append(float(tks[6].strip("]").strip("[")))
        
        # print "n = %d, g = %d, E = %e" % (n, g[-1], E[-1])
        
        line = infile.readline()
        # print line
        tks = line.split()
        
        while tks[0]=="|" and tks[2]!="|":
            if "/" in tks[2]:
                i_slash = tks[2].index("/")
                g.append(int(tks[2][:i_slash]) + 1)
            else:
                g.append(int(tks[2])*2 + 1)
            if tks[4]=="|":
                # missing energy level, record as -1
                E.append(-1)
            else:
                E.append(float(tks[4].strip("]").strip("[")))
            # print " [sub] n = %d, g = %d, E = %e" % (n, g[-1], E[-1])
            line = infile.readline()
            # print line
            tks = line.split()
        
        # now have finished that level, calculate composite g and E
        
        g_sum = sum(g)
        E_av = 0.0; g_sum_av = 0
        
        for i in range(len(g)):
            if E[i]<0:
                # assume missing levels have average energy
                continue
            g_sum_av += g[i]
            E_av += g[i] * E[i]
            
        E_av /= g_sum_av

        # check if this level is above the detected ionization limit
        if E_av > E_ionized: n_ionized_levels+=1        

        # print "level = %d, g_sum = %d, E_av = %e, n = %d" % ( count, g_sum, E_av, n )
        levels.append(Level(n,E_av,g_sum,l,L,S,parity,conf,term))

    infile.close()
    
    # remove all levels with energy higher than the user-defined limit
    n_levels = len(levels)
    loop=True
    while loop:
        loop=False
        for ilev,lev in enumerate(levels):
             if lev.E > E_lev_limit or lev.n > n_lev_limit:
                 # print "Removing level with energy E = ", lev.E
                 levels.pop(ilev)
                 loop=True
                 break                 
    
    # now use individual level data (without fine splitting) to derive further condensed level groupings
    groups = []
    
    level_indices = []
    for i in range(len(levels)): level_indices.append(i)
    
    while len(level_indices)>0:
        i_base = level_indices[0]
        group = [ i_base ]
        # print "base: i = %d, n = %d" % ( i_base, levels[i_base].n )
        for index,item in enumerate(level_indices):
             if item==i_base: continue
             # print "test index = ", index 
             # print "test: i = %d, n = %d" % ( item, levels[item].n )
             if grouping_test( i_base, item ):
                 # print "test success"
                 group.append( item )
        for ilev in group: level_indices.remove(ilev)
        groups.append(group)
        
    # create combined level data, sort in terms of energy
    combined_levels = []
    global ordered_groups
    ordered_groups = []
    
    for ilev,group in enumerate(groups):
        # print "len(combined_levels) = ", len(combined_levels)
        if len(group)==1:
            # single level
            level = levels[group[0]]
        else:
            # combined level - assume n is constant, give -1 for all other Q numbers
            E_av = 0.0; g_sum = 0
            for jlev in group:
                E_av += levels[jlev].E * levels[jlev].g
                g_sum += levels[jlev].g
            E_av /= g_sum
            level = Level( levels[group[0]].n, E_av, g_sum )
        combined_levels.append( level )
        ordered_groups.append(group)
        for i, clevel in enumerate(combined_levels):
            # print "level.E = %f, clevel.E = %f" % ( level.E, clevel.E )
            if level.E < clevel.E:
                combined_levels.insert(i,level)
                ordered_groups.insert(i,group)
                combined_levels.pop(-1)
                ordered_groups.pop(-1)
                break
                
    # write to file
    fout = open('NIST_levels.dat','w')
    if output_format=="librad":
        fout.write("# Levels obtained from NIST ASD:\n\
# http://physics.nist.gov/PhysRefData/ASD/index.html\n\
# Using file: %s\n\
# Col. 1: Level number\n\
# Col. 2: Principal quantum number\n\
# Col. 3: Level energy Ek (cm**-1)\n\
# Col. 4: Level degeneracy g\n\
# Col. 5: Oribital quantum number l\n\
# Col. 6: OAM quantum number L\n\
# Col. 7: Total spin quantum number S\n\
# Col. 8: Parity\n\n" % sys.argv[1])
        tag = X + "_all_levels.levels"
        b1 = "["
        b2 = "]"
        c1 = "'"
        c2 = "'"
        for ilev,level in enumerate(combined_levels):
            fout.write("%s%s%1d%s %s %s%2d %11.2f %5d %4d %4d %4d %4d%s\n" % (tag, b1, ilev, b2, ' = ', c1, level.n,
                                                                                      level.E,
                                                                                      level.g,
                                                                                      level.l,
                                                                                      level.L,
                                                                                      level.S,
                                                                                      level.parity,
                                                                                      c2 ) )
    elif output_format=="libgas":
        fout.write("%s.electronic_levels = {\n" % X)
        fout.write("   n_levels = %d,\n" % len(combined_levels))
        fout.write("   ref = 'NIST ASD: http://physics.nist.gov/PhysRefData/ASD/index.html, %s',\n" % datetime.now())
        fout.write("   -- ===========================================================\n")
        fout.write("   --   No.     n      E(cm-1)      g     l     L     S    parity\n") 
        fout.write("   -- ===========================================================\n")
        for ilev,level in enumerate(combined_levels):
            fout.write("   ilev_%d = { %2d, %11.2f, %5d, %4d, %4d, %4d, %4d },\n" % (ilev, level.n,
                                                                                      level.E,
                                                                                      level.g,
                                                                                      level.l,
                                                                                      level.L,
                                                                                      level.S,
                                                                                      level.parity ) )
        fout.write("   -- ===========================================================\n}\n")
    elif output_format=="tau":
        fout.write("Microstates electronic: ")
        for ilev,level in enumerate(combined_levels):
            fout.write("%5d " % level.g )
        fout.write("\n")
        fout.write("Theta electronic: ")
        for ilev,level in enumerate(combined_levels):
            theta = nu2T(level.E)        # energy as a characteristic temperature in kelvin
            fout.write("%11.2f " % theta )
        fout.write("\n")
    fout.close()

    if 1:
        # plot level energies and write to file
        i_list = []
        E_list = []
        efile = open("level_energies.txt","w")
        efile.write("# Column 1: level number\n")
        efile.write("# Column 2: level energy (cm-1)\n")
        efile.write("# Column 3: level energy (K)\n")
        for ilev,level in enumerate(combined_levels):
            i_list.append( ilev )
            E_list.append( nu2T(level.E) )
            efile.write("%d %e %e\n" % ( ilev, level.E, nu2T(level.E) ) )
        efile.close()
        plt.plot(i_list,E_list, 'o')
        plt.show()
        
    if 1:
        # calculate partition function for a range of temperatures and write to file
        Qfile = open("partition_function.txt","w")
        Qfile.write("# Column 1: Temperature, T (K)\n")
        Qfile.write("# Column 2: Partition function, Q (m**-3)\n")
        T_list = range(200,1100,100) + range(1500,10500,500) + range(11000,21000,1000) + range(25000,55000,5000)
        Q_list = []
        for T in T_list:
            Q = 0.0
            for ilev,level in enumerate(combined_levels):
                Q += level.g * exp( - nu2T(level.E) / T )
            Q_list.append( Q )
            Qfile.write("%e \t %e\n" % ( T, Q ) )
        Qfile.close()
        plt.plot(T_list,Q_list, 'o')
        plt.show()
    
    print "composite levels = %d, total levels = %d, removed levels = %d, ionized levels = %d" % ( len(groups), len(levels), n_levels - len(levels), n_ionized_levels )
    return 0
    
def extract_lines( NIST_line_file, energy_switch ):
    infile = open(NIST_line_file, 'r')
    
    # skip 5 header lines
    for i in range(6):
        line = infile.readline()
        
    # open file for writing
    fout = open('NIST_lines.dat','w')
    fout.write("# Multiplet lines obtained from NIST ASD:\n\
# http://physics.nist.gov/PhysRefData/ASD/index.html\n\
# Using file: %s\n\
# Col. 1: Multiplet number\n\
# Col. 2: Lower state energy Ei (cm**-1)\n\
# Col. 3: Upper state energy Ek (cm**-1)\n\
# Col. 4: Lower state degeneracy gi\n\
# Col. 5: Upper state degeneracy gk\n\
# Col. 6: Einstein coefficient for spontaneous emission Aki\n\
# Col. 7: Lower state energy level group Li\n\
# Col. 8: Upper state energy level group Lk\n\
# Col. 9: Type (0: Allowed, 1: Forbidden)\n\n" % sys.argv[2] )

    line_count = 0
    i_lines = 0
    m_lines = 0
    super_ionized_lines = 0
    
    while line:
        print line
        # get on with some work
        tks = line.split()
        # check if this is the very last line, will be of length 0
        if len(tks)==1:
            # print line
            break
        # this should always be a multiplet header
        multiplet = int(tks[0])
        # extract configurations and terms
        conf_data = True
        if tks[2]!="|":
            conf_i = tks[2]; conf_k = tks[4]
            term_i = tks[6]; term_k = tks[8]
        else:
            # use previous configurations
            term_i = tks[3]; term_k = tks[5]
            conf_data = False
        # search for energy level group
        # print "line = ", line
        Li = find_energy_level_group( conf_i, term_i, ordered_groups )
        Lk = find_energy_level_group( conf_k, term_k, ordered_groups )
        # store multiplet data in case individual lines are missing
        if len(tks)>18:
            # print line
            offset = 0
            if conf_data==False: offset = -3 
            Ei = float(tks[10+offset]); Ek = float(tks[12+offset])
            gi = int(tks[14+offset]); gk = int(tks[16+offset])
            if tks[16+offset]=="|": offset -= 1
            if tks[18+offset]=="|": Aki = 0.0
            else: Aki = float(tks[18+offset])        
            if tks[20+offset]=="|": offset -= 1
            if len(tks)<23+offset: type_str=tks[-1]
            else: type_str = tks[22+offset]
            if type_str=="|": type = 0
            else: type = 1
        # read two lines forward
        line = infile.readline()
        line = infile.readline()
        tks = line.split()
        iml = 0
        gis = []; gks = []; Akis = []; Eis = []; Eks = []; types = []
        while tks[0] == "|":
            if tks[3] == "|":
                line = infile.readline()
                break
            # this is a sublevel with data
            print line
            # store sublevel data
            Eis.append( float( tks[3].replace("[","").replace("]","") ) )
            Eks.append( float( tks[5].replace("[","").replace("]","") ) )
            gis.append( int( tks[7] ) )
            gks.append( int( tks[9] ) )
            if tks[11]=="|":
                # the transition probability is missing
                Akis.append( float( 0.0 ) )
            else:
                Akis.append( float( tks[11] ) )
            if len(tks)==16:
                type_str = tks[15]
            elif len(tks)==15:
                type_str = tks[14]
            else:
                type_str = tks[-1]
            if type_str=="|": type = 0
            else: type = 1
            types.append(type)
            iml += 1
            line = infile.readline()
            tks = line.split()
        # check if the individual line data is empty
        if iml==0:
            print "individual line data missing - using multiplet data!"
            gis.append( gi ); gks.append( gk ); Akis.append(Aki); Eis.append( Ei ); Eks.append( Ek ); types.append(type)
            iml = 1
        # calculate multiplet effectjve line properties
        # copy degeneracy arrays
        nml = iml
        gEi = 0; gi = 0; gEk = 0; gk = 0; gAki = 0;
        gis_copy = copy(gis); gks_copy = copy(gks)
        Eis_copy = copy(gis); Eks_copy = copy(gks)
        # calculate products
        icopies = []; kcopies = []
        for iml in range(nml):
            gAki += gks[iml] * Akis[iml]
            # find unique levels
            for jml in range(iml+1, nml):
                if Eis[iml]==Eis[jml] and iml!=jml and iml not in icopies:
                    icopies.append(jml)
                if Eks[iml]==Eks[jml] and iml!=jml and iml not in kcopies:
                    kcopies.append(jml)
        # i and kcopies should hold list of repeated states
        for iml in range(nml):
            if iml not in icopies:
                gi += gis[iml]
                gEi += gis[iml] * Eis[iml]
            if iml not in kcopies:
                gk += gks[iml]
                gEk += gks[iml] * Eks[iml]
        # check that type is consistent - if not, return -1
        for iml in range(nml):
            if types[0]!=types[iml]: type=-1
        # decode effective line properties
        Ei = gEi / gi; Ek = gEk / gk
        Aki = gAki / gk
        # check that the transition probability is above 0
        if not Aki>0.0:
            continue
        # check if this line begins or ends from above the set ionization limit (but keep in the set)
        if Ek>E_ionized or Ei>E_ionized:
            super_ionized_lines += 1
        # apply energy switch for mulitplets/individuals
        tag = X + "_all_lines.lines"
        b1 = "["
        b2 = "]"
        c1 = "'"
        c2 = "'"
        if (Ek-Ei)<energy_switch:
            # multiplet line
            fout.write("%s%s%1d%s %3s %s%9.2f %11.2f %4d %4d %1.2e %4d %4d %4d%s\n" % (tag, b1, line_count, b2, ' = ', c1, Ei, Ek, gi, gk, Aki, Li, Lk, type, c2))
            line_count += 1
        else:
            # individual lines
            m_lines += 1
            for iml in range(nml):
                fout.write("%s%s%1d%s %3s %s%9.2f %11.2f %4d %4d %1.2e %4d %4d %4d%s\n" % (tag, b1, line_count, b2, ' = ', c1, Eis[iml], Eks[iml], gis[iml], gks[iml], Akis[iml], Li, Lk, types[iml], c2))
                line_count += 1
                i_lines += 1
        # for i in range(len(tks)):
        # print "tks[" + str(i) + "] = " + tks[i] + "\n"

    print "line_count = ", line_count
    print "i_lines = ", i_lines
    print "m_lines = ", m_lines
    print "super_ionized_lines = ", super_ionized_lines
    fout.close()
    infile.close()
    return 0

def printUsage():
    print "extract_NIST_ASD.py"
    print "Work with txt output from http://physics.nist.gov/PhysRefData/ASD/index.html."
    print "Usage(s):"
    print "extract_NIST_ASD.py <level.file>"
    print "extract_NIST_ASD.py <level.file> <line.file>"
    sys.exit(1)

def main():
    if len(sys.argv)!=2 and len(sys.argv)!=3:
        printUsage()
    elif len(sys.argv)==2:
        just_levels = True
    else:
        just_levels = False
        
    # set some parameters
    global output_format
    output_format = raw_input("Select: output format [ 'tau', 'libgas', 'librad' ]: ")
    more_options = int(raw_input("Select: maximum resolution calculation (0), or proceed to grouping and energy options (1): "))
    global test_type, fE, individual_level_limit
    global E_lev_limit, n_lev_limit
    if more_options:
        test_type = int(raw_input("Energy level grouping method [ (0) proximity, (1) Johnston N, (2) Johnston O, (3) by shell, or (4) from groupings.txt ]: "))
        if test_type<0 or test_type>4:
            print "test_type = %d is not allowed" % test_type
            sys.exit()
        elif test_type==0:
            individual_level_limit = float(raw_input("Energy of last individual level:"))
            fE = float(raw_input("Energy level proximity fraction:"))
        elif test_type==4:
            ifile_name = raw_input("Energy grouping filename: ")
            ifile = open(ifile_name,"r")
            lines = ifile.readlines()
            ifile.close()
            global energy_group_limits
            energy_group_limits = []
            for line in lines:
                tks = line.split()
                energy_group_limits.append( float(tks[0]) )
        E_lev_limit = float(raw_input("Energy limit for levels [ eg. E_ionize, N = 117214.0, O = 109836.7, no-limit = 1.0e99 ]: "))
        n_lev_limit = float(raw_input("shell limit for levels [ eg. n>5 use Park data for Oxygen ]: "))
    else:
        test_type = -1
        fE = 0.0
        E_lev_limit = 1.0e99
        n_lev_limit = 9999

    # get the species string
    global X
    X = sys.argv[1][:sys.argv[1].find("-")]
    if         sys.argv[1][sys.argv[1].find("-")+2] == "I":
        X += "_plus"

    # firstly do electronic levels
    NIST_level_file = sys.argv[1]
    extract_and_combine_levels( NIST_level_file )
    
    if not just_levels:
        # now do lines...
        NIST_line_file = sys.argv[2]
        energy_switch = float(raw_input("Multiplet-to-individual-line energy switch [eg. 6.0 eV used by CJ]: "))
        energy_switch *= 8065.6            # convert eV to cm-1
        extract_lines( NIST_line_file, energy_switch )
    
    print "Done."
        
if __name__ == '__main__':
    main()
