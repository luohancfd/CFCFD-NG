#! /usr/bin/env python

## \file    chemkin2libgas.py
## \author  Brendan O'Flaherty
## 
## \version 22 Jan 2008 -- File created
## 
## ---------------------------------------------------------

import ConfigParser
from getopt import getopt
from string import split,atof,atoi,replace,lower

import sys

from math import exp
from numpy import array, arange, log, zeros
from parse_reaction import *
from cea_curves import *

#-----------------------------------------------------------
# 

def is_pd(scanner, parser, r_str):
    tokens = scanner.tokenize(r_str)
    r_dict = parser.parse(tokens)
    for r in r_dict['reactants']:
        if 'pd_term' in r:
            return True
    for r in r_dict['products']:
        if 'pd_term' in r:
            return True
        
    # if we make it this far, assume:
    return False

def has_third_body(scanner, parser, r_str):
    tokens = scanner.tokenize(r_str)
    r_dict = parser.parse(tokens)
    
    for r in r_dict['reactants']:
        for sp in r['sp']:
            if sp == 'M':
                return True
    for r in r_dict['products']:
        for sp in r['sp']:
            if sp == 'M':
                return True

    # if we make it this far, assume:
    return False

def n_species(scanner, parser, r_str):
    tokens = scanner.tokenize(r_str)
    r_dict = parser.parse(tokens)
    nf = 0
    for i in range(len(r_dict['reactants'])):
        nf += r_dict['reactants'][i]['coeff']

    nr = 0
    for i in range(len(r_dict['products'])):
        nr += r_dict['products'][i]['coeff']
        
    return nf, nr

def format_reaction(line_by_space):
    for i in range(len(line_by_space)):
        if line_by_space[i] == '=':
            line_by_space[i] = ' <=> '

    rxn = [''.join(line_by_space[:-3])]
    for n in line_by_space[-3:]:
        rxn.append(n)
    return rxn

def format_line(str_list):
    # do some formatting on this line
    for index in range(len(str_list)):
        str_list[index] = format_word(str_list[index])
    return str_list

def format_word(word):
    # do some formatting on this word
    word = replace(word, 'HCO', 'CHO')
    #if (luaOutput): word = replace(word, 'CH2(S)', 'CH2_S')
    word = replace(word, 'AR', 'Ar')
    word = replace(word, 'HE', 'He')
    word = replace(word, 'E+', 'e')
    word = replace(word, '+', ' + ')
    word = replace(word, '( + M)', ' ( + M )')
    word = replace(word, '>', '> ')
    word = replace(word, '<', ' <')
        
    return word

def get_eff_token(line, eff_token, luaOutput=False):
    # remove spaces
    for i in range(len(line)):
        line[i] = replace(line[i], ' ', '')
    line.pop() # remove newline character
    
    eff_list_token = ''
    if luaOutput:
        for i in range(0,len(line),2):
            eff_list_token += "%s=%.2f, " % (line[i], atof(line[i+1]))
        eff_token = "   efficiencies={%s},\n" % (eff_list_token[:-2])

    else:
        for i in range(0,len(line),2):
            eff_list_token += "(\'%s\', %.2f), " % (line[i], atof(line[i+1]))
        eff_token = ",\n\t\tefficiencies=[%s]" % (eff_list_token[:-2])
    
#   if eff_token == '': # efficiencies
#       eff_token = ",\n\t\tefficiencies=[%s]" % eff_list_token[:-2]
#   else: # line-wrapped efficiencies
#       eff_token = "%s, %s]"% (eff_token[:-1], eff_list_token[:-2])
    
    return eff_token

def get_inf_token(line, luaOutput=False):
    if luaOutput:
        inf_token = "A=%.5e, n=%.5e, T_a=%.5e*S" % (atof(line[0]),
                                                    atof(line[1]),
                                                    atof(line[2]))
    else:
        inf_token = "%.5e, %.5e, %.5e*S" % (atof(line[0]),
                                            atof(line[1]),
                                            atof(line[2]))
    return inf_token

# -----------------------------------------------------------------------------

def extract_lua_thermo_data(i0, fi_lines, newSpecies):
    no_segments = 2

    for i in range(i0, len(fi_lines), 4):
        lineInput = []
        species = fi_lines[i].split()[0]

        if species == 'END':
            return 0

        species = format_word(species)
        species_var = replace(species, '(', '')
        species_var = replace(species_var, ')', '')

        newSpecies.append((species, species_var))

        Trange = map(atof, fi_lines[i].split()[-4:-1])
        for j in range(i+1, i+4): # lines in this block
            for k in range(0, 75, 15): # letters in this line
                try:
                    n = fi_lines[j][k:k+15] # cannot split() chemkin format
                    lineInput.append(n)
                except:
                    continue

        species_mass = 0.0
        if species == 'Ar':
            species_mass = element_mass[species]
	elif species == 'He':
	    species_mass = element_mass[species]
        elif species == 'Rh(s)':
            species_mass = element_mass[species]
        else:
            for element in species:
                try:
                    e = element_mass[element]
                    species_mass += e
                except:
                    try:
                        species_mass += e*(atof(element)-1)
                    except:
                        if element == '(' or element == ')' or lower(element) == 's':
                            continue
                        else:
                            print "# balking on element '%s' in species '%s'" % (element, species)
                            sys.exit(1)
        
        fcurve = []
        gamma = []
        Tmin = []
        Tmax = []
        n = 0
        for m in range(no_segments-1, -1, -1):
            fcurve.append(zeros(9, float))
            Tmin.append(min(Trange[n], Trange[-1]))
            Tmax.append(max(Trange[n], Trange[-1]))
            for i in range(7):
                fcurve[-1][i+2] = lineInput[m*7+i];
            cpor = eval_cp_curve(300.0, fcurve[-1])
            gamma.append(cpor/(cpor-1.0))
            n = n+1

        thermo_source = "The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)"
        M_source = "molecular weight from CEA2"
        fname = "%s.lua" % (species)
        fout = open(fname, 'w')

        fout.write("-- Collater: Brendan T. O'Flaherty\n")
        fout.write("-- Date: 07 Aug 2009\n")
        fout.write("\n")
        fout.write("%s = {}\n" % (species))
        fout.write("%s.M = {\n" % (species))
        fout.write("   value = %.7e,\n" % (species_mass))
        fout.write("   units = 'kg/mol',\n")
        fout.write("   description = 'molecular mass',\n")
        fout.write("   reference = '%s'\n" % (M_source))
        fout.write("}\n")
        fout.write("%s.gamma = {\n" % (species))
        fout.write("   value = %.4e,\n" % (gamma[0]))
        fout.write("   units = 'non-dimensional',\n")
        fout.write("   description = 'ratio of specific heats at 300.0K',\n")
        fout.write("   reference = 'evaluated using Cp/(Cp - R) from Chemkin-II coefficients'\n")
        fout.write("}\n")
        fout.write("%s.CEA_coeffs = {\n" % (species))
        for i in range(n):
            fout.write("   { T_low  = %.1f,\n" % (Tmin[i]))
            fout.write("     T_high = %.1f,\n" % (Tmax[i]))
            fout.write("     coeffs = {")
            for j in range(len(fcurve[i])):
                   fout.write("%10.9e, " % (fcurve[i][j]))
            fout.write("}\n")
            fout.write("   },\n")
        fout.write("   ref='%s'\n" % (thermo_source))
        fout.write("}\n")
    return

def extract_thermo_data(i0, fi_lines, newSpecies):
    # thermodynamic data source
    source_name = "The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)"
    no_segments = 2

    for i in range(i0, len(fi_lines), 4):
        lineInput = []
        species = fi_lines[i].split()[0]

        if species == 'END':
            return 0

        species = format_word(species)
        species_var = replace(species, '(', '')
        species_var = replace(species_var, ')', '')

        newSpecies.append((species, species_var))

        Trange = map(atof, fi_lines[i].split()[-4:-1])
        for j in range(i+1, i+4): # lines in this block
            for k in range(0, 75, 15): # letters in this line
                try:
                    n = fi_lines[j][k:k+15] # cannot split() chemkin format
                    lineInput.append(n)
                except:
                    continue

        species_mass = 0.0
        if species == 'Ar':
            species_mass = element_mass[species]
	elif species == 'He':
	    species_mass = element_mass[species]
        elif species == 'Rh(s)':
            species_mass = element_mass[species]
        else:
            for element in species:
                try:
                    e = element_mass[element]
                    species_mass += e
                except:
                    try:
                        species_mass += e*(atof(element)-1)
                    except:
                        if element == '(' or element == ')' or lower(element) == 's':
                            continue
                        else:
                            print "# balking on element '%s' in species '%s'" % (element, species)
                            sys.exit(1)

        print ("%-5s = {'type'           : 'atom'," % species_var)
        print ("\t 'mol_weight'     : %g," % species_mass)
        print ("\t 'sigma'          : 3.711e-10,")
        print ("\t 'epsilon'        : 1.08468e-21,")
        print ("\t 'warning'        : \"\"\"No data for %s for: sigma, epsilon, polarizability or dipole_moment, viscosity, thermal conductivity. Using Air values.\"\"\"," % species)

        print ("")

        print ("\t#")
        print ("\t# %s" % (source_name))
        print ("\t#")

        print ("\t'no-s-segments'  : %i," % no_segments)
        n = 0
        heat_curves = ''
        for m in range(no_segments-1, -1, -1):
            heat_curves += ("\t'range-s-%i'      : (%.1f, %.1f),\n" % (n, min(Trange[n], Trange[-1]),
                                                                              max(Trange[n], Trange[-1])))
            heat_curves += ("\t'coeffs-s-%i'     : [0.0, 0.0,%s,%s,\n" % (n, lineInput[m*7+0], lineInput[m*7+1]))
            heat_curves += ("\t                    %s,%s,%s,%s,\n" % (lineInput[m*7+2], lineInput[m*7+3], lineInput[m*7+4], lineInput[m*7+5]))
            heat_curves += ("\t                    %s ],\n" % (lineInput[m*7+6]))
            n = n+1

        print "%s," % (heat_curves[0:-2])

        # the rest of the data defaults to that of air
        print "\t'no-mu-segments' : 2,\n\t'range-mu-0'     : (200.0, 1000.0),\n\t'coeffs-mu-0'     : [ 6.51082150e-01, -1.95095944e+01, -1.65482348e+03, 1.60332824e+00 ],\n\t'range-mu-1'     : (1000.0, 8000.0),\n\t'coeffs-mu-1'     : [ 8.18769927e-01,  4.53839453e+02, -1.73904090e+05, 1.39947428e-01 ],\n\t'no-k-segments'  : 3,\n\t'range-k-0'      : (200.0,  1000.0),\n\t'coeffs-k-0'     : [  8.63732593e-01,  1.00133164e+02, -1.11342781e+04, 4.38933582e-01 ],\n\t'range-k-1'      : (1000.0, 5100.0),\n\t'coeffs-k-1'     : [  8.50172611e-01,  3.34676528e+01,  1.49365413e+04, 5.73104994e-01 ],\n\t'range-k-2'      : (5100.0, 8000.0),\n\t'coeffs-k-2'     : [  1.59493808e+00,  6.19368726e+02,  1.10458108e+07,-6.32373477e+00 ]"
        print "}"

    return
# -----------------------------------------------------------------------------

def get_reaction_token(i, fi_lines, speciesList, nc, rate_token, rate_name, 
                       inf_token, pd_reaction, third_body_reaction, luaOutput=False):
    eff_token = ''

    if not pd_reaction:
        if third_body_reaction: 
            # check for efficiencies
            while True:
                i += 1
                try:
                    line_by_slash = fi_lines[i].split('/')
                except:
                    print "Error in input file."
                    print "Reached end of file without 'END' statement"
                    sys.exit()
                    
                line_by_slash = format_line(line_by_slash)
                first_word = replace(line_by_slash[0], ' ', '')

                if '!' in first_word: continue
                elif first_word in speciesList:
                    # these are efficiencies
                        eff_token = get_eff_token(line_by_slash, eff_token, luaOutput)
                        
                else: # end of third body reaction
                    i = i - 1 # revert counter
                    break

        if luaOutput:
            rate_token += "'Arrhenius', %s},\n" % (inf_token)
        else:
            nc_token = "nc=%i" % (nc)
            rate_token += "GA_%s(%s, %s)" % (rate_name, inf_token, nc_token)
    else:
        # PD reaction
        # we expect additional lines of input
        troe_token = ''
        if luaOutput:
            rate_token += "'pressure dependent',\n       k_inf={%s},\n" % (inf_token)
        else:
            nc_token = "nc=%i" % (nc-1)
            rate_token += "PD_%s(inf=(%s), \n" % (rate_name, inf_token)
        while True:
            i +=1
            line_by_slash = fi_lines[i].split('/')
            line_by_slash = format_line(line_by_slash)
            first_word = replace(line_by_slash[0], ' ', '')
            
            if '!' in first_word: continue
            elif first_word == 'LOW':
                # this is low-p reaction data
                low = line_by_slash[1].split()
                if luaOutput:
                    low_token = "A=%.5e, n=%.5e, T_a=%.5e*S" % (atof(low[0]),
                                                                atof(low[1]),
                                                                atof(low[2]))
                else:
                    low_token = "%.5e, %.5e, %.5e*S" % (atof(low[0]),
                                                        atof(low[1]),
                                                        atof(low[2]))
            elif first_word == 'TROE':
                # this is troe data
                troe = line_by_slash[1].split()
                if luaOutput:
                    troe_token += "a=%s, T3=%s, T1=%s, T2=%s" % \
                        (troe[0], troe[1], troe[2], troe[3])
                else:
                    print troe_token
                    for tt in troe:
                        troe_token += "%.4e, " % (atof(tt))
                    troe_token = ",\n\t\ttroe=[%s]" % (troe_token[:-2])
                    
            elif first_word in speciesList:
                # these are efficiencies
                line_by_slash = fi_lines[i].split('/')
                line_by_slash = format_line(line_by_slash)
                eff_token = get_eff_token(line_by_slash, eff_token, luaOutput)

            else:
                # end of pd-reaction, revert counter
                i = i-1
                break
        if luaOutput:
            rate_token += "       k_0={%s},\n" % (low_token)
            if troe_token != "":
                rate_token += "       Troe={%s},\n" % (troe_token)
            rate_token += "   },\n"
        else:
            rate_token += "\t\t              low=(%s), %s)%s" % (low_token, nc_token, troe_token)
    return i, rate_token, eff_token

#-----------------------------------------------------------
# elemental masses

element_mass = {'H' : 1.00794e-3,
                'C' : 12.0107e-3,
                'O' : 15.99940e-3,
                'N' : 14.0067e-3, 
                'Ar': 39.948e-3,
		'He' : 4.003e-3,
                'Rh(s)': 102.90550e-3}

#-----------------------------------------------------------
# command line options

shortOptions = 'h'
longOptions = ['help', 'therm=', 'rate=', 'lua']
    
def printUsage():
    print "Converts files of chemkin2 format to lib/gas or libgas2 format."
    print "Writes --rate to STDOUT, --therm to file."
    print "Note: remember to rename CH2(S).lua to CH2_S.lua"
    print
    print "Usage: python chemkin2libgas.py [options]"
    print " -h, --help             Prints this help"
    print " --rate=<baseFileName>  Expects a file containing a kinetic mechanism in Chemkin-II format."
    print " --therm=<baseFileName> Expects a file containing thermodynamic data in Chemkin-II format."
    print " --lua                  Output written in lib/gas (lua) format. "
    
    return

def main():
    """
    
    """
    
    userOptions, progArgs = getopt(sys.argv[1:], shortOptions, longOptions)

    if len(userOptions) == 0:
        printUsage()
        sys.exit(0)
    else:
        uoDict = dict(userOptions)

    if uoDict.has_key('--help') or uoDict.has_key('-h'):
        printUsage()
        sys.exit(0)

    rateData = False
    thermData = False
    luaOutput = False

    if uoDict.has_key('--rate'):
        rateData = True
        rateInput = uoDict.get('--rate', None)
    if uoDict.has_key('--therm'):
        thermData = True
        thermInput =  uoDict.get('--therm', None)
    if uoDict.has_key('--lua'):
        luaOutput = True
    
    if (thermData):
        # load all data into memory first
        fi = open(thermInput+'.therm', 'r')
        fi_lines = fi.readlines()
        fi.close()

        reactionList = []
        newSpecies = []
        scanner = ReactionScanner()
        parser = ReactionParser()
        i=0 # this is our counter
        if not luaOutput:
            print "# this file has been automatically generated by chemkin2libgas.py\n"
        while True:
            data = fi_lines[i].split()
            if data[0] == '!':
                i = i+1
                continue

            elif data[0] == 'THERMO':
                i = i+2 # skip the temperature range
                continue
            
            else:
                if luaOutput:
                    if not extract_lua_thermo_data(i, fi_lines, newSpecies):
                        break
                else:
                    if not extract_thermo_data(i, fi_lines, newSpecies):
                        break

        if not luaOutput:
            print
            availableSpecies = ''
            for s0, s1 in newSpecies:
                availableSpecies += ("'%s' : %s, " % (s0,s1))
            print ('available_species = {%s}' % availableSpecies[:-2])        

    if (rateData):
        # load all data into memory first
        print rateInput
        fi = open(rateInput+'/'+rateInput+'.dat', 'r')
        fi_lines = fi.readlines()
        fi.close()

        reactionList = []
        speciesList  = []
        scanner = ReactionScanner()
        parser = ReactionParser()
        i=0 # this is our counter

        if luaOutput:
            print "-- this file has been automatically generated by chemkin2libgas.py\n"
        else:
            print "# this file has been automatically generated by chemkin2libgas.py\n"

        for i in range(len(fi_lines)):
            data = fi_lines[i].split()
            # empty line
            if len(data) == 0: 
                continue
            # comment
            if data[0][0] == '!':
                continue
            
            if data[0] == 'ELEMENTS':
                elements = []

                while True:
                    i = i+1
                    element_line = fi_lines[i].split()
                    if element_line[0] == 'END':
                        break
                    for j in element_line:
                        elements.append(j)

                elements = format_line(elements)
                temp = ''
                for e in elements:
                    temp += "\'%s\', " % e

                if luaOutput: continue
                print "elements_list = [%s]" % temp[0:-2]
                print "ndata.declare_elements(elements_list)\n"

            if data[0] == 'SPECIES':
                species = []

                while True:
                    i = i+1
                    species_line = fi_lines[i].split()
                    if species_line[0][0] == '!': 
                        continue
                    if species_line[0] == 'END':
                        break
                    for j in species_line:
                        if j[0] == '!': 
                            break
                        species.append(j)

                species = format_line(species)
                sl = ''
                for s in species:
                    if s[0]=='!': continue
                    sl += "\'%s\', " % s.split()[0]

                speciesList = sl[0:-2]
                if luaOutput: 
                    print "model = \"thermally perfect gas\""
                    print "species = {%s}" % speciesList
                else:
                    print "species_list = [%s]" % speciesList
                    print "ndata.declare_species(species_list)\n"

            if data[0] == 'REACTIONS':
                if luaOutput:
                    print "-- scaling factor, 1/R"
                else:
                    print "# scaling factor, 1/R"
                    
                if len(data) > 1:
                    if data[-1] == 'KJOULES/MOLE':
                        #kj/mol
                        print "S = 1.0/8.314472e-3\n" 
                    elif data[-1] == 'KCAL/MOLE':
                        #kj/mol
                        print "S = 1.0/1.987e-3\n" 
                    else:
                        print "Unknown units. Please modify this script accordingly."
                        sys.exit(1)
                else:
                    # cal/mol
                    print "S = 1.0/1.987\n"

                j = 0
                while True:
                    i = i+1
                    
                    line_by_slash = fi_lines[i].split('/')
                    line_by_space = fi_lines[i].split()

                    line_by_slash = format_line(line_by_slash)
                    line_by_space = format_line(line_by_space)

                    # handle any special terms here...
                    if line_by_slash[0][0]=='!': continue
                    if line_by_space[0][0]=='!': continue

                    if line_by_space[0] == 'DUPLICATE':
                        if luaOutput:
                            print "-- duplicate reaction found, including..."
                        else:
                            print "# duplicate reaction found, including..."
                        continue

                    if line_by_space[0] == 'END':
                        print
                        break
                    
                    # the rest should be part of a reaction
                    line_by_space = format_reaction(line_by_space)
                    br_token = ''

                    if len(line_by_slash) == 1:
                        # here we start a forward reaction
                        r_str = line_by_space[0]
                        if luaOutput:
                            new_reaction_token = "reaction{'%s',\n" % (r_str)
                        else:
                            new_reaction_token = "r%i = make_reaction(\"%s\",\n" % (j, r_str)
                        
                        # we always need these
                        ncf, ncb = n_species(scanner, parser, r_str)
                        inf_token = get_inf_token(line_by_space[1:], luaOutput)
                            
                        # test what kind of reaction this will be
                        if luaOutput:
                            fr_token = "   fr={"
                        else:
                            fr_token = "\t\tfr="
                        fr_name = "forward"
                        pd_reaction = is_pd(scanner, parser, r_str)
                        third_body_flag = has_third_body(scanner, parser, r_str)
                        i, fr_token, eff_token = get_reaction_token(i, fi_lines, speciesList, ncf, 
                                                                    fr_token, fr_name, inf_token, 
                                                                    pd_reaction, third_body_flag, 
                                                                    luaOutput)
                        i += 1
                        
                        # here we start the reverse reaction
                        line_by_slash = fi_lines[i].split('/')
                        line_by_slash = format_line(line_by_slash)
                        first_word = replace(line_by_slash[0], ' ', '')
                        if first_word == 'REV':
                            # this has a reverse rate
                            rev = line_by_slash[1].split()
                            inf_token = get_inf_token(rev, luaOutput)
                            
                            if luaOutput:
                                br_token += "   br={"
                            else:
                                br_token += ",\n\t\tbr="
                            br_name = "backward"
                            i, br_token, eff_token = get_reaction_token(i, fi_lines, speciesList, ncb, 
                                                                        br_token, br_name, inf_token, 
                                                                        pd_reaction, third_body_flag, 
                                                                        luaOutput)
                            
                        else:
                            # no reverse rate, revert counter
                            i = i - 1
                    else:
                        print "What's going on here???"
                        print " line %d" % i
                        print " %s" % line_by_slash
                        sys.exit(1)
                        
                    if luaOutput:
                        label = "   label='r%i'\n" % (j)
                        print new_reaction_token+fr_token+br_token+eff_token+label+'}'
                    else:
                        print new_reaction_token+fr_token+br_token+eff_token+')'
                    reactionList.append(j)
                    j = j+1
        rl = ''
        for rxn in reactionList:
            rl += "r%i, " % rxn
        if not luaOutput:
            print "reactions_list = [%s]" % rl[0:-2]
            print "ndata.declare_reactions(reactions_list)"

    return

if __name__ == '__main__':
    main()



