#!/usr/bin/env python
# nenzfr_format_sensitivity_data_for_LaTeX.py
#
# This scripts loads in sensitivity data that was produced
# by nenzfr_sensitivity.dat and re-writes the data out to a
# tabular format suitable for use with LaTeX and inclusion
# in ones' thesis. 
#
# Although there are LaTeX packages that can do this, I wanted
# to present the data in the same tabular format as used in 
# previous thesis' and David Mee's T4 Uncertainty Departmental
# report. 
#
# This requires a significant reformatting of the data which
# is easier to accomplish using python.
#
# Currently the writing out of the data is quite specific to
# my needs.
#
# Luke Doherty
# 27-Feb-2013

VERSION_STRING = "27-Feb-2013"

import string
import sys, os
E3BIN = os.path.expandvars("$HOME/e3bin")
sys.path.append(E3BIN)

#------------------------------------------------------------
def read_data_file(FileToRead):
    """
    """
    fp = open(FileToRead,'r')

    fp.readline() # Don't need very first line
    # Get names of perturbed variables
    header = fp.readline().strip().split(" ")
    perturbed_vars = [k.strip() for k in header if k!="" and k!="variable"]
    #print perturbed_vars
    fp.readline() # Don't need this line of ------
    # Get the rest of the data
    data = fp.readlines()
    
    # Now assemble a dictionary...
    DataDict = {}
    exitVars = []
    for line in data:
        # Extract line data, removing white/empty space...
        lineData = [k for k in line.strip().split(" ") if k!=""]
        # Assemble exitVars list...
        exit_var = lineData[0]
        exitVars.append(exit_var) 
        # Now put the line data into a dictionary...
        DataDict[exit_var] = {}
        for k in range(len(lineData)-1):
            DataDict[exit_var][perturbed_vars[k]] = lineData[k+1]
     
    return perturbed_vars, exitVars, DataDict

def write_reformated_data_file(FileToWrite,perturbed_vars,exitVars,\
        allDataDicts,separator,EOL):
    """
    """
    sep = separator
    N1 = len(exitVars)
    N2 = len(perturbed_vars)
    
    fp = open(FileToWrite,'w')
    count1 = 0 
    for e_var in exitVars:
        fp.write(("{0:}"+sep).format(e_var))
        count1 += 1
        count2 = 0
        
        if e_var in ['supply_T','supply_h']:
            for p_var in ['p1','T1','Vs','pe']:
                count2 += 1
                fp.write(("{0:}"+sep).format(p_var))
                for k in range(len(allDataDicts)):
                    if k==len(allDataDicts)-1:
                        fp.write("{0:}".format( allDataDicts[k][e_var][p_var] ))
                    else:
                        fp.write(("{0:}"+sep).format( allDataDicts[k][e_var][p_var] ))

                fp.write(EOL)
                fp.write('\n')
                if count2 != 4:
                    fp.write(sep)
                #else:
                #    fp.write('\midrule\n')
        else:
            for p_var in perturbed_vars:
                count2 += 1
                fp.write(("{0:}"+sep).format(p_var))
                for k in range(len(allDataDicts)):
                    if k==len(allDataDicts)-1:
                        fp.write("{0:}".format( allDataDicts[k][e_var][p_var] ))
                    else:
                        fp.write(("{0:}"+sep).format( allDataDicts[k][e_var][p_var] ))

                fp.write(EOL)
                fp.write('\n')
                # If we have reached the end of the perturbed variables the
                # next line will be for the next exit variable and so we
                # don't need a separator at the start of the line
                if count2 != N2:
                    fp.write(sep)
                #else:
                    #if count1 != N1:
                    #    fp.write('\midrule\n')
                    #else:
                    #    fp.write('\\bottomrule\n')
        
        #for p_var in perturbed_vars:
        #    count += 1
        #    if e_var in ['supply_T','supply_h']:
        #        if p_var in ['p1','T1','Vs','pe']:
        #            fp.write(("{0:}"+sep).format(p_var))
        #            for k in range(len(allDataDicts)):
        #                if k==len(allDataDicts)-1:
        #                    fp.write("{0:}".format( allDataDicts[k][e_var][p_var] ))
        #                else:
        #                    fp.write(("{0:}"+sep).format( allDataDicts[k][e_var][p_var] ))
        #            fp.write(EOL)
        #            fp.write('\n')
        #            if count !=
        #    else:
        #        fp.write(("{0:}"+sep).format(p_var))
        #        for k in range(len(allDataDicts)):
        #            if k==len(allDataDicts)-1:
        #                fp.write("{0:}".format( allDataDicts[k][e_var][p_var] ))
        #            else:
        #                fp.write(("{0:}"+sep).format( allDataDicts[k][e_var][p_var] ))
        #    
        #        fp.write(EOL)
        #        fp.write('\n')
        #        # If we have reached the end of the perturbed variables the
        #        # next line will be for the next exit variable and so we 
        #        # don't need a separator at the start of the line
        #        if count != N:
        #            fp.write(sep)
                
        
    fp.close()

def main():
    
    # Am assuming that the perturbed and exit variables
    # are identical in each of the following files...
    
    # Load HP condition RELATIVE sensitivities
    pert_vars, exitVars, hpSensRel = \
            read_data_file("./sensitivity_hp_air/sensitivities.dat")
    # Load HP condition ABSOLUTE sensitivities
    pert_var, exitVars, hpSensAbs = \
            read_data_file("./sensitivity_hp_air/sensitivities_abs.dat")

    # Load LP condition RELATIVE sensitivities
    pert_vars, exitVars, lpSensRel = \
            read_data_file("./sensitivity_lp_air/sensitivities.dat")
    # Load LP condition ABSOLUTE sensitivities
    pert_vars, exitVars, lpSensAbs = \
            read_data_file("./sensitivity_lp_air/sensitivities_abs.dat")

    # I hardcode the exit variables so that the order in which they 
    # get written out is what I want for my thesis.
    exitVars = ['supply_T','supply_h','p','T[0]','rho','vel.x','vel.y',\
                'a','M_local','pitot_p','total_p','total_h','q','m_dot',\
                'Re_u','p/q','mu','k[0]','e[0]','massf[0]-N2','massf[1]-O2',\
                'massf[2]-N','massf[3]-O','massf[4]-NO','mu_t','k_t','tke',\
                'omega','dt_chem']
    
    allDataDicts = [hpSensAbs,hpSensRel,lpSensAbs,lpSensRel]
    #print len(allDataDicts)
    
    write_reformated_data_file("./sensitivities-for-thesis.tex",\
             pert_vars,exitVars,allDataDicts,"&","")

if __name__ == '__main__':
    return_flag = main()
    sys.exit(return_flag)
