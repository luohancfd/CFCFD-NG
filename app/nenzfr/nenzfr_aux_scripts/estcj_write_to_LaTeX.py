#!/usr/bin/env python
"""
This short python script runs estcj with the desired 
inputs and then writes the output to a sexy-looking 
LaTeX table format ready for inclusion in ones' thesis
or report.

Luke Doherty
18-March-2013
"""

VERSION_STRING = "18-March-2013"

import sys, os
E3BIN = os.path.expandvars("$HOME/e3bin")
sys.path.append(E3BIN)

from estcj import reflected_shock_tube_calculation as rstc

def main():
    # Some code to allow me to easily choose between
    # conditions
    condition = 'lp'

    if condition in ['hp',]:
        p1 = 250.e3   # Pa
        T1 = 300.     # K
        Vs = 2194.5   # m/s
        pe = 71.525e6 # MPa
        ar = 1581.165 
        FileToWrite = "hp-nom-estcj-table.tex"
        Label = "estcj_hp_nom"
        Title = "\estcj~output for the nominal high pressure condition."
    elif condition in ['lp',]:
        p1 = 160.e3   # Pa
        T1 = 300.     # K
        Vs = 2299.7   # m/s
        pe = 40.062e6 # MPa
        ar = 1581.165
        FileToWrite = "lp-nom-estcj-table.tex"
        Label = "estcj_lp_nom"
        Title = "\estcj~output for the nominal low pressure condition."
    gas = 'air5species'
    
    # Run estcj
    estc_result = rstc(gas,p1,T1,Vs,pe,None,area_ratio=ar,task='stn')
    
    states = ['state1','state2','state5','state5s','state6','state7']
    #descrip = {'Pre-shock','Post-shock','Reflected-shock','Equilibrium',
    #           'Nozzle Throat','Nozzle Exit'}
    exitVars = ['p','T','rho','e','h','a','s','R','gam','C_p','mu','k','species']

    # Now start writing the LaTeX table file. This is quite specific
    # to the stn calculation
    fp = open(FileToWrite,'w')
    fp.write('% This file requires the packages: siunitx, booktabs\n')
    fp.write('\\begin{table}\n')
    fp.write('\centering\n')
    fp.write('\caption{%s $p_1$: \SI{%g}{\pascal}, $T_1$: \SI{%g}{\kelvin}, $u_{ss}$: \SI{%g}{\metre\per\second}, area ratio: \\num{%g}}\n' % (Title, p1, T1, Vs, ar) )
    fp.write('\sisetup{detect-none,mode=math,table-align-exponent=false,table-format=4.5e+2}\n')
    fp.write('\\begin{tabular}{lS[table-format=6.4e+2]S[table-format=4.4e+2]SSSS[table-format=+6.4e+2]}\n')
    fp.write('\\toprule\n')
    fp.write(' & {State $1$} & {State $2$} & {State $5$} & {State $5$s} & {State $6$} & {State $7$} \\\\\n')
    fp.write(' & {Pre-shock} & {Post-shock} & {Reflected-shock} & {Equilibrium} & {Nozzle Throat} & {Nozzle Exit} \\\\\n')
    fp.write('\\midrule\n')
    # Write all data except species fractions (will do these separately)
    for eV in exitVars:
        # Write first column
        if eV in ['rho']:
            fp.write('$\\rho$ & ')
        elif eV in ['gam']:
            fp.write('$\\gamma$ & ')
        elif eV in ['mu']:
            fp.write('$\\mu$ & ')
        elif eV in ['species']:
            continue
        else:
            fp.write('${0:s}$ & '.format(eV))
        
        # Write data for each state to remaining columns
        for s in states:
            data = vars(estc_result[s])
            if eV not in ['species']:
                if eV in ['mu']:
                    fp.write('{0:1.4e}'.format(data[eV]))
                else:
                    fp.write('{0:g}'.format(data[eV])) 
            if s not in ['state7']:
                fp.write(' & ')
            else:
                fp.write(' \\\\\n')

    fp.write('\\midrule\n')
    
    # Now we deal with the species fractions
    speciesList = ['N2','O2','N','O','NO']
    for sp in speciesList:
        # Write the first column
        fp.write('$Y_{%s}$ & ' % (sp))
        
        # Write data for each state to remaining columns
        for s in states:
            speciesData = vars(estc_result[s])['species']
            if sp in speciesData.keys():
                if sp in ['O']:
                    fp.write('{0:1.4e}'.format(speciesData[sp]))
                else:
                    fp.write('{0:g}'.format(speciesData[sp]))
            else:
                fp.write('0')
            if s not in ['state7']:
                fp.write(' & ')
            else: 
                fp.write(' \\\\\n')

    fp.write('\\midrule\n')
    
    # Now write out velocity information
    #print estc_result.keys()
    # Flow velocities
    fp.write('$u$ & & {0:g} & & & {1:g} & {2:g}\\\\\n'.format(estc_result['V2'],\
               estc_result['V6'], estc_result['V7']))
    # Lab reference frame
    fp.write('$u_g$ & & {0:g} & & & & \\\\\n'.format(estc_result['Vg']))
    # Reflected shock speed
    fp.write('$u_{ss,R}$ & & & %g & & & \\\\\n'% estc_result['Vr'])
    # Pitot
    fp.write('$Pitot$ & & & & & & {0:g} \\\\\n'.format(estc_result['pitot7']))
    fp.write('\\bottomrule\n')
    fp.write('\\end{tabular}\n')
    fp.write('\\label{%s}\n' % Label)
    fp.write('\\end{table}\n')
    fp.close()



    
if __name__ == '__main__':
    return_flag = main()
    sys.exit(return_flag)
