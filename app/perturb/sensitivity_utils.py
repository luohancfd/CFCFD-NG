def perturb_CoreRadiusFraction(var, perturbedVariables,\
                               DictOfCases, levels):
    """
    Perturbation of the CoreRadiusFraction may be completed
    as a post-processing step using the nominal condition
    flow solution. We don't need to compute separate solutions.

    This function copies the nominal condition solution to a 
    new case file and then runs nenzfr with the --just-stats
    option and the perturbed CoreRadiusFraction value. 
    """

    for kk in range(levels):
        if kk != 0:
            caseString = 'case'+"{0:02}".format(\
                 perturbedVariables.index(var))+\
                 "{0:01}".format(kk)
            # Only if required to we proceed with copying the
            # nominal case data and calling 'nenzfr.py'
            if not os.path.exists(caseString):
                print "Perturbing "+var
                #command_text = 'rm -r '+caseString
                #run_command(command_text)

                # Copy nominal case data into a new case
                # directory
                command_text = 'cp -r case000 '+caseString
                run_command(command_text)
                # Change into new directory and call nenzfr --just-stats
                os.chdir(caseString)
                command_text = 'nenzfr.py --just-stats --CoreRadiusFraction='+\
                  str(DictOfCases[caseString][perturbedVariables.index(var)])
                run_command(command_text)
                os.chdir('../')

def get_values(dict, propertyList):
    """
    Adapted from Peter Jacob's function "get_fractions" which
    may be found within "cea2_gas.py".
    """
    valueList = []
    #print propertyList
    for s in propertyList:
        #print s
        if s in dict.keys():
            #print s
            #print dict[s]
            valueList.append(dict[s])
        else:
            valueList.append(0.0)
            print "WARNING: "+s+"was not found in the current case dictionary."
    return valueList

def add_extra_variables(data, var):
    """
    We add some additional variables that of are interest to 
    the input dictionary and list.

    data: dictionary of data
    var: list of variables
    """
    U_squared = data['vel.x']*data['vel.x'] + data['vel.y']*data['vel.y']
    U = U_squared**0.5
    # Dynamic Pressure
    data['q'] = 0.5 * data['rho'] * U_squared
    var.append('q')
    # Mass flow rate per unit area
    data['m_dot'] = data['rho'] * U
    var.append('m_dot')
    # Unit Reynolds number
    data['Re_u'] = data['m_dot'] / data['mu']
    var.append('Re_u')
    # Pressure coefficient
    data['p/q'] = data['p'] / data['q']
    var.append('p/q')

    return data, var

def write_sensitivity_summary(sensitivity, perturbedVariables, exitVar, type):
    """
    Write out a file summarising the sensitivity of each exit flow 
    parameter to each of the input parameters.
    """
    titleFormatDict = {'p1':'{0:{fill}>13}', 'T1':'{0:{fill}>13}',
                       'Vs':'{0:{fill}>13}', 'pe':'{0:{fill}>13}',
                       'Tw':'{0:{fill}>13}', 'BLTrans':'{0:{fill}>13}',
                       'TurbVisRatio':'{0:{fill}>14}',
                       'TurbInten':'{0:{fill}>13}',
                       'CoreRadiusFraction':'{0:{fill}>20}'}
    formatDict = {'p1':'{0:13.5g}', 'T1':'{0:>13.5g}',
                  'Vs':'{0:>13.5g}', 'pe':'{0:>13.5g}',
                  'Tw':'{0:>13.5g}', 'BLTrans':'{0:>13.5g}',
                  'TurbVisRatio':'{0:>14.5g}',
                  'TurbInten':'{0:>13.5g}',
                  'CoreRadiusFraction':'{0:>20.5g}'}
    if type in ['relative']:
        fout = open('sensitivities.dat','w')
    elif type in ['absolute']:
        fout = open('sensitivities_abs.dat','w')

    # Write header information
    fout.write('{0:}\n'.format(type+' sensitivities'))
    fout.write('{0:>13}'.format('variable'))
    for k in perturbedVariables:
        fout.write(titleFormatDict[k].format(k,fill=''))
    fout.write('\n')
    for k in perturbedVariables:
        fout.write(titleFormatDict[k].format('-',fill='-'))
    fout.write('{0:->13}'.format('-'))
    fout.write('\n')

    # Now write out all the data
    for k in exitVar:
        fout.write('{0:>13}'.format(k))
        for kk in perturbedVariables:
            X_kk = sensitivity[kk][exitVar.index(k)]
            fout.write(formatDict[kk].format(X_kk))
        fout.write('\n')

    fout.close()

def write_uncertainty_summary(uncertainty, perturbedVariables, \
                              exitVar, inputUncertainties):
    """
    Write out a file summarising the sensitivity of each exit flow
    parameter to each of the input parameters.
    """
    titleFormatDict = {'p1':'{0:{fill}>10}', 'T1':'{0:{fill}>10}',
                       'Vs':'{0:{fill}>10}', 'pe':'{0:{fill}>10}',
                       'Tw':'{0:{fill}>10}', 'BLTrans':'{0:{fill}>10}',
                       'TurbVisRatio':'{0:{fill}>14}',
                       'TurbInten':'{0:{fill}>11}',
                       'CoreRadiusFraction':'{0:{fill}>20}'}
    formatDict = {'p1':'{0:10.2%}', 'T1':'{0:>10.2%}',
                  'Vs':'{0:>10.2%}', 'pe':'{0:>10.2%}',
                  'Tw':'{0:>10.2%}', 'BLTrans':'{0:>10.2%}',
                  'TurbVisRatio':'{0:>14.2%}',
                  'TurbInten':'{0:>11.2%}',
                  'CoreRadiusFraction':'{0:>20.2%}'}

    fout = open('uncertainties.dat','w')

    # Write header information
    fout.write('{0:>19}'.format('variable'))
    for k in perturbedVariables:
        fout.write(titleFormatDict[k].format(k,fill=''))
    fout.write('{0:>10}'.format('total'))
    fout.write('\n')
    # Write out the uncertainty for each perturbed variable
    fout.write('{0:>19}'.format('(input uncertainty)'))
    for k in perturbedVariables:
        fout.write(formatDict[k].format(inputUncertainties[k]))
    fout.write('{0:>10}'.format(''))
    fout.write('\n')
    for k in perturbedVariables:
        fout.write(titleFormatDict[k].format('-',fill='-'))
    fout.write('{0:->29}'.format('-'))
    fout.write('\n')
    
    # Now write out all the data. We also calculate and write out the 
    # total uncertainty in each exit flow variable due to the contributions
    # from all input values.
    for k in exitVar:
        fout.write('{0:>19}'.format(k))
        X_total = 0.0
        for kk in perturbedVariables:
            X_kk = uncertainty[kk][exitVar.index(k)]
            fout.write(formatDict[kk].format(X_kk))
            X_total += X_kk**2
        fout.write('{0:10.2%}'.format(X_total**0.5))
        fout.write('\n')

    fout.close()


