#! /usr/bin/env python
"""
poshax_python.py:

This is a Python program to interface with Dan Potter's shock
relaxation code poshax.
I wrote this for people in the expansion tube laboratory to use
as the elevated temperature in the freestream conditions proposed
by the expansion tubes often required a bit of messing around to
run CEA manually and then get those equilibrium mass fractions
into poshax. It made messing almost impossible for us.

This code takes the inflow pressure, temperature, and velocity or Mach number,
and starting composition, as well as the information required for poshax,
runs CEA to get an equilibrium inflow, generates the poshax gas file
and cfg file, runs poshax, and then plots the results along with
comparison to the CEA equilibrium post shock condition.

Examples generally assume that this program has been placed in the e3bin folder,
and obviously, poshax needs to be installed.

Chris James (c.james4@uq.edu.au) - 25/04/18

"""

import sys, os, math, shutil, subprocess, copy
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
sys.path.append("") # so that we can find user's scripts in current directory
# We base our calculation of gas properties upon calls to the NASA Glenn CEA code.
from cfpylib.gasdyn.cea2_gas import Gas
from cfpylib.gasdyn.gas_flow import normal_shock

VERSION_STRING = '30-Apr-2018'

def poshax_python(p_inf, T_inf, reactants, inputUnits,
                  species_list, no_of_temperatures,
                  reaction_file, outputUnits = 'massf', u_inf = None, M_inf = None, energy_exchange_file = None,
                  gas_model_file = 'gas_model', cfg_file = 'cfg_file.cfg', output_file = 'output_file.data',
                  source_term_coupling = 'loose', dx = 1.0e-16, adaptive_dx = 'false',
                  final_x = 5.0e-3, plot = True, save_figures = True,
                  temp_result_fig_name = 'temp_result', species_result_fig_name = 'species_result',
                  freestream_normalised_result_fig_name = 'freestream_normalised_result',
                  eq_normalised_result_fig_name = 'equilibrium_normalised_result',
                  title = None, log_x_axis = False, plot_x_limits = None):
    """

    :param p_inf: Freestream pressure (Pa)
    :param T_inf: Freestream temperature (K)
    :param reactants: Reactants dictionary which will be sent to cfypylib's cea2 Gas object.
        So it must conform with the required syntax!! an example is {'N2':0.79,'O2':0.21}
        for 'cfd air' specified by moles, but check the cea2_gas info for more info about
        the required format. inputUnits for this dictionary can be either in 'moles' or 'massf'
        with the standard specified in the next input
    :param inputUnits: inputUnits string to be again set to the cea2 Gas object. Should be either
        'moles' or 'massf'
    :param species_list: Species list for POSHAX!! NOT CEA!! so ions and electons will be N2_plus, e_minus
        This data will be parsed to cea as the input state will only use these species, but the code
        will automatically perform the conversion need to send this list to CEA. This list is for
        sending to the program whcih creates the gas file for poshax
    :param no_of_temperatures: Number of temperatures. Integer input of 1,2, or 3. More than one
        temperature models require the specification of an energy exchange .lua file. (See below.)
    :param reaction_file: String pointing to the location of the reaction scheme.
    :param u_inf: Freestream velocity (m/s). Defaults to None, but either u_inf or M_inf must be specified.
    :param M_inf: Freestream Mach number. Defaults to None, but either u_inf or M_inf must be specified.
    :param energy_exchange_file: String pointing to the location of the energy exchange scheme for
        multi temperature models. Defaults to None, but the code will not run with a multi temperature
        model if this file is not found.
    :param gas_model_file: String name of the gas model file which the code will create as a .inp
        file which the program gasfile uses to make the .lua gasfile for poshax. Defaults to 'gas_model'.
    :param cfg_file: String name of the poshax cfg file. Defaults to 'cfg_file.cfg'.
    :param output_file: String name of the poshax output file. Defaults to 'output_file.data'.
    :param source_term_coupling: poshax source term coupling file. Defaults to 'loose'.
    :param dx: poshax dx (in m) defaults to 1.0e-16 m.
    :param adaptive_dx: Whether to use adaptive_dx with poshax. Defaults to 'true'.
    :param final_x: poshax final x value (in m). Defaults to 5.0e-3.
    :param plot: Whether to plot the result or not. Defaults to True. Will make four plots.
        A temperature plot, mass fraction plot of all species, result normalised by freestream values,
        and a result normalised by equilibrium values.
    :param save_figures: Whether to save figures or not. Defaults to True. results are saved in eps, png,
        and pdf formats.
    :param temp_result_fig_name: Temperature result plot file name. Defaults to 'temp_result'.
    :param species_result_fig_name: Species result plot file name. Defaults to 'species_result'.
    :param freestream_normalised_result_fig_name: Freestream normalised results plot file name.
        Defaults to 'freestream_normalised_result'.
    :param eq_normalised_result_fig_name: Equilibrium normalised results plot file name.
        Defaults to 'equilibrium_normalised_result'.
    :param title: allows a title to be added to the plots. Defaults to None.
        Is a single title so, for example, the condition being simulated can be added to the title.
    :param log_x_axis: Plots the figure with the x axis on a log scale. Defaults to False.
plot_x_limits = None
    :return: Returns a Python dictionary version of the poshax results file.
    """

    print '-'*60
    print "Running poshax python version {0}".format(VERSION_STRING)

    # do any input checking which is required
    # probably don't need to do everything in the world here,
    # but just some important things...
    # probably want to check that p, T, and u values,

    if not isinstance(p_inf, (float, int)):
        raise Exception, "poshax_python(): 'p_inf' input ({0}) is not a float or an int.".format(p_inf)
    if not isinstance(T_inf, (float, int)):
        raise Exception, "poshax_python(): 'T_inf' input ({0}) is not a float or an int.".format(T_inf)
    if not u_inf and not M_inf:
        raise Exception, "poshax_python(): Simulation must have either a 'u_inf' or 'M_inf' value."
    if u_inf and M_inf:
        raise Exception, "poshax_python(): Simulation cannot have have both a 'u_inf' and an 'M_inf' value."
    if u_inf and not isinstance(u_inf, (float, int)):
        raise Exception, "poshax_python(): 'u_inf' input ({0}) is not a float or an int.".format(u_inf)
    if M_inf and not isinstance(M_inf, (float, int)):
        raise Exception, "poshax_python(): 'M_inf' input ({0}) is not a float or an int.".format(M_inf)

    if not isinstance(reactants, dict):
        raise Exception, "poshax_python(): 'reactants' input ({0}) is not a dictionary.".format(reactants)
    if inputUnits not in ['moles', 'massf']:
        raise Exception, "poshax_python(): 'inputUits' input ({0}) is not either 'moles' or 'massf'.".format(inputUnits)

    if not isinstance(reaction_file, str):
        raise Exception, "poshax_python(): 'reaction_file' input ({0}) is not a string.".format(reaction_file)

    if not os.path.exists(reaction_file):
        raise Exception, "poshax_python(): 'reaction_file' input ({0}) does not appear to exist in the specified location.".format(reaction_file)

    if no_of_temperatures > 1 and not energy_exchange_file:
        raise Exception, "poshax_python(): Multi temperature model has been specified without an energy exchange file."

    if no_of_temperatures > 1:
        if not isinstance(energy_exchange_file, str):
            raise Exception, "poshax_python(): 'energy_exchange_file:' input ({0}) is not a string.".format(energy_exchange_file)

        if not os.path.exists(energy_exchange_file):
            raise Exception, "poshax_python(): 'reaction_file' input ({0}) does not appear to exist in the specified location.".format(energy_exchange_file)

    if not isinstance(dx, float):
        raise Exception, "poshax_python(): 'dx' input ({0}) is not a float.".format(dx)
    if not isinstance(final_x, float):
        raise Exception, "poshax_python(): 'final_x' input ({0}) is not a float.".format(final_x)


    print '-'*60
    print "User specified input values are:"

    if no_of_temperatures == 1:
        model = "thermally perfect gas"
    elif no_of_temperatures == 2:
        model = "two temperature gas"
    elif no_of_temperatures == 3:
        model = "three temperature gas"

    print "Gas state input settings:"
    print "p_inf = {0} Pa, T_inf = {1} K, u_inf = {2} m/s".format(p_inf, T_inf, u_inf)
    print "Reactants for start of CEA inflow calculation are:"
    print '{0} (by {1})'.format(reactants, inputUnits)
    print "Species in the calculation are:"
    print species_list
    print "Note: the CEA equilibrium inflow calculation will use only these species."
    print "This means that the inflow may be unrealistic if the species do not match "
    print "the real state of the gas at the user specified temperature and pressure"
    print '-'*60
    print "poshax settings:"
    print "dx = {0} m, final_x = {1} m ({2} mm), with adaptive_dx set to {3}".format(dx, final_x,
                                                                                     final_x*1000,
                                                                                     adaptive_dx)
    print "calculation uses a {0} model, with {1} source term coupling.".format(model, source_term_coupling)
    print "Reaction scheme file is {0}".format(reaction_file)
    if no_of_temperatures > 1:
        print "Energy exchange file is {0}.".format(energy_exchange_file)
    print "Gas model file to be created will be called {0}.lua.".format(gas_model_file)
    print "poshax input file to be created will be called {0}.".format(cfg_file)
    print "poshax output / results file to be created will be called {0}".format(output_file)

    # as cea wants a slightly different form from poshax...
    species_list_for_cea = []

    for species in species_list:
        if '_plus' in species:
            # copy the original species here so we don't mess with the old one!!
            species_copy = copy.copy(species)
            species_list_for_cea.append(species_copy.replace('_plus', '+'))
        elif '_minus' in species:
            # copy the original species here so we don't mess with the old one!!
            species_copy = copy.copy(species)
            species_list_for_cea.append(species_copy.replace('_minus', '-'))
        else:
            # just append the value...
            species_list_for_cea.append(species)

    print '-'*60
    print "Creating user specified equilibrium gas inflow using CEA:"

    # get mass fractions from CEA as poshax needs massf inflow
    state1 = Gas(reactants = reactants, with_ions = True,
                 inputUnits = inputUnits, outputUnits=outputUnits,
                 onlyList=species_list_for_cea)
    state1.set_pT(p_inf, T_inf) #Pa, K

    # print state to the screen...
    state1.write_state(sys.stdout)

    print '-'*60
    print "Preparing the gas .inp file from user specifications..."

    with open('{0}.inp'.format(gas_model_file), 'w') as gas_file:
        gas_file.write('model = "{0}"'.format(model) + '\n')

        species_line = 'species = {'

        for species in species_list:
            if species != species_list[-1]:
                species_line += "'{0}', ".format(species)
            else:
                species_line += "'{0}'".format(species)

        # now add the end bracket...
        species_line += '}'

        gas_file.write(species_line + '\n')

        gas_file.close()

    print "Running gasfile program to generate .lua gasfile."

    subprocess.check_call(["gasfile", "{0}.inp".format(gas_model_file), "{0}.lua".format(gas_model_file)])

    print "Gas file prepared successfully."

    # now we need to build the .cfg file...

    print '-'*60
    print "Building poshax .cfg file from user specifications"

    with open(cfg_file, 'w') as cfg_file_ojbect:
        # add in all of the controls...
        cfg_file_ojbect.write('[controls]' + '\n')
        cfg_file_ojbect.write('source_term_coupling = {0}'.format(source_term_coupling) + '\n')
        cfg_file_ojbect.write('dx = {0}'.format(dx) + '\n')
        if adaptive_dx: # remember that a value of 'false' would still be a string...
            cfg_file_ojbect.write('adaptive_dx = {0}'.format(adaptive_dx) + '\n')
        cfg_file_ojbect.write('final_x = {0}'.format(final_x) + '\n')
        cfg_file_ojbect.write('output_file = {0}'.format(output_file ) + '\n')
        cfg_file_ojbect.write('species_output = {0}'.format(outputUnits) + '\n')

        # the models...
        cfg_file_ojbect.write('[models]' + '\n')
        # a version of this without the .lua was defined above...
        cfg_file_ojbect.write('gas_model_file = {0}.lua'.format(gas_model_file) + '\n')
        cfg_file_ojbect.write('reaction_file = {0}'.format(reaction_file ) + '\n')
        if no_of_temperatures > 1:
            cfg_file_ojbect.write('energy_exchange_file = {0}'.format(energy_exchange_file) + '\n')

        # initial conditions...
        cfg_file_ojbect.write('[initial-conditions]' + '\n')
        cfg_file_ojbect.write('p_inf = {0}'.format(p_inf) + '\n')
        if no_of_temperatures == 1:
            cfg_file_ojbect.write('T_inf = {0}'.format(T_inf) + '\n')
        elif no_of_temperatures == 2:
            cfg_file_ojbect.write('T_inf = {0}, {0}'.format(T_inf) + '\n')
        elif no_of_temperatures == 3:
            cfg_file_ojbect.write('T_inf = {0}, {0}, {0}'.format(T_inf) + '\n')
        if u_inf:
            cfg_file_ojbect.write('u_inf = {0}'.format(u_inf) + '\n')
        if M_inf:
            cfg_file_ojbect.write('M_inf = {0}'.format(M_inf) + '\n')

        # now go through species and we'll make a species line

        if outputUnits == 'moles':
            species_line = 'molef_inf = '
        elif outputUnits == 'massf':
            species_line = 'massf_inf = '

        for species in species_list:
            if species in state1.species:
                if species != species_list[-1]:
                    species_line += '{0}, '.format(state1.species[species])
                else:
                    species_line += '{0}'.format(state1.species[species])
            else:
                if species != species_list[-1]:
                    species_line += '0.0, '
                else:
                    species_line += '0.0'

        cfg_file_ojbect.write(species_line + '\n')

        cfg_file_ojbect.close()

    print "Cfg file building successfully"

    print '-'*60
    print "Now running poshax program..."

    exit_code = subprocess.check_call(["poshax3.x", cfg_file])

    if exit_code == 0:
        print "poshax ran successfully"
    else:
        print "poshax appears to have had some issues. check your inputs!"

    print '-'*60
    print "Performing equilibrium normal shock calculation to compare with poshax result."

    state2 = state1.clone()

    if not u_inf:
        u_inf = M_inf*state1.a

    V2, V2g = normal_shock(state1, u_inf, state2)

    print "Post-shock equilibrium state is :"
    print "post-shock velocity = {0} m/s".format(V2)
    state2.write_state(sys.stdout)

    # now opening the result so it can be returned...

    no_of_species = len(species_list)
    species_start_line = 6 + no_of_temperatures # 2 intro lines, x, p, rho, u

    lines_to_skip = 6 + no_of_temperatures + no_of_species # 2 intro lines, x, p, rho, u

    # I wanted to call my variable output_file (I nromally called the string output_filename)
    # but wanred to try to keep my variables consistent with poshax... so now I have results_file!
    with open(output_file,"r") as results_file:

        results_dict = {}
        columns = []
        output_species_list = []

        # start by grabbing columns names which wek now
        columns.append('x')
        if no_of_temperatures == 1:
            T_list = ['T']
        elif no_of_temperatures == 2:
            T_list = ['T_trans_rotational', 'T_vibrational_electronic']
        elif no_of_temperatures == 3:
            T_list = ['T_trans_rotational', 'T_vibrational_electronic', 'T_electron']

        columns += T_list

        columns += ['p', 'rho', 'u']

        # we'll get the rest of the columns from species...

        for i, row in enumerate(results_file):

            if i >= species_start_line and i < lines_to_skip:
                # we could just get the species here, but good to keep it this way
                # so we can abstract the plotting into a function in the future...
                # (i.e. so it could run with no knowledge of the species beforehand, just the number of them...
                # here we want to split the line about the '-' so we can get the final value
                # need to use .strip() to get rid of the new line character on teh end too...
                split_row = row.strip().split('-')
                # last should be what we want...
                columns.append(split_row[-1])
                output_species_list.append(split_row[-1])

            elif i >= lines_to_skip: # to skip the header and know when the data starts...

                if i == lines_to_skip:
                    # we need to add our columns to the dictionary!
                    for column in columns:
                        results_dict[column] = []

                # now we start filling in data...

                split_row = row.split()

                # go through the split row and fill the data away in the right columns..
                for column, value in zip(columns, split_row):
                    results_dict[column].append(float(value))

    if plot:
        print '-'*60
        print "Now plotting the result"
        import matplotlib.pyplot as plt
        import numpy as np

        #---------------------------------------------------------------------------
        # temp plot
        fig, ax = plt.subplots()

        x_list = np.array(results_dict['x']) * 1000.0 # convert to mm

        for T in T_list:
            ax.plot(x_list, results_dict[T], label = T)

        # add eq temp from CEA...
        ax.plot([x_list[0], x_list[-1]], [state2.T, state2.T], label = 'eq T from CEA')

        ax.set_xlabel('post-shock distance (mm)')
        ax.set_ylabel('temperature (K)')

        if log_x_axis:
            ax.set_xscale('log')

        if plot_x_limits:
            ax.set_xlim(plot_x_limits)

        plt.legend(loc = 'best')

        if title:
            plt.title(title)

        if save_figures:
            plt.savefig(temp_result_fig_name + '.png', dpi = 300)
            plt.savefig(temp_result_fig_name + '.eps', dpi = 300)
            plt.savefig(temp_result_fig_name + '.pdf', dpi = 300)

        plt.show()

        #-------------------------------------------------------------------------
        # now do massf plot
        species_fig, species_ax = plt.subplots()

        for i, species in enumerate(output_species_list):
            # use solid lines if we are line 7 or less
            # change over for lines past this
            # this gives maximum 12 species, I guess, before overlap
            # works with earlier version of matplotlib will only 7 colours too..
            if i < 7:
                linestyle = '-'
            else:
                linestyle = '--'

            # change the _plus or _minus to + or - here to make the labels nicer...
            if '_plus' in species:
                species_label = species.replace('_plus', '+')
            elif '_minus' in species:
                species_label = species.replace('_minus', '-')
            else:
                species_label = species
            species_ax.plot(x_list, results_dict[species], linestyle = linestyle, label = species_label)

        species_ax.set_xlabel('post-shock distance (mm)')
        if outputUnits == 'moles':
            species_ax.set_ylabel('moles per original mole')
        elif outputUnits == 'massf':
            species_ax.set_ylabel('mass fraction of species')

        if log_x_axis:
            species_ax.set_xscale('log')

        if plot_x_limits:
            species_ax.set_xlim(plot_x_limits)

        plt.legend(loc = 'best')

        if title:
            plt.title(title)

        if save_figures:
            plt.savefig(species_result_fig_name + '.png', dpi = 300)
            plt.savefig(species_result_fig_name + '.eps', dpi = 300)
            plt.savefig(species_result_fig_name + '.pdf', dpi = 300)

        plt.show()

        # -------------------------------------------------------------------
        # now do the freestream normalised plot
        freestream_normalised_fig, freestream_normalised_ax = plt.subplots()

        for T in T_list:
            freestream_normalised_ax.plot(x_list, np.array(results_dict[T])/T_inf, label = T)

        # now do p, rho, and u
        # I commented out pressure here as it changes by too much and messes with the plot
        #freestream_normalised_ax.plot(x_list, np.array(results_dict['p'])/p_inf, label = 'p')
        # have to get rho from state 1 as we didn't need to specify it for inflow calculation...
        freestream_normalised_ax.plot(x_list, np.array(results_dict['rho'])/state1.rho, label = 'rho')
        freestream_normalised_ax.plot(x_list, np.array(results_dict['u'])/u_inf, label = 'u')

        # add eq values from CEA
        freestream_normalised_ax.plot([x_list[0], x_list[-1]], [state2.T/T_inf, state2.T/T_inf], label = 'eq T from CEA')
        freestream_normalised_ax.plot([x_list[0], x_list[-1]], [state2.rho/state1.rho, state2.rho/state1.rho], label = 'eq rho from CEA')
        freestream_normalised_ax.plot([x_list[0], x_list[-1]], [V2/u_inf, V2/u_inf], label = 'eq u from CEA')


        freestream_normalised_ax.set_xlabel('post-shock distance (mm)')
        freestream_normalised_ax.set_ylabel('quantity normalised by freestream value')

        if log_x_axis:
            freestream_normalised_ax.set_xscale('log')

        if plot_x_limits:
            freestream_normalised_ax.set_xlim(plot_x_limits)

        plt.legend(loc = 'best')

        if title:
            plt.title(title)

        if save_figures:
            plt.savefig(freestream_normalised_result_fig_name + '.png', dpi = 300)
            plt.savefig(freestream_normalised_result_fig_name + '.eps', dpi = 300)
            plt.savefig(freestream_normalised_result_fig_name + '.pdf', dpi = 300)

        plt.show()

        #------------------------------------------------------------------------
        # now do the equilibrium normalised plot
        eq_normalised_fig, eq_normalised_ax = plt.subplots()

        for T in T_list:
            eq_normalised_ax.plot(x_list, np.array(results_dict[T])/state2.T, label = T)

        # now do p, rho, and u
        eq_normalised_ax.plot(x_list, np.array(results_dict['p'])/state2.p, label = 'p')
        eq_normalised_ax.plot(x_list, np.array(results_dict['rho'])/state2.rho, label = 'rho')
        eq_normalised_ax.plot(x_list, np.array(results_dict['u'])/V2, label = 'u')

        eq_normalised_ax.set_xlabel('post-shock distance (mm)')
        eq_normalised_ax.set_ylabel('quantity normalised by equilibrium value from CEA')

        if log_x_axis:
            eq_normalised_ax.set_xscale('log')

        if plot_x_limits:
            eq_normalised_ax.set_xlim(plot_x_limits)

        plt.legend(loc = 'best')

        if title:
            plt.title(title)

        if save_figures:
            plt.savefig(eq_normalised_result_fig_name + '.png', dpi = 300)
            plt.savefig(eq_normalised_result_fig_name + '.eps', dpi = 300)
            plt.savefig(eq_normalised_result_fig_name + '.pdf', dpi = 300)

        plt.show()

    return results_dict

if __name__=='__main__':

    # This is Rowan's Air Mach12-3 single temperature example...
    title = 'Air Mach 12.28'
    p_inf = 133.322  # Pa
    T_inf = 300.0  # K
    M_inf = 12.28

    reactants = {'N2': 0.79, 'O2': 0.21}  # reactant dict for cea calc
    inputUnits = 'moles'  # for cea
    outputUnits = 'moles' # for poshax output

    species_list = ['O2', 'N2', 'NO', 'O', 'N', 'NO+', 'e-']

    no_of_temperatures = 1

    # name of the .inp gas file which we will make...
    # and the .lua gas file which the gasfile program will make
    # don't add the .lua here as we have two different files to deal with here
    gas_model_file = 'air-model'

    # these are both lua files so add the .lua
    # reaction file is always needed
    reaction_file = 'gupta_etal_with_Keq.lua'
    # energy exchange file is only needed for multiple temperatures
    energy_exchange_file = None

    # cfg file is the file whcih we'll send to poshax to run the simulation
    cfg_file = 'cfg_file.cfg'

    # output_file is the result table with distance from the shock
    output_file = 'output_file.data'

    # simulation values for poshax
    source_term_coupling = 'loose'
    dx = 1.0e-10
    adaptive_dx = None
    final_x = 30.0e-2

    plot = True
    save_figures = True

    poshax_python(p_inf = p_inf, T_inf = T_inf, M_inf = M_inf,
                  reactants = reactants, inputUnits = inputUnits, outputUnits = outputUnits,
                  species_list = species_list, no_of_temperatures = no_of_temperatures,
                  reaction_file = reaction_file, energy_exchange_file = energy_exchange_file,
                  gas_model_file = gas_model_file, cfg_file = cfg_file,
                  source_term_coupling = source_term_coupling, dx = dx, adaptive_dx = adaptive_dx,
                  final_x = final_x, plot = plot, save_figures = save_figures, title = title,
                  log_x_axis = True, plot_x_limits= [0.01,100])