#! /usr/bin/env python
"""
pitot_multiple.py: program to easily run a few different runs of pitot

This file will take a normal pitot input file where the pressure inputs
are lists of floats instead of single floats. The program will automatically 
populate various folders with all of the tests requested.

I was getting sick of doing this kind of stuff manually, so I made this/

Chris James (c.james4@uq.edu.au) - 07/04/14

"""

VERSION_STRING = "07-Apr-2014"


import sys

from pitot import run_pitot
from pitot_input_utils import *

def run_pitot_multiple(cfg = {}, config_file = None):
    """
    
    Chris James (c.james4@uq.edu.au) 07/04/14
    
    run_pitot_multiple(dict, bool) - > None
    
    """
    
    import copy
    import os
   
    #---------------------- get the inputs set up --------------------------
    
    if config_file:
        cfg = config_loader(config_file)
        
    #----------------- check inputs ----------------------------------------
    
    # add new value to the cfg file so it knows to not bail when it sees
    # that the pressure inputs are lists not floats
    cfg['pitot_multiple'] = True
    cfg = input_checker(cfg)
    
    #turn clean up on here if they haven't, as we don't want to be leaving
    # cea temp files in a heap of folders
    cfg['cleanup'] = True
    
    # store the pressure input lists
    
    if cfg['secondary']:
        psd1_list = copy.deepcopy(cfg['psd1'])
    if cfg['tunnel_mode'] == 'expansion-tube':
        p1_list = copy.deepcopy(cfg['p1'])
        p5_list = copy.deepcopy(cfg['p5'])
    elif cfg['tunnel_mode'] == 'nr-shock-tunnel' or cfg['tunnel_mode'] == 'reflected-shock-tunnel':
        p1_list = copy.deepcopy(cfg['p1'])
        
    # check they're all the same length here too, or it could cause the user
    # some issues further down the track!
    
    print "Checking all the input lists are the same length."
    
    if cfg['secondary'] and cfg['tunnel_mode'] == 'expansion-tube':
        if len(psd1_list) == len(p1_list) == len(p5_list):
            print "All lists are the same length. Continuing..."
            print '-'*60
        else:
            print "One list is not the same length as the others. Bailing out."
            raise Exception, "pitot.run_pitot_multiple() Pressure lists are different lengths. Check your input."
    elif not cfg['secondary'] and cfg['tunnel_mode'] == 'expansion-tube':
        if len(p1_list) == len(p5_list):
            print "Both lists are the same length. Continuing..."
            print '-'*60
        else:
            print "One list is not the same length as the other. Bailing out."
            raise Exception, "pitot.run_pitot_multiple() Pressure lists are different lengths. Check your input."
    elif cfg['secondary'] and cfg['tunnel_mode'] == 'nr-shock-tunnel' or cfg['secondary'] and cfg['tunnel_mode'] == 'reflected-shock-tunnel':
        if len(psd1_list) == len(p1_list):
            print "Both lists are the same length. Continuing..."
            print '-'*60
        else:
            print "One list is not the same length as the other. Bailing out."
            raise Exception, "pitot.run_pitot_multiple() Pressure lists are different lengths. Check your input."
        
    # take a copy of the initial cfg
    
    initial_cfg = copy.deepcopy(cfg)
    
    # now get running
    
    for i in range(0, len(p1_list)):
        # make a new cfg
        cfg = copy.deepcopy(initial_cfg)
        
        # set the correct pressure values and set our correct filename
        filename = ''
        if cfg['secondary']:
            cfg['psd1'] = psd1_list[i]
            filename = filename + 'psd1-{0}-'.format(int(cfg['psd1']))
        if cfg['tunnel_mode'] == 'expansion-tube':
            cfg['p1'] = p1_list[i]
            cfg['p5'] = p5_list[i]
            filename = filename + 'p1-{0}-p5-{1}'.format(int(cfg['p1']), int(cfg['p5']))
        elif cfg['tunnel_mode'] == 'nr-shock-tunnel' or cfg['tunnel_mode'] == 'reflected-shock-tunnel':
            cfg['p1'] = p1_list[i]
            filename = filename + 'p1-{0}'.format(int(cfg['p1']))
            
        cfg['filename'] = filename
        
        print "Now running test {0} of {1}.".format(i+1, len(p1_list))
        
        # now move ourselves into a new folder and work from there
        try:
            os.mkdir(cfg['filename'])
            os.chdir(cfg['filename'])
            change_directory_back = True
            run_test = True
        except Exception as e:
            print '-'*60
            print "Error {0}".format(str(e))
            print "Cannot run this test as the folder is already there."
            print "If you want to re-run this test. Delete the folder and try again."
            print '-'*60
            change_directory_back = False
            run_test = False
           
        try:
            if run_test:
                cfg, states, V, M = run_pitot(cfg)
        except Exception as e:
            print "Error {0}".format(str(e))
            print "Test failed. Moving onto the next one..."
        
        # now move back into the original directory ready to start again
        
        if change_directory_back:
            os.chdir('../')
    
    return
                                
#----------------------------------------------------------------------------

def main():
    
    import optparse  
    op = optparse.OptionParser(version=VERSION_STRING)   
    op.add_option('-c', '--config_file', dest='config_file',
                  help=("filename where the configuration file is located"))    

    opt, args = op.parse_args()
    config_file = opt.config_file
           
    run_pitot_multiple(cfg = {}, config_file = config_file)
    
    return
    
#----------------------------------------------------------------------------

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print "pitot_multiple.py - Pitot Equilibrium expansion tube simulator easy multi run tool"
        print "start with --help for help with inputs"
        
    else:
        main()
