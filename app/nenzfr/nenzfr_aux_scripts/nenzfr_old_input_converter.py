#!/usr/bin/env python
"""
Program that accepts the old nenzfr inputs and then writes a cfg file for the 
new input.

Chris James (c.james4@uq.edu.au) - 03-May-13
"""

VERSION_STRING = "03-May-2013"

from string import upper
import sys, os
import optparse
E3BIN = os.path.expandvars("$HOME/e3bin")
sys.path.append(E3BIN)

def main():
    """
    Examine the command-line options to decide the what to do
    and then coordinate the calculations done by estcj and Eilmer3.
    """

    op = optparse.OptionParser(version=VERSION_STRING)
    op.add_option('--facility', dest='facility', default='reflected-shock-tunnel',
                  choices=['reflected-shock-tunnel','expansion-tube','gun-tunnel'],
                  help=("type of facility: "
                        "reflected-shock-tunnel ; " 
                        "gun-tunnel ;"
                        "expansion-tube ; "
                        " [default: %default]"))   
    op.add_option('--gas', dest='gasName', default='air',
                  choices=['air', 'air5species', 'n2', 'co2', 'h2ne'],
                  help=("name of gas model: "
                        "air; " "air5species; " "n2; " "co2; " "h2ne"
                        " [default: %default]"))
    op.add_option('--p1', dest='p1', type='float', default=None,
                  help=("shock tube fill pressure or static pressure, in Pa;"
                        "used for reflected shock tunnel facilities"))
    op.add_option('--T1', dest='T1', type='float', default=None,
                  help=("shock tube fill temperature, in degrees K;"
                        "used for reflected shock tunnel facilities"))
    op.add_option('--Vs', dest='Vs', type='float', default=None,
                  help=("incident shock speed, in m/s;"
                        "used for reflected shock tunnel facilities"))
    op.add_option('--pe', dest='pe', type='float', default=None,
                  help=("equilibrium pressure (after shock reflection), in Pa;"
                        "used for reflected shock tunnel facilities"))
    op.add_option('--p7', dest='p7', type='float', default=None,
                  help=("unsteadily expanded test gas pressure, in Pa;"
                        "used for expansion tube facilities"))
    op.add_option('--T7', dest='T7', type='float', default=None,
                  help=("unsteadility expanded test gas temperature, in degrees K;"
                        "used for expansion tube facilities"))
    op.add_option('--V7', dest='V7', type='float', default=None,
                  help=("unsteadily expanded test gas velocity, in m/s;"
                        "used for expansion tube facilities"))
    op.add_option('--p0', dest='p0', type='float', default=None,
                  help=("stagnation pressure, in Pa"
                        "used for gun tunnel facilities"))
    op.add_option('--T0', dest='T0', type='float', default=None,
                  help=("stagnation temperature, in degrees K"
                        "used for gun tunnel facilities"))                        
    op.add_option('--chem', dest='chemModel', default='eq',
                  choices=['eq', 'neq', 'frz', 'frz2'], 
                  help=("chemistry model: " "eq=equilibrium; " 
                        "neq=non-equilibrium; " "frz=frozen " 
                        "[default: %default]"))
    op.add_option('--area', dest='areaRatio', default=1581.165,
                  help=("nozzle area ratio. only used for estcj calc. "
                        "use when --cfile(--gfile) are "
                        "specified. [default: %default]"))
    op.add_option('--job', dest='jobName', default='nozzle',
                  help="base name for Eilmer3 files [default: %default]")
    op.add_option('--cfile', dest='contourFileName', 
                  default='Bezier-control-pts-t4-m10.data',
                  help="file containing Bezier control points "
                       "for nozzle contour [default: %default]")
    op.add_option('--gfile', dest='gridFileName', default='None',
                  help="file containing nozzle grid. "
                  "overrides --cfile if both are given "
                  "[default: %default]")
    op.add_option('--exitfile', dest='exitSliceFileName', 
                  default='nozzle-exit.data',
                  help="file for holding the nozzle-exit data [default: %default]")
    op.add_option('--just-stats', dest='justStats', action='store_true', 
                  default=False,
                  help="skip the detailed calculations and "
                  "just retrieve exit-flow statistics")
    op.add_option('--block-marching', dest='blockMarching', action='store_true', 
                  default=False, help="run nenzfr in block-marching mode")
    # The following defaults suit Luke's Mach 10 calculations.
    op.add_option('--nni', dest='nni', type='int', default=1800,
                  help=("number of axial cells [default: %default]"))
    op.add_option('--nnj', dest='nnj', type='int', default=100,
                  help=("number of radial cells [default: %default]"))
    op.add_option('--nbi', dest='nbi', type='int', default=180,
                  help=("number of axial blocks for the divergence section (nozzle_blk) "
                        "[default: %default]"))
    op.add_option('--nbj', dest='nbj', type='int', default=1,
                  help=("number of radial blocks [default: %default]"))
    op.add_option('--bx', dest='bx', type='float', default=1.10,
                  help=("clustering in the axial direction [default: %default]"))
    op.add_option('--by', dest='by', type='float', default=1.002,
                  help=("clustering in the radial direction [default: %default]"))
    op.add_option('--max-time', dest='max_time', type='float', default=6.0e-3,
                  help=("overall simulation time for nozzle flow [default: %default]"))
    op.add_option('--max-step', dest='max_step', type='int', default=800000,
                  help=("maximum simulation steps allowed [default: %default]"))
    #
    op.add_option('--Twall', dest='Tw', type='float', default=300.0,
                  help=("Nozzle wall temperature, in K "
                        "[default: %default]"))
    op.add_option('--BLTrans', dest='BLTrans', default="x_c[-1]*1.1",
                  help=("Transition location for the Boundary layer. Used "
                        "to define the turbulent portion of the nozzle. "
                        "[default: >nozzle length i.e. laminar nozzle]"))
    op.add_option('--TurbVisRatio', dest='TurbVisRatio', type='float',
                  default=100.0, help=("Turbulent to Laminar Viscosity Ratio "
                  "[default: %default]"))
    op.add_option('--TurbIntensity', dest='TurbInten', type='float', 
                  default=0.05, help=("Turbulence intensity at the throat "
                  "[default: %default]"))
    op.add_option('--CoreRadiusFraction', dest='coreRfraction', type='float',
                  default=2.0/3.0, help=("Radius of core flow as a fraction of "
                  "the nozzle exit radius [default: %default]"))
    opt, args = op.parse_args()
    
    config_file = open('nenzfr_inputs.cfg',"w")  #txt_output file creation
    print "Opening file 'nenzfr_inputs.cfg'."
    
    print "Starting the process of adding inputs to the config file."

    facility_string = "facility = '{0}'".format(opt.facility)
    config_file.write(facility_string+'\n')
    
    empty_line_string = " \n"
    
    config_file.write(empty_line_string)
    
    gasName_string = "gasName = '{0}'".format(opt.gasName)
    config_file.write(gasName_string+'\n')
    
    config_file.write(empty_line_string)

    chemModel_string = "chemModel = '{0}'".format(opt.chemModel)
    config_file.write(chemModel_string+'\n')
    
    config_file.write(empty_line_string)    

    p1_string = "p1 = {0}".format(opt.p1)
    config_file.write(p1_string+'\n')    

    T1_string = "T1 = {0}".format(opt.T1)
    config_file.write(T1_string+'\n')     

    Vs_string = "Vs = {0}".format(opt.Vs)
    config_file.write(Vs_string+'\n')

    pe_string = "pe = {0}".format(opt.pe)
    config_file.write(pe_string+'\n')
    
    config_file.write(empty_line_string)   

    areaRatio_string = "areaRatio = {0}".format(opt.areaRatio)
    config_file.write(areaRatio_string+'\n')

    config_file.write(empty_line_string)
    
    jobName_string = "jobName = '{0}'".format(opt.jobName)
    config_file.write(jobName_string+'\n')

    config_file.write(empty_line_string)
    
    contourFileName_string = "contourFileName = '{0}'".format(opt.contourFileName)
    config_file.write(contourFileName_string+'\n')
    
    config_file.write(empty_line_string)

    gridFileName_string = "gridFileName = '{0}'".format(opt.gridFileName)
    config_file.write(gridFileName_string+'\n')      

    config_file.write(empty_line_string)

    exitSliceFileName_string = "exitSliceFileName = '{0}'".format(opt.exitSliceFileName)
    config_file.write(exitSliceFileName_string+'\n') 

    config_file.write(empty_line_string)

    justStats_string = "justStats = {0}".format(opt.justStats)
    config_file.write(justStats_string+'\n')
    
    config_file.write(empty_line_string)

    blockMarching_string = "blockMarching = {0}".format(opt.blockMarching)
    config_file.write(blockMarching_string+'\n')
    
    config_file.write(empty_line_string)

    nni_string = "nni = {0}".format(opt.nni)
    config_file.write(nni_string+'\n')
    
    nnj_string = "nnj = {0}".format(opt.nnj)
    config_file.write(nnj_string+'\n')

    config_file.write(empty_line_string)    

    nbi_string = "nbi = {0}".format(opt.nbi)
    config_file.write(nbi_string+'\n')
    
    nbj_string = "nbj = {0}".format(opt.nbj)
    config_file.write(nbj_string+'\n')
    
    config_file.write(empty_line_string)
    
    bx_string = "bx = {0}".format(opt.bx)
    config_file.write(bx_string+'\n')
    
    by_string = "by = {0}".format(opt.by)
    config_file.write(by_string+'\n')

    config_file.write(empty_line_string)
    
    max_time_string = "max_time = {0}".format(opt.max_time)
    config_file.write(max_time_string+'\n')  
    
    config_file.write(empty_line_string)
    
    max_step_string = "max_step = {0}".format(opt.max_step)
    config_file.write(max_step_string+'\n')
    
    config_file.write(empty_line_string)
    
    Tw_string = "Tw = {0}".format(opt.Tw)
    config_file.write(Tw_string+'\n')

    config_file.write(empty_line_string)
    
    BLTrans_string = "BLTrans = '{0}'".format(opt.BLTrans)
    config_file.write(BLTrans_string+'\n')   
    
    config_file.write(empty_line_string)
    
    TurbVisRatio_string = "TurbVisRatio = {0}".format(opt.TurbVisRatio)
    config_file.write(TurbVisRatio_string+'\n')  
    
    config_file.write(empty_line_string)
    
    TurbInten_string = "TurbInten = {0}".format(opt.TurbInten)
    config_file.write(TurbInten_string+'\n')  
    
    config_file.write(empty_line_string)
    
    coreRfraction_string = "coreRfraction = {0}".format(opt.coreRfraction)
    config_file.write(coreRfraction_string+'\n')  

    config_file.close()                         
    
if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print "Program to change old nenzfr inputs to the new one"
        print "   Version:", VERSION_STRING
        print "   To get some useful hints, invoke the program with option --help."
        sys.exit(0)
    return_flag = main()
    sys.exit(return_flag)
