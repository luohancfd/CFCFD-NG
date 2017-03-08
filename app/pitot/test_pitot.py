#!/usr/bin/env python
"""
test_pitot.py -- test script for the pitot program

MUST BE RAN FROM THE SAME FOLDER AS THE PROGRAM

Chris James - 22-Apr-2013

Version: 
    08-Mar-2017: I decided to start again with this and build it back up.
        I want it to run in not too much time, but also test everything.

"""

import sys, os
sys.path.append(os.path.expandvars("$HOME/e3bin"))
import unittest
from pitot import run_pitot


class TestPitot(unittest.TestCase):
 
    def test_air_pressures_specified_eq_example(self):
        print " "
        cfg = {'solver':'eq', 'test':'fulltheory-pressure', 
               'facility':'x2', 'nozzle':True, 'secondary': False, 'piston':'lwp-2mm-new-paper',
               'driver_gas':'He:0.8,Ar:0.2', 'test_gas':'air',
               'p1':3000.0,'p5':10.0, 'area_ratio':5.64, 'expand_to':'shock-speed',
               'mode':'return','shock_over_model':True,'conehead':True, 'cone_angle':15.0}
        print '-'*60
        print "Running equilibrium air test with pressures specified."
        cfg, states, V, M = run_pitot(cfg=cfg)
        #check fill pressures
        self.assertAlmostEqual(states['s1'].p, 3000.0, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 10.0, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(cfg['Vs1'], 4196.99223582, delta=1.0)        
        self.assertAlmostEqual(cfg['Vs2'], 8437.66523168, delta=1.0)
        #check expanded driver gas condition
        self.assertAlmostEqual(states['s3'].p, 551780.00 , delta=1.0)  
        self.assertAlmostEqual(states['s3'].T, 608.94, delta=1.0)  
        self.assertAlmostEqual(V['s3'], 3778.14, delta=1.0)
        #check shocked test gas condition
        self.assertAlmostEqual(states['s2'].p, 558830.00, delta=1.0)  
        self.assertAlmostEqual(states['s2'].T, 4751.53, delta=1.0)  
        self.assertAlmostEqual(V['s2'], 3778.14, delta=1.0)
        # check shocked test gas condition
        self.assertAlmostEqual(states['s6'].p, 7869.00, delta=1.0)  
        self.assertAlmostEqual(states['s6'].T, 6864.85, delta=1.0)  
        self.assertAlmostEqual(V['s6'], 7972.63, delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p, 4468.00, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T, 2647.63, delta=1.0)  
        self.assertAlmostEqual(V['s7'], 8437.67, delta=1.0)
        #check nozzle exit conditions      
        self.assertAlmostEqual(states['s8'].p, 562.07, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 1966.46, delta=1.0)  
        self.assertAlmostEqual(V['s8'], 8601.79, delta=1.0)
        self.assertAlmostEqual(cfg['stagnation_enthalpy'], 38956333.8, delta=1.0)
        #check pg shock over model   
        self.assertAlmostEqual(states['s10f'].p, 65451.20, delta=1.0)  
        self.assertAlmostEqual(states['s10f'].T, 26956.80, delta=1.0)  
        self.assertAlmostEqual(V['s10f'], 7589.48, delta=1.0)
        #check eq shock over model        
        self.assertAlmostEqual(states['s10e'].p, 69449.00, delta=1.0)  
        self.assertAlmostEqual(states['s10e'].T, 8031.63, delta=1.0)  
        self.assertAlmostEqual(V['s10e'], 8057.78, delta=1.0)
        #check 15 degree cone calc 
        self.assertAlmostEqual(states['s10c'].p, 5799.00, delta=1.0)  
        self.assertAlmostEqual(states['s10c'].T, 3023.75, delta=1.0)  
        self.assertAlmostEqual(V['s10c'], 8256.20, delta=1.0) 
        return 
        
    def test_air_pressures_specified_pg_example(self):
        print " "
        cfg = {'solver':'pg', 'test':'fulltheory-pressure', 
               'facility':'x2', 'nozzle':True, 'secondary': False, 'piston':'lwp-2mm-new-paper',
               'driver_gas':'He:0.8,Ar:0.2', 'test_gas':'air',
               'p1':3000.0,'p5':10.0, 'area_ratio':5.64, 'expand_to':'shock-speed',
               'mode':'return','shock_over_model':True,'conehead':True, 'cone_angle':15.0}
        print '-'*60
        print "Running perfect gas air test with pressures specified."
        cfg, states, V, M = run_pitot(cfg=cfg)
        #check fill pressures
        self.assertAlmostEqual(states['s1'].p, 3000.0, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 10.0, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(cfg['Vs1'], 4519.34147214, delta=1.0)        
        self.assertAlmostEqual(cfg['Vs2'], 9911.38713043, delta=1.0)
        #check expanded driver gas condition
        self.assertAlmostEqual(states['s3'].p, 588992.33, delta=1.0)  
        self.assertAlmostEqual(states['s3'].T, 624.49, delta=1.0)  
        self.assertAlmostEqual(V['s3'], 3744.03, delta=1.0)
        #check shocked test gas condition
        self.assertAlmostEqual(states['s2'].p, 596149.14, delta=1.0)  
        self.assertAlmostEqual(states['s2'].T, 10164.17, delta=1.0)  
        self.assertAlmostEqual(V['s2'], 3744.03, delta=1.0)
        # check shocked test gas condition
        self.assertAlmostEqual(states['s6'].p, 9564.02, delta=1.0)  
        self.assertAlmostEqual(states['s6'].T, 47815.05, delta=1.0)  
        self.assertAlmostEqual(V['s6'], 8249.42, delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p, 788.33, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T, 1529.68, delta=1.0)  
        self.assertAlmostEqual(V['s7'], 9911.39, delta=1.0)
        #check nozzle exit conditions   
        self.assertAlmostEqual(states['s8'].p, 69.21, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 763.37, delta=1.0)  
        self.assertAlmostEqual(V['s8'], 9988.76, delta=1.0)
        self.assertAlmostEqual(cfg['stagnation_enthalpy'], 50354796.50955628, delta=1.0)
        #check pg shock over model   
        self.assertAlmostEqual(states['s10f'].p, 26251.49, delta=1.0)  
        self.assertAlmostEqual(states['s10f'].T, 48999.17, delta=1.0)  
        self.assertAlmostEqual(V['s10f'], 8298.37, delta=1.0)
        #check eq shock over model        
        self.assertAlmostEqual(states['s10e'].p, 29562.00, delta=1.0)  
        self.assertAlmostEqual(states['s10e'].T, 10040.72, delta=1.0)  
        self.assertAlmostEqual(V['s10e'], 9348.68, delta=1.0)
        #check 15 degree cone calc 
        self.assertAlmostEqual(states['s10c'].p, 2329.95, delta=1.0)  
        self.assertAlmostEqual(states['s10c'].T, 4822.44, delta=1.0)  
        self.assertAlmostEqual(V['s10c'], 9576.30, delta=1.0) 
        
        return
        
    def test_air_pressures_specified_eq_example_with_rs_out_of_st(self):
        print " "
        cfg = {'solver':'eq', 'test':'fulltheory-pressure', 
               'facility':'x2', 'nozzle':True, 'secondary': False, 'piston':'lwp-2mm-new-paper',
               'driver_gas':'He:0.8,Ar:0.2', 'test_gas':'air-cfd-with-ions', 'rs_out_of_st':True, 'Mr_st':'max',
               'p1':3000.0,'p5':10.0, 'area_ratio':5.64, 'expand_to':'shock-speed',
               'mode':'return','shock_over_model':True,'conehead':True, 'cone_angle':15.0}
        print '-'*60
        print "Running equilibrium cfd air with ions test with pressures specified and reflected shock at secondary diaphragm."
        cfg, states, V, M = run_pitot(cfg=cfg)
        #check fill pressures
        self.assertAlmostEqual(states['s1'].p, 3000.0, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 10.0, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(cfg['Vs1'], 4197.88706885, delta=1.0)        
        self.assertAlmostEqual(cfg['Vs2'], 8802.594999483183, delta=1.0)
        #check expanded driver gas condition
        self.assertAlmostEqual(states['s3'].p, 549930.00, delta=1.0)  
        self.assertAlmostEqual(states['s3'].T, 608.12, delta=1.0)  
        self.assertAlmostEqual(V['s3'], 3779.88, delta=1.0)
        #check shocked test gas condition
        self.assertAlmostEqual(states['s2'].p, 557000.00, delta=1.0)  
        self.assertAlmostEqual(states['s2'].T, 4730.77, delta=1.0)  
        self.assertAlmostEqual(V['s2'], 3779.88, delta=1.0)
        #check shocked test gas condition
        self.assertAlmostEqual(states['s2r'].p, 6541900.00, delta=1.0)  
        self.assertAlmostEqual(states['s2r'].T, 7878.00, delta=1.0)  
        self.assertAlmostEqual(V['s2r'], 0.0, delta=1.0)
        self.assertAlmostEqual(cfg['Mr-st'], 3.2308761949, delta=1.0)
        # check shocked test gas condition
        self.assertAlmostEqual(states['s6'].p, 8557.00, delta=1.0)  
        self.assertAlmostEqual(states['s6'].T, 7315.34, delta=1.0)  
        self.assertAlmostEqual(V['s6'], 8310.58, delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p, 5031.00, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T, 3099.11, delta=1.0)  
        self.assertAlmostEqual(V['s7'], 8802.59, delta=1.0)
        #check nozzle exit conditions      
        self.assertAlmostEqual(states['s8'].p, 687.68, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 2577.38, delta=1.0)  
        self.assertAlmostEqual(V['s8'], 8998.86, delta=1.0)
        self.assertAlmostEqual(cfg['stagnation_enthalpy'], 44148900.0, delta=1.0)
        #check pg shock over model   
        self.assertAlmostEqual(states['s10f'].p, 66945.86, delta=1.0)  
        self.assertAlmostEqual(states['s10f'].T, 18015.86, delta=1.0)  
        self.assertAlmostEqual(V['s10f'], 8352.87, delta=1.0)
        #check eq shock over model        
        self.assertAlmostEqual(states['s10e'].p, 67249.00, delta=1.0)  
        self.assertAlmostEqual(states['s10e'].T, 8941.42, delta=1.0)  
        self.assertAlmostEqual(V['s10e'], 8390.80, delta=1.0)
        #check 15 degree cone calc 
        self.assertAlmostEqual(states['s10c'].p, 5793.00, delta=1.0)  
        self.assertAlmostEqual(states['s10c'].T, 3457.33, delta=1.0)  
        self.assertAlmostEqual(V['s10c'], 8631.26, delta=1.0) 
        return
        
if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestPitot)
    unittest.TextTestRunner(verbosity=2).run(suite)
