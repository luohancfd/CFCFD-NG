#!/usr/bin/env python
"""
test_pitot.py -- test script for the pitot program

MUST BE RAN FROM THE SAME FOLDER AS THE PROGRAM

Chris James - 22-Apr-2013
"""

import sys, os
sys.path.append(os.path.expandvars("$HOME/e3bin"))
import unittest
from pitot import run_pitot


class TestPitot(unittest.TestCase):
    def test_air_shocks_specified_eq_example(self):
        print " "
        cfg = {'solver':'eq', 'test':'fulltheory-shock', 
               'facility':'x2', 'nozzle':True, 'secondary': False,
               'driver_gas':'He:1.0', 'test_gas':'air',
               'Vs1':5645.0,'Vs2':11600.0, 'mode':'return'}
        cfg, states, V, M = run_pitot(cfg=cfg)
        #check fill pressures
        self.assertAlmostEqual(states['s1'].p, 3002.78, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 10.75, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(cfg['Vs1'], 5645.0, delta=1.0)        
        self.assertAlmostEqual(cfg['Vs2'], 11600.0, delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p, 15031.0, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T, 3866.2, delta=1.0)  
        self.assertAlmostEqual(V['s7'], 10471.7, delta=1.0)
        #check nozzle exit conditions         
        self.assertAlmostEqual(states['s8'].p, 4948.0, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 3313.6, delta=1.0)  
        self.assertAlmostEqual(V['s8'], 10597.7, delta=1.0)
        return
        
    def test_air_shocks_specified_pg_example(self):
        print " "
        cfg = {'solver':'pg', 'test':'fulltheory-shock', 
               'facility':'x2', 'nozzle':True, 'secondary': False,
               'driver_gas':'He:1.0', 'test_gas':'air',
               'Vs1':6120.0,'Vs2':13500.0, 'mode':'return'}
        cfg, states, V, M = run_pitot(cfg=cfg)
        #check fill pressures
        self.assertAlmostEqual(states['s1'].p, 3003.6, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 9.05, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(cfg['Vs1'], 6120.0, delta=1.0)        
        self.assertAlmostEqual(cfg['Vs2'], 13500.0, delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p, 15790.9, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T, 5492.0, delta=1.0)  
        self.assertAlmostEqual(V['s7'], 11243, delta=1.0)
        #check nozzle exit conditions         
        self.assertAlmostEqual(states['s8'].p, 4297.2, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 3787.0, delta=1.0)  
        self.assertAlmostEqual(V['s8'], 11394.4, delta=1.0)
        return
        
    def test_air_pressures_specified_eq_example(self):
        print " "
        cfg = {'solver':'eq', 'test':'fulltheory-pressure', 
               'facility':'x2', 'nozzle':True, 'secondary': False,
               'driver_gas':'He:1.0', 'test_gas':'air',
               'p1':3000.0,'p5':10.0, 'mode':'return'}
        cfg, states, V, M = run_pitot(cfg=cfg)
        #check fill pressures
        self.assertAlmostEqual(states['s1'].p, 3000.0, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 10.0, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(cfg['Vs1'], 5645.76, delta=1.0)        
        self.assertAlmostEqual(cfg['Vs2'], 11680.41, delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p, 14148.0, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T, 3829.1, delta=1.0)  
        self.assertAlmostEqual(V['s7'], 10534.7, delta=1.0)
        #check nozzle exit conditions         
        self.assertAlmostEqual(states['s8'].p, 4672.0, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 3292.5, delta=1.0)  
        self.assertAlmostEqual(V['s8'], 10658.4, delta=1.0)
        return 
        
    def test_air_pressures_specified_pg_example(self):
        print " "
        cfg = {'solver':'pg', 'test':'fulltheory-pressure', 
               'facility':'x2', 'nozzle':True, 'secondary': False,
               'driver_gas':'He:1.0', 'test_gas':'air',
               'p1':3000.0,'p5':10.0, 'mode':'return'}
        cfg, states, V, M = run_pitot(cfg=cfg)
        #check fill pressures
        self.assertAlmostEqual(states['s1'].p, 3000.0, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 10.0, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(cfg['Vs1'], 6119.87, delta=1.0)        
        self.assertAlmostEqual(cfg['Vs2'], 13390.60 , delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p, 17180.56, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T, 5628.0, delta=1.0)  
        self.assertAlmostEqual(V['s7'], 11151.8, delta=1.0)
        #check nozzle exit conditions         
        self.assertAlmostEqual(states['s8'].p, 4671.8, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 3879.7 , delta=1.0)  
        self.assertAlmostEqual(V['s8'], 11308.2, delta=1.0)
        return
        
    def test_titan_pressures_specified_eq_example(self):
        print " "
        cfg = {'solver':'eq', 'test':'fulltheory-pressure',
               'mode':'return','conehead':True, 'area_ratio':3.0,
               'facility':'x2', 'nozzle':True, 'secondary': False,
               'driver_gas':'He:0.80,Ar:0.20', 'test_gas':'titan',
               'p1':3200.0,'p5':10.0}
        cfg, states, V, M = run_pitot(cfg=cfg)
        #check fill pressures
        self.assertAlmostEqual(states['s1'].p, 3200.0, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 10.0, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(cfg['Vs1'], 4254.32, delta=1.0)        
        self.assertAlmostEqual(cfg['Vs2'],  8564.29 , delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p, 7973.0, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T, 2638.6, delta=1.0)  
        self.assertAlmostEqual(V['s7'], 8089.9, delta=1.0)
        #check nozzle exit conditions         
        self.assertAlmostEqual(states['s8'].p, 2142.0, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 2191.9, delta=1.0)  
        self.assertAlmostEqual(V['s8'], 8215.3, delta=1.0)
        #now use this test to check the conehead stuff
        self.assertAlmostEqual(states['s10c'].p, 16969.0, delta=1.0)  
        self.assertAlmostEqual(states['s10c'].T, 3210.6, delta=1.0)  
        self.assertAlmostEqual(V['s10c'], 7872.9, delta=1.0)
        return
        
    def test_titan_pressures_specified_pg_example(self):
        print " "
        cfg = {'solver':'pg', 'test':'fulltheory-pressure',
               'mode':'return','conehead':True, 'area_ratio':3.0,
               'facility':'x2', 'nozzle':True, 'secondary': False,
               'driver_gas':'He:0.80,Ar:0.20', 'test_gas':'titan',
               'p1':3200.0,'p5':10.0}
        cfg, states, V, M = run_pitot(cfg=cfg)
        #check fill pressures
        self.assertAlmostEqual(states['s1'].p, 3200.0, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 10.0, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(cfg['Vs1'], 4548.94, delta=1.0)        
        self.assertAlmostEqual(cfg['Vs2'], 10002.00, delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p, 9584.76, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T, 2993.8, delta=1.0)  
        self.assertAlmostEqual(V['s7'], 8325.3, delta=1.0)
        #check nozzle exit conditions         
        self.assertAlmostEqual(states['s8'].p, 2030.84, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 1934.1, delta=1.0)  
        self.assertAlmostEqual(V['s8'], 8461.3, delta=1.0)
        #now use this test to check the conehead stuff
        self.assertAlmostEqual(states['s10c'].p, 20239.4, delta=1.0)  
        self.assertAlmostEqual(states['s10c'].T, 4853.8, delta=1.0)  
        self.assertAlmostEqual(V['s10c'], 8086.8, delta=1.0)
        return
        
    def test_titan_experiment_eq_example(self):
        print " "
        cfg = {'solver':'eq', 'test':'experiment',
               'mode':'return','conehead':True, 'area_ratio':3.0,
               'facility':'x2', 'nozzle':True, 'secondary': False,
               'driver_gas':'He:0.80,Ar:0.20', 'test_gas':'titan',
               'p1':3200.0,'p5':10.0, 'Vs1': 4100.0, 'Vs2': 8620.0}
        cfg, states, V, M = run_pitot(cfg=cfg)
        #check fill pressures
        self.assertAlmostEqual(states['s1'].p, 3200.0, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 10.0, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(cfg['Vs1'], 4100.00, delta=1.0)        
        self.assertAlmostEqual(cfg['Vs2'], 8620.00, delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p, 5056.0, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T, 2380.1, delta=1.0)  
        self.assertAlmostEqual(V['s7'], 8141.4, delta=1.0)
        #check nozzle exit conditions         
        self.assertAlmostEqual(states['s8'].p, 1305.0, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 1884.4, delta=1.0)  
        self.assertAlmostEqual(V['s8'], 8254.1, delta=1.0)
        return                        
                               
    def test_titan_experiment_pg_example(self):
        print " "
        cfg = {'solver':'pg', 'test':'experiment',
               'mode':'return','conehead':True, 'area_ratio':3.0,
               'facility':'x2', 'nozzle':True, 'secondary': False,
               'driver_gas':'He:0.80,Ar:0.20', 'test_gas':'titan',
               'p1':3200.0,'p5':10.0, 'Vs1': 4100.0, 'Vs2': 8620.0}
        cfg, states, V, M = run_pitot(cfg=cfg)
        #check fill pressures
        self.assertAlmostEqual(states['s1'].p, 3200.0, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 10.0, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(cfg['Vs1'], 4100.00, delta=1.0)        
        self.assertAlmostEqual(cfg['Vs2'], 8620.00, delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p, 12252.11, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T, 2782.1, delta=1.0)  
        self.assertAlmostEqual(V['s7'], 7172.0, delta=1.0)
        #check nozzle exit conditions         
        self.assertAlmostEqual(states['s8'].p, 2581.37, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 1794.5, delta=1.0)  
        self.assertAlmostEqual(V['s8'], 7318.8, delta=1.0)
        return 
    
    def test_h2he_pressures_specified_eq_example(self):
        print " "
        cfg = {'solver':'eq', 'test':'fulltheory-pressure',
               'mode':'return','conehead':False, 'area_ratio':2.5,
               'facility':'x2', 'nozzle':True, 'secondary': True,
               'driver_gas':'He:1.0', 'test_gas':'gasgiant_h215he',
               'psd1':17500.0, 'p1':4700.0, 'p5':6.37, 
               'shock_over_model':True}
        cfg, states, V, M = run_pitot(cfg=cfg)
        #check fill pressures
        self.assertAlmostEqual(states['sd1'].p, 17500, delta=1.0)
        self.assertAlmostEqual(states['s1'].p, 4700.0, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 6.37, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(cfg['Vsd'], 6999.99, delta=1.0)        
        self.assertAlmostEqual(cfg['Vs1'], 9076.89, delta=1.0)
        self.assertAlmostEqual(cfg['Vs2'], 16428.12, delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p, 18095.0, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T, 1638.2, delta=1.0)  
        self.assertAlmostEqual(V['s7'], 14953.6, delta=1.0)
        #check nozzle exit conditions         
        self.assertAlmostEqual(states['s8'].p, 4917.0, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 1144.7, delta=1.0)  
        self.assertAlmostEqual(V['s8'], 15379.0, delta=1.0)
        #should check the shock over model stuff here too
        self.assertAlmostEqual(states['s10f'].p, 235100.6, delta=1.0)  
        self.assertAlmostEqual(states['s10f'].T, 10177.5, delta=1.0)  
        self.assertAlmostEqual(V['s10f'], 12519.6, delta=1.0)
        self.assertAlmostEqual(states['s10e'].p, 262180.0, delta=1.0)  
        self.assertAlmostEqual(states['s10e'].T, 3978.6, delta=1.0)  
        self.assertAlmostEqual(V['s10e'], 13993.7, delta=1.0)
        return  
        
    def test_h2he_pressures_specified_pg_example(self):
        print " "
        cfg = {'solver':'pg', 'test':'fulltheory-pressure',
               'mode':'return','conehead':False, 'area_ratio':2.5,
               'facility':'x2', 'nozzle':True, 'secondary': True,
               'driver_gas':'He:1.0', 'test_gas':'gasgiant_h215he',
               'psd1':17500.0, 'p1':4700.0, 'p5':6.37, 
               'shock_over_model':True}
        cfg, states, V, M = run_pitot(cfg=cfg)
        #check fill pressures
        self.assertAlmostEqual(states['sd1'].p, 17500, delta=1.0)
        self.assertAlmostEqual(states['s1'].p, 4700.0, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 6.37, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(cfg['Vsd'], 6999.43, delta=1.0)        
        self.assertAlmostEqual(cfg['Vs1'], 9458.49, delta=1.0)
        self.assertAlmostEqual(cfg['Vs2'], 17777.01, delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p, 19400.56, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T, 1681.5, delta=1.0)  
        self.assertAlmostEqual(V['s7'], 14809.1, delta=1.0)
        #check nozzle exit conditions         
        self.assertAlmostEqual(states['s8'].p, 5014.24, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 1119.3, delta=1.0)  
        self.assertAlmostEqual(V['s8'], 15255.9, delta=1.0)
        #should check the shock over model stuff here too
        self.assertAlmostEqual(states['s10f'].p, 237934.0, delta=1.0)  
        self.assertAlmostEqual(states['s10f'].T, 10482.3, delta=1.0)  
        self.assertAlmostEqual(V['s10f'], 12245.0, delta=1.0)
        self.assertAlmostEqual(states['s10e'].p, 268850.0, delta=1.0)  
        self.assertAlmostEqual(states['s10e'].T,  3966.6, delta=1.0)  
        self.assertAlmostEqual(V['s10e'], 13872.7, delta=1.0)
        return 
        
    def test_h2ne_pressures_specified_eq_example(self):
        print " "
        cfg = {'solver':'eq', 'test':'fulltheory-pressure',
               'mode':'return','conehead':False, 'area_ratio':2.5,
               'facility':'x2', 'nozzle':True, 'secondary': True,
               'driver_gas':'He:1.0', 'test_gas':'gasgiant_h285ne',
               'psd1':30000.0, 'p1':2400.0, 'p5':2.5, 
               'shock_over_model':True}
        cfg, states, V, M = run_pitot(cfg=cfg)
        #check fill pressures
        self.assertAlmostEqual(states['sd1'].p, 30000, delta=1.0)
        self.assertAlmostEqual(states['s1'].p, 2400.0, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 2.5, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(cfg['Vsd'], 6498.78, delta=1.0)        
        self.assertAlmostEqual(cfg['Vs1'], 7432.99, delta=1.0)
        self.assertAlmostEqual(cfg['Vs2'], 15581.87, delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p, 6285.0, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T, 3225.1, delta=1.0)  
        self.assertAlmostEqual(V['s7'], 14080.3, delta=1.0)
        #check nozzle exit conditions         
        self.assertAlmostEqual(states['s8'].p, 2061.0, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 2771.9, delta=1.0)  
        self.assertAlmostEqual(V['s8'], 14202.9, delta=1.0)
        #should check the shock over model stuff here too
        self.assertAlmostEqual(states['s10f'].p, 268153.0, delta=1.0)  
        self.assertAlmostEqual(states['s10f'].T, 33370.4, delta=1.0)  
        self.assertAlmostEqual(V['s10f'], 12889.1, delta=1.0)
        self.assertAlmostEqual(states['s10e'].p, 269270.0, delta=1.0)  
        self.assertAlmostEqual(states['s10e'].T, 19669.0, delta=1.0)  
        self.assertAlmostEqual(V['s10e'], 12942.9, delta=1.0)
        return   
        
    def test_h2ne_pressures_specified_pg_example(self):
        print " "
        cfg = {'solver':'pg', 'test':'fulltheory-pressure',
               'mode':'return','conehead':False, 'area_ratio':2.5,
               'facility':'x2', 'nozzle':True, 'secondary': True,
               'driver_gas':'He:1.0', 'test_gas':'gasgiant_h285ne',
               'psd1':30000.0, 'p1':2400.0, 'p5':2.5, 
               'shock_over_model':True}
        cfg, states, V, M = run_pitot(cfg=cfg)
        #check fill pressures
        self.assertAlmostEqual(states['sd1'].p, 30000.0, delta=1.0)
        self.assertAlmostEqual(states['s1'].p, 2400.0, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 2.5, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(cfg['Vsd'], 6498.56, delta=1.0)        
        self.assertAlmostEqual(cfg['Vs1'], 7947.77, delta=1.0)
        self.assertAlmostEqual(cfg['Vs2'], 17283.97, delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p,  7125.08, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T, 3999.4, delta=1.0)  
        self.assertAlmostEqual(V['s7'], 14398.1, delta=1.0)
        #check nozzle exit conditions         
        self.assertAlmostEqual(states['s8'].p, 1606.16, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 2277.4, delta=1.0)  
        self.assertAlmostEqual(V['s8'], 14548.0, delta=1.0)
        #should check the shock over model stuff here too
        self.assertAlmostEqual(states['s10f'].p, 239974.4, delta=1.0)  
        self.assertAlmostEqual(states['s10f'].T, 81445.3, delta=1.0)  
        self.assertAlmostEqual(V['s10f'], 11065.8, delta=1.0)
        self.assertAlmostEqual(states['s10e'].p, 284370.0, delta=1.0)  
        self.assertAlmostEqual(states['s10e'].T, 19860.3, delta=1.0)  
        self.assertAlmostEqual(V['s10e'], 13295.4, delta=1.0)
        return  
        
    def test_scramjet_pressures_specified_eq_example(self):
        print " "
        cfg = {'solver':'eq', 'test':'fulltheory-pressure',
               'mode':'return', 'area_ratio':2.5,
               'facility':'x2', 'nozzle':True, 'secondary': True,
               'driver_gas':'He:0.80,Ar:0.20', 'test_gas':'air',
               'psd1':100000.0, 'p1':486000.0, 'p5':1500, 'shock_switch':True}                 
        cfg, states, V, M = run_pitot(cfg=cfg)
        #check fill pressures
        self.assertAlmostEqual(states['sd1'].p, 100000.0, delta=1.0)
        self.assertAlmostEqual(states['s1'].p, 486000.0, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 1500.0, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(cfg['Vsd'], 4290.03, delta=1.0)        
        self.assertAlmostEqual(cfg['Vs1'], 1588.02, delta=1.0)
        self.assertAlmostEqual(cfg['Vs2'], 3427.67, delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p, 91285.0, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T, 472.3, delta=1.0)  
        self.assertAlmostEqual(V['s7'], 3065.3, delta=1.0)
        #check nozzle exit conditions         
        self.assertAlmostEqual(states['s8'].p, 24898.0, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 327.1, delta=1.0)  
        self.assertAlmostEqual(V['s8'], 3112.9, delta=1.0)
        return
        
    def test_scramjet_shocks_specified_eq_example(self):
        print " "
        cfg = {'solver':'eq', 'test':'fulltheory-shock',
               'mode':'return', 'area_ratio':2.5,
               'facility':'x2', 'nozzle':True, 'secondary': True,
               'driver_gas':'He:0.80,Ar:0.20', 'test_gas':'air',
               'Vsd':4290.0, 'Vs1':1588.0, 'Vs2':3424.0}                    
        cfg, states, V, M = run_pitot(cfg=cfg)
        #check fill pressures
        self.assertAlmostEqual(states['sd1'].p, 100020.0, delta=1.0)
        self.assertAlmostEqual(states['s1'].p, 486080.0, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 1520.993, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(cfg['Vsd'], 4290, delta=1.0)        
        self.assertAlmostEqual(cfg['Vs1'], 1588.00, delta=1.0)
        self.assertAlmostEqual(cfg['Vs2'], 3424.00, delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p, 92394.0, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T, 473.9, delta=1.0)  
        self.assertAlmostEqual(V['s7'], 3061.5, delta=1.0)
        #check nozzle exit conditions         
        self.assertAlmostEqual(states['s8'].p, 25199.0, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 328.1, delta=1.0)  
        self.assertAlmostEqual(V['s8'], 3109.4, delta=1.0)
        return

if __name__ == '__main__':
    pitot_input = raw_input("This is the test program for Pitot. Type 'yes' to begin or anything else to quit: ")
    if pitot_input == 'yes':
        suite = unittest.TestLoader().loadTestsFromTestCase(TestPitot)
        unittest.TextTestRunner(verbosity=2).run(suite)
