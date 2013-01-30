#!/usr/bin/env python
"""
test_pitot.py -- test script for the pitot program

MUST BE RAN FROM THE SAME FOLDER AS THE PROGRAM

Chris James - 22-Jan-2013
"""

import sys, os
sys.path.append(os.path.expandvars("$HOME/e3bin"))

import unittest
import pitot


class TestPitot(unittest.TestCase):
    def test_air_shocks_specified_eq_example(self):
        print " "
        sys.argv = ['pitot.py', '--mode', 'return', '--solver','eq','--test','fulltheory-shock', 
                    '--Vs1','5645.0', '--Vs2','11600.0']
        states, V, M, Vs1, Vs2 = pitot.main()
        #check fill pressures
        self.assertAlmostEqual(states['s1'].p, 3002.0, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 10.58, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(Vs1, 5645.0, delta=1.0)        
        self.assertAlmostEqual(Vs2, 11600.0, delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p, 13698.0, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T, 3807.9, delta=1.0)  
        self.assertAlmostEqual(V['s7'], 10471.8, delta=1.0)
        #check nozzle exit conditions         
        self.assertAlmostEqual(states['s8'].p, 4531.0, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 3280.4, delta=1.0)  
        self.assertAlmostEqual(V['s8'], 10595.4, delta=1.0)
        return
        
    def test_air_shocks_specified_pg_example(self):
        print " "
        sys.argv = ['pitot.py', '--mode', 'return', '--solver','pg','--test','fulltheory-shock', 
                    '--Vs1','6120.0', '--Vs2','13500.0']
        states, V, M, Vs1, Vs2 = pitot.main()
        #check fill pressures
        self.assertAlmostEqual(states['s1'].p, 3003.6, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 8.91, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(Vs1, 6120.0, delta=1.0)        
        self.assertAlmostEqual(Vs2, 13500.0, delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p, 14310.0, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T, 5340.0, delta=1.0)  
        self.assertAlmostEqual(V['s7'], 11243, delta=1.0)
        #check nozzle exit conditions         
        self.assertAlmostEqual(states['s8'].p, 3896.3, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 3682.5, delta=1.0)  
        self.assertAlmostEqual(V['s8'], 11390.2, delta=1.0)
        return
        
    def test_air_pressures_specified_eq_example(self):
        print " "
        sys.argv = ['pitot.py', '--mode','return', '--solver','eq', '--p1','3000.0','--p5','10.0']
        states, V, M, Vs1, Vs2 = pitot.main()
        #check fill pressures
        self.assertAlmostEqual(states['s1'].p, 3000.0, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 10.0, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(Vs1, 5645.69, delta=1.0)        
        self.assertAlmostEqual(Vs2, 11660.78, delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p, 13062.0, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T, 3779.9, delta=1.0)  
        self.assertAlmostEqual(V['s7'], 10519.3, delta=1.0)
        #check nozzle exit conditions         
        self.assertAlmostEqual(states['s8'].p, 4331.0, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 3264.4, delta=1.0)  
        self.assertAlmostEqual(V['s8'], 10641.2, delta=1.0)
        return 
        
    def test_air_pressures_specified_pg_example(self):
        print " "
        sys.argv = ['pitot.py', '--mode','return', '--solver','pg', '--p1','3000.0','--p5','10.0']
        states, V, M, Vs1, Vs2 = pitot.main()
        #check fill pressures
        self.assertAlmostEqual(states['s1'].p, 3000.0, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 10.0, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(Vs1, 6119.87, delta=1.0)        
        self.assertAlmostEqual(Vs2, 13376.76 , delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p, 15821.07, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T, 5497.0, delta=1.0)  
        self.assertAlmostEqual(V['s7'], 11140.3, delta=1.0)
        #check nozzle exit conditions         
        self.assertAlmostEqual(states['s8'].p, 4303.906, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 3789.9 , delta=1.0)  
        self.assertAlmostEqual(V['s8'], 11293.2, delta=1.0)
        return
        
    def test_titan_pressures_specified_eq_example(self):
        print " "
        sys.argv = ['pitot.py', '--driver_gas','He:0.80,Ar:0.20','--test_gas','titan', 
                        '--p1','3200.0','--p5','10.0', '--ar','3.0', '--conehead','yes',
                        '--mode','return',  '--solver', 'eq']
        states, V, M, Vs1, Vs2 = pitot.main()
        #check fill pressures
        self.assertAlmostEqual(states['s1'].p, 3200.0, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 10.0, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(Vs1, 4253.75, delta=1.0)        
        self.assertAlmostEqual(Vs2,  8549.44 , delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p, 7386.0, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T, 2612.6 , delta=1.0)  
        self.assertAlmostEqual(V['s7'], 8076.1, delta=1.0)
        #check nozzle exit conditions         
        self.assertAlmostEqual(states['s8'].p, 1979.0, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 2163.4, delta=1.0)  
        self.assertAlmostEqual(V['s8'], 8200.5, delta=1.0)
        #now use this test to check the conehead stuff
        self.assertAlmostEqual(states['s10c'].p, 15813.0, delta=1.0)  
        self.assertAlmostEqual(states['s10c'].T,3184.4, delta=1.0)  
        self.assertAlmostEqual(V['s10c'], 7859.0, delta=1.0)
        return
        
    def test_titan_pressures_specified_pg_example(self):
        print " "
        sys.argv = ['pitot.py', '--driver_gas','He:0.80,Ar:0.20','--test_gas','titan', 
                        '--p1','3200.0','--p5','10.0', '--ar','3.0', '--conehead','yes',
                        '--mode','return', '--solver', 'pg']
        states, V, M, Vs1, Vs2 = pitot.main()
        #check fill pressures
        self.assertAlmostEqual(states['s1'].p, 3200.0, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 10.0, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(Vs1, 4548.41, delta=1.0)        
        self.assertAlmostEqual(Vs2,  9990.38 , delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p, 8825.643, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T, 2924.6, delta=1.0)  
        self.assertAlmostEqual(V['s7'], 8315.6 , delta=1.0)
        #check nozzle exit conditions         
        self.assertAlmostEqual(states['s8'].p, 1870.876, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 1889.6, delta=1.0)  
        self.assertAlmostEqual(V['s8'], 8448.6, delta=1.0)
        #now use this test to check the conehead stuff
        self.assertAlmostEqual(states['s10c'].p, 18981.3, delta=1.0)  
        self.assertAlmostEqual(states['s10c'].T, 4795.5, delta=1.0)  
        self.assertAlmostEqual(V['s10c'], 8075.3, delta=1.0)
        return
        
    def test_titan_experiment_eq_example(self):
        print " "
        sys.argv = ['pitot.py', '--test','experiment', '--driver_gas','He:0.80,Ar:0.20',
                        '--test_gas','titan', '--p1','3200.0','--p5','10.0', 
                        '--Vs1','4100.0','--Vs2','8620.0','--ar','3.0', 
                        '--mode','return', '--solver', 'eq']
        states, V, M, Vs1, Vs2 = pitot.main()
        #check fill pressures
        self.assertAlmostEqual(states['s1'].p, 3200.0, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 10.0, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(Vs1, 4100.00, delta=1.0)        
        self.assertAlmostEqual(Vs2, 8620.00, delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p, 4520.0, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T,2340.8, delta=1.0)  
        self.assertAlmostEqual(V['s7'], 8141.4, delta=1.0)
        #check nozzle exit conditions         
        self.assertAlmostEqual(states['s8'].p, 1160.0, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 1839.5, delta=1.0)  
        self.assertAlmostEqual(V['s8'], 8252.3, delta=1.0)
        return                        
                               
    def test_titan_experiment_pg_example(self):
        print " "
        sys.argv = ['pitot.py', '--test','experiment', '--driver_gas','He:0.80,Ar:0.20',
                        '--test_gas','titan', '--p1','3200.0','--p5','10.0', 
                        '--Vs1','4100.0','--Vs2','8620.0','--ar','3.0', 
                        '--mode','return', '--solver', 'pg']
        states, V, M, Vs1, Vs2 = pitot.main()
        #check fill pressures
        self.assertAlmostEqual(states['s1'].p, 3200.0, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 10.0, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(Vs1, 4100.00, delta=1.0)        
        self.assertAlmostEqual(Vs2, 8620.00, delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p, 11407.81, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T, 2726.7, delta=1.0)  
        self.assertAlmostEqual(V['s7'], 7172.0, delta=1.0)
        #check nozzle exit conditions         
        self.assertAlmostEqual(states['s8'].p, 2404.83, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 1759.0, delta=1.0)  
        self.assertAlmostEqual(V['s8'], 7315.9, delta=1.0)
        return 
    
    def test_h2he_pressures_specified_eq_example(self):
        print " "
        sys.argv = ['pitot.py', '--test','fulltheory-pressure', '--config',
                    'sec-nozzle','--driver_gas','He:1.0', '--test_gas',
                    'gasgiant_h215he','--psd1','17500','--p1','4700',
                    '--p5','6.37','--mode','return','--solver','eq',
                    '--shock_over_model','True']
        states, V, M, Vsd, Vs1, Vs2 = pitot.main()
        #check fill pressures
        self.assertAlmostEqual(states['sd1'].p, 17500, delta=1.0)
        self.assertAlmostEqual(states['s1'].p, 4700.0, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 6.37, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(Vsd, 6999.99, delta=1.0)        
        self.assertAlmostEqual(Vs1, 9076.83, delta=1.0)
        self.assertAlmostEqual(Vs2, 16426.43, delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p, 17405.0, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T, 1621.4, delta=1.0)  
        self.assertAlmostEqual(V['s7'], 14951.9 , delta=1.0)
        #check nozzle exit conditions         
        self.assertAlmostEqual(states['s8'].p, 4726.0, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 1131.9, delta=1.0)  
        self.assertAlmostEqual(V['s8'], 15373.1, delta=1.0)
        #should check the shock over model stuff here too
        self.assertAlmostEqual(states['s10f'].p, 228266.3, delta=1.0)  
        self.assertAlmostEqual(states['s10f'].T, 10168.8, delta=1.0)  
        self.assertAlmostEqual(V['s10f'], 12513.6, delta=1.0)
        self.assertAlmostEqual(states['s10e'].p, 254660.0, delta=1.0)  
        self.assertAlmostEqual(states['s10e'].T, 3969.4, delta=1.0)  
        self.assertAlmostEqual(V['s10e'], 13991.0, delta=1.0)
        return  
        
    def test_h2he_pressures_specified_pg_example(self):
        print " "
        sys.argv = ['pitot.py', '--test','fulltheory-pressure', '--config',
                    'sec-nozzle','--driver_gas','He:1.0', '--test_gas',
                    'gasgiant_h215he','--psd1','17500','--p1','4700',
                    '--p5','6.37','--mode','return','--solver','pg',
                    '--shock_over_model','True']
        states, V, M, Vsd, Vs1, Vs2 = pitot.main()
        #check fill pressures
        self.assertAlmostEqual(states['sd1'].p, 17500, delta=1.0)
        self.assertAlmostEqual(states['s1'].p, 4700.0, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 6.37, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(Vsd, 6999.42, delta=1.0)        
        self.assertAlmostEqual(Vs1, 9458.47, delta=1.0)
        self.assertAlmostEqual(Vs2, 17775.37, delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p, 18686.78, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T, 1662.6, delta=1.0)  
        self.assertAlmostEqual(V['s7'], 14807.8, delta=1.0)
        #check nozzle exit conditions         
        self.assertAlmostEqual(states['s8'].p, 4832.016, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 1106.9, delta=1.0)  
        self.assertAlmostEqual(V['s8'], 15249.5, delta=1.0)
        #should check the shock over model stuff here too
        self.assertAlmostEqual(states['s10f'].p, 231669.9, delta=1.0)  
        self.assertAlmostEqual(states['s10f'].T, 10462.9, delta=1.0)  
        self.assertAlmostEqual(V['s10f'], 12243.0, delta=1.0)
        self.assertAlmostEqual(states['s10e'].p, 261780.0, delta=1.0)  
        self.assertAlmostEqual(states['s10e'].T,  3958.0, delta=1.0)  
        self.assertAlmostEqual(V['s10e'], 13869.3, delta=1.0)
        return 
        
    def test_h2ne_pressures_specified_eq_example(self):
        print " "
        sys.argv = ['pitot.py', '--test','fulltheory-pressure', '--config',
                    'sec-nozzle','--driver_gas','He:1.0', '--test_gas',
                    'gasgiant_h285ne','--psd1','30000','--p1','2400',
                    '--p5','2.5','--mode','return','--solver','eq',
                    '--shock_over_model','True']
        states, V, M, Vsd, Vs1, Vs2 = pitot.main()
        #check fill pressures
        self.assertAlmostEqual(states['sd1'].p, 30000, delta=1.0)
        self.assertAlmostEqual(states['s1'].p, 2400.0, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 2.5, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(Vsd, 6498.78, delta=1.0)        
        self.assertAlmostEqual(Vs1, 7433.07, delta=1.0)
        self.assertAlmostEqual(Vs2, 15533.92, delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p, 5778.0, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T, 3182.6, delta=1.0)  
        self.assertAlmostEqual(V['s7'], 14031.1, delta=1.0)
        #check nozzle exit conditions         
        self.assertAlmostEqual(states['s8'].p, 1901.0, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 2744.8, delta=1.0)  
        self.assertAlmostEqual(V['s8'], 14152.0, delta=1.0)
        #should check the shock over model stuff here too
        self.assertAlmostEqual(states['s10f'].p, 248841.4, delta=1.0)  
        self.assertAlmostEqual(states['s10f'].T, 33005.4, delta=1.0)  
        self.assertAlmostEqual(V['s10f'], 12851.8, delta=1.0)
        self.assertAlmostEqual(states['s10e'].p, 249730.0, delta=1.0)  
        self.assertAlmostEqual(states['s10e'].T, 19541.0, delta=1.0)  
        self.assertAlmostEqual(V['s10e'], 12898.5, delta=1.0)
        return   
        
    def test_h2ne_pressures_specified_pg_example(self):
        print " "
        sys.argv = ['pitot.py', '--test','fulltheory-pressure', '--config',
                    'sec-nozzle','--driver_gas','He:1.0', '--test_gas',
                    'gasgiant_h285ne','--psd1','30000','--p1','2400',
                    '--p5','2.5','--mode','return','--solver','pg',
                    '--shock_over_model','True']
        states, V, M, Vsd, Vs1, Vs2 = pitot.main()
        #check fill pressures
        self.assertAlmostEqual(states['sd1'].p, 30000.0, delta=1.0)
        self.assertAlmostEqual(states['s1'].p, 2400.0, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 2.5, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(Vsd, 6498.56, delta=1.0)        
        self.assertAlmostEqual(Vs1, 7947.77, delta=1.0)
        self.assertAlmostEqual(Vs2, 17242.44, delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p, 6415.755, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T, 3844.0, delta=1.0)  
        self.assertAlmostEqual(V['s7'], 14363.5, delta=1.0)
        #check nozzle exit conditions         
        self.assertAlmostEqual(states['s8'].p, 1447.083, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 2189.4, delta=1.0)  
        self.assertAlmostEqual(V['s8'], 14507.9, delta=1.0)
        #should check the shock over model stuff here too
        self.assertAlmostEqual(states['s10f'].p, 223673.8, delta=1.0)  
        self.assertAlmostEqual(states['s10f'].T, 80929.6, delta=1.0)  
        self.assertAlmostEqual(V['s10f'], 11038.3, delta=1.0)
        self.assertAlmostEqual(states['s10e'].p, 266190.0, delta=1.0)  
        self.assertAlmostEqual(states['s10e'].T, 19743.0, delta=1.0)  
        self.assertAlmostEqual(V['s10e'], 13261.9, delta=1.0)
        return  
        
    def test_scramjet_pressures_specified_eq_example(self):
        print " "
        sys.argv = ['pitot.py', '--test','fulltheory-pressure', '--config',
                    'sec-nozzle','--driver_gas','He:0.80,Ar:0.20', '--test_gas',
                    'air','--psd1','100000','--p1','486000', '--p5','1500.0',
                    '--shock_switch','True', '--mode','return','--solver','eq']
                    
        states, V, M, Vsd, Vs1, Vs2 = pitot.main()
        #check fill pressures
        self.assertAlmostEqual(states['sd1'].p, 100000.0, delta=1.0)
        self.assertAlmostEqual(states['s1'].p, 486000.0, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 1500.0, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(Vsd, 4290.03, delta=1.0)        
        self.assertAlmostEqual(Vs1, 1588.00, delta=1.0)
        self.assertAlmostEqual(Vs2, 3423.95, delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p, 83945.0, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T, 461.3, delta=1.0)  
        self.assertAlmostEqual(V['s7'], 3061.6, delta=1.0)
        #check nozzle exit conditions         
        self.assertAlmostEqual(states['s8'].p, 22894.0, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 319.3, delta=1.0)  
        self.assertAlmostEqual(V['s8'], 3108.2, delta=1.0)
        return
        
    def test_scramjet_shocks_specified_eq_example(self):
        print " "
        sys.argv = ['pitot.py', '--test','fulltheory-shock', '--config',
                    'sec-nozzle','--driver_gas','He:0.80,Ar:0.20', '--test_gas',
                    'air', '--Vsd','4290.0','--Vs1','1588.0', '--Vs2','3424.0',
                    '--mode','return','--solver','eq']
                    
        states, V, M, Vsd, Vs1, Vs2 = pitot.main()
        #check fill pressures
        self.assertAlmostEqual(states['sd1'].p, 100010.0, delta=1.0)
        self.assertAlmostEqual(states['s1'].p, 486000.0, delta=1.0)
        self.assertAlmostEqual(states['s5'].p, 1498.673, delta=1.0)
        #check shock speeds      
        self.assertAlmostEqual(Vsd, 4290, delta=1.0)        
        self.assertAlmostEqual(Vs1, 1588.00, delta=1.0)
        self.assertAlmostEqual(Vs2, 3424.00, delta=1.0)
        #check nozzle entry conditions
        self.assertAlmostEqual(states['s7'].p, 83931.0, delta=1.0)  
        self.assertAlmostEqual(states['s7'].T, 461.3, delta=1.0)  
        self.assertAlmostEqual(V['s7'], 3061.7, delta=1.0)
        #check nozzle exit conditions         
        self.assertAlmostEqual(states['s8'].p, 22890.0, delta=1.0)  
        self.assertAlmostEqual(states['s8'].T, 319.3, delta=1.0)  
        self.assertAlmostEqual(V['s8'], 3108.3, delta=1.0)
        return

if __name__ == '__main__':
    pitot_input = raw_input("This is the test program for Pitot. Type 'yes' to begin or anything else to quit: ")
    if pitot_input == 'yes':
        suite = unittest.TestLoader().loadTestsFromTestCase(TestPitot)
        unittest.TextTestRunner(verbosity=2).run(suite)
