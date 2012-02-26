#!/usr/bin/env python
"""
test_ideal_gas_flow.py -- test script only
"""

import sys, os
sys.path.append(os.path.expandvars("$HOME/e3bin"))

from cfpylib.gasdyn.ideal_gas_flow import demo

demo()
