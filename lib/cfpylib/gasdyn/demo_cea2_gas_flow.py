#!/usr/bin/env python
"""
test_cea2_gas_flow.py -- test script
"""

import sys, os
sys.path.append(os.path.expandvars("$HOME/e3bin"))
from cfpylib.gasdyn.cea2_gas_flow import demo

demo()
