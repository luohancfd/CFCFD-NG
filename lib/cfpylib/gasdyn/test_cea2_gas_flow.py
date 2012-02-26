#!/usr/bin/env python
"""
test_cea2_gas_flow.py -- test script only
"""

import sys, os
sys.path.append(os.path.expandvars("$HOME/e3bin"))

from cfpylib.gasdyn.cea2_gas import Gas
from cfpylib.gasdyn.cea2_gas_flow import demo

demo()
