#!/bin/bash

radmodel.py -i argon-radiators-NIST-TB-EQ.py -L rad-model.lua
./LOS.py
