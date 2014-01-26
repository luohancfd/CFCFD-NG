#!/bin/bash

radmodel.py -i argon-radiators-NIST-TB-EQ.py -L radmodel.lua
./LOS.py
