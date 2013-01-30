#!/bin/bash
radmodel.py -i all-radiators.py -L rad-model.lua
gasfile gmodel.inp gas-model.lua
radGUI.py
