#!/bin/bash
time ./octvce.exe -gas ov_gas.par -bc ov_BC.par -ic ov_IC.par -par ov.par -geom bodies.pvtp -2D=axisymmetric -ncpus 1 -count=energy 
