#!/usr/bin/env python
import os

ifile = open("hayabusa.control")
lines = ifile.readlines()
ifile.close()
for line in lines:
    tks = line.split()
    if len(tks)<2: continue
    if tks[0]=="max_step": initial_max_step = int(tks[2])

prev_max_step = initial_max_step
max_steps = [ 100, 100, 1000, 1000, 1000, 1000, 1000 ]

for iteration,max_step in enumerate(max_steps):
    print "Performing flowfield-radiation coupling iteration %d" % iteration
    os.system("e3rad.exe -f hayabusa -r > E3RAD_LOG_%d" % iteration)
    os.system("sed -ie 's/max_step = %d/max_step = %d/g' hayabusa.control" % ( prev_max_step, max_step ) ) 
    os.system("mpirun -np 4 e3mpi.exe -f hayabusa -t %d -r > E3MPI_LOG_%d" % (iteration+1,iteration))
    prev_max_step = max_step
