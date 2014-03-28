#!/usr/bin/env python

import sys

qrad = {}
qrad["measured"] = float(sys.argv[2])

ifile = open(sys.argv[1],"r")
lines = ifile.readlines()
ifile.close()

for line in lines:
    tks = line.split()
    if len(tks)==0: continue
    if tks[0]=="#": continue
    elif float(tks[0])==0.0:
        qrad["calculated"] = float(tks[3])
        break

error = abs(qrad["measured"] - qrad["calculated"])/qrad["measured"] * 100.0

print "qrad error = %0.1f percent" % ( error )
