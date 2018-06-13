#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  6 19:40:14 2017

@author: Han Luo

Convert poshax3 result to tecplot format
"""
import sys
import os
import numpy as np
import configparser
from scipy.constants import N_A

# load configuration file
cfgfile = sys.argv[1]
cfg = configparser.ConfigParser()
if not os.path.exists(cfgfile):
    raise FileExistsError
cfg.read(cfgfile)

# get shock velocity
datfile = cfg['controls']['output_file']
datfile = os.path.join(os.path.dirname(cfgfile), datfile)

if 'u_inf' in cfg['initial-conditions']:
    uinf = float(cfg['initial-conditions']['u_inf'])
else:
    print('Shock velocity not found in %s' % (cfgfile,))
    uinf = input('Enter the velocity in m/s:')

varnames = list()
data = list()
with open(datfile) as f:
    line = f.readline()
    while line:
        if 'Columns' in line:
            line = f.readline()
            while '#' in line:
                varnames.append(line.split(':')[1][:-1])
                line = f.readline()
            while line:
                data.append(line)
                line = f.readline()
        line = f.readline()

for i in range(len(data)):
    data[i] = [float(j) for j in data[i].split()]
data = np.array(data)
m = data.shape[0]
time = data[:, 1] / uinf   # time in sec
varnames.append('time (s)')
data = np.concatenate((data, np.reshape(time, (-1, 1))), 1)
for k0, vname in enumerate(varnames):
    if 'moles' in vname:
        break

for kf in range(k0, len(varnames)):
    if 'moles' not in varnames[kf]:
        break
if 'moles' in varnames[kf]:
    kf = len(varnames)

number_density = np.zeros((m, kf-k0))
for i in range(k0, kf):
    number_density[:, i-k0] = data[:, i] * N_A
    varnames.append(varnames[i].replace('moles', 'number_density'))
data = np.concatenate((data, number_density), 1)

filename = os.path.splitext(os.path.basename(datfile))[0]
tecfile = os.path.join(os.path.dirname(cfgfile), '%s.dat' % (filename,))
print(tecfile)
with open(tecfile, 'w') as f:
    f.write('TITLE = "poshax3 cfg:%s"\n' % (cfgfile,))
    f.write('VARIABLES = ')
    for i in varnames:
        f.write('"'+i.strip()+'" ')
    f.write('\n')
    f.write('ZONE I=%d T="u=%8.5f km/s"\n' % (m, uinf/1000,))
    f.write('AUXDATA u_inf = "%9.5f"\n' % (uinf,))
    f.write('AUXDATA p_inf = "%s"\n' % (cfg['initial-conditions']['p_inf']))
    f.write('AUXDATA T_inf = "%s"\n' % (cfg['initial-conditions']['T_inf']))
    for i in data:
        for j in i:
            f.write('%e ' % (j,))
        f.write('\n')
