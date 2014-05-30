#!/usr/bin/env python

import numpy as np

qrad, x, y = np.loadtxt("part3-viscous-with-radiation/hf_profile.data", delimiter=' ', usecols=(3, 9, 10), unpack=True)
theta = np.arctan( y / -x ) * 180 / np.pi

qrad *= 1.0e-7
negtheta = -1.0 * theta[::-1]
negqrad = qrad[::-1]

result = np.ndarray( ( len(theta)*2, 2 ) )
result[0:len(theta),0] = negtheta
result[len(theta):len(theta)*2,0] = theta
result[0:len(qrad),1] = negqrad
result[len(qrad):len(qrad)*2,1] = qrad


np.savetxt("calculated/qrad.txt", result, fmt='%.18e', delimiter=' ', newline='\n', header='theta qrad', footer='', comments='')

