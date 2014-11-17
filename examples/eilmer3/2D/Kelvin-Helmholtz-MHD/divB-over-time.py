# -*- coding: utf-8 -*-
"""
Created on Wed Nov 12 13:23:08 2014

@author: daryl
"""

import os
import gzip
import numpy as np
import matplotlib.pyplot as plt

src = ['./clean/flow',
       './un-clean/flow']

full_data = []

for pid, path in enumerate(src):
    
    data = []
    
    for fid, folder in enumerate(os.listdir(path)):
        time = 0.0
        max_divB = 0.0
        min_divB = 1.0
        for gid, gz_name in  enumerate(os.listdir(os.path.join(path,folder))):
            name = os.path.join(path,folder,gz_name)
            gz = gzip.open(name, 'r')
            file_content = gz.readlines()
            
            time = float(file_content[0])
            
            address = file_content[1].split()
            idx = address.index('"divB"')
            
            divB = []

            for line in file_content[3::]:
                line_data = line.split()
                divB.append(float(line_data[idx]))
                
            max_divB = max(max_divB, np.max(divB))
            min_divB = min(min_divB, np.min(divB))
        
        data.append([time, min_divB, max_divB])

    data = np.array(data)
    
    I = np.argsort(data[:,0])
    
    data = data[I,:]
    
    full_data.append(data)
    
data_0 = full_data[0]
data_1 = full_data[1]

plt.semilogy(data_0[:,0]/data_0[-1,0], data_0[:,2],'k-', label='clean')
plt.semilogy(data_1[:,0]/data_1[-1,0], data_1[:,2],'k--', label='not-clean')
plt.legend(loc=2)
plt.xlabel('t/t_final')
plt.ylabel('max divB')
plt.savefig("divB-trace.pdf")
plt.show()
    


        
