from pylab import *
import numpy as np
import ctrans


def makeGaussian(size):
    Qmax = array([1.0, 1.0, 1.0])
    Qmin = array([-1.0, -1.0, -1.0])
    dQN = array([size,size, size])
    sigma = 0.1
    grid = np.mgrid[0:dQN[0], 0:dQN[1], 0:dQN[2]]
    r = (Qmax - Qmin) / dQN

    X = grid[0] * r[0] + Qmin[0]
    Y = grid[1] * r[1] + Qmin[1]
    Z = grid[2] * r[2] + Qmin[2]
    
    out = np.exp(-(X**2 + Y**2 + Z**2) / (2 * sigma**2))

    out = array([np.ravel(X),
                 np.ravel(Y),
                 np.ravel(Z),
                 np.ravel(out)])
    
    return out.T

def griddata(data, size):

    gridData, gridOccu, gridOut = ctrans.grid3d(data,
                                                array([-1.0, -1.0, -1.0]),
                                                array([1.0, 1.0, 1.0]),
                                                array([size, size, size]), norm = 1)
    print "gridout = ", gridOut,
    print "gridOccu sum = ", gridOccu.sum(),
    boxsize = size / 8
    midp = size / 2
    box = gridData[midp-boxsize:midp+boxsize+1,
                   midp-boxsize:midp+boxsize+1,
                   midp-boxsize:midp+boxsize+1]
    boxOccu = gridOccu[midp-boxsize:midp+boxsize+1,
                       midp-boxsize:midp+boxsize+1,
                       midp-boxsize:midp+boxsize+1]
    print (boxOccu == 0).sum(), 
    print gridData.sum(), box.sum() / (size**3),
    print boxOccu.size / (boxOccu != 0).sum(),
    #figure()
    #subplot(1, 2, 1)
    #imshow(gridData.sum(0))
    #subplot(1, 2, 2)
    #imshow(box.sum(0))

if __name__ == "__main__":
    for j in range(10, 110, 10):
        data = makeGaussian(j)
        for i in range(10, 60, 10):
            griddata(data, i)     
        

