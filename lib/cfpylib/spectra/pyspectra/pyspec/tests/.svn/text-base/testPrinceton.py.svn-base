from pyspec.ccd.PrincetonSPE import PrincetonSPEFile
import princeton
import time
import numpy as np
import pylab

def test():
    t1 = time.time()
    conIm = PrincetonSPEFile('testimage2.spe')
    conIm = conIm.getData().astype(np.float32)
    t2 = time.time()
    print "Read took %.3f msec" % ((t2 - t1) * 1000)

    pylab.figure()
    pylab.imshow(conIm.sum(0))

    t1 = time.time()
    conIm = princeton.readSPE('testimage2.spe').astype(np.float32)
    t2 = time.time()
    print "Read took %.3f msec" % ((t2 - t1) * 1000)

    pylab.figure()
    pylab.imshow(conIm.sum(0))
    pylab.show()
    
if __name__ == "__main__":
    test()
