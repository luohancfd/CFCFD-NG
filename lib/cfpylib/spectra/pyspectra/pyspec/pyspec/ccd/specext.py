#
# specext.py (c) Stuart B. Wilkins 2011
#
# $Id: fit.py 176 2011-02-13 19:08:16Z stuwilkins $
# $HeadURL: https://pyspec.svn.sourceforge.net/svnroot/pyspec/trunk/pyspec/fit.py $
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Part of the "pyspec" package
#

import os
import numpy
from pyspec.spec import SpecExtension

class CCDSpecExtension(SpecExtension):
    def __init__(self):
        self.ccddark = "-DARK"
        self.ccdpath = None
        self.ccdtail = ".spe"
    def initSpec(self, object):
        object.ccdpath = self.ccdpath
        object.ccddark = self.ccddark
        object.ccdtail = self.ccdtail
    def initSpecScan(self, object):
        object.ccdAcquireTime = 0.0
        object.ccdAcquirePeriod = 0.0
        object.ccdNumExposures = 1
        object.ccdNumImages = 1
        object.ccdNumAcquisitions = 1
    def getName(self):
        return "SPEC / EPICS CCD Extension"
    def parseSpecScanHeader(self, object, line):
        if line[0:6] == "#UCCD2":
            print "---- Reading CCD Header information"
            try:
                pos = line[6:].strip().split()
                pos = map(float, pos)
            except:
                print "**** Unable to parse CCD data (UCCD2)"
                return

            object.ccdAcquireTime = pos[0]
            object.ccdAcquirePeriod = pos[1]
            object.ccdNumExposures = int(pos[2])
            object.ccdNumImages = int(pos[3])
            if hasattr(object, 'ccdSubImages'):
                print "---- CCD Sub Images set to %d" % object.ccdSubImages
                object.ccdNumAcquisitions = object.ccdSubImages
            else:
                object.ccdNumAcquisitions = int(pos[4])
            
        return

    def concatenateSpecScan(self, object, a):
        object.ccdAcquireTime = numpy.concatenate((object.ccdAcquireTime, a.ccdAcquireTime))
        object.ccdAcquirePeriod = numpy.concatenate((object.ccdAcquirePeriod, a.ccdAcquirePeriod))
        object.ccdNumExposures = numpy.concatenate((object.ccdNumExposures, a.ccdNumExposures))
        object.ccdNumImages = numpy.concatenate((object.ccdNumImages, a.ccdNumImages))
        object.ccdNumAcquisitions = numpy.concatenate((object.ccdNumAcquisitions, a.ccdNumAcquisitions))
        object.ccdFilenames = object.ccdFilenames + a.ccdFilenames
        object.ccdDarkFilenames = object.ccdDarkFilenames + a.ccdDarkFilenames
        return 
    def postProcessSpecScanHeader(self, object):

        # Define parameters from header
        object.ccdAcquireTime = numpy.ones(object.data.shape[0]) * object.ccdAcquireTime
        object.ccdAcquirePeriod = numpy.ones(object.data.shape[0]) * object.ccdAcquirePeriod
        object.ccdNumExposures = numpy.ones(object.data.shape[0], dtype = numpy.int) * object.ccdNumExposures 
        object.ccdNumImages = numpy.ones(object.data.shape[0], dtype = numpy.int) * object.ccdNumImages
        object.ccdNumAcquisitions = numpy.ones(object.data.shape[0], dtype = numpy.int) * object.ccdNumAcquisitions

        #Define CCD Filenames

        if object.datafile.ccdpath is not None:
            _path = object.datafile.ccdpath + os.sep
        else:
            _path = ""

        _datafile = object.datafile.filename.split(os.sep)

        if object.data.ndim == 1:
            ndps = 1
        else:
            ndps = object.data.shape[0]

        allfilenames = []
        for dark in ['', object.datafile.ccddark]:
            filenames = []
            for (i, scan, cna) in zip(object.scandatum, object.scanno, object.ccdNumAcquisitions):
                _fnames = []
                for j in range(cna):
                    _f = "%s_%04d-%04d%s_%04d%s" % (_datafile[-1], 
                                                    scan, 
                                                    i, 
                                                    dark, 
                                                    j,
                                                    object.datafile.ccdtail)
                    _fnames.append("%s%s" % (_path, _f))
                filenames.append(_fnames)
            allfilenames.append(filenames)
        object.ccdFilenames = allfilenames[0]
        object.ccdDarkFilenames = allfilenames[1]
    

