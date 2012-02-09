
#
# transformations.py (c) Stuart B. Wilkins 2010 and (c) Sven Partzsch 2010
#
# $Id$
# $HeadURL$
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
import gc
#gc.set_debug(gc.DEBUG_LEAK | gc.DEBUG_STATS)
import time
import sys
import struct
import numpy as np
from   scipy.optimize import leastsq
import matplotlib.pyplot as plt
from   pyspec import fit, fitfuncs
from   pyspec.ccd import utils as ccdutils
from   pyspec.diffractometer import Diffractometer

try:
    from   pyspec.ccd.files import *
except:
    pass

try:
    from   pyspec.ccd.plotter import PlotImages, PlotGrid, PlotGrid2, PlotWindow
except:
    pass

try:
    import pyspec.ccd.ctrans as ctrans
except:
    try:
        import ctrans
    except:
        pass

from ConfigParser import ConfigParser

#gc.set_debug(gc.DEBUG_STATS | gc.DEBUG_LEAK)
gc.enable()

__version__   = "$Revision$"
__author__    = "Stuart B. Wilkins <stuwilkins@mac.com>, Sven Partzsch <SvenPartzsch@gmx.de>"
__date__      = "$LastChangedDate$"
__id__        = "$Id$"

#
# FileProcessor Class
#

class FileProcessor():
    """FileProcessor Class

    This class processes CCD files and returns a numpy array of 
    float values for all the images. Images can be saved to disk
    for faster processing in the future.

    The basic building block of this concept is a list of all filenames
    to process. This should be provided through the helper functions or
    can be set through a specfile object.

    """

    def __init__(self, filenames = None, 
                 darkfilenames = None, 
                 norm = None, format = 'SPE',
                 spec = None, process = False):
        """Initialize the class

        filenames      : list of strings.
           List of filenames for all images (2d list will bin images).
        darkfilenames  : list of strings.
           List of filenames fror all dark images. (2d list will bin images).
        spec           : SpecScan object.
           SpecScan object from which to obtain image file names.
        norm           : ndarray.
           Normalize individual images to this numpy array.
        format         : string.
           Data file format. 'SPE' princeton format.
        process        : Bool.
           If True, then process the images once class is initialized.
        """
        self._format = format
        self.filenames = filenames
        self.darkfilenames = darkfilenames
        self.normData = None
        self.meanMonitor = False

        self._cropOnRead = None
        self._binOnRead = None
        
        if spec is not None:
            self.setFromSpec(spec)
        if norm is not None:
            self.normData = norm
        self.bgndParams = np.array([])

        if process:
            self.process()

    def setCropOnRead(self, roi = None):
        """Sets the portion of the image to crop on loading each image.

        This function sets the region of interest which will be taken from
        each image after reading from disk and before placing in the
        image processor object. This can be used to reduce the size of an image
        stack.

        roi : None or tuple
            The ROI. If None use whole image. The tuple should be formated
            as (xstart, ystart, ystop, ystop). Please note, this uses python
            standard indexing. The stop values are NOT included in the crop.

        """

        self._cropOnRead = roi

    def setCropOnRead(self,bins = None):
        """Sets the portion of the image to bin on loading each image.

        This function will bin the images before loading them into the image processor
        
        bins : None or tuple
               The number of bins in the x and y directions. If None use (1x1) binning.
               The tuple should be formated as (xbins, ybins).
        """

        self._binOnRead = bins

    def setFilenames(self, filenames = None, darkfilenames = None):
        """Set the list of filenames and darkfilenames

        filenames     : list of strings
           Set the list of filenames to process for the light images
        darkfilenames : list of strings
           Set the list of filenames to process the dark images
        """
        self.filenames = filenames
        self.darkfilenames = darkfilenames

    def setMeanMonitor(self, b):
        """Set if the images are normalized by the mean of the monitor

        b : bool
           If True, normalize to mean of monitor counts."""
        
        self.meanMonitor = b

    def setFromSpec(self, scan, mon = 'Monitor'):
        """Set the filenames from a SpecScan instance

        This function sets the following from a SpecScan instance:
           
           CCD Filenames
           CCD Dark Filenames
           Monitor normalization values.

        scan    : SpecScan instance
           SpecScan to use for data
        mon     : String
           Name of monitor in SpecScan instance"""
        
        self.filenames     = scan.ccdFilenames
        self.darkfilenames = scan.ccdDarkFilenames
        self.normData      = scan.values[mon]

    def process(self, dark = True, norm = True, 
                dtype = np.float, quiet = False,
                crop = False, BG = False):
        """Read in images and process into 3D array.
        
        dark  : bool
           If True, subtract the dark images from the data
        norm  : bool
           If True, normalize by monitor.
        dtype : datatype 
           numpy datatype of processed array.
        quiet : bool
           If True, dont write to screen when reading images."""
        
        images = []
        darkimages = []

        if norm:
            if self.normData is None:
                normData = np.ones(len(self.filenames))
                print "XXXX No normalization data found"
            else:
                if self.meanMonitor:
                    normData = np.ones(len(self.filenames)) * self.normData.mean()
                    print "---- Normalizing data (using mean value)."
                else:
                    normData = self.normData
                    print "---- Normalizing data."
        else:
            normData = np.ones(len(self.filenames))

        if dark:
            print "---- Correcting for dark current"

        if quiet:
            print "---- Reading Images"

        if len(self.darkfilenames) == 0:
            self.darkfilenames = self.filenames

        for i, (iname, diname, normVal) in enumerate(zip(self.filenames, self.darkfilenames, normData)):
            if type(iname) == list:
                _images = None
                _darkimages = None 
                #Start reading the light images
                for j, _in in enumerate(iname): 
                    if os.path.exists(_in):
                        image = self._getRawImage(_in).astype(dtype)
                        if _images is not None:
                            _images = _images + image
                        else:
                            _images = image    
                        #_images.append(image)
                        if not quiet:
                            print "---- Reading image %-3d of %-3d (sub image %-3d of %-3d)     \r" % (i + 1, len(self.filenames), j + 1, len(iname)),
                            sys.stdout.flush()
                    else:
                        if not quiet:
                            print "---- Missing image %-3d of %-3d (sub image %-3d of %-3d)\r" % (i + 1, len(self.filenames), j + 1, len(iname)),
                            sys.stdout.flush()

                if dark:
                    for j, _din in enumerate(diname):
                        if os.path.exists(_din):
                            darkimage =  self._getRawImage(_din).astype(dtype)
                            if _darkimages is not None:
                                _darkimages = _darkimages + darkimage
                            else:
                                _darkimages = darkimage
                        #_darkimages.append(darkimage)
                            if not quiet:
                                print "---- Reading dark image %-3d of %-3d (sub image %-3d of %-3d)\r" % (i + 1, len(self.darkfilenames), j + 1, len(diname)),
                                sys.stdout.flush()
                            else:
                                if not quiet:
                                    print "---- Missing dark image %-3d of %-3d (sub image %-3d of %-3d)\r" % (i + 1, len(self.darkfilenames), j + 1, len(diname)),
                                    sys.stdout.flush()

                image = _images
                if _darkimages is not None:
                    darkimage = _darkimages
                    darkimages.append(darkimage)

            else:
                # Process only single image pair
                image = self._getRawImage(iname).astype(dtype)

                #if os.path.exists(diname):
                #    darkimage =  self._getRawImage(diname)#.astype(dtype)
                #    darkimages.append(darkimage)
                if not quiet:
                    print "---- Reading image %-3d of %-3d\r" % (i, len(self.filenames)),

            if dark:
                if len(darkimages):
                    image = image - darkimages[-1]
                else:
                    print "XXXX Unable to dark currect correct. No Image found"
            
            if norm:
                image = image / normVal


            if crop:
                print "---- Cropping Image"
                #image=image[0:len(image)/2]
                print "---- Correcting for background"
                image.resize(250,570)
                TempIm=image.transpose()
                BG = sum(TempIm[:,30:50], axis=1)/20
                BG2= sum(TempIm[:,231:250], axis=1)/20
                BGprofile = []
                for i in range(0,570):
                    for j in range(0,250):
                            BGprofile.append(BG[i]+((BG2[i]-BG[i])/550*(j-10)))
                    TempIm[i,:] = TempIm[i,:]-BGprofile
                    BGprofile=[] 
                image=TempIm#.transpose()

            images.append(image)

        print "\n---- Processed %d images (%d dark images)" % (len(images), len(darkimages))

        self.images = np.array(images)

        print "---- Done. Array size %s (%.3f Mb)." % (str(self.images.shape), 
                                                       self.images.nbytes / 1024**2) 
            
    def _getRawImage(self, iname):
        """Read raw image"""
        
        if self._format == 'SPE':
            img = PrincetonSPEFile(iname).getBinnedData()
        elif self._format == 'LCLS':
            img = LCLSdataformat(iname)
        else:
            raise Exception("Unknown file format \"%s\"" % self._format)

        if self._cropOnRead is not None:
            img = img[self._cropOnRead[0]:self._cropOnRead[2],
                      self._cropOnRead[1]:self._cropOnRead[3]]
        if self._binOnRead is not None:
            img = ccdutils.binArray(img, self._binOnRead)

        return img
    
    def _computeMeanImage(self):
        N = self.images.shape[0]
        self.mean = self.images.sum(0) / N
        self.stderr = (self.images - self.mean)**2

    def getMeanImage(self):
        """Return the mean image after processing.

        Calculates the mean of all images processed. Returns a tuple
        of the mean image and the standard error"""
        _computeMeanImage()
        return self.mean, self.stderr

    def getImage(self, n = None):
        """Return the image data

        n : int or None
           Image number to return. If None, return all images."""

        if n is None:
            return self.images
        else:
            return self.images[n]

    def saveImage(self, filename, inum = None, dtype = np.float32):
        """Save image to binary file of specified datatype

        filename : string
           Filename to save image stack to.
        inum     : Integer
           Image number to save
        dtype    : numpy data type
           Datatype for images

        """
        if inum is None:
            self.images.astype(dtype).tofile(filename)
        else:
            self.images[n].astype(dtype).tofile(filename)
                               

    def save(self, filename, compressed = False):
        """Save an image sequance to a numpy binary file
        
        filename   : string
           Filename to save to.
        compressed : bool
           If True save as compressed file"""
        
        obj = {'images'   : self.images,
               'normData' : self.normData}

        np.savez(filename, **obj)
        print "**** Image saved to %s" % filename

    def load(self, filename):
        """Load images from previous save operation

        filename : string
           Filename of images previously stored by save()"""
        
        print "**** Loading image from %s" % filename
        obj = np.load(filename)
        self.images = obj['images']
        self.normData = obj['normData']

    def __str__(self):
        """Return text string of FileProcessor stats"""
        

#
# image processor class
#

class ImageProcessor():
    """Image Processor class

    This class provides the processing of single and sets of CCD-images
    in more detail: each pixel is transformed into reciprocal space
    the set of reciprocal vectors and intensities is gridded on a regular cuboid
    the needed informations can be provided by a spec scan from the pyspec package"""

    def __init__(self, fP = None, configfile = "ccd.cfg", ccdname = None, spec = None):
        """Initialize the image processor

        fP         : pyspec.ccd.transformations.FileProcessor instance
                     File processor object for getting the images
        configfile : string
                     File defining config of CCD
        ccdname    : string
                     Name of CCD as defined in config file. If none, then 
                     do not process config file
        spec       : pyspec.spec.SpecScan instance
                     Spec Scan used to obtain information on current set"""
        
        # file processor to provied the needed images
        if fP is not None:
            self.setFileProcessor(fP)
             
        # set parameters to configure the CCD setup
        # detector distance 30cm and detector pixel size 20um
        self.setDetectorProp(0.020, 0.020, 1300, 1340, 650, 670)
        self.setDetectorPos(300, 0.0)
        self.setDetectorMask()
        
        # Overide if a config file is used.
        #self.ccdname = ccdname
        #self.ccdConfigFile = configfile
        #if ccdname is not None:
        #    self.readConfigFile(configfile)
        
        # Default to frame mode 1
        self.setFrameMode(1)
        
        self.totSet   = None
        
        self.gridData   = None
        self.gridOccu   = None
        self.gridOut    = None
        self.gridStdErr = None

        self.maskData   = None
        self.maskOccu   = None

        self.Qmin  = None
        self.Qmax  = None
        self.dQN   = None

    def readConfigFile(self, filename):
        """Read config file and set detector parameters

        This routine will read the detector vital statistics from 
        a config file of name 'filename'. If the config file does not
        exitst in the path, then standard locations are searched.

        filename : string
           filename of config file.

        """
        config = CCDParamsConfigParser()
        fname = config.readAllLocations(filename)
        print "---- Reading config '%s' from %s" % (self.ccdName, fname)
        xsize = config.getInt(self.ccdName, 'xsize', 1300)
        ysize = config.getInt(self.ccdName, 'ysize', 1340)
        pixx  = config.getFloat(self.ccdName, 'pixel_xsize', 0.020)
        pixy  = config.getFloat(self.ccdName, 'pixel_ysize', 0.020)
        xcen  = config.getFloat(self.ccdName, 'xcen', 650.0)
        ycen  = config.getFloat(self.ccdName, 'ycen', 570.0)
        dist  = config.getFloat(self.ccdName, 'sample_det_distance', 300.0)
        rot   = config.getFloat(self.ccdName, 'detector_rotation', 0.0)
        self.setDetectorProp(pixx, pixy, xsize, ysize, xcen, ycen)
        self.setDetectorPos(dist, rot)
        
    #
    # set and get part
    #

    def getGrid(self):
        """Convinience function to return useful grid values

        This function returns the X, Y and Z coordinates of the grid
        along with the intensity and standard error grids. For example::

            >>>X, Y, Z, I, E, N = ip.getGrid()

        """
        X, Y, Z = self.getGridMesh()
        I = self.getGridData()
        E = self.getGridStdErr()
        N = self.getGridOccu()
        return X, Y, Z, I, E, N

    def getImageData(self):
        """Get the totat image data, transformed into Q"""
        return self.totSet

    def getGridData(self):
        """Return the intensity grid"""
        return self.gridData

    def getGridStdErr(self):
        """Return the standard error grid"""
        return self.gridStdErr

    def getGridStdDev(self):
        """Return the standard deviation grid"""
        return self.gridStdErr * sqrt(self.gridOccu)
    
    def getGridOccu(self):
        """Return the occupation of the grid"""
        return self.gridOccu

    def setDetectorProp(self, detPixSizeX, detPixSizeY, detSizeX, detSizeY, detX0, detY0):
        """Set properties of used detector

        detPixSizeX : float (mm)
             Detector pixel size in detector X-direction.
        detPixSizeY : float (mm)
             Detector pixel size in detector Y-direction.
        detSizeX    : integer
             Detector no. of pixels (size) in detector X-direction.
        detSizeY    : integer
             Detector no. of pixels (size) in detector Y-direction.
        detX0       : float
             Detector X-coordinate of center for reference.
        detY0       : float
             Detector Y-coordinate of center for reference."""
        
        self.detPixSizeX = detPixSizeX  
        self.detPixSizeY = detPixSizeY
        self.detSizeX    = detSizeX
        self.detSizeY    = detSizeY
        self.detX0       = detX0
        self.detY0       = detY0

    def getDetectorProp(self):
        """Get properties of used detector returned as a tuple

        detPixSizeX : float (mm)
             Detector pixel size in detector X-direction.
        detPixSizeY : float (mm)
             Detector pixel size in detector Y-direction.
        detSizeX    : integer
             Detector no. of pixels (size) in detector X-direction.
        detSizeY    : integer
             Detector no. of pixels (size) in detector Y-direction.
        detX0       : float
             Detector X-coordinate of center for reference.
        detY0       : float
             Detector Y-coordinate of center for reference."""

        return self.detPixSizeX, self.detPixSizeY, self.detSizeX, self.detSizeY, self.detX0, self.detY0

    def setDetectorPos(self, detDis = 300.0, detAng = 0.0):
        """Set the detector position

        detDis : float (mm)
           Detector distance from sample.
        detAng : float (deg)
           Detector miss alignement angle (rotation around incident beam)"""
        
        self.detDis = detDis
        self.detAngle = detAng

    def getDetectorPos(self, detDis = 30.0, detAng = 0.0):
        """Get the detector position

        detDis : float (mm)
           Detector distance from sample.
        detAng : float (deg)
           Detector miss alignement angle (rotation around incident beam)"""
        
        return self.detDis, self.detAngle

    def setDetectorMask(self, mask = None):
        """Set the detector mask, to exclude data points

        mask : ndarray
           Array of same form as CCD images, with 1's for valid data points

        """
        self.detMask = mask

    def setBins(self, binX = 1, binY = 1):
        """Set detector binning.

        This function sets the detector binning, and modifies the detector paremeters accordingly.

        binX : integer
           Binning in x-direction.
        binY : integer
           Binning in y-direction.
           
        """
    
        # try to get old bins
        try:
            oldBinX = self.binX
        except:
            oldBinX = 1
        try:
            oldBinY = self.binY
        except:
            oldBinY = 1

        # scaling ratios
        ratX = 1.0*binX/oldBinX
        ratY = 1.0*binY/oldBinY

        # apply changes to detector probperties
        self.detPixSizeX *= ratX
        self.detPixSizeY *= ratY
        self.detSizeX     = int(self.detSizeX / ratX)
        self.detSizeY     = int(self.detSizeY / ratY)
        self.detX0       /= ratX
        self.detY0       /= ratY

        self.binX = binX
        self.binY = binY
        
    def getBins(self, binX, binY):
        """Set no. of bins. Takes them into acount for pixel size, detector size and detector center
        
        binX : no. of pixels along detector X-direction which are bined
        binY : no. of pixels along detector Y-direction which are bined"""

        return self.binX, self.binY

    def setSetSettings(self, waveLen, settingAngles, UBmat):
        """Set the settings for the set 

        The set settings are:
        waveLen       : float
           Wavelength of the X-rays (Angstroms).
        settingAngles : setting angles of the diffractometer at each image (data point)
        UBmat         : UB matrix (orientation matrix) to transform the HKL-values into the sample-frame (phi-frame)
        """

        self.waveLen       = waveLen
        self.settingAngles = settingAngles
        self.UBmat         = UBmat

    def setSpecScan(self, conScan):
        """Set the settings for the set from the considered pyspec scan object

        conScan : pyspec scan object which contains the needed information

        The set settings are:
        waveLen       : wavelength of the X-rays  (float in Angstrom)
        energy        : photon energy             (float in eV)
        imFileNames   : filenames for each image
        darkFileNames : filenames of the dark images
        settingAngles : setting angles of the diffractometer at each image (data point)

        The file Processor processes the corresponding images"""

        self.conScan = conScan
        self._setFromSpecScan()
   
    def getSpecScan(self):
        """Get the pyspec scan object which was used for the set settings
        
        returns pyspec scan object which contains the needed information"""

        return self.conScan

    def setFrameMode(self, mode):
        """Set the mode of the output frame for (Qx, Qy, Qz)

        mode : Integer
           Mode (see below)

        The image processor uses a number of modes which defile the
        coordinates of the final grid. These are:

        mode 1 : 'theta'    : Theta axis frame.  
        mode 2 : 'phi'      : Phi axis frame.
        mode 3 : 'cart'     : Crystal cartesian frame.
        mode 4 : 'hkl'      : Reciproal lattice units frame."""

        if mode == 'theta':
            self.frameMode = 1
        elif mode == 'phi':
            self.frameMode = 2
        elif mode == 'cart':
            self.frameMode = 3
        elif mode == 'hkl':
            self.frameMode = 4
        else:
            self.frameMode = mode

    def getFrameMode(self):
        """Get the mode of the output frame for (Qx, Qy, Qz)

        mode 1 : 'theta'    : Theta axis frame.  
        mode 2 : 'phi'      : Phi axis frame.
        mode 3 : 'cart'     : Crystal cartesian frame.
        mode 4 : 'hkl'      : Reciproal lattice units frame."""

        return self.frameMode

    def setGridOptions(self, *args, **kwargs):
        """Depreciated function

        .. warning::
           This function is depreciated. See :func:setGridSize()
        """
        raise Exception("Depriciated Function, use : setGridSize()")

    def getGridOptions(self, *args, **kwargs):
        """Depreciated function

        .. warning::
           This function is depreciated. See :func:getGridSize()
        """
        raise Exception("Depriciated Function, use : getGridSize()")

    def setGridSize(self, Qmin = None, Qmax = None, dQN = None):
        """Set the options for the gridding of the dataset

        Qmin : ndarray
           Minimum values of the cuboid [Qx, Qy, Qz]_min
        Qmax : ndarray
           Maximum values of the cuboid [Qx, Qy, Qz]_max
        dQN  : ndarray
           No. of grid parts (bins)     [Nqx, Nqy, Nqz]"""

        if Qmin is not None:
            self.Qmin = np.array(Qmin)
        if Qmax is not None:
            self.Qmax = np.array(Qmax)
        if dQN is not None:
            self.dQN  = np.array(dQN)

    def getGridSize(self):
        """Get the options for the gridding of the dataset

        returns tuple of (Qmin, Qmax, dQN)"""

        return self.Qmin, self.Qmax, self.dQN

    def getGridMesh(self):
        """Return the grid vectors as a mesh.

        This function returns the X, Y and Z coordinates of the grid as 3d
        arrays.

        Example:

        X, Y, Z = ip.getGridMesh()

        These values can be used for obtaining the coordinates of each voxel.
        For instance, the position of the (0,0,0) voxel is given by

        x = X[0,0,0]
        y = Y[0,0,0]
        z = Z[0,0,0]"""
        
        grid = np.mgrid[0:self.dQN[0], 0:self.dQN[1], 0:self.dQN[2]]
        r = (self.Qmax - self.Qmin) / self.dQN

        X = grid[0] * r[0] + self.Qmin[0]
        Y = grid[1] * r[1] + self.Qmin[1]
        Z = grid[2] * r[2] + self.Qmin[2]
        
        return X, Y, Z

    #
    # get set functions for input output
    #

    def setFileProcessor(self, fp = None):
        """Set the FileProcessor object for treating the CCD images

        fp : FileProcessor instance
           FileProcessor to use to obtain CCD images"""
        
        self.fileProcessor = fp

    def getFileProcessor(self):
        """Get the FileProcessor object for treating the CCD images

        fp : FileProcessor object with the CCD images"""
        
        return self.fileProcessor

    #
    # help function part
    #

    def _calcBMatrix(self, angles):
        """Calculate the B matrix from reciprocal space angles

        anges: ndarray
           Array of real space values, [a, b, c, alpha, beta, gamma]

        returns B matrix as (3x3) ndarray.
        """
        B = np.ones((3,3))
        B[0,0] = angles[0]
        B[1,0] = 0
        B[2,0] = 0
        B[0,1] = angles[1] * cos(angles[5])
        B[1,1] = angles[1] * sin(angles[5])
        B[2,1] = 0
        B[0,2] = angles[2] * cos(angles[4])
        B[1,2] = -1.0 * angles[2] * sin(angles[4]) * cos(angles[3])
        B[2,2] = 2 * np.pi / angles[2]

        return B
    
    def _setFromSpecScan(self):
        """Set the settings for the set from the considered pyspec scan object

        The set settings are:
        waveLen       : wavelength of the X-rays  (float in Angstrom)
        energy        : photon energy             (float in eV)
        imFileNames   : filenames for each image
        darkFileNames : filenames of the dark images
        settingAngles : setting angles of the diffractometer at each image (data point)
        intentNorm    : normalization factor (division) for each image
        UBmat         : UB matrix (orientation matrix) to transform the HKL-values into the sample-frame (phi-frame)
        setName       : name of the considered set, e.g. 'Scan #' in the spec case
        setNum        : no. to determine the set, e.g. 244 in the spec case
        setSize       : no. of images in the set, e.g. 81"""

        self.waveLen       = self.conScan.wavelength  # in Angstrom
        self.energy        = Diffractometer.hc_over_e / self.conScan.wavelength # in eV
        self.imFileNames   = self.conScan.ccdFilenames
        self.darkFileNames = self.conScan.ccdDarkFilenames
        self.settingAngles = self.conScan.getSIXCAngles()
        self.intentNorm    = self.conScan.Ring
        self.UBmat         = self.conScan.UB
        self.setName       = 'Scan #'
        self.setNum        = self.conScan.scan
        self.setSize       = self.settingAngles.shape[0]
        self.procImSelect  = range(self.setSize)
        
    
    def _XYCorrect(self, xVal, yVal):
        """Correct the miss alignement of the CCD camera

        xVal : measured X-values of the detector
        yVal : measured Y-values of the detector

        return
        xNew : corrected X-values of the detector
        yNew : corrected Y-values of the detector"""

        # detetoc angle in rad
        detAn = self.detAngle /180.0 * np.pi

        xNew = np.cos(detAn) * xVal - np.sin(detAn) * yVal
        yNew = np.sin(detAn) * xVal + np.cos(detAn) * yVal

        return xNew, yNew

    #
    # help functions for input / output
    #

    def __str__(self):
        """Create the information about the grid process"""
 
        _s = "Detector Parameters:\n\n"

        _s += "Detector size        :%d x %d [pixels]\n" % (self.detSizeX, self.detSizeY)
        _s += "Detector pixel size  :%f x %f [mm]\n" % (self.detPixSizeX, self.detPixSizeY)
        _s += "Zero (center) of CCD :(%f, %f) [pixels]\n" % (self.detX0, self.detY0)
        _s += "Sample Detector dist.:%f [mm]\n" % self.detDis
        _s += "Detector rot. ang.   :%f [deg]\n" % self.detAngle

        if self.totSet is not None:
            _s += "\n\nData Set:\n" 

            _s += "Number of Pixels : %.2e\n" % self.totSet.shape[0]
            _s += 'Q_min            : [%.3e %.3e %.3e]\n' % (self.totSet[:,0].min(), 
                                                             self.totSet[:,1].min(), 
                                                             self.totSet[:,2].min())
            _s += 'Q_max            : [%.3e %.3e %.3e]\n' % (self.totSet[:,0].max(), 
                                                             self.totSet[:,1].max(), 
                                                             self.totSet[:,2].max())
            _s += 'I_min            : %.3e\n' % self.totSet[:,3].min()
            _s += 'I_max            : %.3e\n' % self.totSet[:,3].max()
            _s += 'I_ave            : %.3e\n' % self.totSet[:,3].mean()

        _s += "\n\nGrid Parameters:\n\n"
        mode = self.frameMode
        if mode == 1:
            _s += 'Frame Mode 1: (Qx, Qy, Qz) in theta-frame and (1/Angstrom)\n' 
        elif mode == 2:
            _s += 'Frame Mode 2: (Qx, Qy, Qz) in phi-frame and (1/Angstrom)\n'
        elif mode == 3:
            _s += 'Frame Mode 3: (Qx, Qy, Qz) in cartesian-frame and (1/Angstrom)\n'
        elif mode == 4:
            _s += 'Frame Mode 4: (H, K, L) in hkl-frame and (reciprocal lattice units)\n'

        _s += '\n\nGrid Vectors:\n\n'
        if self.Qmin is not None:
            _s += 'Q_min     : [%.3e %.3e %.3e]\n' % (self.Qmin[0], self.Qmin[1], self.Qmin[2])
        if self.Qmax is not None:
            _s += 'Q_max     : [%.3e %.3e %.3e]\n' % (self.Qmax[0], self.Qmax[1], self.Qmax[2])
        if self.dQN is not None:
            _s += 'Grid Size : [%.3e %.3e %.3e]\n' % (self.dQN[0], self.dQN[1], self.dQN[2])

        if self.gridData is not None:
            _s += '\n\nGrid Statistics:\n\n'
            _s += 'No. of voxels in grid           : \t %.2e\n' % (self.gridData.size)
            _s += 'No. of data points outside grid : \t %.2e\n' % (self.gridOut)
            _s += 'No. of empty voxels             : \t %.2e\n' % ((self.gridOccu == 0).sum())
            
        return _s

    #
    # process part
    #
    
    def processToQ(self):
        """Process selcted images of the full set into (Qx, Qy, Qz, I)

        This function takes the images provided by a FileProcessor instance, and
        converts them to a set of (Qx, Qy, Qz, I) values in the frame mode which
        is set prevously in this class.

        """
        ccdToQkwArgs = {}
        if self.totSet is not None:
            del self.totSet
            gc.collect()

        if not self.fileProcessor:
            raise Exception("No FileProcessor specified.")

        # if images not yet processed, do it
        if getattr(self.fileProcessor, 'images', None) is None:
            self.fileProcessor.process()

        if self.settingAngles is None:
            raise Exception("No setting angles specified.")

        print "\n"
        print "---- Setting angle size :", self.settingAngles.shape
        print "---- CCD Size :", (self.detSizeX, self.detSizeY)
        print "**** Converting to Q"
        t1 = time.time()
        self.totSet = ctrans.ccdToQ(angles      = self.settingAngles * np.pi / 180.0, 
                                    mode        = self.frameMode,
                                    ccd_size    = (self.detSizeX, self.detSizeY),
                                    ccd_pixsize = (self.detPixSizeX, self.detPixSizeY),
                                    ccd_cen     = (self.detX0, self.detY0),
                                    dist        = self.detDis,
                                    wavelength  = self.waveLen,
                                    UBinv       = np.matrix(self.UBmat).I,
                                    **ccdToQkwArgs)
        t2 = time.time()
        print "---- DONE (Processed in %f seconds)" % (t2 - t1)
        print "---- Setsize is %d" % self.totSet.shape[0]
        self.totSet[:,3] = np.ravel(self.fileProcessor.getImage())  
       
        if self.detMask is not None:
            print "---- Masking data"
            totMask = self.detMask.copy().ravel()
            for i in range(1, self.settingAngles.shape[0]):
                totMask = concatenate((totMask, self.detMask.ravel()))

            print "---- Mask size ", totMask.shape
            print "---- totSet size ", self.totSet.shape
            
            self.totSet =  self.totSet[(totMask != 0),:]   
            print "---- totSet size ", self.totSet.shape

        # for info file
        #self.opProcInfo += '\n\n**** Image Set processed to %.2e %s sets (Took %f seconds)' % (self.totSet.shape[0], self.setEntLabel, (t2 - t1))

    def process(self):
        """Convinience function to perform all processing"""
        self.processToQ()
        self.processGrid()

    def processGrid(self):
        """Process imageset of (Qx, Qy, Qz, I) into grid 

        This function, process the set of (Qx, Qy, Qz, I) values and grid the data."""
        print "---- Total data is %f MBytes\n" % (self.totSet.nbytes / 1024.0**2)
        
        if self.totSet is None:
            raise Exception("No set of (Qx, Qy, Qz, I). Cannot process grid.")

        # prepare min, max,... from defaults if not set
        if self.Qmin is None:
            self.Qmin = np.array([ self.totSet[:,0].min(), self.totSet[:,1].min(), self.totSet[:,2].min() ])
        if self.Qmax is None:
            self.Qmax = np.array([ self.totSet[:,0].max(), self.totSet[:,1].max(), self.totSet[:,2].max() ])
        if self.dQN is None:
            self.dQN = [100, 100, 100]

        # 3D grid of the data set 
        print "**** Gridding Data."
        t1 = time.time()
        gridData, gridOccu, gridStdErr, gridOut = ctrans.grid3d(self.totSet, self.Qmin, self.Qmax, self.dQN, norm = 1)
        t2 = time.time()
        print "---- DONE (Processed in %f seconds)" % (t2 - t1)
        emptNb = (gridOccu == 0).sum()
        if gridOut != 0:
            print "---- Warning : There are %.2e points outside the grid (%.2e bins in the grid)" % (gridOut, gridData.size)
        if emptNb:
            print "---- Warning : There are %.2e values zero in the grid" % emptNb

       # store intensity, occupation and no. of outside data points of the grid
        self.gridData   = gridData
        self.gridOccu   = gridOccu
        self.gridOut    = gridOut
        self.gridStdErr = gridStdErr

        # masks for the data and occupation no.
        self.maskData = (gridData == 0)
        self.maskOccu = (gridOccu == 0)

    def makeGridData(self, *args, **kwargs):
        """Convinience to call makeGrid"""
        self.process(*args, **kwargs)


##################################################################################################################
##                                                                                                              ##
##  GridProcessorClass to perform cutting of grid                                                               ##
##                                                                                                              ##
##################################################################################################################

class GridProcessor():
    """Class to process grid data

    This class is used to process grid data, and perform line and sum integration from the 
    gridded data. 

    This class assumes that you have processed the data using an ImageProcessor and requires 
    the use of that when initializing the class."""

    def __init__(ip = None):
        if ip is None:
            raise Exception("This class must be initialized with a valid ImageProcessor")
        else:
            self.ip = ip

    def get1DCut(pos, axis):
        """Make a 1D cut of the grid along a principal axis."""
        ax = self._processAxis()

    def _processAxis(axis):
        if type(axis) == int:
            return axis
        elif type(axis) == str:
            if axis.upper() == "X":
                return 0
            elif axis.upper() == "Y":
                return 1
            elif axis.upper() == "Z":
                return 2
            else:
                raise Exception("Invalid string. Axis must be 'X', 'Y', or 'Z'.")
        else:
            raise Exception("Invalid type in axis. Must be integer or string")
        return Null
    
##################################################################################################################
##                                                                                                              ##
##  CCD Config Parser, to read config file from disk                                                            ##
##                                                                                                              ##
##################################################################################################################        
        
    
class CCDParamsConfigParser(ConfigParser):
    """Class to read config file which defines all CCD parameters"""
    def readAllLocations(self, filename):
        
        if not os.path.isfile(filename):
            locations = []
            if os.name is 'posix':
                if os.environ.has_key('HOME'):
                    locations.append(os.environ['HOME'] + os.path.sep + ".pyspec")
                locations.append('/usr/local/pyspec/etc')
                locations.append('/etc/pyspec')
            
            for l in locations:
                _f = l + os.path.sep + filename
                if os.path.isfile(_f):
                    self.read(_f)
                    return _f
        else:
            self.read(filename)
            return filename

        return None
            
    def getWithDefault(self,section, option, default):
        try:
            r = self.get(section, option)
        except ConfigParser.NoOptionError:
            r = default
            print "**** Using default value for %s." % option
        return r
    
    def _getWithConvert(self,_conv_fcn, section, option, default):
        try:
            val = self.getWithDefault(section, option, default)
        except:
            raise Exception("Unable to read option %s from config file." % (option))
        try:
            val = _conv_fcn(val)
        except:
            raise Exception("Unable to convert option %s to correct datatype." % (option))
        return val

    def getFloat(self, *args, **kwargs):
        return self._getWithConvert(float, *args, **kwargs)
    def getInt(self, *args, **kwargs):
        return self._getWithConvert(int, *args, **kwargs)
    def getFloat(self, *args, **kwargs):
        return self._getWithConvert(float, *args, **kwargs)

        
####################################
#
# main program, for test
#
####################################


if __name__ == "__main__":

    from pyspec import spec
    
    sf   = spec.SpecDataFile('/home/tardis/spartzsch/data/ymn2o5_oct10/ymn2o5_oct10_1', 
    			     ccdpath = '/mounts/davros/nasshare/images/oct10')
    scan = sf[360]

    fp = FileProcessor(spec = scan)
    fp.process()

    testData = ImageProcessor(fp)
    #testData.setDetectorPos(detAng = -1.24)
    testData.setBins(4, 4)
    testData.setSpecScan(scan)
    testData.setFrameMode(4)

    testData.process()

    print testData

    
