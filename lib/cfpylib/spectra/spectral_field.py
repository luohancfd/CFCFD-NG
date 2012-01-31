#    filename: read_ASCII_data.py
# description: A set of classes to read in X2 and EAST ASCII spectrometer data

import sys
from numpy import *
from random import randint
from YvX import YvX
from caldata import calData

class spectralField(object):
    """base class for holding spectral field data"""
    # define the datatypes that this class can own
    __slots__ = 'wavelengths', 'distances', 'spatial_strips', 'specfile_name', \
                'iframe', 'nframe', 'calibrated'
    # -----------------------------------------------------------------------
    def __init__(self, wavelengths=[], spatial_strips=[], specfile_name="none", iframe=0, nframes=0):
        self.wavelengths = array( wavelengths )
        self.distances = array(distances)
        self.spatial_strips = array( spatial_strips )
        self.specfile_name = specfile_name
        self.iframe = iframe
        self.nframes = nframes
        self.calibrated = False
    # -----------------------------------------------------------------------
    def spatially_average( self, pixel_range=[0,-1], x_range=None ):
        if x_range!=None:
            # search distances array for corresponding strips
            istrip = searchsorted( self.distances, x_range[0] )
            fstrip = searchsorted( self.distances, x_range[1] ) - 1
        else:
            # assume a pixel range has been given
            istrip = pixel_range[0]
            fstrip = pixel_range[1]
        if fstrip<0: fstrip = len(self.spatial_strips)+1+fstrip
        nstrips = fstrip - istrip
        av_counts = zeros( len(self.wavelengths) )
        for strip in self.spatial_strips[istrip:fstrip]:
            av_counts += strip
        av_counts /= nstrips
        CvL = YvX(self.wavelengths,av_counts)
        return CvL
    # -----------------------------------------------------------------------
    def get_spatial_strip( self, istrip ):
        return self.spatial_strips[istrip]
    # -----------------------------------------------------------------------
    def get_random_spatial_strip( self ):
        rstrip = randint( 0, len(self.spatial_strips)-1 )
        return self.spatial_strips[rstrip]
    # -----------------------------------------------------------------------
    def spectrally_integrate( self, wavel_i=0.0, wavel_f=99999.0 ):
        if self.calibrated == False:
            print "Spectral field has not been calibrated - will not attempt to spectrally integrate!"
            sys.exit()
        # 0. create an array to store spatial intensity data
        x_intensities = zeros(len(self.spatial_strips))
        # 1. loop over spatial strips
        for ix,xstrip in enumerate(self.spatial_strips):
            # 2. loop over wavelengths
            for iw,wavel in enumerate(self.wavelengths):
                # 3. check if wavelength is in-range
                I = xstrip[iw]
                if wavel>=wavel_i and wavel_i<=wavel_f and iw!=0:
                    # this wavelength is in-range, add integrated contribution
                    x_intensities[ix] += 0.5*( I + I_prev ) * ( wavel - wavel_prev ) * 1.0e-3   # wavel[nm] -> wavel[um]
                I_prev = I; wavel_prev = wavel
        # 4. Create a YvX() class for the returned data
        IvX = YvX(self.distances, x_intensities)
        return IvX
        
class X2spectralField( spectralField ):
    """derived class describing X2 spectral field data"""
    def __init__(self,specfile_name,aspr=[0,999999],awpr=[0,999999],avbr="none",frame="all",bkgfile_name="none"):
        """specfile_name: spectrometer data file name
           aspr: active spatial pixel range
           awpr: active wavelength pixel range
           frame: requested frame index
           avbr: average background wavelength pixel range
           bkgfile_name: background counts file name"""
        # read in the background counts file if one has been given
        if bkgfile_name!="none":
            if avbr!="none":
                print "Both a background file and strip range have been given - choose one!"
                sys.exit()
            bkgField = X2spectralField( bkgfile_name, aspr, awpr )
        # read in and process the spectrometer data file
        specFile = open( specfile_name, "r" )
        # discard "Pixel" line
        line = specFile.readline()
        if line.find('Pixel') < 0:
            print "unexpected entry in data file: ", line
            sys.exit()
        tks = line.split()
        # wavelength line
        line = specFile.readline()
        tks = line.split()
        if line.find('Wavelength') < 0:
            print "unexpected entry in data file: ", line
            sys.exit()
        # at what index do the wavelengths start?
        for i,tk in enumerate(tks):
            if tk=='"Wavelength"':
                tk_start = i+1
                break
        wavelengths_list = []
        for i,tk in enumerate(tks[tk_start:]):
            if i<awpr[0] or i>awpr[1]: continue
            wavelengths_list.append( float(tk) )      # wavelength in nm
        self.wavelengths = array(wavelengths_list)
        nwavels = len(wavelengths_list)
        # "Strip" [counts] lines
        line = specFile.readline()
        tks = line.split()
        if line.find('Strip') < 0:
            print "unexpected entry in data file: ", line
            sys.exit()
        # possibility of multiple frames in one file
        if line.find('Frame') < 0:
            # single frame data
            self.nframes = 1
            self.iframe = 1
            self.spatial_strips = []
            while line:
                line = specFile.readline()
                tks = line.split()
                if len(tks)==0: break
                # check if strip is inside requested range
                istrip = int(tks[0])-1
                if istrip<aspr[0] or istrip>aspr[1]: continue
                counts_list = []
                for i,tk in enumerate(tks[1:]): counts_list.append(float(tk))
                # make background counts array
                if avbr!="none":
                    bkg_counts = array([average(counts_list[avbr[0]:avbr[1]+1])]*nwavels)
                elif bkgfile_name!="none":
                    bkg_counts = bkgField.get_spatial_strip( istrip )
                else: bkg_counts = zeros( nwavels )
                counts = array(counts_list[awpr[0]:awpr[1]+1]) - bkg_counts
                # check that the counts array is the same size as the wavelength array
                if len(counts)!=len(self.wavelengths):
                    print "Dimension mismatch: len(counts) = %d, len(wavelengths) = %d" % ( len(counts), len(self.wavelengths) )
                    sys.exit()
                self.spatial_strips.append(counts)
        else:
            # multi frame data
            if frame!="all":
                print "Only averaging of multiple-frame data is currently supported."
                sys.exit()
            self.iframe = frame
            self.spatial_strips = []
            first_frame = True; prev_frame = 0; self.nframes = 0
            while line:
                line = specFile.readline()
                tks = line.split()
                if len(tks)==0: break
                # check if strip is inside requested range
                istrip = int(tks[1])-1
                if istrip<aspr[0] or istrip>aspr[1]: continue
                counts_list = []
                for i,tk in enumerate(tks[2:]): counts_list.append(float(tk))
                # make background counts array
                if avbr!="none":
                    bkg_counts = array([average(counts_list[avbr[0]:avbr[1]+1])]*nwavels)
                if bkgfile_name!="none":
                    bkg_counts = bkgField.get_spatial_strip( istrip )
                counts = array(counts_list[awpr[0]:awpr[1]+1]) - bkg_counts
                # check that the counts array is the same size as the wavelength array
                if len(counts)!=len(self.wavelengths):
                    print "Dimension mismatch: len(counts) = %d, len(wavelengths) = %d" % ( len(counts), len(self.wavelengths) )
                    sys.exit()
                this_frame = int(tks[0])
                if this_frame!=prev_frame:
                    # print "processing frame: ", self.n_frames
                    self.nframes += 1
                if self.nframes==1:
                    self.spatial_strips.append(counts)
                else:
                    # print "strip = %d, len(spatial_strips) = %d" % ( istrip, len(self.spatial_strips) )
                    self.spatial_strips[istrip] += counts
                prev_frame = this_frame
            # complete averaging by dividing by number of frames
            for i in range(len(self.spatial_strips)):
                self.spatial_strips[i] /= self.nframes
        # finished with specFile
        specFile.close()
        # convert spatial_strips list to an array
        self.spatial_strips = array(self.spatial_strips)
        # initialize the distances array simply as (original) pixel indices
        self.distances = arange( float(aspr[0]), float(aspr[0]+len(self.spatial_strips)), 1.0 )
        # set calibration flag
        self.calibrated = False
        print "len(self.distances) = %d, len(self.spatial_strips) = %d, len(self.wavelengths) = %d, len(self.spatial_strips[0]) = %d" % ( len(self.distances), len(self.spatial_strips), len(self.wavelengths), len(self.spatial_strips[0]) )
        print "Successfully created an X2spectralField class of dimensions: %d x %d from file: %s averaging over %d frames\n" % ( nwavels, len(self.spatial_strips), specfile_name, self.nframes )
    # -----------------------------------------------------------------------
    def apply_calibration_vector( self, calfile_name ):
        """Apply a calibration vector (one calibration factor for each wavelength strip) from file"""
        if self.calibrated:
           print "Data has already been calibrated!"
           sys.exit()
        fvL = calData( calfile_name )
        # apply spectral calibration
        f_array = fvL.y_from_x( self.wavelengths )
        I_peak = -1.0e99
        for strip in self.spatial_strips:
            strip *= f_array            # counts -> W/cm2-um-sr
            if strip.argmax() > I_peak: I_peak = strip.argmax()
        # apply spatial calibration
        self.distances /= fvL.pixels_per_mm * 1000.0
        
        print "Successfully applied calibration vector from file: ", calfile_name
        print "Summary: I_av = %f W/cm2-um-sr, I_peak = %d W/cm2-um-sr, x_range = [%e,%e] m\n" % ( average( self.spatial_strips ), I_peak, self.distances[0], self.distances[-1] )
        self.calibrated = True
