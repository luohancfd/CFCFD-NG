#    filename: caldata.py
# description: A class derived from YvX for handling calibration vectors (one factor per wavelength)

from YvX import YvX
from numpy import average

class calData(YvX):
    """A class derived from YvX for handling calibration vectors (one factor per wavelength)"""
    def __init__(self,infile_name):
        # grab solid angle data from calibration file then use YvX constructor for the data
        calFile = open( infile_name, 'r' )
        line = calFile.readline()
        self.f_optics = -99; self.pixels_per_mm = -99
        while line:
            tks = line.split()
            if len(tks) < 2:
                print "unexpected format for calibration file"
                sys.exit()        
            elif tks[0]=="#" and tks[1]=="f_optics":
                self.f_optics = float(tks[3])
                line = calFile.readline()
                continue
            elif tks[0]=="#" and tks[1]=="pixels_per_mm":
                self.pixels_per_mm = float(tks[3])
                line = calFile.readline()
                continue
            elif tks[0]=="#":
                line = calFile.readline()
                continue
            if self.f_optics<0:
                print "failed to find '# f_optics = ...' line in the calibration file."
                sys.exit()
            if self.pixels_per_mm<0:
                print "pixels_per_mm = ", pixels_per_mm
                print "failed to find '# pixels_per_mm = ...' line in the calibration file."
                sys.exit()
            break
        calFile.close()
        # parse infile name to the YvX constructor to do its work
        self.create_from_file( infile_name )
        # apply optical factor to counts
        self.y_array *= self.f_optics   # W/m**2-nm -> W/cm**2-micron-sr
        self.recompute_spline()
        print "Successfully created a calibData class from file: ", infile_name
        print "Summary: f_optics = %f m2-nm/cm2-um-sr, pixels_per_mm = %f p/mm, f_cal_av = %f W/cm2-um-sr\n" % ( self.f_optics, self.pixels_per_mm, average(self.y_array) )