#! spectra_example.py
# Example usage for FZ's spectra stuff
#
# Fabian Zander
# February 9th 2012
# Dungeon, Hawken Bldg

# Import the required functions
from pyspec.ccd.files import PrincetonSPEFile as speFile
import spectral_functions as sf

# Read in the example spectra file supplied
# This is an image taken of the flouros in the x-labs optics lab
# using the visible near infrared spectrometer
spe = speFile('fluoro_spectra.spe')

# Extract the wavelengths of the recorded image
wave = spe.getWavelengths()

# Extract the data into an array
data = spe.getData()

'''
The functions in spectral_fucntions.py are all written for the 
UQ specific system and may not translate directly to other 
systems.
'''

# Plot a basic spectral image
sf.plot_2d_spectra(wave, data)
# This can also be done directly using
# sf.plot_2d_spectra(spe.getWavelengths(), spe.getData())
# however I prefer to extract the data first. I'll show the
# same for the other plotting examples as well.

# Plot a single wavelength against slit position
sf.plot_wave(546, wave, data)
#sf.plot_wave(546, spe.getWavelengths(), spe.getData())
# As expected for an image of a flouro this is a 'straight'
# line (this is not a great spectra....).

# Plot a single slit position against wavelength
sf.plot_slit_pos(128, wave, data)
#sf.plot_slit_pos(128, spe.getWavelengths(), spe.getData())
