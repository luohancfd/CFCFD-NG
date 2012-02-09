#
# fitfuncs.py (c) Stuart B. Wilkins 2008
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

import numpy
import time as _time

try:
    import wx as _wx
except ImportError:
    pass
try:
    import pylab as _pylab
except ImportError:
    pass

__version__   = "$Revision$"
__author__    = "Stuart B. Wilkins <stuwilkins@mac.com>"
__date__      = "$LastChangedDate$"
__id__        = "$Id$"

def peakguess(x,y):
    """Guess the "vital statistics" of a peak.
    
    This function guesses the peak's center
    width, integral, height and linear background
    (m, c) from the x and y data passed to it.
    """

    # Calculate linear background

    m = (y[-1] - y[0]) / (x[-1] - x[0])
    c = y[0] - x[0] * m

    # Now subtract this background

    ty = y - (m * x) - c

    # Calculate integral of background
    # subtracted data

    integral = abs(x[1] - x[0]) / 2
    dsum = ty[1:-1].sum()

    integral = integral * (ty[0] + ty[-1] + (2 * dsum))

    xmax = 0
    ymax = 0
    for i in range(len(x)):
        # find max
        if ty[i] > ymax:
            ymax = ty[i]
            xmax = i

    centre = x[xmax]

    lxmax = 0
    rxmax = len(x)

    i = xmax - 1
    while (i >= 0):
        if ty[i] < (ymax / 2):
            lxmax = i
            break
        i -= 1

    i = xmax + 1
    while (i <= len(x)):
        if ty[i] < (ymax / 2):
            rxmax = i
            break
        i += 1

    width = abs(x[rxmax] - x[lxmax])

    height = ty.max()

    out = [centre, width, integral, height, m, c]

    return out

def lor2a(x, p, mode='eval'):
    """Lorentzian squared function defined by area"""
    if mode == 'eval':
        out = ((numpy.sqrt(2)*p[2])/(numpy.pi*p[1]) / (1 + 0.5*((x - p[0])/p[1])**2)**2);
    elif mode == 'params':
        out = ['cent', 'width', 'area']
    elif mode == 'name':
        out = "Lorentzian Squared"
    elif mode == 'guess':
        g = peakguess(x, p)
        out = g[0:3]
    else:
        out = []

    return out
def lor2(x, p, mode='eval'):
    """Lorentzian squared defined by amplitude

    Function:
      :math:`f(x) = p_2\left(\\frac{1}{1 + \left(\\frac{x - p_0}{p_1}\\right)^2}\\right)^2`

    The HWHM is related to the parameter :math:`\Gamma` by the relation:
      :math:`\kappa = \sqrt{\sqrt{2} - 1}\Gamma`

    """
    if mode == 'eval':
        out = p[2] * (1 / (1 + ((x - p[0]) / p[1])**2))**2
    elif mode == 'params':
        out = ['Cent', 'Gamma', 'Height']
    elif mode == 'name':
        out = "Lorentzian Squared"
    elif mode == 'guess':
        g = peakguess(x, p)
        out = g[0:3]
    elif mode == 'graphguess':
        x = GetGraphInput()
        print "Click on (in order)"
        print "\tPeak center"
        print "\tPeak width (left)"
        print "\tPeak width (right)"
        print "\tBackground at center"
        c = x(4, 0)
        width = abs(c[1][0] - c[2][0])
        area = abs(c[0][1] - c[3][1])
        out = [c[0][0], width, area]
    else:
        out = []

    return out

def linear(x, p, mode='eval'):
    """Linear (strait line)

    Function:
       :math:`f(x) = p_0 x + p_1`

    """
    if mode == 'eval':
        out = (p[0] * x) + p[1]
    elif mode == 'params':
        out = ['grad','offset']
    elif mode == 'name':
        out = "Linear"
    elif mode == 'guess':
        g = peakguess(x, p)
        out = g[4:6]
    else:
        out = []

    return out

def constant(x, p, mode='eval'):
    """Single constant value

    Function:
       :math:`f(x) = p_0`

    """
    if mode == 'eval':
        out = p[0]
    elif mode == 'params':
        out = ['value']
    elif mode == 'name':
        out = "Constant"
    elif mode == 'guess':
        g = peakguess(x, p)
        out = [(g[4] * (x[0] + x[-1]) + 2*g[5]) / 2]
    else:
        out = []

    return out

def lor(x, p, mode='eval'):
    """Lorentzian defined by amplitude

    Function:
       :math:`f(x) = p_2\left(\\frac{1}{1 + \\left(\\frac{x - p_0}{p_1}\\right)^2}\\right)`

    """
    if mode == 'eval':
        out = p[2] / (1 + ((x - p[0]) / p[1])**2)
    elif mode == 'params':
        out = ['Center', 'Gamma', 'Height']
    elif mode == 'name':
        out = "Lorentzian"
    elif mode == 'guess':
        g = peakguess(x, p)
        out = g[0:3]
    else:
        out = []

    return out

def lorr(x, p, mode='eval'):
    """Lorentzian defined by area

    Function:
       :math:`f(x) = \\frac{p_2}{p_1\pi} \left(\\frac{1}{1 + 4\left(\\frac{x - p_0}{p_1}\\right)^2}\\right)`

    """
    if mode == 'eval':
        out = (2*p[2]/p[1]/numpy.pi) / (1 + 4*((x-p[0])/p[1])**2)
    elif mode == 'params':
        out = ['Center', 'width', 'area']
    elif mode == 'name':
        out = "Lorentzian"
    elif mode == 'guess':
        g = peakguess(x, p)
        out = g[0:3]
    else:
        out = []

    return out

def pvoight(x, p, mode='eval'):
    """Pseudo Voight function

    Function:
       :math:'f(x) = '

    """
    if mode == 'eval':
        cent=p[0];wid=p[1];area=p[2];lfrac=p[3];
        out = area / wid / ( lfrac*numpy.pi/2 + (1-lfrac)*numpy.sqrt(numpy.pi/4/numpy.log(2)) ) * ( lfrac / (1 + 4*((x-cent)/wid)**2) + (1-lfrac)*numpy.exp(-4*numpy.log(2)*((x-cent)/wid)**2) );
    elif mode == 'params':
        out = ['cent', 'FWHM', 'area', 'Lorr. Fac.']
    elif mode == 'name':
        out = "Pseudo Voight"
    elif mode == 'guess':
        g = peakguess(x, p)
        out = g[0:3]
        out[3] = 1.0
    else:
        out = []

    return out

def gauss(x, p, mode='eval'):
    """Gaussian defined by amplitide

    Function:
       :math:`f(x) = p_2 \exp\left(\\frac{-(x - p_0)^2}{2p_1^2}\\right)`

    """
    if mode == 'eval':
        cent=p[0];wid=p[1];amp=p[2];
        out = amp * numpy.exp(-1.0 * (x - cent)**2 / (2 * wid**2))
    elif mode == 'params':
        out = ['cent', 'sigma', 'amp']
    elif mode == 'name':
        out = "Gaussian"
    elif mode == 'guess':
        g = peakguess(x, p)
        out = [g[0], g[1] / (4 * numpy.log(2)), g[3]]
    else:
        out = []

    return out

def power(x, p, mode = 'eval'):
    """Power function

    Function:
       :math:`f(x) = p_0 p_1^x`

    """
    if mode == 'eval':
        out = p[0] * pow(p[1], x)
    elif mode == 'params':
        out = ['A', 'B']
    elif mode == 'name':
        out = "Power"
    elif mode == 'guess':
        out = [1, 1]
    else:
        out = []

    return out


def stokes(x, p, mode='eval'):
	if mode == 'eval':
		out = p[0] * 0.5 * (1 + (p[1]*numpy.cos(x*numpy.pi/90)) + (p[2]*numpy.sin(x*numpy.pi/90)))
	elif mode == 'params':
		out = ['int', 'P1', 'P2']
	elif mode == 'name':
		out = "Stokes"
	elif mode == 'guess':
		g = peakguess(x, p)
		out = g[0:3]
	else:
		out = []

	return out	


def stokes_rot(x, p, mode='eval'):
	# ASSUMES COMPLETELY LINEAR LIGHT! (P1^2+P2^2 = 1)
	if mode == 'eval':
		out = p[1]*(1+numpy.cos((x-p[0])*numpy.pi/90))
	elif mode == 'params':
		out = ['Rotation','S0']
	elif mode == 'name':
		out = "Stokes (assuming entirely linear)"
	elif mode == 'guess':
		out = [1,2]
	else:
		out = []
	return out
