*********************************
Module to process SPEC data files
*********************************

============
Introduction
============

This module defines classes, objects and routines for reading SPEC datafiles
into python classes. There are three main classes for the data, SpecDataFile,
SpecScan and SpecData and are arranged as follows:

SpecDataFile, contains a single data file. The datafile is indexed upon
initializing the class to allow for fast searching of large data files.

To retrieve a scan from the datafile, an item is requested from the object. For
example::

    >>>sd = SpecDataFile('data.01')
    >>>scan = sd[100]

The construct sd[100] will return scan 100 as a SpecScan object. By passing a range
of integers, for example::

    >>>scans = sd[[100, 101]]

these data will be concatenated into a single object. This is useful, for example,
when two scans are used to take the data (one wide, one close) and you need to fit
the data.

The SpecScan object contains members for all motors and counters for the scan in
question. These can be accessed in two main ways. For example, to obtain the
detector counts for the scan as an array::

   >>>Det = sd[1000].Detector

or::

   >>>scan = sd[1000]
   >>>Det = scan.Detector

can be used. For situations where you need machine adjustable code, the construct::

   >>>Det = scan.values['Detector']

can be used, where the member 'values' is a python dictionary of all scan variables.
In addition to the counters and motors, values defined in the header are also included.
See the SpecScan class documentation for further information.

For example, to plot (using matplotlib) the *theta* values against the
*detector* normalized to the *monitor* one could run the code::

   >>>sd = SpecDataFile('data.01')
   >>>scan = sd[100]
   >>>plot(scan.Theta, scan.Detector / scan.Monitor)

and error-bars can be plotted using::

   >>>x = scan.Theta
   >>>y = scan.Detector / scan.Monitor
   >>>e = sqrt(scan.Detetor) / scan.Monitor
   >>>errorbar(x, y, yerr = e)

One can easily perform calculations on the data with this setup. For
example, if measurements are taken with a photodiode, often the x-ray
flux is the required value. The ``pyspec.utils`` contains the function
:func:`itophotons()` which can be used to calculate the flux. For
example:: 

   >>>from pyspec.utils import itophotons
   >>>x = scan.Theta
   >>>y = itophotons(scan.pdiode, scan.energy)
   >>>plot(x, y, 'ro')

Finally, the ``SpecScan`` class defines some helper routines for
plotting, fitting etc. which can be used to help in data analysis. See
the documentation on the ``SpecScan`` class for more
information. These classes form the basis for the SPEC datafile
handling. There are three objects, the ``SpecDataFile`` the
``SpecScan`` and the ``SpecData``. 

=============
Plotting Data
=============

There are a number of helper routines which allow for quick plotting
of scans. These try to guess the best they can what the last scan
was viewing. This is done like CPLOT [#cplot]_ where the first column
in the data file is assumed to be the scan axis and the last column
the counter. Data is normalized by a beam intensity monitor which is
assumed to be the second from last column. To plot a scan from a spec
file::

   >>>sd = SpecDataFile('data.01')
   >>>scan = sd[100]
   >>>scan.plot()

The default behavior can be overridden by a number of keyword
arguments passed to the plot function. These allow for the user to
easily plot different counters, or to normalize to a different
monitor. This is accomplished by the use of the ``xcol``, ``ycol`` and
``mcol`` keyword arguments. These arguments can be either a column
number (negative numbers count from the last column) or a column
(counter) name. For example, to plot the *Theta* axis against *mdet1*
with *Ring* as the monitor::

   >>>scan.plot(xcol='Theta', ycol='mdet1', mcol = 'Ring')

The ``scan.plot()`` function is really the overloaded
``SpecPlot.show()`` function, and returns a ``SpecPlot`` instance. The
possible keyword arguments are defined by this function and are listed below:

.. autoclass:: pyspec.spec.SpecPlot
   :members: show

=============================
Obtaining Information on Data
=============================

You can obtain data on the scan, to find out which variables are
defined in the class for a give scan. This is accomplished using the
``str`` method. Whenever a string is requested of the class, a string
containing a formatted output of the information is provided. For example:: 

   >>>print sd[100]
   Scan:

	500

   Datafile:

	/home/tardis/swilkins/Dropbox/Data/CDI_2011_02_01/CDI_2011_02_01

   Scan Command:

	round_roi_scan npbx -0.0075 0.0075 npbz -0.0075 0.0075 0.0012 5 0.5

   Scan Constants:

   ccdAcquireTime      Lattice             RLattice            ccdAcquirePeriod    
   scandatum           energy              or0                 or1                 
   ccdNumExposures     or1Lambda           or1Angles           wavelength          
   scan_command        comments            scan_type           or0Angles           
   UB                  header              data                ccdNumImages        
   alphabeta           Qvec                scanno              scan                
   or0Lambda           ccdNumAcquisitions  omega               azimuth             


   Motors:

   Chi                 MIR_ROLL            ExitSlit            NPBY                
   XDiff1              Theta               YS                  XDiff2              
   EXSyh               XDiff               MuT                 Phi                 
   AZS                 MuR                 ZDiff               NPTX                
   EXSh                TMM                 NPTZ                RDiff               
   NPTY                SGM                 Gamma               Delta               
   ExitSlitY           PDiff               BPM_Y               BPM_X                
   Mu                  DiffSlitVert        MIR_X               XS                  
   EntSlit             YDiff               AYS                 ZS                  
   SExitSlit           MIR_YAW             NIX                 ZDiff3              
   ZDiff2              ZDiff1              

   Scan Variables:

   kpts                Detector            DELPos              Current             
   NPBX                Epoch               CCPres              NPBZ                
   Monitor             SGME                ccdnet3             Seconds             
   cnptx               sdet2               IONPress            ccdtot2             
   Monitor2            K                   ccdnet2             mdet1               
   mdet2               H                   L                   Treg                
   Tsam                ccdnet5             Ring                sdet1               
   SGMENC              Shutter             ccdnet1             ccdnet4             
   THPos               ccdtot3             ccdtot1             ccdtot4             
   ccdtot5             



====================
Spec Data File Class
====================

.. autoclass:: pyspec.spec.SpecDataFile
   :members:

===============
Spec Scan Class
===============

.. autoclass:: pyspec.spec.SpecScan
   :members:

.. toctree::
   :maxdepth: 2

   extensions

.. rubric:: Footnotes

.. [#cplot] CPLOT Certified Scientific Software <http://www.certif.com/>.



