*************************
CCD Image Format Routines
*************************

Introduction
------------

This section governs the various routines to read CCD (2D Image file
formats) Currently the following formats are implemented.

* Princeton Instruments SPE Files (Winview / WinSpec)

Princeton SPE Files
-------------------

Princeton files can be read using the following example code::

   >>> from pyspec.ccd.files import PrincetonSPEFile
   >>> f = PrincetonSPEFile('myfile.spe')
   >>> f.getData()
   array([[[2855, 3125, 3041, ..., 1366, 1371, 1325],
        [2853, 3074, 3073, ..., 1384, 1362, 1371],
        [2830, 3098, 3105, ..., 1373, 1386, 1358],
        ..., 
        [3086, 3376, 3375, ..., 1442, 1424, 1405],
        [3126, 3422, 3492, ..., 1450, 1409, 1390],
        [3232, 3592, 3561, ..., 1501, 1486, 1467]]], dtype=uint16)
   >>> f.getData().shape
   (1, 1340, 1300)
   >>> f.getSize()
   (1, 1340, 1300)

The class documentation is provided below for convinience:
   

.. autoclass:: pyspec.ccd.files.PrincetonSPEFile
   :members:

