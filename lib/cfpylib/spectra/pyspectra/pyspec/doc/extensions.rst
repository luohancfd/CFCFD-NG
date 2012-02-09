=============================================
Defining Extensions to the Spec Data Routines
=============================================

**************
The base class
**************

The base class for each spec extension is defined from the ``SpecExtension`` base class::

    class SpecExtension:
        """Class to define extensions to SpecDataFile"""
        def __init__(self):
            return
        def getName(self):
            """Return string of name of extension"""
            return "Dummy"
        def initSpec(self, object):
            """Initialize SpecDataFile class"""
            return
        def initSpecScan(self, object):
            """Initialize SpecScan class"""
            return
        def parseSpecHeader(self, object, line):
            """Parse a line of a spec header and modify object with results"""
            return
        def parseSpecScanHeader(self, object, line):
            """Parse a line of a spec scan header and modify object with results"""
            return
        def postProcessSpecHeader(self, object):
            """Post process spec header"""
        return
        def postProcessSpecScanHeader(self, object):
            """Post process spec scan header"""
            return
        def concatenateSpecScan(self, object, a):
	    """Concatenate"""
            return


