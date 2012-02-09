"""Interactive user module for the pyspec package"""

import IPython.ipapi
import pyspec.spec as spec
import tempfile
import matplotlib.pyplot as pyplot

__version__   = "$Revision: 176 $"
__author__    = "Stuart B. Wilkins <stuwilkins@mac.com>"
__date__      = "$LastChangedDate: 2011-02-14 05:08:16 +1000 (Mon, 14 Feb 2011) $"
__id__        = "$Id: ipython.py 176 2011-02-13 19:08:16Z stuwilkins $"

def magic_printfig(self, args):
	"""Print current figure"""
	api = self.api
	fo = tempfile.NamedTemporaryFile(delete = False)
	print fo.name
	pyplot.savefig(fo)
	fo.close()
	ip = IPython.ipapi.get()
	ip.system('ls -l %s' % tempfile)

def magic_loadspec(self, args):
	api = self.api
	api.ex("SPECFILE = spec.SpecDataFile(%s)" % args)

def magic_getspec(self,args):
	api = self.api
	api.ex("SPECSCAN = SPECFILE[%s]" % args)
	#api.ex("SPECPLOT = SPECSCAN.plot()")
	
	# Now load variables
	
	specscan = ip.user_ns['SPECSCAN']
	for i in specscan.scandata.values.keys():
		foo = specscan.scandata.values[i]
		eval("ip.user_ns.update(dict(%s=foo))" % i)
	
def magic_plotspec(self, args):
	api = self.api
	api.ex("SPECPLOT = SPECSCAN.plot()")
	
def magic_prints(self,args):
	api = self.api
	api.ex("SPECPLOT.prints()")

def magic_fitto(self,args):
	api = self.api
	
	func = args.split()
	
	f = "[%s" % func[0]
	
	for i in func[1:]:
		f = f + ", %s" % i 
		
	f = f + "]"
	
	print "---- Fitting to %s" % f
	api.ex("SPECPLOT.fit(%s)" % f)
	
def magic_reload(self, args):
	api = self.api
	api.ex("SPECFILE.index()")
	
def magic_header(self, args):
	self.api.ex("print SPECSCAN.header")
	
def install():
	ip = IPython.ipapi.get()
	ip.expose_magic('getspec', magic_getspec)	
	ip.expose_magic('loadspec', magic_loadspec)
	ip.expose_magic('prints', magic_prints)
	ip.expose_magic('fitto', magic_fitto)
	ip.expose_magic('updatespec', magic_reload)
	ip.expose_magic('header', magic_header)
	ip.expose_magic('plotspec', magic_plotspec)

	ip.expose_magic('printfig', magic_printfig)
	
	ip.ex('import pyspec.spec as spec')
	ip.ex('import pyspec.fitfuncs as fitfuncs')

install()
