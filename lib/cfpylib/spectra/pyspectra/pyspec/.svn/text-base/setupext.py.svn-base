import os
import sys
import ConfigParser
import numpy as np
import copy
from distutils.core import Extension

__version__   = "$Revision$"
__author__    = "Stuart B. Wilkins <stuwilkins@mac.com>"
__date__      = "$LastChangedDate$"
__id__        = "$Id$"

options = {'build_levmar'    : False ,
           'build_ctrans'   : False ,
           'build_sginfo'    : False }

ext_default  = {'include_dirs' : [np.get_include()],
                'library_dirs' : [],
                'libraries'    : [],
                'define_macros': []}

setup_files = ['setup.cfg.%s' % sys.platform, 'setup.cfg']

def detectCPUs():
    # Linux, Unix and MacOS:
    if hasattr(os, "sysconf"):
        if os.sysconf_names.has_key("SC_NPROCESSORS_ONLN"):
            # Linux & Unix:
            ncpus = os.sysconf("SC_NPROCESSORS_ONLN")
            if isinstance(ncpus, int) and ncpus > 0:
                return ncpus
            else: # OSX:
                return int(os.popen2("sysctl -n hw.ncpu")[1].read())
    # Windows:
    if os.environ.has_key("NUMBER_OF_PROCESSORS"):
        ncpus = int(os.environ["NUMBER_OF_PROCESSORS"]);
    if ncpus > 0:
        return ncpus
    return 1 # Default

def parseExtensionSetup(name, config, default):
    default = copy.deepcopy(default)
    try: default['include_dirs'] = config.get(name, "include_dirs").split(os.pathsep)
    except: pass
    try: default['library_dirs'] = config.get(name, "library_dirs").split(os.pathsep)
    except: pass
    try: default['libraries'] = config.get(name, "libraries").split(",")
    except: pass

    return default

setupfile = None
for f in setup_files:
    if os.path.exists(f):
        setupfile = f
        break

if setupfile is not None:
    print 'Reading config file %s' % setupfile
    config = ConfigParser.SafeConfigParser()
    config.read(setupfile)

    try: options['build_levmar'] = config.getboolean("levmar","build")
    except: pass
    try: options['build_ctrans'] = config.getboolean("ctrans","build")
    except: pass
    try: options['build_sginfo'] = config.getboolean("sginfo","build")
    except: pass

    levmar = parseExtensionSetup('levmar', config, ext_default)
    levmar['libraries'].append('levmar')

    ctrans = parseExtensionSetup('ctrans', config, ext_default)
    threads = False
    try: threads = config.getboolean("ctrans", "usethreads")
    except: pass
    nthreads = detectCPUs() * 2
    try: nthreads = config.getint("ctrans", "max_threads")
    except: pass
    
    if threads:
        ctrans['define_macros'].append(('USE_THREADS', None))
        ctrans['define_macros'].append(('NTHREADS', nthreads))
    
    sginfo = parseExtensionSetup('sginfo', config, ext_default)


ext_modules = []
if options['build_levmar']:
    ext_modules.append(Extension('pylevmar', ['src/pylevmar.c'],
                                 extra_compile_args = ['-g'],
                                 depends = ['src/pylevmar.h'],
                                 **levmar))
if options['build_ctrans']:
    ext_modules.append(Extension('ctrans', ['src/ctrans.c'],
                                 depends = ['src/ctrans.h'],
                                 **ctrans))

if options['build_sginfo']:
    ext_modules.append(Extension('pyspec.calcs.sginfo', 
                                 ['src/sginfomodule.c',
                                  'src/sgclib.c',
                                  'src/sgio.c',
                                  'src/sgfind.c',
                                  'src/sghkl.c',
                                  'src/sgsi.c'],
                                 **sginfo))
