Getting the codes and preparing to run them
===========================================

The code repository
-------------------
The codes are available for download from a Mercurial repository.
To make a clone of the repository::

  $ cd $HOME
  $ hg clone https://cfcfdlocal@triton.pselab.uq.edu.au/cfcfd3-hg/cfcfd3-hg/ cfcfd3

This takes about 40 seconds on campus at UQ.  
It may take much longer, depending on your internet connection.

To see what's changed::

  $ cd cfcfd3
  $ hg incoming https://cfcfdlocal@triton.pselab.uq.edu.au/cfcfd3-hg/cfcfd3-hg/
  ...
  $ hg pull -u https://cfcfdlocal@triton.pselab.uq.edu.au/cfcfd3-hg/cfcfd3-hg/

Notes

#. You will need a password for any access.  Please ask.
#. You can read but not write with the "cfcfdlocal" username.
#. Some usernames (by negotiation) may push changesets back to the repository.


Licence
-------
CFCFD program collection is a set of flow simulation tools for compressible fluids.
Copyright (C) 1991-2012 Peter Jacobs, Rowan Gollan, Daniel Potter,
Brendan O'Flaherty, Fabian Zander, Wilson Chan, Peter Blyton and
other members of the CFCFD group.

This collection is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or any later version.

This program collection is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
See the GNU General Public License for more details.

You should have received a copy of the GNU-General-Public-License_
along with this program.  If not, see <http://www.gnu.org/licenses/>.

.. _GNU-General-Public-License: ./_static/gpl.txt


Your computational environment
------------------------------
The code collection comes as source code only so,
to use any of them, you will need to compile and install them.

To build and run the newer codes, you will need the following:

* a Unix-like system with GNU-make, C and C++ compilers
* popt (command-line parser) library and development files
* readline library (including the header files, libreadline-dev on Ubuntu)
* Python + the numpy, matplotlib and scipy extensions
* SWIG
* Tcl/Tk + the BWidget library (to run the GUI program e3console.tcl)

We have been able to get the programs to build on Linux, MacOS-X 
(with a recent Xcode development environment) and Cygwin 1.7 (on MS-Windows).

On MS-Windows, install the full kit of Cygwin (Python, X-Windows and all)
and be careful not to have another Python installed outside of Cygwin.
The multiple installations of Python seem not to play well together.

Some other things that are useful:

* awk
* MetaPost (mpost) or, more recently, InkScape (for looking at and editing svg files)
* GNUplot
* Paraview or MayaVi or VisIt

To a basic Fedora 16 installation, you should add the following packages:

#. mercurial
#. gcc
#. gcc-c++
#. m4
#. openmpi
#. openmpi-devel (to use openmpi on Fedora, :ref:`the module must be loaded <label-openmpi-fedora>`)
#. gcc-gfortran
#. libgfortran.i686 and glibc-devel.i686 (to compile the 32-bit CEA code on 64-bit Fedora)
#. swig
#. python-devel
#. numpy
#. python-matplotlib
#. scipy
#. readline-devel (for Lua)
#. popt-devel

To a basic Ubuntu 10.04 (or any recent Debian derivative) installation, 
you should add the following packages and their dependencies:

#. mercurial
#. g++
#. m4
#. mpi-default-dev
#. mpi-default-bin
#. gfortran
#. gfortran-multilib (for compiling 32-bit CEA2 on a 64-bit system)
#. swig
#. python-dev
#. python-numpy
#. python-matplotlib
#. python-scipy
#. libreadline-dev
#. libpopt-dev
#. tk
#. bwidget
#. gnuplot
#. tcl-dev (if you want to build IMOC)
#. maxima (to run the Method-of-Manufactured-Solutions test case for Eilmer3)

Using the codes on MS-Windows
-----------------------------
The codes assemble most conveniently on a Linux/Unix-like environment.
They should also build and run within Cygwin (http://cygwin.com/), however,
it may be convenient to run a full linux installation within 
VirtualBox (https://www.virtualbox.org/), on your MS-Windows computer.

Using the codes on Apple OSX
----------------------------
The codes can be compiled and run on OSX as this is a Unix based OS.
The Xcode development environment (https://developer.apple.com/xcode/) 
should be downloaded and installed to provide Apple's versions of the 
GNU Compiler Collection, Python and the make utility, amongst other
development tools.
popt, readline, SWIG and Tcl/Tk can either be installed from source
or via a package manager such as MacPorts (http://www.macports.org/) or 
Fink (http://www.finkproject.org/).

Notes:

#. If possible, it is recommended to install these dependencies from source.
#. The required Python packages (numpy, scipy and matplotlib) are all available
   as pre-packaged binaries for OSX on sourceforge.net, although they can also
   be installed from source if necessary.
#. Ingo has had a good experience installing binary packages from MacPorts,
   the only subtly being the need to install swig and swig-python.

SSH access to the repository for developers
-------------------------------------------
Alternative access to the Mercurial repository for developers is possible via ssh.
You will need the password or your public key installed for any access.  Please ask.

::

  $ cd ~
  $ hg clone ssh://cfcfd3@triton/cfcfd3-hg cfcfd3
  $ cd cfcfd3/extern/
  $ hg clone ssh://cfcfd3@triton/cea2-hg cea2
  $ hg clone ssh://geothermal@triton/refprop-hg refprop



Notes about Mercurial and https certificate warnings
----------------------------------------------------
For versions of Mercurial greater than 1.7.3, a warning will be issued
about the certificate not being verified when accessing the repository
over https. To satisy Mercurial's complaints, you will need to configure
the Certificate Authorities (CAs) which it uses. There are two ways to
do this:

1. configure HTTPS certificate authorities; or
2. verify ``triton.pselab.uq.edu.au`` individually using its fingerprint.

In either case, you will need to edit your hg configuration file which
can be a repository-specific file ``.hg/hgrc`` or set globally in
``~/.hgrc``.

To configure the certificate authorities, the value for ``cacerts`` need to
be set correctly for your system. For example, a Fedora (or Fedora-like) linux system,
this can be done by adding the following to the ``hgrc`` file::

  [web]
  cacerts = /etc/ssl/certs/ca-bundle.crt

Examples for other linux systems can be found at MercurialCAs_.

The alternative is to configure the host fingerprint for
``triton.pselab.uq.edu.au`` explicitly. To do this, add
the following to your hg config file::

  [hostfingerprints]
  triton.pselab.uq.edu.au = 1d:33:32:b0:6c:e2:5c:13:67:35:ba:e6:60:cc:4e:c1:03:63:5a:2e

More information about configuring Mercurial to use your system's certificate
authorities is available at MercurialCAs_.


.. _MercurialCAs: http://mercurial.selenic.com/wiki/CACertificates
