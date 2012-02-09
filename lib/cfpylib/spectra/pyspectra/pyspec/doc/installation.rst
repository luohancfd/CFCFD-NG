======================
Installation of pyspec
======================

*****************************
Installing using easy_install
*****************************

``pyspec`` is hosted on PyPi_ under the name `pyspec`. This version is never "bleading edge" but is provided to allow easy installation. To install pyspec using ``easy_install``, from the command line run::

   easy_install pyspec

To force an upgrade run (if ``pyspec`` is already installed)::

   easy_install -U pyspec

Note: you can install a user copy by running::

   easy_install -U --user pyspec

just remember to set your `PYTHONPATH` environmental variable!

.. _PyPi: http://pypi.python.org/pypi/pyspec


**************************
Installing from SVN server
**************************

To obtain the latest bleading edge version of pyspec, it is best to check out a working copy using the subversion control system. For an introduction to SVN please see the SVN book here_.

.. _here: http://svnbook.red-bean.com/en/1.1/ch01.html

The SVN source is hosted on sourceforge and can be obtained by checking out a working copy, and building in the standard python way::

   svn co https://pyspec.svn.sourceforge.net/svnroot/pyspec/trunk pyspec
   cd pyspec
   python setup.py build
   sudo python setup.py install

which is of course for unix type systems. Building of extension modules is handled by the setup file `setup.cfg.<patform>` where `<platform>` is the type of computer you are building for. For example the `setup.cfg.darwin` (for Mac OS X systems) is shown below.

.. literalinclude:: ../setup.cfg.darwin

For each subheading the `build =` directive is used to switch on and off the building of each extension module. By default on windows systems, none of the extension modules are built to enable easy installation of the ``pyspec`` package.
