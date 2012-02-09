UQ Python Spectral Functions
Fabian Zander
February 9th 2012.

This file describes the basic usage of the spectral analysis files as modified and developed by FZ.
This includes the installation of the pyspec package (© Copyright 2010, Stuart B. Wilkins; Sven Partzsch,
http://packages.python.org/pyspec/index.html), the modifications by FZ and the basic additional
functions written by FZ.

This file is current as of February 9th 2012 and is a work in progress, report problems or ask questions
to f.zander@uq.edu.au


These instructions assume that you are using a linux system and have a working understanding of linux systems.

1. Install the pyspec package from this directory.
	cd pyspec
	python setup.py build
	sudo python setup.py install

2. In order to read the winspec files used by the UQ spectrometer systems the pyspec installed files needed 
to be modified to read the required data. Copy the modified 'files.py' file into the correct directory.
The example shown below is the location for FZ's system.
	cd ..
	sudo cp files.py /usr/local/lib/python2.7/dist-packages/pyspec-0.2_r223-py2.7-linux-x86_64.egg/pyspec/ccd/.

3. Copy the spectral_functions.py file into your pythonpath if you wish to use FZ's functions that have 
been developed for specific applications. It is recommended you change these as required for your work. It is
anticipated that these will be developed into something more significant in the future. The plotting parts of
these functions require matplotlib.

4. Run the example for a demonstration of the output
	python spectra_example.py
