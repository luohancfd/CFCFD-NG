imoc
====

Interactive method-of-characteristics code for irrotational flow of an ideal-gas.
It sounds limited but MOC has the nice design capability of being able to specify 
what upstream flow that you have, what downstream flow you would like, 
and then be able to compute the flow field in between.

Report-3-00_ describes the user interface to the Interactive Method-of-Characteristics program.
Note that, although the current web copy is a bit broken, the code and
documentation is in the code repository.

.. _Report-3-00: ./imoc/index.html

Basic install instructions
--------------------------
These instructions are for a linux system that is already running other cfcfd3 codes (ie. Eilmer3).

Ensure that you have tcl-dev installed, then::

 $ cd cfcfd3/app/imoc/unix
 $ make
 $ make install

This should create a folder in your home directory called imoc_bin. IMOC can then be run from within
this folder using the command::

 $ cd ~/imoc_bin
 $ ./moc_startup.tcl
