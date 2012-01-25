The CFCFD Code Collection
=========================

Current Codes
-------------
For local and external collaborators, these are available via a Mercurial repository.

#. *eilmer3*: The current compresible-flow Navier-Stokes simulation code
   for 2D and 3D flows. 
   This has the best bits of mbcns2 and elmer2, together with a rebuilt
   thermochemistry module.
   We have changed the name to Eilmer, the correct spelling of the name of the
   patron monk of CFD and aeronautics.
   See the wikipedia article http://en.wikipedia.org/wiki/Eilmer_of_Malmesbury
#. *L1d3*: Lagrangian simulation of quasi-one-dimensional flow.
#. *imoc*: interactive method-of-characteristics code, written in C+Tcl.

.. toctree::
   :maxdepth: 2

   getting-started
   eilmer3
   l1d3

Maintained Codes
----------------
Not actively developed but available to local collaborators through the Bazaar-ng repository.
 
#. *bos*: Background-oriented schlieren (current but needs maintenance; new numpy)
#. *elmer*: The first version of our time-dependent 3D flow simulation code.
   C version built for LES simulations. (obsolete)
#. *elmer2*: Second version using libgas2 thermochemsitry and written in C++. (obsolete)
#. *L1d*: Lagrangian simulation of quasi-one-dimensional flow. (obsolete)
#. *L1d2*: second version of this code build with libgas2 thermochemistry.
   Written in C++. (current)
#. *mb_cns*: Multiple-block compressible Navier-Stokes flow simulator.
   This is the code corresponding to the 1996 report. (obsolete)
#. *mbcns2*: second version using libgas2 and written in C++. (obsolete)
#. *newt*: Newtonian flow approximation for hypersonics flows.
#. *octvce*: Joseph's fast Euler code with hierarchical adaptive grid for 3D flows.
#. *poshax*: one-dimensional flow with complex thermochemistry.
#. *sdas*: firmware for simple data acquisition system for student Zuni projects.
#. *sf3d_v4.1*: Ian Johnston's shock-fitting code for 3D flow.
#. *sm_3d_plus*: space-marching code for 3D flow.  
   (obsolete now that block-sequencing is part of eilmer3)
#. *tds*: tunnel data server, software for both the server and the client


