Semi-retired (Old) Codes
========================

Not actively developed but available to local collaborators 
through the Bazaar-ng repository.
 
#. *bos*: Background-oriented schlieren (current but needs maintenance; new numpy)
#. *elmer*: The first version of our time-dependent 3D flow simulation code.
   C version built for LES simulations. (obsolete)
#. *elmer2*: Second version using libgas2 thermochemsitry and written in C++.
   (obsolete)
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


Getting the old codes
---------------------
The older cfcfd2 code collection is still available using bazaar-ng.
To get (read-only) access via http, use the command:

| cd $HOME
| bzr branch http://cfcfdlocal@triton.pselab.uq.edu.au/bzr-cfcfd2 cfcfd2

Notes

#. You will need a password for bzr access.  Please ask.
#. To push code changes, you will need access via ssh.
#. Even if you have ssh access, we'd like to keep this repository
   static at revision 1212.  All new work should be going into the 
   cfcfd3 repository.


Other Archived Codes
--------------------
These are in various states of repair and are not being updated. 
However, they continue to be useful for the odd shock tunnel job within the group.
They can be found in the alt_app directory of the old repository.

#. *MONC Version 4.8* Allan Paull's data acquisition and display program.
#. *SM_3D* Space Marching Code as per Report 2/95.
#. *MOC_87* Interactive Turbo-Pascal-3 program for the Method of Characteristics (1987)
#. *MOC_95* Peter Eang's (incomplete) work on MOC (1995). 
   I think that we'll finish this job with one of the newer scripting languages 
   such as Python or Tcl/Tk. See IMOC above.
#. *ConeFlow* Simon Sanderson's cone flow code. 
   Includes FORTRAN source code with DOS executable from 1989.
#. *ESTC* Malcolm McIntosh's Equilibrium Shock Tube Conditions code. 
   Includes FORTRAN source code with DOS executable from some date between 1968 and 1990. 
   (Use ESTCj instead if you have a Python interpreter and the CEA software.)
#. *NENZF* Nonequilibrium nozzle flow code from Lordi et.al. 
   Includes FORTRAN source code with DOS executable from 1988. 
   This is another old code which has been updated to run on PC's.
#. *SHOCK* A replacement for ESTC at low temperatures. 
   Originally written by Richard Morgan and Rob Casey in BASIC. 
   Updated for fortran on an IBM 3081 by Craig Brescianini, 1988.
#. *STN* A replacement for ESTC and NENZF for equilibrium air and nitrogen. 
   Written in FORTRAN-77 by Rob Krek and PJ, 1993. Updated 1995.
#. *EQSTATE/CREK* Generalised one dimensional flow with finite-rate chemical reactions.
   FORTRAN code that used to be run on an IBM 3081. There seems to be a DOS version also. 
   Bob Bakos might know some more about this program.
#. *SHARC* Craig Brescianini's finite difference code; 
   developed for his PhD work on scramjet simulations.
#. *NASA MOC* This neat Method of Characteristics program was written by J. M. Forster.
   Updated for DOS compilers by PJ, May 1997. 
   See the beginning of the source code for more information.
#. *Gary Allen's Tube code* A program (in C) for shock tube simulations; 
   derived from Michael Macrossan's FORTRAN program.


