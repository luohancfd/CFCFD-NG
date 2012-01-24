<?php include("layout-top.php"); ?>

<p align="center" class="headline">
<!-- Insert page header here -->
Code Collection for the Hypersonics Group.
</p>

<!-- Insert text underneath -->
<P>
This archive contains source code, example files and documentation for 
a number of compressible-flow simulation codes.
It is intended for use by members of the UQ Hypersonics Group and the
Australian Centre for Hypersonics.
If you are browsing this archive via our web server, not all of the 
directories will be complete.
</P>


<A NAME="research_codes"><H3>Research Codes</H3></A>
<P>
<OL>
   <LI> Compressible Flow Simulation - One-dimensional
      <UL>
         <LI> <A HREF="./code/L1d/doc/">L1d</A> 
              Lagrangian 1-D flow simulation program. </LI>
      </UL></LI>
   <LI> Compressible Flow Simulation - Two-dimensional
      <UL>
         <LI> <A HREF="./code/imoc/doc/">IMOC</A> 
              Method-of-Characteristics for two-dimensional supersonic flow. </LI>
         <LI> <A HREF="./code/mb_cns/doc/">MB_CNS</A>
              Multi-block Navier-Stokes solver for 2D/axisymmetric geometries.
	      There is an MPI parallel flavour, also. </LI>
         <LI> <A HREF="./code/pamela_who/">PAMELA</A>
              Andrew's parallel Navier-Stokes program. </LI>
         <LI> <A HREF="./code/u2de/">U2DE</A>
              Paul's unstructured grid code for 2-D flows. </LI>
      </UL></LI>
   <LI> Compressible Flow Simulation - Three-dimensional
      <UL>
         <LI> <A HREF="./code/newt/">Newt</A>
              Kevin's Newtonian flow program. </LI>
         <LI> <A HREF="./code/elmer/">Elmer</A>
              Navier-Stokes solver for 3D geometries. </LI>
         <LI> <A HREF="./code/sm_3d_plus/">SM_3D+</A>
              Chris' space-marching program.
              This supercedes SM_3D listed below. </LI>
         <LI> <A HREF="./code/sf3d_v4.1/">sf3d_v4.1</A>
              Ian's 3-D shock-fitting code. </LI>
      </UL></LI>
   <LI> Compressible Flow Simulation - Utility
      <UL>
         <LI> <A HREF="./code/tds/">TDS</A>
              Shock Tunnel Data Server and Browser</LI>
         <LI> <A HREF="./code/gas_models/source/estcj.py">ESTCj</A> 
	      A Python program for the calculation of shock-tube conditions
	      assuming chemical equilibrium. 
	      (Needs the NASA Lewis CEA software as an executable code.)</LI>
      </UL></LI>
   <LI> Common components for the Research codes.
      <UL>
         <LI> <A HREF="./code/gas_models/">Code for thermodynamic models</A> </LI>
         <LI> <A HREF="./code/plot/">Graphics and Plotting Codes</A> </LI>
         <LI> <A HREF="./code/util/">Utility code</A> </LI>
         <LI> <A HREF="./code/nm/">Numerical Methods</A> </LI>
      </UL>
   <LI> Tar-balls for some of the code.
      <UL> 
         <LI> <A HREF="./download/mb_cns.tar.gz">mb_cns.tar.gz</A> </LI>
         <LI> <A HREF="./download/cfd_util.tar.gz">cfd_util.tar.gz</A> 
	 (Includes the gas-model code that is needed by L1d and mb_cns.)</LI>
         <LI> <A HREF="./download/cfd_plot.tar.gz">cfd_plot.tar.gz</A> </LI>
         <LI> <A HREF="./download/L1d.tar.gz">L1d.tar.gz</A> </LI>
         <LI> <A HREF="./download/imoc.tar.gz">imoc.tar.gz</A> </LI>
         <LI> <A HREF="./download/imoc.zip">imoc.zip</A> </LI>
         <LI> <A HREF="./download/newt.tar.gz">newt.tar.gz</A> </LI>
         <LI> <A HREF="./download2/tds.tar.gz">tds.tar.gz</A>
	 (Current version; Works on both Linux/UNIX and MS-Windows.) </LI>
         <LI> <A HREF="./download2/tds_w32.exe">tds_w32.exe</A> 
	 (Old version; Self-extracting archive for Windows.)</LI>
      </UL></LI>
</OL>
</P>


<A NAME="archived_codes"><H3>Other Archived Codes</H3></A>
<P>
   These are in various states of repair and are not being updated.
   However, they continue to be useful for the odd shock tunnel job.
   <UL>
      <LI> <A HREF="./code/others/monc_v4.8/">MONC Version 4.8</A>
           Allan Paull's data acquisition and display program. </LI>
      <LI> <A HREF="./code/sm_3d_95/">SM_3D</A>
           Space Marching Code as per Report 2/95. </LI>
      <LI> <A HREF="./code/others/moc_87/">MOC_87</A>
           Interactive Turbo-Pascal-3 program for 
           the Method of Characteristics (1987) </LI>
      <LI> <A HREF="./code/others/moc_eang_95/">MOC_95</A> 
           Peter Eang's (incomplete) work on MOC (1995).
           I think that we'll finish this job with one of
           the newer scripting languages such as Python or Tcl/Tk. 
           See IMOC above. </LI>
      <LI> <A HREF="./code/others/maccoll_89/">ConeFlow</A>
           Simon Sanderson's cone flow code.
           Includes FORTRAN source code with DOS executable from 1989. </LI>
      <LI> <A HREF="./code/others/estc_90/">ESTC</A>
           Malcolm McIntosh's Equilibrium Shock Tube Conditions code.
           Includes FORTRAN source code with DOS executable 
           from some date between 1968 and 1990. 
	   (Use ESTCj instead if you have a Python interpreter 
	   and the CEA software.)</LI>
      <LI> <A HREF="./code/others/nenzf_88/">NENZF</A>
           Nonequilibrium nozzle flow code from Lordi et.al. 
           Includes FORTRAN source code with DOS executable from 1988. 
           This is another old code which has been updated to run on PC's. </LI>
      <LI> <A HREF="./code/others/shock_88/">SHOCK</A>
           A replacement for ESTC at low temperatures.
           Originally written by Richard Morgan and Rob Casey in BASIC.
           Updated for fortran on an IBM 3081 by Craig Brescianini, 1988. </LI>
      <LI> <A HREF="./code/others/stn_95/">STN</A>
           A replacement for ESTC and NENZF for equilibrium air and nitrogen.
           Written in FORTRAN-77 by Rob Krek and PJ, 1993.
           Updated 1995. </LI>
      <LI> <A HREF="./code/others/crek_92/">EQSTATE/CREK</A>
           Generalised one dimensional flow with finite-rate chemical reactions.
           FORTRAN code that used to be run on an IBM 3081.
           There seems to be a DOS version also.  
           Bob Bakos might know some more about this program. </LI>
      <LI> <A HREF="./code/others/sharc_91/">SHARC</A>
           Craig Brescianini's finite difference code; developed for his PhD 
           work on scramjet simulations. </LI>
      <LI> <A HREF="./code/others/nasa_moc/">NASA MOC</A>
           This neat Method of Characteristics program 
           was written by J. M. Forster.
           Updated for DOS compilers by PJ, May 1997.
           See the beginning of the source code for more information. </LI>
      <LI> <A HREF="./code/others/tube_gary_91">Gary's Tube code</A> 
           A program (in C) for shock tube simulations; derived from MNM's
           FORTRAN program. </LI>
   </UL>
</P>

<p><small>
   Last Updated by PJ 30-Aug-2005 <BR>
</small></p>
<!-- End text -->

<?php include("layout-bottom.php"); ?>
