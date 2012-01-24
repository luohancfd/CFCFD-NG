// full-cyl-schematic.asy
// Asymptote description of the cylinder and symmetry plane.

import three;
import solids;
import bsp; // need this for hidden-line removal

size(350);
unitsize(1mm);
defaultpen(1.0);
currentprojection = perspective(-30,30,40, up=Y);

// parameters defining the cylinder, from the Python description
real D = 15.0;     // Diameter of cylinder, millimetres
real L2 = D;       // (axial) half-length of cylinder
real Rc = D/2.0;   // cylinder radius
    
// Define a few key nodes.
triple a = (-Rc,0,0);  // stagnation point on the cylinder
triple b = (0,Rc,0);   // top of cylinder
triple c = (0,0,0);    // centre of curvature

// Coordinates of shock from Billig's correlation,
// as computed by the job input script.
triple[] d = new triple[6];
d[0] = 0.8*(-10.5685325882,0,0);
d[1] = 0.8*(-10.1060297773,0.5*Rc,0);
d[2] = 0.8*(-8.71958558487,1.0*Rc,0);
d[3] = 0.8*(-6.4123805245,1.5*Rc,0);
d[4] = 0.8*(-3.18967508786,2.0*Rc,0);
d[5] = 0.8*(+0.941249679354,2.5*Rc,0);

// Draw first so that we can draw over the top...
guide3 symmetry_plane = (-D,-D,0)--(-D,D,0)--(D,D,0)--(D,-D,0)--cycle;
draw(symmetry_plane,grey+2.0);

// A set of labelled axes, so that we know which way is up.
draw((0,0,0)--(17,0,0), red, Arrow3); label(XY()*"x", (17,0,0), N);
draw((0,0,0)--(0,17,0), red, Arrow3); label(XY()*"y", (0,17,0), E);
draw((0,0,0)--(0,0,21), red, Arrow3); label(ZY()*"z", (0,0,21), N);

// Cylinder as a solid in two pieces so that one can be drawn dashed.
revolution cyl0 = cylinder((0,0,-L2),Rc,L2,Z);
revolution cyl1 = cylinder((0,0,0),Rc,L2,Z);
draw(cyl0,black+dashed+1.0);
draw(cyl1,black+1.0);
label(XY()*"symmetry plane", (D,D,0), SW);

// Scribe the expected shock onto the symmetry plane.
guide3 shock_on_plane = (xpart(d[5]),-ypart(d[5]),zpart(d[5]));
for (int i=4; i > 0; --i) 
  shock_on_plane = shock_on_plane .. (xpart(d[i]),-ypart(d[i]),zpart(d[i]));
for (int i=0; i < 6; ++i) 
  shock_on_plane = shock_on_plane .. d[i];
draw(shock_on_plane,red);

// free-stream indicator
label(XY()*"$U_\infty$", (-13,0,0), W);
draw((-13,0,0)--(-10,0,0), blue, Arrow3);

