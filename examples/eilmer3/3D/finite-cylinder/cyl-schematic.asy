// cyl-schematic.asy
// Asymptote description of the cylinder and some of the gas domain boundaries.

import three;

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
d[0] = 0.87*(-10.5685325882,0,0);
d[1] = 0.87*(-10.1060297773,0.5*Rc,0);
d[2] = 0.87*(-8.71958558487,1.0*Rc,0);
d[3] = 0.87*(-6.4123805245,1.5*Rc,0);
d[4] = 0.87*(-3.18967508786,2.0*Rc,0);
d[5] = 0.87*(+0.941249679354,2.5*Rc,0);

triple e = shift(Z*L2)*d[0];
triple f = (-Rc,0,L2);
triple g = (-Rc/2,0,L2);
triple c2 = (0,0,L2);
triple h = (0,Rc/2,L2);
triple i = (0,Rc,L2);
triple j = shift(Rc*Z)*e;
triple k = shift(Rc*Z)*f;

// A set of labelled axes, so that we know which way is up.
draw((0,0,0)--(15,0,0), red, Arrow3); label(XY()*"x", (15,0,0), N);
draw((0,0,0)--(0,20,0), red, Arrow3); label(XY()*"y", (0,20,0), E);
draw((0,0,0)--(0,0,25), red, Arrow3); label(ZY()*"z", (0,0,25), N);

// Paths on symmetry plane...
guide3 inlet = d[0];
for (int i=0; i < 6; ++i) {
  inlet = inlet .. d[i];
}
guide3 cyl_edge = arc(c, b, a, normal=Z, direction=CCW);
guide3 symmetry_plane = d[0]--inlet--cyl_edge--cycle;
draw(symmetry_plane,dashed+1.0);
// Cylinder surfaces.
surface cyl_surf = extrude(cyl_edge, Z*L2);
draw(cyl_surf,nu=15,nv=2,surfacepen=nullpen,meshpen=currentpen+1.0);
guide3 arc_gh = arc(c2,g,h,normal=Z,direction=CW);
guide3 end_face = g--arc_gh--arc(c2,i,f,normal=Z,direction=CCW)--cycle;
draw(surface(end_face),nu=15,nv=1,surfacepen=nullpen,meshpen=currentpen+1.0);
// Draw some outlines of block edges so we can see the extent of the domain.
draw(d[0]--e--j--k--shift(Rc*Z)*g--g,dashed);
draw(shift(Rc*Z)*arc_gh,dashed);
draw(h--shift(Rc*Z)*h,dashed);
draw(i--shift(Rc*Z)*i,dashed);
draw(d[5]--shift(L2*Z)*d[5]--shift((L2+Rc)*Z)*d[5]--shift(Rc*Z)*i,dashed);
draw(shift(Rc*Z)*h--shift(Rc*Z)*i,dashed);

// Some dotted points with labels to match input script Nodes.
dot(inlet, red);
label("\small d[0]", d[0], W);
label("\small d[5]", d[5], NE);
dot(a,red); label("\small a", a, NW);
dot(b,red); label("\small b", b, NE);
dot(c,red); label("\small c", c, NW);
dot(e,red); label("\small e", e, SW);
dot(f,red); label("\small f", f, SE);
dot(g,red); label("\small g", g, NE);
dot(h,red); label("\small h", h, NE);
dot(i,red); label("\small i", i, NE);
dot(j,red); label("\small j", j, S);
dot(k,red); label("\small k", k, S);

// free-stream indicator
label(XY()*"$U_\infty$", (-13,0,L2/2), W);
draw((-13,0,L2/2)--(-10,0,L2/2), blue, Arrow3);

