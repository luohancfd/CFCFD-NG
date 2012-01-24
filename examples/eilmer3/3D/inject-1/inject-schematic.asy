// inject-schematic.asy
// Asymptote description of the duct and injector.

import three;

size(0,300);
currentprojection = perspective(-20,-10,3.5);

// parameters (in cm)  defining the duct, as per the ptyhon description
real L0 = 20.0;
real L1 = 5.0;
real L2 = 1.0;
real Whalf0 = 5.0;
real Whalf1 = 0.5;
real H = 5.0;

draw((0,0,0)--(L0,0,0)--(L0,Whalf0,0)--(0,Whalf0,0)--cycle3); // Bottom
draw((0,0,H)--(L0,0,H)--(L0,Whalf0,H)--(0,Whalf0,H)--cycle3); // Top

draw((0,0,0)--(0,0,H)--(0,Whalf0,H)--(0,Whalf0,0)--cycle3); // West
label("inflow", (0,Whalf0/2,0), S);

draw((L0,0,0)--(L0,0,H)--(L0,Whalf0,H)--(L0,Whalf0,0)--cycle3); // East
label("outflow", (L0,Whalf0/2,0), N);

draw((0,0,0)--(L0,0,0)--(L0,0,H)--(0,0,H)--cycle3); // South
draw((0,Whalf0,0)--(L0,Whalf0,0)--(L0,Whalf0,H)--(0,Whalf0,H)--cycle3); // North

guide3 port = (L1,0,0)--(L1+L2,0,0)--(L1+L2,Whalf1,0)--(L1,Whalf1,0)--cycle3;
filldraw(port,lightgrey);
label("port", (L1+L2/2,0,0), SE);

draw((L0+0.5,0,0)--(L0+2.5,0,0), EndArrow);
label("x", (L0+2.5,0,0), SE);
draw((0,0,H+0.1)--(0,0,H+0.6), EndArrow);
label("z", (0,0,H+0.6), E);
draw((0,Whalf0+0.1,0)--(0,Whalf0+0.6,0), EndArrow);
label("y", (0,Whalf0+0.6,0), S);
