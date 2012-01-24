// block3d-defn.asy
// Asymptote description of the generic 3D block.

import three;

size(0,250);
currentprojection = perspective(-3,-4,2);
guide3 bottom = (0,0,0)--(1,0,0)--(1,1,0)--(0,1,0)--cycle3;
guide3 top = (0,0,1)--(1,0,1)--(1,1,1)--(0,1,1)--cycle3;
guide3 west = (0,0,0)--(0,0,1)--(0,1,1)--(0,1,0)--cycle3;
guide3 east = (1,0,0)--(1,0,1)--(1,1,1)--(1,1,0)--cycle3;
guide3 south = (0,0,0)--(1,0,0)--(1,0,1)--(0,0,1)--cycle3;
guide3 north = (0,1,0)--(1,1,0)--(1,1,1)--(0,1,1)--cycle3;

filldraw(bottom, lightgrey);
filldraw(east, lightgrey);
filldraw(north, lightgrey);
draw(top);
draw(west);
draw(south);
label("BOTTOM", (0.5,0.5,0));
label("TOP", (0.5,0.5,1));
// label("WEST", (0,0.5,0.5));
label("EAST", (1,0.5,0.5));
// label("SOUTH", (0.5,0,0.5));
label("NORTH", (0.5,1,0.5));

draw((1.05,0,0)--(1.2,0,0), EndArrow);
label("$i$", (1.2,0,0), S);
draw((0,1.05,0)--(0,1.2,0), EndArrow);
label("$j$", (0,1.2,0), S);
draw((0,0,1.1)--(0,0,1.25), EndArrow);
label("$k$", (0,0,1.25), E);

dot(unitcube);
label("0", (0,0,0), S);
label("1", (1,0,0), S);
label("2", (1,1,0), S);
label("3", (0,1,0), S);
label("4", (0,0,1), N);
label("5", (1,0,1), N);
label("6", (1,1,1), N);
label("7", (0,1,1), N);
