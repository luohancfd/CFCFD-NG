#!/usr/bin/env python
# euler_source.py
from math import sin, cos, pow, pi
from blockgrid2d import *

L = 1.0;
gam = 1.4;

rho0 = 1.0;
rhox = 0.15;
rhoy = -0.1;

uvel0 = 800.0;
uvelx = 50.0;
uvely = -30.0;

vvel0 = 800.0;
vvelx = -75.0;
vvely = 40.0;

wvel0 = 0.0;

press0 = 1.0e5;
pressx = 0.2e5;
pressy = 0.5e5;


def rho_source(x, y):

    f_m = (3*pi*uvelx*cos((3*pi*x)/(2.*L))*(rho0 + rhoy*cos((pi*y)/(2.*L)) + rhox*sin((pi*x)/L)))/(2.*L) + \
   (2*pi*vvely*cos((2*pi*y)/(3.*L))*(rho0 + rhoy*cos((pi*y)/(2.*L)) + rhox*sin((pi*x)/L)))/(3.*L) + \
   (pi*rhox*cos((pi*x)/L)*(uvel0 + uvely*cos((3*pi*y)/(5.*L)) + uvelx*sin((3*pi*x)/(2.*L))))/L - \
	(pi*rhoy*sin((pi*y)/(2.*L))*(vvel0 + vvelx*cos((pi*x)/(2.*L)) + vvely*sin((2*pi*y)/(3.*L))))/(2.*L)

    return f_m


def xmom_source(x, y):
    f_x = (3*pi*uvelx*cos((3*pi*x)/(2.*L))*(rho0 + rhoy*cos((pi*y)/(2.*L)) + rhox*sin((pi*x)/L))* \
           (uvel0 + uvely*cos((3*pi*y)/(5.*L)) + uvelx*sin((3*pi*x)/(2.*L))))/L + \
           (2*pi*vvely*cos((2*pi*y)/(3.*L))* \
            (rho0 + rhoy*cos((pi*y)/(2.*L)) + rhox*sin((pi*x)/L))* \
            (uvel0 + uvely*cos((3*pi*y)/(5.*L)) + uvelx*sin((3*pi*x)/(2.*L))))/(3.*L) + \
            (pi*rhox*cos((pi*x)/L)*pow(uvel0 + uvely*cos((3*pi*y)/(5.*L)) + \
                                         uvelx*sin((3*pi*x)/(2.*L)),2))/L - (2*pi*pressx*sin((2*pi*x)/L))/L - \
                                         (pi*rhoy*(uvel0 + uvely*cos((3*pi*y)/(5.*L)) + uvelx*sin((3*pi*x)/(2.*L)))* \
                                          sin((pi*y)/(2.*L))*(vvel0 + vvelx*cos((pi*x)/(2.*L)) + vvely*sin((2*pi*y)/(3.*L)))) \
                                          /(2.*L) - (3*pi*uvely*(rho0 + rhoy*cos((pi*y)/(2.*L)) + rhox*sin((pi*x)/L))* \
                                                     sin((3*pi*y)/(5.*L))*(vvel0 + vvelx*cos((pi*x)/(2.*L)) + \
                                                                           vvely*sin((2*pi*y)/(3.*L))))/(5.*L)

    return f_x


def ymom_source(x, y):
    f_y = (pi*pressy*cos((pi*y)/L))/L - (pi*vvelx*sin((pi*x)/(2.*L))* \
                                         (rho0 + rhoy*cos((pi*y)/(2.*L)) + rhox*sin((pi*x)/L))* \
                                         (uvel0 + uvely*cos((3*pi*y)/(5.*L)) + uvelx*sin((3*pi*x)/(2.*L))))/(2.*L) +  \
                                         (3*pi*uvelx*cos((3*pi*x)/(2.*L))* \
                                          (rho0 + rhoy*cos((pi*y)/(2.*L)) + rhox*sin((pi*x)/L))* \
                                          (vvel0 + vvelx*cos((pi*x)/(2.*L)) + vvely*sin((2*pi*y)/(3.*L))))/(2.*L) +  \
                                          (4*pi*vvely*cos((2*pi*y)/(3.*L))* \
                                           (rho0 + rhoy*cos((pi*y)/(2.*L)) + rhox*sin((pi*x)/L))* \
                                           (vvel0 + vvelx*cos((pi*x)/(2.*L)) + vvely*sin((2*pi*y)/(3.*L))))/(3.*L) +  \
                                           (pi*rhox*cos((pi*x)/L)*(uvel0 + uvely*cos((3*pi*y)/(5.*L)) +  \
                                                                   uvelx*sin((3*pi*x)/(2.*L)))* \
                                            (vvel0 + vvelx*cos((pi*x)/(2.*L)) + vvely*sin((2*pi*y)/(3.*L))))/L -  \
                                            (pi*rhoy*sin((pi*y)/(2.*L))*pow(vvel0 + vvelx*cos((pi*x)/(2.*L)) + \
                                                                              vvely*sin((2*pi*y)/(3.*L)),2))/(2.*L)

    return f_y
    


def energy_source(x, y):

    f_e = (uvel0 + uvely*cos((3*pi*y)/(5.*L)) + uvelx*sin((3*pi*x)/(2.*L)))* \
          ((-2*pi*pressx*sin((2*pi*x)/L))/L + (rho0 + rhoy*cos((pi*y)/(2.*L)) + rhox*sin((pi*x)/L))* \
           ((-2*pi*pressx*sin((2*pi*x)/L))/((-1 + gam)*L*(rho0 + rhoy*cos((pi*y)/(2.*L)) + rhox*sin((pi*x)/L))) + \
            ((3*pi*uvelx*cos((3*pi*x)/(2.*L))*(uvel0 + uvely*cos((3*pi*y)/(5.*L)) + uvelx*sin((3*pi*x)/(2.*L))))/L - \
             (pi*vvelx*sin((pi*x)/(2.*L))*(vvel0 + vvelx*cos((pi*x)/(2.*L)) + vvely*sin((2*pi*y)/(3.*L))))/L)/2. - \
            (pi*rhox*cos((pi*x)/L)*(press0 + pressx*cos((2*pi*x)/L) + pressy*sin((pi*y)/L)))/ \
            ((-1 + gam)*L*pow(rho0 + rhoy*cos((pi*y)/(2.*L)) + rhox*sin((pi*x)/L),2))) + \
           (pi*rhox*cos((pi*x)/L)*((pow(wvel0,2) + pow(uvel0 + uvely*cos((3*pi*y)/(5.*L)) + uvelx*sin((3*pi*x)/(2.*L)),2) + \
                                        pow(vvel0 + vvelx*cos((pi*x)/(2.*L)) + vvely*sin((2*pi*y)/(3.*L)),2))/2. + \
                                       (press0 + pressx*cos((2*pi*x)/L) + pressy*sin((pi*y)/L))/((-1 + gam)*(rho0 + rhoy*cos((pi*y)/(2.*L)) + rhox*sin((pi*x)/L)))) \
            )/L) + (3*pi*uvelx*cos((3*pi*x)/(2.*L))*(press0 + pressx*cos((2*pi*x)/L) + pressy*sin((pi*y)/L) + \
                                                         (rho0 + rhoy*cos((pi*y)/(2.*L)) + rhox*sin((pi*x)/L))* \
                                                         ((pow(wvel0,2) + pow(uvel0 + uvely*cos((3*pi*y)/(5.*L)) + uvelx*sin((3*pi*x)/(2.*L)),2) + \
                                                           pow(vvel0 + vvelx*cos((pi*x)/(2.*L)) + vvely*sin((2*pi*y)/(3.*L)),2))/2. + \
                                                          (press0 + pressx*cos((2*pi*x)/L) + pressy*sin((pi*y)/L))/((-1 + gam)*(rho0 + rhoy*cos((pi*y)/(2.*L)) + rhox*sin((pi*x)/L)))) \
                                                         ))/(2.*L) + (2*pi*vvely*cos((2*pi*y)/(3.*L))*(press0 + pressx*cos((2*pi*x)/L) + pressy*sin((pi*y)/L) + \
                                                                                                           (rho0 + rhoy*cos((pi*y)/(2.*L)) + rhox*sin((pi*x)/L))* \
                                                                                                           ((pow(wvel0,2) + pow(uvel0 + uvely*cos((3*pi*y)/(5.*L)) + uvelx*sin((3*pi*x)/(2.*L)),2) + \
                                                                                                             pow(vvel0 + vvelx*cos((pi*x)/(2.*L)) + vvely*sin((2*pi*y)/(3.*L)),2))/2. + \
                                                                                                            (press0 + pressx*cos((2*pi*x)/L) + pressy*sin((pi*y)/L))/((-1 + gam)*(rho0 + rhoy*cos((pi*y)/(2.*L)) + rhox*sin((pi*x)/L)))) \
                                                                                                           ))/(3.*L) + (vvel0 + vvelx*cos((pi*x)/(2.*L)) + vvely*sin((2*pi*y)/(3.*L)))* \
                                                                                                           ((pi*pressy*cos((pi*y)/L))/L - (pi*rhoy*sin((pi*y)/(2.*L))* \
                                                                                                                                               ((pow(wvel0,2) + pow(uvel0 + uvely*cos((3*pi*y)/(5.*L)) + uvelx*sin((3*pi*x)/(2.*L)),2) + \
                                                                                                                                                 pow(vvel0 + vvelx*cos((pi*x)/(2.*L)) + vvely*sin((2*pi*y)/(3.*L)),2))/2. + \
                                                                                                                                                (press0 + pressx*cos((2*pi*x)/L) + pressy*sin((pi*y)/L))/((-1 + gam)*(rho0 + rhoy*cos((pi*y)/(2.*L)) + rhox*sin((pi*x)/L)))) \
                                                                                                                                               )/(2.*L) + (rho0 + rhoy*cos((pi*y)/(2.*L)) + rhox*sin((pi*x)/L))* \
                                                                                                            ((pi*pressy*cos((pi*y)/L))/((-1 + gam)*L*(rho0 + rhoy*cos((pi*y)/(2.*L)) + rhox*sin((pi*x)/L))) + \
                                                                                                             ((-6*pi*uvely*(uvel0 + uvely*cos((3*pi*y)/(5.*L)) + uvelx*sin((3*pi*x)/(2.*L)))*sin((3*pi*y)/(5.*L)))/(5.*L) + \
                                                                                                              (4*pi*vvely*cos((2*pi*y)/(3.*L))*(vvel0 + vvelx*cos((pi*x)/(2.*L)) + vvely*sin((2*pi*y)/(3.*L))))/(3.*L))/2. + \
                                                                                                             (pi*rhoy*sin((pi*y)/(2.*L))*(press0 + pressx*cos((2*pi*x)/L) + pressy*sin((pi*y)/L)))/ \
                                                                                                             (2.*(-1 + gam)*L*pow(rho0 + rhoy*cos((pi*y)/(2.*L)) + rhox*sin((pi*x)/L),2))))
    return f_e

def main():
    import sys

    grid_name = sys.argv[1]
    gfp = open(grid_name, "r")
    while(1):
        line = gfp.readline()
        if not line:
            break
        tks = line.split()
        ni = int(tks[0]) + 1
        nj = int(tks[1]) + 1
        grid = BlockGrid2D(ni, nj)
        grid.read_block_in_classic_mbcns_format(gfp)
    gfp.close()

    solution = []
    for i in range(grid.ni - 1):
        for j in range(grid.nj - 1):
            xSE = grid.x[i+1][j];   ySE = grid.y[i+1][j]
            xNE = grid.x[i+1][j+1]; yNE = grid.y[i+1][j+1]
            xNW = grid.x[i][j+1];   yNW = grid.y[i][j+1]
            xSW = grid.x[i][j];     ySW = grid.y[i][j]
            # Cell area in the (x,y)-plane.
            xyarea = 0.5 * ((xNE + xSE) * (yNE - ySE) + (xNW + xNE) * (yNW - yNE) +
                            (xSW + xNW) * (ySW - yNW) + (xSE + xSW) * (ySE - ySW))
            # Cell Centroid.
            centre_x = 1.0 / (xyarea * 6.0) * \
                       ((yNE - ySE) * (xSE * xSE + xSE * xNE + xNE * xNE) + 
                        (yNW - yNE) * (xNE * xNE + xNE * xNW + xNW * xNW) +
                        (ySW - yNW) * (xNW * xNW + xNW * xSW + xSW * xSW) + 
                        (ySE - ySW) * (xSW * xSW + xSW * xSE + xSE * xSE))
            centre_y = -1.0 / (xyarea * 6.0) * \
                       ((xNE - xSE) * (ySE * ySE + ySE * yNE + yNE * yNE) + 
                        (xNW - xNE) * (yNE * yNE + yNE * yNW + yNW * yNW) +
                        (xSW - xNW) * (yNW * yNW + yNW * ySW + ySW * ySW) + 
                        (xSE - xSW) * (ySW * ySW + ySW * ySE + ySE * ySE));

            rho_s = rho_source(centre_x, centre_y)
            u_s = xmom_source(centre_x, centre_y)
            v_s = ymom_source(centre_x, centre_y)
            E_s = energy_source(centre_x, centre_y)

            solution.append([rho_s, u_s, v_s, E_s])
    print "Begin writing analytical solution..."
    fp = open("Euler_source_terms.vtk", "w")
    fp.write("# vtk DataFile Version 2.0\n")
    fp.write("Analytical solution for verification purposes\n")
    fp.write("ASCII\n")
    fp.write("\n")
    fp.write("DATASET UNSTRUCTURED_GRID\n")
    no_points = grid.ni*grid.nj
    fp.write("POINTS %d float\n" % (no_points))
    for i in range(grid.ni):
        for j in range(grid.nj):
            fp.write("%15.6e%15.6e%15.6e\n" % (grid.x[i][j], grid.y[i][j], 0.0))
    no_cells = (grid.ni-1)*(grid.nj-1)
    fp.write("CELLS %d %d\n" % (no_cells, no_cells*5))
                
    for i in range(grid.ni-1):
        for j in range(grid.nj-1):
            SW = i*grid.nj + j
            NW = SW + 1
            SE = SW + grid.nj
            NE = SE + 1
            fp.write("%d %d %d %d %d\n" % (4, SW, SE, NE, NW))

    fp.write("CELL_TYPES %d\n" % (no_cells))
    for i in range(grid.ni-1):
        for j in range(grid.nj-1):
            fp.write("%d\n" % (9)) # VTK_QUAD == 9
    fp.write("CELL_DATA %d\n" % (no_cells))
    fp.write("SCALARS rho_source float 1\n")
    fp.write("LOOKUP_TABLE default\n")
    for sol in solution:
        fp.write("%g\n" % (sol[0]))
    fp.write("SCALARS xmom_source float 1\n")
    fp.write("LOOKUP_TABLE default\n")
    for sol in solution:
        fp.write("%g\n" % (sol[1]))
    fp.write("SCALARS ymom_source float 1\n")
    fp.write("LOOKUP_TABLE default\n")
    for sol in solution:
        fp.write("%g\n" % (sol[2]))
    fp.write("SCALARS energy_source float 1\n")
    fp.write("LOOKUP_TABLE default\n")
    for sol in solution:
        fp.write("%g\n" % (sol[3]))
    fp.close()                    
            
if __name__ == '__main__': main()
    
