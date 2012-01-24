#!/usr/bin/env python
# Author: Rowan J. Gollan
# Date: 04-Sep-2009
# Place: NASA Langley, Hampton, Virginia, USA
#
# p2e.py - This program can be used to convert
# a Plot3D file format to an Eilmer3 file format.
#

import sys

from libprep3 import *
from e3prep import gdata, write_times_file
from e3_grid import *
from e3_block import *
from e3_flow import *

from optparse import OptionParser
try:
    from numpy import array, zeros
except:
    try:
        from Numeric import array, zeros
    except:
        print "Could import neither numpy nor Numeric."

parser = OptionParser()

parser.add_option("-d", "--dimensions", dest="dim", type="int",
                  help="dimensions of problem, 2 or 3")

parser.add_option("-j", "--job", dest="jobname", type="string",
                  help="write eilmer3 files with jobname JOB",
                  metavar="JOB")
parser.add_option("-p", "--plot3d", dest="p3dname",
                  help="read from Plot3D files with basename BASE",
                  metavar="BASE")
parser.add_option("-m", "--gas-model", dest="gmodel_file", metavar="GFILE",
                  help="gas model filename")
parser.add_option("-g", "--grid-only", action="store_true", dest="grid_only",
                  default=False,
                  help="only produce a grid, no solution files")

# Map Vulcan variable names to local symbols
var_map = {'static density' : 'rho',
           'u-velocity;velocity' : 'u',
           'v-velocity' : 'v',
           'w-velocity' : 'w',
           'static pressure' : 'p',
           'total pressure' : 'p_t',
           'static temperature' : 'T',
           'total temperature' : 'T_t',
           'entropy' : 's',
           'Mach no.' : 'M'
           }

def read_plot3d_function_file(bfn, dim):
    """Reads a Plot3D function file (and name file), returns dict of data.
    
    Assumptions:
    - formatted
    - multi-block
    """
    pfile = bfn + ".f"
    nfile = bfn + ".nam"
    nf = open(nfile, "r")
    var_names = []
    while 1:
        line = nf.readline()
        if not line:
            break
        var_names.append(line[:-1]) # strip off "\n"

    vars = []
    for v in var_names:
        if v in var_map:
            vars.append(var_map[v])
        else:
            print "Variable '%s' is unknown.  It will be ignored." % v

    nf.close()

    f = open(pfile, "r")
    tks = f.read().split()
    f.close()
    nblocks = int(tks[0])

    data = []
    pos = 1
    for i in range(nblocks):
        ni = int(tks[pos]); pos += 1
        nj = int(tks[pos]); pos += 1
        nk = 1
        if dim == 3:
            nk = int(tks[pos]); pos += 1

        # Ignore number indicating how many fields of data
        # We learnt this from the .nam file (so long as they
        # are consistent).
        pos += 1

        data.append({})
        for v in vars:
            data[-1][v] = zeros((ni, nj, nk), 'd')
            for k in range(nk):
                for j in range(nj):
                    for i in range(ni):
                        data[-1][v][i,j,k] = float(tks[pos]); pos += 1

    return data

def prepare_data_for_eilmer3(data, grids, dim):
    """Transform Plot3D node-centred data to finite-volume data.

    This function interpolates node-centres data to
    cell-centered data.
    """

    # 1. interpolation
    data2 = []
    for ib, db in enumerate(data):
        data2.append({})
        for v in db.keys():
            (ni, nj, nk) = db[v].shape
            data2[-1][v] = zeros((ni-1, nj-1, nk-1), 'd')
        data2[-1]['pos.x'] = zeros((ni-1, nj-1, nk-1), 'd')
        data2[-1]['pos.y'] = zeros((ni-1, nj-1, nk-1), 'd')
        data2[-1]['pos.z'] = zeros((ni-1, nj-1, nk-1), 'd')
        data2[-1]['vol'] = zeros((ni-1, nj-1, nk-1), 'd')
        
        for k in range(nk-1):
            for j in range(nj-1):
                for i in range(ni-1):
                    if dim == 3:
                        p0 = Vector(grids[ib].x[i,j,k],       grids[ib].y[i,j,k],       grids[ib].z[i,j,k])
                        p1 = Vector(grids[ib].x[i+1,j,k],     grids[ib].y[i+1,j,k],     grids[ib].z[i+1,j,k])
                        p2 = Vector(grids[ib].x[i+1,j+1,k],   grids[ib].y[i+1,j+1,k],   grids[ib].z[i+1,j+1,k])
                        p3 = Vector(grids[ib].x[i,j+1,k],     grids[ib].y[i,j+1,k],     grids[ib].z[i,j+1,k])
                        p4 = Vector(grids[ib].x[i,j,k+1],     grids[ib].y[i,j,k+1],     grids[ib].z[i,j,k+1])
                        p5 = Vector(grids[ib].x[i+1,j,k+1],   grids[ib].y[i+1,j,k+1],   grids[ib].z[i+1,j,k+1])
                        p6 = Vector(grids[ib].x[i+1,j+1,k+1], grids[ib].y[i+1,j+1,k+1], grids[ib].z[i+1,j+1,k+1])
                        p7 = Vector(grids[ib].x[i,j+1,k+1],   grids[ib].y[i,j+1,k+1],   grids[ib].z[i,j+1,k+1])
                        centre = hexahedron_centroid(p0, p1, p2, p3, p4, p5, p6, p7)
                        vol = hexahedron_volume(p0, p1, p2, p3, p4, p5, p6, p7)

                        a = quad_centroid(p0, p1, p2, p3)
                        b = quad_centroid(p1, p5, p6, p2)
                        c = quad_centroid(p4, p5, p6, p7)
                        d = quad_centroid(p0, p4, p7, p3)
                        e = quad_centroid(p0, p1, p5, p4)
                        f = quad_centroid(p3, p2, p6, p7)

                        p01 = 0.5*(p0 + p1)
                        p12 = 0.5*(p1 + p2)
                        p23 = 0.5*(p2 + p3)
                        p30 = 0.5*(p3 + p0)
                        p15 = 0.5*(p1 + p5)
                        p56 = 0.5*(p5 + p6)
                        p62 = 0.5*(p6 + p2)
                        p67 = 0.5*(p6 + p7)
                        p74 = 0.5*(p7 + p4)
                        p45 = 0.5*(p4 + p5)
                        p73 = 0.5*(p7 + p3)
                        p40 = 0.5*(p4 + p0)

                        w0 = hexahedron_volume(centre, b, p62, f, c, p56, p6, p67) / vol
                        w1 = hexahedron_volume(d, centre, f, p73, p74, c, p67, p7) / vol
                        w2 = hexahedron_volume(p40, e, centre, d, p4, p45, c, p74) / vol
                        w3 = hexahedron_volume(e, p15, b, centre, p45, p5, p56, c) / vol
                        w4 = hexahedron_volume(a, p12, p2, p23, centre, b, p62, f) / vol
                        w5 = hexahedron_volume(p30, a, p23, p3, d, centre, f, p73) / vol
                        w6 = hexahedron_volume(p0, p01, a, p30, p40, e, centre, d) / vol
                        w7 = hexahedron_volume(p01, p1, p12, a, e, p15, b, centre) / vol

                        for v in db.keys():
                            data2[-1][v][i,j,k] = w0*db[v][i,j,k] + \
                                                  w1*db[v][i+1,j,k] + \
                                                  w2*db[v][i+1,j+1,k] + \
                                                  w3*db[v][i,j+1,k] + \
                                                  w4*db[v][i,j,k+1] + \
                                                  w5*db[v][i+1,j,k+1] + \
                                                  w6*db[v][i+1,j+1,k+1] + \
                                                  w7*db[v][i,j+1,k+1]
                        data2[-1]['pos.x'][i,j,k] = centre.x
                        data2[-1]['pos.y'][i,j,k] = centre.y
                        data2[-1]['pos.z'][i,j,k] = centre.z
                        data2[-1]['vol'][i,j,k] = vol
                    else:
                        print "2D version not yet implemented."
                        print "Bailing out!"
                        sys.exit(1)

    return data2

def write_flow_file(f, d):
    f.write("%20.12e\n" % 0.0)
    f.write("%s\n" % variable_list_for_cell(gdata))
    (nni, nnj, nnk) = d['pos.x'].shape
    f.write("%d %d %d\n" % (nni, nnj, nnk)) # number of cells in each dir
    for k in range(nnk):
        for j in range(nnj):
            for i in range(nni):
                fc = FlowCondition(p=d['p'][i,j,k],
                                   u=d['u'][i,j,k],
                                   v=d['v'][i,j,k],
                                   w=d['w'][i,j,k],
                                   T=d['T'][i,j,k])
                cell_properties = fc.to_dict()
                cell_properties['pos.x'] = d['pos.x'][i,j,k]
                cell_properties['pos.y'] = d['pos.y'][i,j,k]
                cell_properties['pos.z'] = d['pos.z'][i,j,k]
                cell_properties['volume'] = d['vol'][i,j,k]
                write_cell_data(f, cell_properties, gdata)
    return
    
def main():
    (options, arg) = parser.parse_args()

    if options.dim is None:
        print "The dimensions for the grid must be specified with -d DIM"
        print ""
        parser.print_help()
        sys.exit(1)

    if options.jobname is None:
        print "A job filename must be specified with --job=JOB"
        print ""
        parser.print_help()
        sys.exit(1)

    if options.p3dname is None:
        print "A plot3d base name must be specified with --plot3d=BASE"
        print ""
        parser.print_help()
        sys.exit(1)

    if options.gmodel_file is None:
        print "A gas model file must be specified with --gas-model=GFILE"
        print ""
        parser.print_help()
        sys.exit(1)

    bfn, ext = os.path.splitext(options.jobname)

    # 1. Work on grid.

    # 1a. Read in grid.

    gname = options.p3dname + ".g"
    print "Reading Plot3D grid(s) from: ", gname
    grids = read_plot3d_grid(gname, options.dim)

    # 1b. Write out grid.
    
    print "Writing Eilmer3 grid(s)"
    gridPath = os.path.join("grid", "t0000")
    if not os.access(gridPath, os.F_OK):
        os.makedirs(gridPath)
    for i, g in enumerate(grids):
        fname = bfn+(".grid.b%04d.t0000" % i)
        fname = os.path.join(gridPath, fname)
        f = open(fname, 'w')
        g.write(f)
        f.close()
        os.system("gzip -f "+fname)

    if options.grid_only:
        print "Done."
        return

    # 2. Work on data
    print "Reading Plot3D data from: ", options.p3dname
    data = read_plot3d_function_file(options.p3dname, options.dim)
    data2 = prepare_data_for_eilmer3(data, grids, options.dim)

    set_gas_model(options.gmodel_file)

    print "Writing Eilmer3 flow data file(s)"
    flowPath = os.path.join("flow", "t0000")
    if not os.access(flowPath, os.F_OK):
        os.makedirs(flowPath)
    for i,d in enumerate(data2):
        fname = bfn+(".flow.b%04d.t0000" % i)
        fname = os.path.join(flowPath, fname)
        f = open(fname, 'w')
        write_flow_file(f, d)
        f.close()
        os.system("gzip -f "+fname)

    print "Writing dummy .times file."
    write_times_file(bfn)

    print "Writing dummy block_lebels.llst."
    bll = open('block_labels.list', 'w')
    bll.write('# indx label\n')
    for i in range(len(data)):
        bll.write('%d block-%d\n' % (i, i))
    bll.close()
    

    print "Done."

if __name__ == '__main__':
    main()

    
    



                  
