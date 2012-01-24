## \file grid.py
## \ingroup geom
## \brief Block-structured 3D grid functions.
## \author PJ
##
## \version 02-Dec-2004 was import_grid.py originally
## \version 05-Apr-2005 put into this (grid) module
## 

import sys
import Numeric
from roberts import distribute_points_1
from geom_elmer import Edge3D, ParametricVolume

class BlockGrid3D(object):
    """
    Storage and service functions for a mesh of points defining
    the cell vertices within a block-structured grid.
    """
    def __init__(self, ni=None, nj=None, nk=None, label=None, qlist=[]):
        """
        @param ni: number of points in the i-index direction
        @type ni: int
        @param nj: number of points in the j-index direction
        @type nj: int
        @param nk: number of points in the j-index direction
        @type nk: int
        @param label: string label for the block
        @type label: string
        @param qlist: list of property names to be stored at the
            mesh points
        @type qlist: list of strings
        
        @note: If the number of vertices are specified in each direction,
            we actually create the storage arrays now.
        """
        self.ni = ni
        self.nj = nj
        self.nk = nk
        if ni != None and nj != None and nk != None:
            print "New BlockGrid3D: number of vertices ni=", ni, "nj=", nj, "nk=", nk
            self.init_arrays(qlist)
        else:
            print "New BlockGrid3D: unknown size"
            self.x = self.y = self.z = None
        self.label = label
        return

    def init_arrays(self, qlist=[]):
        """
        Create the storage arrays at the previously-specified sizes.
        """
        self.x = Numeric.zeros((self.ni, self.nj, self.nk), 'd')
        self.y = Numeric.zeros((self.ni, self.nj, self.nk), 'd')
        self.z = Numeric.zeros((self.ni, self.nj, self.nk), 'd')
        self.iblank = Numeric.zeros((self.ni, self.nj, self.nk), 'i')
        self.q = []
        for qname in qlist:
            self.q.append(Numeric.zeros((self.ni, self.nj, self.nk), 'd'))
        self.v = []
        velocity_index = 0
        for velocity_index in range(3):
            self.v.append(Numeric.zeros((self.ni, self.nj, self.nk), 'd'))
        return
    
    def read_from_plot3d_whole_grid(self, f):
        """
        Read one block from plot3D whole-grid format.

        Note that
        (1) the file is already opened.
        (2) the size of the grid has already been read.
        """
        print "Start reading plot3D block in whole-grid format..."
        # print "Read x-coordinates"
        np = self.ni * self.nj * self.nk
        numbers = []
        while len(numbers) < np:
            lineContent = f.readline()
            words = lineContent.split()
            for word in words:
                numbers.append(float(word))
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    self.x[i][j][k] = numbers[i+self.ni*(j+self.nj*k)]
        # print "Read y-coordinates"
        numbers = []
        while len(numbers) < np:
            lineContent = f.readline()
            words = lineContent.split()
            for word in words:
                numbers.append(float(word))
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    self.y[i][j][k] = numbers[i+self.ni*(j+self.nj*k)]
        # print "Read z-coordinates"
        numbers = []
        while len(numbers) < np:
            lineContent = f.readline()
            words = lineContent.split()
            for word in words:
                numbers.append(float(word))
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    self.z[i][j][k] = numbers[i+self.ni*(j+self.nj*k)]
        # print "Read blanking"
        numbers = []
        while len(numbers) < np:
            lineContent = f.readline()
            words = lineContent.split()
            for word in words:
                numbers.append(int(word))
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    self.iblank[i][j][k] = numbers[i+self.ni*(j+self.nj*k)]
        print "Finished reading block from plot3d whole-grid format."
        return

    def read_from_plot3d_in_planes(self, f):
        """
        Read one block from plot3D in-planes format...

        Note that
        (1) the file is already opened.
        (2) the size of the grid has already been read.

        This seems to be the format written by ICEM software
        as used by Bianca and Rowan at EPFL, Lausanne.
        
        The only place that this format seemed to be documented is the
        original Plot3D manual.
        """
        print "Start reading plot3D block in planes format..."
        np = self.ni * self.nj
        for k in range(self.nk):
            # print "Plane:", k, "Read x-coordinates"
            numbers = []
            while len(numbers) < np:
                lineContent = f.readline()
                words = lineContent.split()
                for word in words:
                    numbers.append(float(word))
            for j in range(self.nj):
                for i in range(self.ni):
                    self.x[i][j][k] = numbers[i+self.ni*j]
            # print "Plane:", k, "Read y-coordinates"
            numbers = []
            while len(numbers) < np:
                lineContent = f.readline()
                words = lineContent.split()
                for word in words:
                    numbers.append(float(word))
            for j in range(self.nj):
                for i in range(self.ni):
                    self.y[i][j][k] = numbers[i+self.ni*j]
            # print "Plane:", k, "Read z-coordinates"
            numbers = []
            while len(numbers) < np:
                lineContent = f.readline()
                words = lineContent.split()
                for word in words:
                    numbers.append(float(word))
            for j in range(self.nj):
                for i in range(self.ni):
                    self.z[i][j][k] = numbers[i+self.ni*j]
            # print "Plane:", k, "Read iblanking"
            numbers = []
            while len(numbers) < np:
                lineContent = f.readline()
                words = lineContent.split()
                for word in words:
                    try:
                        numbers.append(int(word))
                    except ValueError:
                        print "Line:", lineContent
                        raise Exception, "Cannot proceed"
            for j in range(self.nj):
                for i in range(self.ni):
                    self.iblank[i][j][k] = numbers[i+self.ni*j]
            # print "End plane", k
        print "Finished reading block from plot3d in planes format."
        return

    def write_block_in_VTK_format(self, f):
        """
        Writes the grid to an already open file.

        Note that this function writes all of the VTK header lines so
        that it is implicit that only one block grid goes into the file.
        """
        print "Begin write block: label=", self.label
        f.write("# vtk DataFile Version 2.0\n")
        f.write("%s\n" % self.label)
        f.write("ASCII\n")
        f.write("\n")
        f.write("DATASET STRUCTURED_GRID\n")
        f.write("DIMENSIONS %d %d %d\n" % (self.ni, self.nj, self.nk))
        f.write("POINTS %d float\n" % (self.ni * self.nj * self.nk))
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    f.write("%e %e %e\n" % (self.x[i][j][k], self.y[i][j][k],
                                            self.z[i][j][k]))
        print "End write block."
        return

    def write_regrid_block_in_VTK_format(self, f):
        """
        Writes the grid to an already open file.

        Note that this function writes all of the VTK header lines so
        that it is implicit that only one block grid goes into the file.
        """
        print "Begin write block: label=", self.label
        f.write("# vtk DataFile Version 2.0\n")
        f.write("%s\n" % self.label)
        f.write("ASCII\n")
        f.write("\n")
        f.write("DATASET STRUCTURED_GRID\n")
        f.write("DIMENSIONS %d %d %d\n" % (self.ni, self.nj, self.nk))
        f.write("POINTS %d float\n" % (self.ni * self.nj * self.nk))
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    f.write("%e %e %e\n" % (self.x[i][j][k], self.y[i][j][k],
                                            self.z[i][j][k]))
        f.write("POINT_DATA %d\n" % (self.ni * self.nj * self.nk))
        f.write("VECTORS velocity float\n")                
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    f.write("%e %e %e\n" % (self.v[0][i][j][k], self.v[1][i][j][k],
                                            self.v[2][i][j][k]))
                    
        f.write("SCALARS pressure float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    f.write("%e\n" % (self.q[1][i][j][k]))
                    
        f.write("SCALARS temperature float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    f.write("%e\n" % (self.q[2][i][j][k]))        
                    
        f.write("SCALARS density float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    f.write("%e\n" % (self.q[3][i][j][k]))
                    
        f.write("SCALARS internal-energy float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    f.write("%e\n" % (self.q[4][i][j][k]))  
                  
        f.write("SCALARS sound-speed float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    f.write("%e\n" % (self.q[5][i][j][k]))  
                  
        f.write("SCALARS mu float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    f.write("%e\n" % (self.q[6][i][j][k]))

                  
        f.write("SCALARS mu_t float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    f.write("%e\n" % (self.q[7][i][j][k]))  

        
                  
        f.write("SCALARS Pr float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    f.write("%e\n" % (self.q[8][i][j][k]))  
                  
        f.write("SCALARS mass-fraction[0] float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    f.write("%e\n" % (self.q[9][i][j][k]))  

        
        print "End write block."
        return

    def write_stats_block_in_VTK_format(self, f):
        """
        Writes the grid to an already open file.

        Note that this function writes all of the VTK header lines so
        that it is implicit that only one block grid goes into the file.
        """
        print "Begin write block: label=", self.label
        f.write("# vtk DataFile Version 2.0\n")
        f.write("%s\n" % self.label)
        f.write("ASCII\n")
        f.write("\n")
        f.write("DATASET STRUCTURED_GRID\n")
        f.write("DIMENSIONS %d %d %d\n" % (self.ni, self.nj, self.nk))
        f.write("POINTS %d float\n" % (self.ni * self.nj * self.nk))
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    f.write("%e %e %e\n" % (self.x[i][j][k], self.y[i][j][k],
                                            self.z[i][j][k]))

        f.write("POINT_DATA %d\n" % (self.ni * self.nj * self.nk))
        f.write("SCALARS u_dash float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    f.write("%e\n" % (self.q[0][i][j][k]))
                    
        f.write("SCALARS v_dash float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    f.write("%e\n" % (self.q[1][i][j][k]))  

                    
        f.write("SCALARS w_dash float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    f.write("%e\n" % (self.q[2][i][j][k])) 

        
        f.write("SCALARS tke float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    f.write("%e\n" % (self.q[3][i][j][k]))
                    
        f.write("SCALARS Reynolds_stress_XX float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    f.write("%e\n" % (self.q[4][i][j][k]))  

                    
        f.write("SCALARS Reynolds_stress_XY float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    f.write("%e\n" % (self.q[5][i][j][k])) 
                    
        f.write("SCALARS Vorticity_magnitude float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    f.write("%e\n" % (self.q[6][i][j][k]))
                    
        print "End write block."
        return

    def write_single_block_in_VTK_format(self, f, qlist=[]):
        """
        Writes the grid to an already open file.

        Note that this function writes all of the VTK header lines so
        that it is implicit that only one block grid goes into the file.
        """
        print "Begin write block: label=", self.label
        f.write("# vtk DataFile Version 2.0\n")
        f.write("%s\n" % self.label)
        f.write("ASCII\n")
        f.write("\n")
        f.write("DATASET STRUCTURED_GRID\n")
        f.write("DIMENSIONS %d %d %d\n" % (self.ni, self.nj, self.nk))
        f.write("POINTS %d float\n" % (self.ni * self.nj * self.nk))
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    f.write("%e %e %e\n" % (self.x[i][j][k], self.y[i][j][k],
                                            self.z[i][j][k]))
        f.write("POINT_DATA %d\n" % (self.ni * self.nj * self.nk))

        if len(qlist) == 0: return
        #
        print "write flow property data"

        for q_index in range(len(qlist)):
            if (qlist[q_index] == 'velocity'):
                f.write("VECTORS velocity float\n")
                for k in range(self.nk):
                    for j in range(self.nj):
                        for i in range(self.ni):
                            f.write("%e %e %e\n" % (self.v[0][i][j][k], self.v[1][i][j][k],
                                                    self.v[2][i][j][k]))
            else:
                f.write("SCALARS %s float 1\n" % qlist[q_index])
                f.write("LOOKUP_TABLE default\n")
                for k in range(self.nk):
                    for j in range(self.nj):
                        for i in range(self.ni):
                            f.write("%e\n" % (self.q[q_index][i][j][k]))                

        print "End write block."
        return
    
    def read_block_in_VTK_format(self, f, qlist=[]):
        """
        Reads the grid from an already open file.

        @param f     : file handle of the already-opened file
        @param qlist : list of flow quantities to be read from the file

        Note that
        (1) it is assumed that there is only one block grid in the file
        (2) the dimensions of the grid are obtained from within the file.
        """
        print "Begin read block:", 
        lineText = f.readline()  # should be "# vtk DataFile Version 2.0\n"
        self.label = f.readline()[:-1]  # drop the newline char
        print "label=", self.label
        lineText = f.readline()  # should be "ASCII\n"
        lineText = f.readline()  # should be "\n"
        lineText = f.readline()  # should be "DATASET STRUCTURED_GRID\n"
        lineText = f.readline()  # should be "DIMENSIONS %d %d %d\n"
        words = lineText.split()
        self.ni = int(words[1])
        self.nj = int(words[2])
        self.nk = int(words[3])
        print "expecting to read points for:",
        print "  ni=", self.ni, "nj=", self.nj, "nk=", self.nk
        self.init_arrays(qlist)
        #
        print "read point locations"
        lineText = f.readline()  # should be "POINTS %d float\n"
        words = lineText.split()
        nprod = float(words[1])
        assert nprod == self.ni * self.nj * self.nk, "Mismatch in number of points"
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    lineText = f.readline()  # should be "%e %e %e\n"
                    words = lineText.split()
                    self.x[i][j][k] = float(words[0])
                    self.y[i][j][k] = float(words[1])
                    self.z[i][j][k] = float(words[2])
            print ".",
        print
        if len(qlist) == 0: return
        #
        print "read flow property data"
        lineText = f.readline()  # should be "POINT_DATA %d\n"
        words = lineText.split()
        nprod = float(words[1])
        assert nprod == self.ni * self.nj * self.nk, "Mismatch in number of points"
        lineText = f.readline()
        words = lineText.split() # should be velocity
        if (qlist[0] == 'velocity'):
            which_q = qlist.index(words[1])
            print "reading %s into v[%d]" % (words[1], which_q)
            for k in range(self.nk):
                for j in range(self.nj):
                    for i in range(self.ni):
                        lineText = f.readline()  # should be "%e %e %e\n"
                        words = lineText.split()
                        self.v[0][i][j][k] = float(words[0])
                        self.v[1][i][j][k] = float(words[1])
                        self.v[2][i][j][k] = float(words[2])
            # After velocities,
            # assume qlist is in order of appearance in vtk file
            lineText = f.readline()
            words = lineText.split()
        else:
            print "Skip over velocity vectors."
            while words[0] != "SCALARS":
                lineText = f.readline()
                words = lineText.split()
        # Scalar data should be:
        # pressure, temperature, density, internal energy,
        # sound speed, mu, mu_t, Pr, f[0]
        while 1:
            print "Start of indefinite loop..."
            print "qlist=", qlist, "words[1]=", words[1]
            try:
                qlist.index(words[1])
            except ValueError:
                print "End of parameter list"
                break
            which_q = qlist.index(words[1])
            print "reading %s into q[%d]" % (words[1], which_q)
            lineText = f.readline()  # should be "LOOKUP_TABLE default"
            for k in range(self.nk):
                for j in range(self.nj):
                    for i in range(self.ni):
                        lineText = f.readline()  # should be "%e\n"
                        words = lineText.split()
                        self.q[which_q][i][j][k] = float(words[0])
            # After reading the data for a particular scalar...
            lineText = f.readline()
            if not lineText:
                # We seem to have reached the end of the file.
                break
            else:
                # Get header of next parameter.
                words = lineText.split()
        #
        print "End read block."
        return

    def make_grid_via_TFI(self, edge_list, surface_list=[None,]*6):
        """
        Given a list of edges and/or surfaces, create the grid via TFI.
        
        A list of 12 edges or a list of 6 surfaces can define a wire-frame volume.
        If both are specified, use each surface to define the shape if available,
        else fall back to using the edges.
        
        The clustering information always comes from the edges and so a full compliment
        of 12 must always be supplied, even if the geometry/path data is not definitive.
        """
        print "Begin make grid, label=", self.label
        # Set up distributions of points along each of the nondimensional edges.
        end0, end1, beta = edge_list[0].cluster_tuple
        r01 = distribute_points_1(0.0, 1.0, self.ni-1, end0, end1, beta)
        end0, end1, beta = edge_list[2].cluster_tuple
        r32 = distribute_points_1(0.0, 1.0, self.ni-1, end0, end1, beta)
        end0, end1, beta = edge_list[1].cluster_tuple
        s12 = distribute_points_1(0.0, 1.0, self.nj-1, end0, end1, beta)
        end0, end1, beta = edge_list[3].cluster_tuple
        s03 = distribute_points_1(0.0, 1.0, self.nj-1, end0, end1, beta)
        #
        end0, end1, beta = edge_list[4].cluster_tuple
        r45 = distribute_points_1(0.0, 1.0, self.ni-1, end0, end1, beta)
        end0, end1, beta = edge_list[6].cluster_tuple
        r76 = distribute_points_1(0.0, 1.0, self.ni-1, end0, end1, beta)
        end0, end1, beta = edge_list[5].cluster_tuple
        s56 = distribute_points_1(0.0, 1.0, self.nj-1, end0, end1, beta)
        end0, end1, beta = edge_list[7].cluster_tuple
        s47 = distribute_points_1(0.0, 1.0, self.nj-1, end0, end1, beta)
        #
        end0, end1, beta = edge_list[8].cluster_tuple
        t04 = distribute_points_1(0.0, 1.0, self.nk-1, end0, end1, beta)
        end0, end1, beta = edge_list[9].cluster_tuple
        t15 = distribute_points_1(0.0, 1.0, self.nk-1, end0, end1, beta)
        end0, end1, beta = edge_list[10].cluster_tuple
        t26 = distribute_points_1(0.0, 1.0, self.nk-1, end0, end1, beta)
        end0, end1, beta = edge_list[11].cluster_tuple
        t37 = distribute_points_1(0.0, 1.0, self.nk-1, end0, end1, beta)
        #
        # Assemble the volume from surfaces and/or edges.
        volume = ParametricVolume(surface_list=surface_list, edge_list=edge_list)
        #
        # Now, work through the mesh, blending the stretched parameter values
        # and creating the actual vertex coordinates in Cartesian space.
        #
        for k in range(self.nk):
            t = float(k) / (self.nk - 1)
            for j in range(self.nj):
                s = float(j) / (self.nj - 1)
                for i in range(self.ni):
                    r = float(i) / (self.ni - 1)
                    tdash = (1.0-r)*(1.0-s)*t04[k] + r*s*t26[k] + \
                            (1.0-s)*r*t15[k] + s*(1.0-r)*t37[k]
                    sdash = (1.0-t)*(1.0-r)*s03[j] + t*r*s56[j] + \
                            (1.0-t)*r*s12[j] + t*(1-r)*s47[j]
                    rdash = (1.0-s)*(1.0-t)*r01[i] + s*t*r76[i] + \
                            (1.0-s)*t*r45[i] + s*(1.0-t)*r32[i]
                    p = volume.eval(rdash, sdash, tdash)
                    self.x[i][j][k], self.y[i][j][k], self.z[i][j][k] = p.x, p.y, p.z
            print ".",
            sys.stdout.flush()
        print 
        print "End make grid."
        return
    
    def get_vertex_coords(self, ivtx):
        """
        Returns a tuple of coordinates for a single vertex.

        The indexing for ivtx corresponds to the VTK convention
        for a hexahedral cell (scaled up to the whole block).
        """
        if ivtx == 0:
            i = 0; j = 0; k = 0
        elif ivtx == 1:
            i = self.ni-1; j = 0; k = 0
        elif ivtx == 2:
            i = self.ni-1; j = self.nj-1; k = 0
        elif ivtx == 3:
            i = 0; j = self.nj-1; k = 0
        elif ivtx == 4:
            i = 0; j = 0; k = self.nk-1
        elif ivtx == 5:
            i = self.ni-1; j = 0; k = self.nk-1
        elif ivtx == 6:
            i = self.ni-1; j = self.nj-1; k = self.nk-1
        elif ivtx == 7:
            i = 0; j = self.nj-1; k = self.nk-1
        else:
            raise ValueError, ("vertex index: %d" % ivtx)
        return self.x[i][j][k], self.y[i][j][k], self.z[i][j][k]
    
#--------------------------------------------------------------------

if __name__ == '__main__':
    print "Begin demo of block3d_grid.py."

    print "Generate a grid in a box."
    b0 = BlockGrid3D(3, 4, 5, label="Some box")
    from simple_boxes import *
    edge_list = simpleBox()
    edge_list[3].cluster_tuple = (0,1,1.1)
    edge_list[5].cluster_tuple = (0,1,1.1)
    b0.make_grid_via_TFI(edge_list)

    print "Write to VTK file."
    fout = open("test.vtk", "w")
    b0.write_block_in_VTK_format(fout)
    fout.close()

    print "Read from VTK file."
    fin = open("test.vtk", "r")
    b1 = BlockGrid3D()
    b1.read_block_in_VTK_format(fin)
    fin.close()

    print "Extract corners of block."
    for ivtx in range(8):
        print "vertex", ivtx, ":", b1.get_vertex_coords(ivtx)
    
    print "Done."
