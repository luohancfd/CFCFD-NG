#!/usr/bin/env python
# tecplot_feblock.py
# Load a 2D unstructured-block of quad cells from a Tecplot file.
# PJ, 24-Apr-2010

import math

class FEBlock2D(object):
    """
    Stores a 2D unstructured-block of quad cells.
    """

    def __init__(self):
        self.title = ""
        self.variables = [] # will contain the field variable names
        self.zone = {'T':None,
                     'N':None,
                     'E':None,
                     'ET':'QUADRILATERAL',
                     'F':'BLOCK'}
        self.data = {}  # will contain lists of data values, one for each variable
        self.vlist = [] # will contain lists of indices specifying cell vertices
        return

    def read_tecplot_file(self, fileName):
        """
        Read from an ASCII Tecplot file.

        For the moment, we make a number of restrictive assumptions.
        (1) There is onely one zone.
        (2) The data is in FBLOCK format.
        (3) The cells are quadrilaterals.
        """
        fp = open(fileName, 'r')
        #
        # Get title string from first line.
        #
        line = fp.readline()
        items = line.strip().split('=')
        keyword = items[0].strip().upper()
        if keyword != 'TITLE':
            print "readTecplotFile: First line, expected title."
            print "    actual line:", line
            return
        self.title = items[-1].strip(' "')
        print 'title=', self.title
        #
        # Get variable names.
        #
        line = fp.readline()
        items = line.strip().split('=')
        keyword = items[0].strip().upper()
        if keyword != 'VARIABLES':
            print "readTecplotFile: Second line, expected variable list."
            print "    actual line:", line
            return
        for name in items[-1].split(','):
            self.variables.append(name.strip(' "'))
        print 'variables=', self.variables
        #
        # Get zone parameters.
        #
        line = fp.readline()
        items = line.strip().split(' ')
        keyword = items[0].strip().upper()
        if keyword != 'ZONE':
            print "readTecplotFile: Third line, expected zone parameters."
            print "    actual line:", line
            return
        # Re-assemble the line so that we can further split into
        # items of the form name=value
        line = ''
        for item in items[1:]: line += item
        items = line.split(',')
        for item in items:
            name, value = item.split('=')
            if name == 'T': self.zone['T'] = value.strip(' "')
            if name == 'N': self.zone['N'] = int(value)
            if name == 'E': self.zone['E'] = int(value)
            if name == 'ET' and value != 'QUADRILATERAL':
                print "readTecplotFile: Third line, expected QUADRILATERAL elements."
                return
            if name == 'F' and value != 'FEBLOCK':
                print "readTecplotFile: Third line, expected FEBLOCK format."
                return
        print "zone=", self.zone
        #
        # Now, accumulate the floating-point data for each field variable.
        #
        def accumulate_floats(fp, N):
            """
            Accumulate the N floating-point data values and return as a list.
            """
            data = []
            while len(data) < N:
                line = fp.readline().strip(' \n')
                for item in line.split(): data.append(float(item))
            return data
        for var in self.variables:
            self.data[var] = accumulate_floats(fp, self.zone['N'])
            print 'found data for', var, 'length=', len(self.data[var])
        #
        # Element vertex definitions.
        # Be aware that Tecplot starts counting from 1, we will start from 0.
        #
        for e in range(self.zone['E']):
            line = fp.readline().strip(' \n')
            items = line.split()
            self.vlist.append([int(items[0])-1, int(items[1])-1,
                               int(items[2])-1, int(items[3])-1])
        print "number of elements found=", len(self.vlist)
        #
        fp.close()
        return

    def find_nearest(self, x, y, z):
        """
        Assuming that X,Y,Z are available, find the vertex closest to given point.
        """
        dx = x - self.data['X'][0]
        dy = y - self.data['Y'][0]
        dz = z - self.data['Z'][0]
        dist = math.sqrt(dx*dx + dy*dy + dz*dz)
        min_dist = dist
        vindx = 0
        for i in range(1,self.zone['N']):
            dx = x - self.data['X'][i]
            dy = y - self.data['Y'][i]
            dz = z - self.data['Z'][i]
            dist = math.sqrt(dx*dx + dy*dy + dz*dz)
            if dist < min_dist:
                min_dist = dist
                vindx = i
        return vindx
    
if __name__ == '__main__':
    print "Read Tecplot unstructured grid. Begin..."
    feblock = FEBlock2D()
    feblock.read_tecplot_file('tecplot_feblock_example_file.tec')
    indx = feblock.find_nearest(0.000,0.025,0.020)
    print "x=", feblock.data['X'][indx], "y=", feblock.data['Y'][indx], \
          "z=", feblock.data['Z'][indx], "p=", feblock.data['pressure'][indx]
    print "Done."
