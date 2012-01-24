#! /usr/bin/env python
# csv_util.py
# Load a cloud of points from a CSV file.

import math

class CSVCloud(object):
    """
    Hold a cloud of points read from a CSV file.
    """

    def __init__(self):
        """
        Create an empty container.
        """
        self.title = ""
        self.variables = [] # will contain the field variable names
        self.N = 0          # number of points in cloud
        self.data = {}      # will contain lists of data values, one for each variable
        return

    def read(self, csvFileName):
        """
        Read from the CSV file, assuming a restricted format.

        The model is the CSV file produced by CFX postprocessor.
        """
        fp = open(csvFileName, 'r')
        #
        # Get title string from third line.
        #
        line = fp.readline()
        line = fp.readline()
        line = fp.readline()
        self.title = line.strip(' "\n')
        print 'title=', self.title
        #
        # Get variable list from sixth line
        #
        line = fp.readline()
        line = fp.readline()
        line = fp.readline()
        items = line.strip(' \n').split(',')
        for var in items: 
            self.variables.append(var.strip())
        print "variables in file=", self.variables
        # Override
        self.variables = ['I', 'X', 'Y', 'Z', 'pressure', 'S']
        print "variables renamed to", self.variables
        #
        # Read the field data points.
        #
        for var in self.variables:
            self.data[var] = []
        line = fp.readline().strip(' \n')
        while len(line) > 0:
            items = line.split(',')
            self.data['I'].append(int(items[0]))
            self.data['X'].append(float(items[1]))
            self.data['Y'].append(float(items[2]))
            self.data['Z'].append(float(items[3]))
            self.data['pressure'].append(float(items[4]))
            self.data['S'].append(float(items[5]))
            try:
                line = fp.readline().strip(' \n')
            except:
                break
        self.N = len(self.data['X'])
        print "N=", self.N
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
        indx = 0
        for i in range(1,self.N):
            dx = x - self.data['X'][i]
            dy = y - self.data['Y'][i]
            dz = z - self.data['Z'][i]
            dist = math.sqrt(dx*dx + dy*dy + dz*dz)
            if dist < min_dist:
                min_dist = dist
                indx = i
        return indx


if __name__ == '__main__':
    print "Demo CSV cloud reader."
    c = CSVCloud()
    c.read("Pressure_xyz_StreamwiseLocation_RotorBlade_CFX_Inviscid.csv")
    print c.data.keys()
    indx = c.find_nearest(0.000,0.025,0.020)
    print "indx=", indx
    print "x=", c.data['X'][indx], "y=", c.data['Y'][indx], \
          "z=", c.data['Z'][indx], "p=", c.data['pressure'][indx]
    print "Done."

