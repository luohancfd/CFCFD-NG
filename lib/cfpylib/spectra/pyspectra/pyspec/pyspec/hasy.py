"""
Python module for handling of HASYLAB data files.

Written by TAW Beale

Written for Spectra data files, as used on BW5, hasylab.   
Here call HASY, rather than SPECTRA to avoid confusion with SPEC.

This is distinctly different from spec.py (in pyspec) as each scan at HASYLAB is 
held in an individual file.    Therefore the name of the file is passed
to the routine.

Example:

	# Load datafile and index
	
	scan = HasyScanFile('fourc.01')		
	
	# Access to variables
	
	scan.data							
	scan.values['TTH']					 
	scan.TTH							

	# access to suplimentary information
	
	scan.header										

"""

import time
import sys
import os
from numpy import *
from scipy import *
from pylab import *

__verbose__ = 0x01


NoFile = "Scan file does not exist"
	
def removeIllegals(key):
	illegal = ['/', '-', ' ']
	for j in illegal:
		key = key.replace(j,'')
	if key[0].isalpha() == False:
		key = "X" + key
		
	return key

		
class HasyScanFile:
	"""
	Datafile class for handling a HASYLAB file
	"""
	def __init__(self, prefix, scan):
		"""
		Read scan data from Data class and set variables
		to all the data contained in the scan
		
		"""
		self.prefix = prefix
		self.suffix = ".fio"
		self.scanno = scan		
		self.dir = dir

		self.constructfilename() 

		self.readheader()
		self.readdata()

		self.file.close()
		return

	def constructfilename(self):
		"""
		Construct the filename from the scan number and prefix
		"""
		sc = "%05d" %(self.scanno)
		self.filename = ''.join([self.prefix,'_',sc,self.suffix])
		
		return

	def readheader(self):
		"""
		 Read the hasylab header from file 
		"""

		self.file = open(self.filename, 'r')
		if __verbose__:
			print "...reading scan",self.filename

		self.motors = {}
 
		self.file.seek(0,0)
		line = self.file.readline()
		
		# find first %c instance
		while line[0:2] != '%c':
			line = self.file.readline()
		
		# move to next line
		line = self.file.readline()
		# all proceeding lines not starting with ! until %p are header...

		self.header=''

		while line[0:2] != '%p':
			if line[0:1]!="!":
				self.header = self.header + line
				# Find sampling time - only appears to be in header
				if line[0:5]==" Name":
					time = line.split("sampling")
					self.time = float(time[1].split("s")[0])
				# Do same for lattice parameters
				if line[0:12]==" lattice-par":
					self.lattice = [float(s) for s in line[14:-1].split("  ")]
					
			line = self.file.readline()

		# move to next line to read motor positions
		line = self.file.readline()
		
	
		while line[0:2] != '%d':
			if line[0:1]!="!":
				#line = removeIllegals(line)
				br = line.split("=")
				# remove \n from value
				br[1] = br[1].replace('\n','')
				br[0]=removeIllegals(br[0])
				self.motors[br[0]]=float(br[1])
			line = self.file.readline()
		return

	def readdata(self):
		"""
		Read the colomn data from the datafile
		"""
		kk=0
		self.columns=[]
		self.data=[]
		line = self.file.readline()
	
		while line[0:5] == " Col ":
			self.columns.append(line)
			line = self.file.readline()

		dataline=zeros((1,len(self.columns)))

		while line != "":
			br = line.split(" ")
			jj=0
			for ii in range(0,len(br)-1):
				if br[ii][0:1]!=" ":
					if br[ii]!="":
						dataline[0,jj]=float(br[ii])
						jj=jj+1
			if self.data==[]:
				self.data=dataline
			self.data=vstack((self.data,dataline))
			if kk==0:
				self.data = self.data[0,:]
				kk=1
			line=self.file.readline()

		# x column is always the first
		self.x = self.data[:,0]
		# counter is always the second column
		self.y = self.data[:,1]
		# find column with DORIS in (normally 5th!)
		for n in range(0,len(self.columns)-1):
			exist = self.columns[n].find("DORIS")
			if exist != -1:
				self.doris = self.data[:,n]
		# find column with C6 in (diode)
		for n in range(0,len(self.columns)-1):
			exist = self.columns[n].find("C6")
			if exist != -1:
				self.diode = self.data[:,n]
		# find column with T_CONTROL in it (t_reg)
		for n in range(0,len(self.columns)-1):
			exist = self.columns[n].find("T_CONTROL")
			if exist != -1:
				self.treg = self.data[:,n]
		# find column with T_SAMPLE in it (tsam)
		for n in range(0,len(self.columns)-1):
			exist = self.columns[n].find("T_SAMPLE")
			if exist != -1:
				self.tsam = self.data[:,n]
			
		return
