import time
import sys
import os
import warnings

try:
	import h5py
except ImportError:
	warnings.warn("There is no h5py package.")
	
from numpy import *
from scipy import *
from pylab import *

class gdanxs:
	def __init__(self, prefix, scannumber):
		"""
		Read scan data and set variables
		from nexus file format
		"""
		
		self.prefix = prefix
		self.scannumber = scannumber
		sc = str(self.scannumber)
		self.suffix = ".nxs"
		self.filename = ''.join([self.prefix,sc,self.suffix])
		self.nexus=h5py.File(self.filename,"r")
		self.readdata()
		
		return
		
	def readdata(self):
		# all this data is in the self.nexus tree, but 
		#extract some for ease of use.
		
		# read scan command
		self.command = self.nexus['entry1/scan_command'][0]
		
		# extract counters
		for each_key in self.nexus['entry1/instrument'].keys():
			if each_key != ('name'):
				if each_key != ('source'):
					for sub_keys in self.nexus[''.join(['entry1/instrument/',each_key])].keys():
						if sub_keys == (each_key):
							setattr(self, each_key, self.nexus[''.join(['entry1/instrument/',each_key,'/',each_key])][:])
						if sub_keys == ('data'):
							setattr(self, each_key, self.nexus[''.join(['entry1/instrument/',each_key,'/data'])][:])
		return

class gdaScan:
	def __init__(self, scannumber):
		"""
		Read scan data and set variables
		to all the data contained in the scan
		
		"""

		# Modified from SpecScan
		
		scan="%05d" %(scannumber)
		self.filename = ''.join([scan,'.dat'])
		
		self.scandata = gdaData()
		self.data = []	
		self.values={}
		self.line=''	

		self.scanfile = open(self.filename)
		self.readheader()
		self.readdata()
		self.scanfile.close()
		
		return

	def readheader(self):

		#
		# Read the spec header and place the data into this class
		#
		
		self.line = self.scanfile.readline()

		# extract metadata
		while self.line[0:17] != '<MetaDataAtStart>':
			self.line = self.scanfile.readline()
			if self.line[0:5]==' &END':
				print "No metadata in this scan file";
				return

		# move to next line
		self.line = self.scanfile.readline()

		while self.line[0:18] != '</MetaDataAtStart>':
			if self.line == "" or self.line=="\n":
				self.line=self.scanfile.readline()
			pos = self.line.split("=")
			if self.is_number(pos[1])==True:
				self.scandata.setValue(pos[0],float(pos[1]))
			else:
				self.scandata.setValue(pos[0],pos[1])
			self.line=self.scanfile.readline()
		return
		

	def readdata(self):
		# find beginnning of data
		while self.line[0:5] != ' &END':
			self.line=self.scanfile.readline()

		# move to column headers
		self.line=self.scanfile.readline()
		
		# extract column names removing last CR 
		self.cols=self.line.split('\t')
		self.cols[-1]=self.cols[-1].rstrip('\n')
		# first scan data line
		self.line=self.scanfile.readline()
		while self.line!="":
			dataline = self.line.split('\t')
			dataline[-1]=dataline[-1].rstrip('\n')
			dataArray=zeros((len(dataline)))
			for ii in range (0,len(dataline)):
				dataArray[ii] = float(dataline[ii])
			if self.data==[]:
				self.data =dataArray
			else:
				self.data=vstack([self.data,dataArray])
			self.line=self.scanfile.readline()
		
		for ii in range (0,len(self.cols)):
			self.scandata.setValue(self.cols[ii],self.data[:,ii])
			
		for i in self.scandata.values.keys():
			setattr(self,i , self.scandata.values[i])

		return
		
	def show(self):
		"""
		Show statistics on scan
		"""
		print "Datafile:"
		print self.filename
		self.scandata.show()
		return

	def is_number(self, s):
	    try:
	        float(s)
	        return True
	    except ValueError:
	        return False
	
			
class gdaData:	
	def __init__(self):
		self.values = {}
	def __call__(self, key):
		print key
		
	def setValue(self, key, data, setdata = True):
		self.values[key] = data
	
	def get(self, key):
		if self.values.has_key(key):
			return self.values[key]
		else:
			return None
			
	def show(self):
		"""
		Show statistics on data (motors, scalars)
		
		example   : SpecData.show()
		"""
		
		print "Header Variables"
		for i in self.values.keys():
			if self.values[i].shape == (1,):
				print "%10s " % i,
		
		print "Scan Variables:"
		for i in self.values.keys():
			if self.values[i].shape != (1,):
				print "%10s " % i,
		
		return None
