"""
	this is an interface module to the program idcoef by Mark Abney.
	idcoef software is copyright Mark Abney 2004, 2005, 2007, 2008 and licensed
	under the GPL 2.0. I may be contacted by email at abney@uchicago.edu. 
	
	This module knows how to:
		*) run the programs given the input files and 
		*) read the output and calculate kinship for a given pair
		*) clean temp files

"""
from numpy import array
import sys
import os
from subprocess import call, check_call
DEBUG=False
PROG_PATH = 'idcoef/Idcoefs2.1.1/idcoefs'


# helper
class condensedIBDcoeffs():
	"""
		a class to manipulate the 9 condensed IBD states
		also provides K, Z0, Z1 and Z2 calculation
		Notes:
		When wishing to use the kinship coeff for heritability calc it
			is advised to use the calc implemented here from the condensed coeffs
		However, if Z1 and Z2 are required the calc here will return the correct sharing
			provided that you do not use Kz=1.5*Z1+Z2. this calc will give the wrong answer
			if inbreeding is present!
	"""
	def __init__(self,state_list=None):
		self.state_list=state_list
		if state_list==None:
			self.state_list=[0]*9
		if len(self.state_list)!=9:
			print '<condensedIBDcoeffs> wrong input for state list: ', state_list
			raise ValueError
			
		delta=dict()
		for i in range(len(self.state_list)):
			delta[i+1]=self.state_list[i]
		self.delta=delta
	
	def __str__(self):
		return str(self.state_list)
		
	def calcK(self):
		""" calc the kinship coeff pi-hat  """
		D=self.delta
		return 2*D[1] + D[3]+D[5]+D[7] + 0.5*D[8]
	
	def calcZ1(self):
		"""
			The output here depends on the definition of Z1
			Note that in the current definition 1/2*Z1+Z2 != K(condensed kinship)
				unless there is no inbreeding
		"""
		D=self.delta
		return D[3]+D[5]+D[8]
	
	def calcZ2(self):
		"""
			The output here depends on the definition of Z1
			Note that in the current definition 1/2*Z1+Z2 != K(condensed kinship)
				unless there is no inbreeding
		"""
		D=self.delta
		return D[1] + D[7]
		
class progUI():
	def __init__(self,prog_path=PROG_PATH):
		self.prog_path = prog_path
	
	def run(self,pedfile,studyfile,outfile, ram=500):
		"""
			ram is in megabytes, i.e. 1Gb is 1000
			example of command line: idcoefs -p ex.pedigree -s ex.study -o new.output -r 200
		"""
		cmd = self.prog_path + ' -p ' + pedfile + ' -s ' + studyfile + ' -o ' + outfile + ' -r ' + str(ram)

		call_result=None
		if DEBUG==False:
			call_result=call(cmd,shell=True)
			#print cmd+'\n'+'idcoef program return with:',call_result
			#sys.stdout.flush()
		return call_result			
	
	def readOutput(self,outfile):
		"""
			reads the output file of the program
			returns the condensed stated for each pair of individuals 
			the returned object is a class that manipulated condensed states
			
		"""
		file_lines = open(outfile).read().split('\n')
		if len(file_lines[-1])<2:
			file_lines = file_lines[:-1]
		# parse the lines and init a hash table for each pair
		# key is a tuple (id1,id2) (id's are strings)
		idcoeff_table=dict()
		uniq_id_pairs=[]
		for line in file_lines:
			splt=line.split()
			# get the ids
			ID1 = splt[0]
			ID2 = splt[1]
			uniq_id_pairs.append((ID1,ID2))
			# get the states
			STATES = array(splt[2:],dtype=float)
			idcoeff_table[(ID1,ID2)] = idcoeff_table[(ID2,ID1)] =  condensedIBDcoeffs(STATES)
			
		return idcoeff_table, uniq_id_pairs

RAM=1000		
def globalrunner(input_files):
	"""	running the program on the input parameters """
	pedfile   = input_files[0]
	studyfile = input_files[1]
	outfile   = input_files[2]
	global RAM
	ram = RAM
	# debug
	# print 'running RAM = ', ram
	launcher = progUI()
	call_result = launcher.run(pedfile,studyfile,outfile,ram)
	if call_result!=0:
		return 1
	return 0


