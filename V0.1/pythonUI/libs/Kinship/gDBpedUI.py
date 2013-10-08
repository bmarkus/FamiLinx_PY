"""

	an interface for a pedigree class

"""

import pedInterface as PEDI
import sys
import traceback
import MySQLdb as mdb
from numpy import *

DEBUG=True
DEBUG1=False

# utils

# -----------
# main class
# -----------
"""
"""
class gPed(PEDI.PedI):
	"""
		Sadly, this was the implementation when I just started ... Could do nicer code
		
		db_connection should be connected already with the database
	"""
	
	def __init__(self,db_connection):
		self.db_connection=db_connection
	
	def getParents(self,ID):
		"""
			I know, this interface should have been somewhere else ... but what can you do :(
			
			return a list of parents IDs 
			return [] if no parents
		"""	
		parentID=[]			
		try:
			with self.db_connection:	
				cur = self.db_connection.cursor()			
				cur.execute("SELECT Parent_id FROM %s where Child_id = '%d'"%('relationship',ID))
				res = cur.fetchall()
				if res!=None:
					parentID=[int(item[0]) for item in res]
				cur.close()							
		except: # some specific exception you want to wrap
			exp_value = sys.exc_info()[1]
			exp_type = sys.exc_info()[0]
			print 'Excpetion_type:',exp_type
			print 'Excpetion_value:',exp_value
			traceback.print_tb(sys.exc_info()[2],file=sys.stdout)		
			parentID=[]
		
		return parentID
	
	def getDescen(self,ID):	
		"""
			I know, this interface should have been somewhere else ... but what can you do :(
			
			return a list of descendants IDs (direct children) 
			return [] is no descendants
		"""
		descID=[]			
		try:
			with self.db_connection:	
				cur = self.db_connection.cursor()			
				cur.execute("SELECT Child_id FROM %s where Parent_id = '%d'"%('relationship',ID))
				res = cur.fetchall()
				if res!=None:
					descID=[int(item[0]) for item in res]
				cur.close()							
		except: # some specific exception you want to wrap
			exp_value = sys.exc_info()[1]
			exp_type = sys.exc_info()[0]
			print 'Excpetion_type:',exp_type
			print 'Excpetion_value:',exp_value
			traceback.print_tb(sys.exc_info()[2],file=sys.stdout)		
			descID=[]
		return descID
	
		
	def getAllSibs(self,ID):
		"""
			return all siblings and half siblings
		"""
		# 1st get the parents
		parents = self.getParents(ID)
		if len(parents)==0:
			return []		
		sibs = []
		for p in parents:
			d = self.getDescen(p)
			sibs.extend(d)
		# make unique
		sibs=unique(sibs).tolist()	
		# take out the test indiviual
		sibs.pop(sibs.index(ID))
		return sibs
				
		
		
	def getIDs(self):
		"""
			return a list of IDs for all entries in the pedigree
			since this is too big for the large database it is disabled!
		"""
		return []

	# -----------------------------------------
	# general methods
	
	def exportLinkage(self,IDs,SORT_PED=True):
		"""
			given a list of IDs, export a sub pedigree to a list with linkage format
			Return 
				a PEDI::pedFixer object with the pedigree info
							
			Input	
				SORT_PED: bool, should the pedigree be fixed and sorted prior to writing
			
			In this version we also fix individuals that have only 1 parent.
			
		""" 
		
		# init a fixer class to sort and fix the pedigree into a valid linkage format
		fixer = PEDI.pedFixer_random_sex()
		for ID in IDs:
			parents = self.getParents(ID)
			if len(parents)==0:
				parents=[0,0]
			elif len(parents)==1:
				# debug 
				# print 'found 1 parent!'
				# add another parent with arbitrary num (-1*other parent id)
				# this is kosher since we map the ids to integers before export!
				parents.append(-parents[0])
			# add an entry to the fixer 
			fixer.addEntry(ID,parents[0],parents[1])
		if SORT_PED:
			fixer.sortPed()			
		return fixer


def init_ped(user_name, password, db_name):
	"""	Global function for database interface connection """
	#print 'initiating a connection to the database!'
	db_connection=None;
	db_connection = mdb.connect('localhost', user_name, password, db_name)
	if db_connection==None:
		print 'Error: could not establish a connection to the database'
		return None
	db_ped = gPed(db_connection)
	return db_ped
	


# -----------------------------------------------
# interface functions
# -----------------------------------------------
class UI():
	"""	
		since we initiate a connection to the database within the function it is safe for multi threading  
		these functions return a list of IDS or an empty list if nothing is found
	"""
	@staticmethod
	def getAllSibs(ID,user_name=None, password=None, db_name=None):	
		db_ped = init_ped(user_name, password, db_name)
		return db_ped.getAllSibs(ID)
	@staticmethod
	def getChildren(ID,user_name=None, password=None, db_name=None):
		db_ped = init_ped(user_name, password, db_name)
		return db_ped.getDescen(ID)		
	@staticmethod
	def getParents(ID,user_name=None, password=None, db_name=None):
		db_ped = init_ped(user_name, password, db_name)
		return db_ped.getParents(ID)
	
		

