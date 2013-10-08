"""
	This Module serves as an interface to the Familinx MySQL database
	It uses python modules from: libs/
	The interface with the database is made with the MySQLdb library
	Some interfaces could be extended and new interface could easily be written with the examples in this module
	
	Configuring this module with the database: 
		Make sure you have a valid user name with privileges to access the MySQL database. 
		The database name, user name and password for accessing the database are hard-coded in the file: libs/__init__.py
	
	Interface classes:
		Getter: provides functions for fetching data by individual Id
		IBDctorUI: provides function for the calculation and update of pw-IBD
		DBconnection: allows a connection to the data base with user specific queries
		tableItor: allows to iterate over a table
	
	Use ? over the classes for help, for example:
	>> ?familinxUI.Getter
	
	Dependencies:
	Python 2.6.5 or above
	MySQLdb (http://mysql-python.sourceforge.net/)	
	NumPy (http://www.numpy.org/)
	
"""

import MySQLdb as mdb
import os
import sys	
import multiprocessing

def updatePathSingle(path_name):
    """ update sys.path with a single pathName """
    sys.path.append(path_name);
    
# ADD LIBS FROM DB-MODULES
#updatePathSingle('..')

# additional libraries
from libs.Getters.DBgetter import UI as DBgetterUI
from libs.Kinship.gDBpedUI import UI as pedUI
import libs.Kinship.pairIBDcalctor as pwIBDctor
import libs.Getters.DBconnect as DBconnect
from libs import DATABASE, USER, PASSWORD



class tableItor():
	"""
		methods for iterating over tables
			tableItor.Itor
			tableItor.tableLen
	"""
	@staticmethod	
	def Itor(table_name,start=0,end=None):
		"""
			return an iterator for a table with optional positioning 
			start - index starting from 0
			end - index (or None to iterate to the last location in the table)
			
			example - iterate through the parent table:
			
			itor = familinxUI.tableItor.Itor('relationship')
			for i in range(10):
				print itor.next()
			
			>>	(1L, 1L, 4L)
				(2L, 1L, 5L)
				(3L, 2L, 4L)
				(4L, 2L, 5L)
				(5L, 3L, 4L)
				(6L, 3L, 5L)
				(7L, 6L, 10L)
				(8L, 6L, 11L)
				(9L, 7L, 6L)
				(10L, 7L, 12L)

			
		"""
		db_obj=DBconnect.init_database_obj(db_name=DATABASE, user_name=USER,  password=PASSWORD)
		# initiate an iterator
		return db_obj.tableItorChunks(table_name,start=start,end=end,jump=10000)
		
	@staticmethod
	def tableLen(table_name):
		"""
			tableLen(table_name)
			returns the length of a table
			Example:
			>> familinxUI.tableItor.tableLen('age')
			>> 1039321L

		"""
		return DBconnect.utils.tableLen(table_name,db_name=DATABASE, user_name=USER,  password=PASSWORD)

class DBconnection():
	"""
		Interface for a direct database connection 
		Available functions:
			getDBconnection
			makeQuery	
	"""
	@staticmethod	
	def getDBconnection():
		"""
			getDBconnection()
	
			connect to the database and returns a connection (MySQLdb connection object)
		"""
		db_connection=None
		user_name=USER
		password=PASSWORD
		db_name=DATABASE
		db_connection = mdb.connect('localhost', user_name, password, db_name);
		return db_connection
	
	@staticmethod		
	def makeQuery(query_string,cursor_type='reg'):	
		"""
			makeQuery(query_string,cursor_type='reg')
			
			return a regular/dictionary cursor with the query
			parameters:
				cursor_type='reg','dict'
			
			example:
			
			# get a cursor of dictionaries for each entry in the query
			cur=familinxUI.DBconnection.makeQuery("select * from location limit 10",cursor_type='dict')
			
			# print the 1st 	
			print cur.fetchone()
			>>	{'Country': 'RU', 'Lon': 31.863, 'Continent': 'EU', 'Res': 'Town', 'Lat': 61.369799999999998, 'Id': 33124684L}


			# go back to the 1st entry	
			cur.scroll(0, mode='absolute')
			
			# make a list of all Lon,Lat coordinates
			coords=[]
			for i in range(cur.rowcount):
				current_entry=cur.fetchone()
				coords.append((current_entry['Lon'],current_entry['Lat']))

			print coords
			>>		[(31.863, 61.369799999999998),
					 (-98.308999999999997, 56.954700000000003),
					 (-93.389899999999997, 41.938299999999998),
					 (-2.4611700000000001, 54.969499999999996),
					 (-83.310199999999995, 42.295900000000003),
					 (-85.690899999999999, 37.822400000000002),
					 (-3.7192500000000002, 50.728099999999998),
					 (-3.2027700000000001, 55.9542),
					 (-82.113, 29.2911),
					 (0.35602, 51.784700000000001)]
			
				
		"""
		db=DBconnection.getDBconnection()
		if cursor_type=='dict':
			cur=db.cursor(mdb.cursors.DictCursor)
		else:
			cur=db.cursor()
		cur.execute(query_string)
		return cur
	

	
class Getter():
	"""
		Interfaces using the database identifier for each individual (Id):	
		------------------------------------------------------------------
		Getter.getAge(Id)
			returns the age from the Age table (int or None)
			note that not all individuals with birth and death years have entry in this table
			for more info please look at the documentation 
			
		Getter.getYears(Id)	
			return birth and death years 
			[None, None] if no entry
			
		Getter.getGender(Id)
			return gender (1=male, 2=female, None)
		

		Getter.getCoord(Id)
			return Longitude and Latitude
			[None, None] if no entry

		Getter.getChildren(Id)	

		Getter.getParents(Id)
			return a list of children/parents
	
	"""
	@staticmethod
	def getAge(Id):
		"""
			Returns age from the age table
			Note: this is pruned data and therefore represents a small portion of the available data
		"""
		return DBgetterUI.getAge(Id,db_name=DATABASE, user_name=USER,  password=PASSWORD)
		
	@staticmethod
	def getYears(Id):	
		"""
			(birth, death) Years
		"""
		return DBgetterUI.getYears(Id,db_name=DATABASE, user_name=USER,  password=PASSWORD)

	@staticmethod	
	def getGender(Id):	
		"""
			1=male 2=female 
		"""
		return DBgetterUI.getGender(Id,db_name=DATABASE, user_name=USER,  password=PASSWORD)
			
	@staticmethod	
	def getCoord(Id):
		"""
			return Longitude and Latitude
		"""
		return DBgetterUI.getCoord(Id,db_name=DATABASE, user_name=USER,  password=PASSWORD)
	
	@staticmethod	
	def getChildren(ID):	
		return pedUI.getChildren(ID,db_name=DATABASE, user_name=USER,  password=PASSWORD)
	
	@staticmethod
	def getParents(ID):
		return pedUI.getParents(ID,db_name=DATABASE, user_name=USER,  password=PASSWORD)
	
		
class IBDctorUI():
	@staticmethod
	def calcPW(id1,id2,G=7,RAM=2000,TMP_DIR='./tmp'):
		"""
			K, CID = calcPW(id1,id2)
			K, CID = calcPW(id1,id2,G=7,RAM=2000,TMP_DIR='.')
			==================================================
			Calculate 9 IBD coefficients and the overall kinship coefficient for the pair (id1,id2) using idcoef (Mark Abney) program. 
			Note that the program relies on an interface class for running idcoef.
			In particular, the interface class holds the path to the program.
			(see: idcoefInterface.py)

			Parameters:
			-----------
				G = num of generations (7 is the default)
				RAM = memory to use (2GB is the default)
				TMP_DIR = directory for temporary files (input/output) for idcoeff program
				
			Returned value:
			--------------
				K = overall kinship coefficient
				CID = a list of 9 condensed id coefficients
		"""
		res = pwIBDctor.UIcaller(id1,id2,G=G,RAM=RAM,TMP_DIR=TMP_DIR,db_name=DATABASE, user_name=USER,  password=PASSWORD)
		K = res.calcK()
		CID = res.state_list
		return K,CID 
	
	@staticmethod
	def calcPWlist(list_of_pairs,Nthreads=5,G=7,RAM=2000,TMP_DIR='./tmp'):
		"""
			outputlist = calcPWlist(list_of_pairs,Nthreads=5,G=7,RAM=2000,TMP_DIR='./tmp')
			
			Calculates PW-IBD coeffs as in calcPW() over a list of pairs with multi-processing option
			list_of_pairs = list(tuples) with pair ids, i.e. [(1,2),(1,3),(2,3) ...]
			Nthreads = number of processors to apply (Note, should be smaller than the available cores :(o) )
			Other parameters as in calcPW()
			
			returns 
				list of tuples (K,(9 id-coeffs))
				item i in the output correspond to pair i
				
			example:
			# get the 1st 10 entries in the age table
			cur=familinxUI.DBconnection.makeQuery('select Id from age limit 10')
			res = cur.fetchall()			
			ids = [int(item[0]) for item in res]
			# make combination of pairs (45 unique pairs)
			pairids=[]
			for i in range(10):
				for j in range(i+1,10):
					pairids.append((ids[i],ids[j]))
			# run IBD calculation for these pairs using 10 threads (fast)
			ibds = familinxUI.IBDctorUI.calcPWlist(pairids,Nthreads=10)
			# show the results for the kinship coefficient
			for i,p in enumerate(pairids):
				print p,': ',ibds[i][0]
			
		"""
		pool=multiprocessing.Pool(Nthreads)
		params = [[pair[0],pair[1],G,RAM,TMP_DIR,DATABASE,USER,PASSWORD] for pair in list_of_pairs]
		#for p in params:
		#	print p
		print 'Performing IBD calculations with ', Nthreads,' threads!'
		res=pool.map(pwIBDctor.parameterWrapper,params)
		pool.close()
		pool.join()
		# parse the output
		output = [(item.calcK(),item.state_list) for item in res]
		return output
		
