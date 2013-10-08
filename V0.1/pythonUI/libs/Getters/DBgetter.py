"""
	exploring the table by key=Id
"""
from DBconnect import DBconnect, DBconnectBase
import sys
import traceback


class DbGetByID():
	
	def __init__(self,db_name=None, user_name=None, password=None):
		self.db_name=db_name
		self.user_name=user_name
		self.password=password
		self.db = DBconnectBase(db_name=self.db_name, user_name=self.user_name, password=self.password)
	
	def getAge(self,Id):
		"""
			return Age from the age table
		"""
		age=None
		db_connection=self.db.db_connection
		try:
			with db_connection:	
				cur = db_connection.cursor()			
				table_name='age'
				cur.execute("SELECT Age FROM %s where Id = '%d'"%(table_name,Id))
				res = cur.fetchone()
				if res!=None:
					#print 'not none: ', res
					res=res[0]
					age = int(res)
				cur.close()				
			

		except: # some specific exception you want to wrap
			exp_value = sys.exc_info()[1]
			exp_type = sys.exc_info()[0]
			print 'Excpetion_type:',exp_type
			print 'Excpetion_value:',exp_value
			traceback.print_tb(sys.exc_info()[2],file=sys.stdout)		
		
		return age
	
	def getGender(self,Id):
		"""
			return gender
		"""
		gender=None
		db_connection=self.db.db_connection
		try:
			with db_connection:	
				cur = db_connection.cursor()			
				table_name='gender'
				cur.execute("SELECT Gender FROM %s where Id = '%d'"%(table_name,Id))
				res = cur.fetchone()
				if res!=None:
					#print 'not none: ', res
					res=res[0]
					gender = int(res)
				cur.close()				
			

		except: # some specific exception you want to wrap
			exp_value = sys.exc_info()[1]
			exp_type = sys.exc_info()[0]
			print 'Excpetion_type:',exp_type
			print 'Excpetion_value:',exp_value
			traceback.print_tb(sys.exc_info()[2],file=sys.stdout)		
		
		return gender

	def getYears(self,Id):
		"""
			return birth and death years
		"""
		years=[None,None]
		db_connection=self.db.db_connection
		try:
			with db_connection:	
				cur = db_connection.cursor()			
				table_name='years'
				cur.execute("SELECT Byear, Dyear FROM %s where Id = '%d'"%(table_name,Id))
				res = cur.fetchall()
				if res!=None and len(res)>0:
					#print 'not none: ', res
					res=res[0]
					years = [int(res[0]),int(res[1])]
				cur.close()				
			

		except: # some specific exception you want to wrap
			exp_value = sys.exc_info()[1]
			exp_type = sys.exc_info()[0]
			print 'Excpetion_type:',exp_type
			print 'Excpetion_value:',exp_value
			traceback.print_tb(sys.exc_info()[2],file=sys.stdout)		
		
		return years
		

	def getCoord(self,Id):
		"""
			return Longitude and Latitude
		"""
		coord=[None,None]
		db_connection=self.db.db_connection
		try:
			with db_connection:	
				cur = db_connection.cursor()			
				table_name='location'
				cur.execute("SELECT Lon, Lat FROM %s where Id = '%d'"%(table_name,Id))
				res = cur.fetchall()
				if res!=None and len(res)>0:
					#print 'not none: ', res
					res=res[0]
					coord = [float(res[0]),float(res[1])]
				cur.close()				
			

		except: # some specific exception you want to wrap
			exp_value = sys.exc_info()[1]
			exp_type = sys.exc_info()[0]
			print 'Exception_type:',exp_type
			print 'Exception_value:',exp_value
			traceback.print_tb(sys.exc_info()[2],file=sys.stdout)		
		
		return coord
		

	
		
		
class UI():
	@staticmethod	
	# get age from the pruned age table (not by birth-death years)
	def getAge(Id,db_name=None, user_name=None, password=None):
		age_expl=DbGetByID(db_name=db_name, user_name=user_name, password=password)
		return age_expl.getAge(Id)
		
	@staticmethod	
	# return birth and death years
	def getYears(Id,db_name=None, user_name=None, password=None):
		expl=DbGetByID(db_name=db_name, user_name=user_name, password=password)	
		return 	expl.getYears(Id)

	@staticmethod	
	# 1=male 2=female 
	def getGender(Id,db_name=None, user_name=None, password=None):
		expl=DbGetByID(db_name=db_name, user_name=user_name, password=password)
		return expl.getGender(Id)
			
	@staticmethod	
	# return Longitude and Latitude
	def getCoord(Id,db_name=None, user_name=None, password=None):
		expl=DbGetByID(db_name=db_name, user_name=user_name, password=password)
		return expl.getCoord(Id)
	

