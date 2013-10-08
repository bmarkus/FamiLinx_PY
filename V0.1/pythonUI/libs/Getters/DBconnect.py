"""
	functions that explore the database
	
"""
import MySQLdb as mdb

class Struct():
	pass

# ======================
class IndexIterator():
	"""
		indexIterator(max_idx,jump)
		a basic iterator to give chunk of indices based on a jump value and a limit
		see the test example below
		Note: it will reach the max_idx in the iterations and stop
	"""
	def __init__(self,start_idx=0,max_idx=10,jump=1):
		self.max_idx=max_idx
		self.jump=jump
		self.current_loc=start_idx
	
	def indexIterator(self):
		while True:
			start = self.current_loc
			if start > self.max_idx:
				raise StopIteration
			end = self.current_loc + self.jump
			if self.current_loc + self.jump > self.max_idx:
				end = self.max_idx+1
			self.current_loc=end
			yield range(start,end)

	def regionIterator(self):
		""" yields the start and end instead of the indices """
		while True:
			start = self.current_loc
			if start > self.max_idx:
				raise StopIteration
			end = self.current_loc + self.jump
			if self.current_loc + self.jump > self.max_idx:
				end = self.max_idx+1
			self.current_loc=end
			yield (start,end)
		
	
# test	indexIterator			
def test_index_iterator(start_idx=0, max_idx=99, jump=10):
	itor1 = IndexIterator(start_idx,max_idx,jump)
	for I in itor1.regionIterator():
		print I


	
class DBconnectBase():
	def __init__(self ,db_name=None, user_name=None, password=None):
		self.db_name=db_name
		self.user_name=user_name 
		self.password=password
		
		# iteraor counter
		self.iter_counter=None # location in a table
		
		# connect to database
		self.db_connection = mdb.connect('localhost', user_name, password, db_name);
		
	def __del__(self):
		if	self.db_connection:
			self.db_connection.close()
			
class DBconnect(DBconnectBase):
	def tableItor(self,table_name=None):
		"""
			a table iterator for entire table
			yields the rows one by one
			Note: on large tables (>1000000) this will take too long
		"""
		with self.db_connection:	
			cur = self.db_connection.cursor()
			cur.execute("SELECT * FROM %s"%table_name)
			for i in range(cur.rowcount):
				yield cur.fetchone()
			
	def tableItorLimited(self,table_name=None,start=1,howmany=100):
		"""
			a table iterator from "start" returning "howmany"
			yields the rows one by one
		"""
		with self.db_connection:	
			cur = self.db_connection.cursor()
			cur.execute("SELECT * FROM %s limit %d , %d"%(table_name,start,howmany))
			#cur.execute("SELECT * FROM %s where %s > %d limit %d"%(table_name,primary_key_name,start,howmany))
			for i in range(cur.rowcount):
				yield cur.fetchone()
	
	def tableItorChunks(self,table_name,start=0,end=None,jump=10000):
		"""
			iterate from start to end (table rows) using limited jump size	each time 
		"""
		# fix the end point
		if end==None:
			max_idx=self.tableLength(table_name)
		else:	
			max_idx = min(end,self.tableLength(table_name))
		# debug
		#print 'end=',end
		#print 'start=',start	
		if max_idx-start+1 < jump:
			jump = max_idx-start+1
		itor1 = IndexIterator(start,max_idx,jump)
		for I in itor1.regionIterator():
			cur_start=I[0];
			for item in self.tableItorLimited(table_name=table_name,start=cur_start,howmany=jump):
				yield item
		
	def tableLength(self,table_name):
		"""
			the length of a table
		"""
		with self.db_connection:	
			# preparing the table
			cur = self.db_connection.cursor()
			cur.execute("select count(*) from %s"%table_name)
			(number_of_rows,)=cur.fetchone()	
			#print 'number_of_rows=',	number_of_rows
			return number_of_rows



			
	

# =============================================================================
# interface
# =============================================================================
		
class utils():
	@staticmethod
	def getTableItor(table_name,start=0,end=None,db_name=None, user_name=None, password=None):		
		db = DBconnect(db_name=db_name, user_name=user_name, password=password);
		itor = db.tableItorChunks(table_name,start=start,end=end)
		return itor
	
	@staticmethod
	def tableLen(table_name,db_name=None, user_name=None, password=None):
		db = DBconnect(db_name=db_name, user_name=user_name, password=password);
		return db.tableLength(table_name)
	
		
def init_database_obj(db_name=None, user_name=None, password=None):
	db = DBconnect(db_name=db_name, user_name=user_name, password=password);
	return db

def test_table_itor(start=0,end=None,db_name=None, user_name=None, password=None):
	db = DBconnect(db_name=db_name, user_name=user_name, password=password);
	table_name='location'
	itor = db.tableItorChunks(table_name,start=start,end=end,jump=10000)
	return itor
	#for i in range(10):
	#	print itor.next()
		

