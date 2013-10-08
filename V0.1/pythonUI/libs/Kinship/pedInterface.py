"""
	an abstract interface class which is required to be implemented by CAsearch.py classes
	
"""
import sys
DEBUG=False

class Struct():
	pass

def printDebug(*args):
	if DEBUG:
		print '# <pedInterface> ',
		for item in args:
			print item,
		print ''


class PedI():
	# -----------------------------------------
	# abstract methods
	
	def getParents(self,ID):
		"""
			return a list of parents IDs 
			return [] if no parents
		"""
		pass
	
	def getDescen(self,ID):	
		"""
			return a list of descendants IDs 
			return [] is no descendants
		"""
		pass
	
	def getIDs(self):
		"""
			return a list of IDs for all entries in the pedigree
		"""
		pass
	
	# -----------------------------------------
	# general methods
	
	def exportLinkage(self,IDs=None,fout=None,SORT_PED=True):
		"""
			given a list of IDs, export a sub pedigree to a file with linkage format
			if fout == None -> write to stdout
				note, there is no check that the pedigree is connected or valid
				we just print out the given IDs with their parents
			SORT_PED: bool, should the pedigree be fixed and sorted prior to writing
							i.e. ids of parents should appear prior to their children
		""" 
		pass




# Helper
class pedFixer():
	"""
		this is a utility class that fixes and sorts pedigrees
		the input should be a list of ids in the following format:
		id father mother
		
		all ids should be INTs ...
		
	"""		
	def __init__(self,ped_name=None):
		self.par_table=dict() # parent hash table
		self.dec_table=dict() # desc hash table
		self.ids=[]
		self.ped_name=ped_name

	def mapIDs(self):
		"""
			Map IDs (strings) to ints and build a pedigree with integers
			Not the most efficient but at least readable to some of us. 
		"""
		# mapping IDs to integers
		#	id_map.new2orig[int] -> original_id
		#	id_map.orig2new[original_id] -> int
		id_map=Struct()		
		id_map.new2orig=dict()
		id_map.orig2new=dict()
		# add the unknown indv to the map
		id_map.new2orig[0]=0
		id_map.orig2new[0]=0
		# start iteration
		intID=1 # increment int as new IDs
		for ID in self.ids:
			if ID in id_map.orig2new:
				print '<pedFixer::mapIDs> warning, found duplicates for ID:',ID
			else:
				id_map.orig2new[ID]=intID
				id_map.new2orig[intID]=ID
				intID+=1
		return id_map
			
	def exportPed(self):
		"""
			ped_lines, id_map = exportPed()
	
			export pedigree as a list(list(integers))
			ped_lines[i]=[id father mothers]
			id_map = mapping IDs (strings) to integers using self.mapIDs
	
		"""
		id_map = self.mapIDs()
		# get the pedigree as a list
		ped_lines_str = self.asList()
		ped_lines_int=[]
		for line_str in ped_lines_str:
			ped_lines_int.append([id_map.orig2new[ID] for ID in line_str])		
		return ped_lines_int, id_map
		
	def addEntry(self,ID,F,M):
		
		self.par_table[ID]=[F,M]
		self.ids=self.par_table.keys()
		
		# if we did not see the children of ID -> init with [] (no children at the moment)
		if not self.dec_table.has_key(ID):
			self.dec_table[ID]=[]
		# add the current ID to the parents as  a child
		if self.dec_table.has_key(F):
			self.dec_table[F].append(ID)
		else:
			self.dec_table[F]=[ID]
		if self.dec_table.has_key(M):
			self.dec_table[M].append(ID)
		else:
			self.dec_table[M]=[ID]

	def getIDs(self):
		"""
			return all ids in the pedigree that have an entry
		"""
		return self.ids
		
	def fixMissingParents(self):
		"""
			if some entries contain parents that do not exist in the pedigree
			-> add these as founders
		"""	
		all_IDs = self.getIDs()
		for ID in all_IDs:
			# are the parents in the list of IDs as well?
			PARENTS=self.par_table[ID]
			for P in PARENTS:
				if P == 0:
					continue
				if P not in all_IDs:
					self.addEntry(P,0,0)
	
	def sortPed(self):
		"""
			sort the pedigree so that parents are before their children
			the sorted pedigree is manifested in the variable self.ids
			(ad hoc implementation and not a pretty sight - maybe could do better with sorting algorithms)
			
			Algorithm
			find all founders -> head of the list: sortedIDs
			for each f in founder (make a list of dec)
				dlist<-Enqueue(f.children()) (dlist is a FIFO queue - 1st in 1st out)
			while dlist not empty
				id<-Dequeue(dlist)
				if both parents are in sortedIDs
					sortedIDs <- id
					dlist<-Enqueue(id.children())
				else
					dlist<-Enqueue(id) // put at the end of the queue
		"""
		# fix the pedigree for founders
		self.fixMissingParents()
		# find the founders		
		Fids=self.findFounders()
		# make a list of descendants 
		dlist=[]
		# this is the sorted ids - starting with the founders
		sorted_ids = list(Fids)
		
		for ID in Fids:
			#dlist.extend(self.dec_table[ID])
			for Did in self.dec_table[ID]:
				if Did not in dlist:
					dlist.append(Did)
		# debug
		printDebug('Fids=',Fids)
		printDebug('dlist[Fids]=',dlist)
			
		counter=0; ped_size = len(self.ids)
		while len(dlist)>0:
			counter+=1
			if counter>ped_size**2:			
				print '\nError in <pedFixer::sortPed> '
				print 'sorting did not converge for pedigree', self.ped_name,'\n'
				raise AssertionError	

			cur_id = dlist.pop(0)
			# get the parents
			cur_id_f, cur_id_m = self.par_table[cur_id]
			# are both parents in the sorted list so far?
			# debug
			printDebug('-'*10)
			printDebug( 'begin counter: ',counter)
			printDebug('cur_id=',cur_id)
			printDebug('cur_id_f=',cur_id_f,' ; cur_id_m',cur_id_m)
			printDebug('sorted_ids=',sorted_ids )
			if (cur_id_f in sorted_ids) and (cur_id_m in sorted_ids):
				printDebug( 'adding this id to sorted list')
				# add the current id and extract it's children
				sorted_ids.append(cur_id)
				for Did in self.dec_table[cur_id]:
					if Did not in dlist:
						dlist.append(Did)
				# else, push this id to the end of the line
			else:
				printDebug( 'Enqueuing back ',cur_id )
				dlist.append(cur_id)

			printDebug('dlist = ',dlist)
		# update the ids				
		self.ids=sorted_ids		
			
	def sortPed_debug(self):
		"""
			sort the pedigree so that parents are before their children
			the sorted pedigree is manifested in the variable self.ids
			
			Algorithm
			find all founders -> head of the list: sortedIDs
			for each f in founder (make a list of dec)
				dlist<-Enqueue(f.children()) (dlist is a queue FIFO - 1st in 1st out)
			while dlist not empty
				id<-Dequeue(dlist)
				if both parents are in sortedIDs
					sortedIDs <- id
					dlist<-Enqueue(id.children())
				else
					dlist<-Enqueue(id) // put at the end of the queue
		"""
		# fix the pedigree for founders
		self.fixMissingParents()
		# find the founders		
		Fids=self.findFounders()
		# make a list of descendants 
		dlist=[]
		# this is the sorted ids - starting with the founders
		sorted_ids = list(Fids)
		
		for ID in Fids:
			print 'founder ID:',ID,' children are:', self.dec_table[ID]
			
			#dlist.extend(self.dec_table[ID])
			for Did in self.dec_table[ID]:
				if Did not in dlist:
					dlist.append(Did)
		print 'dlist=',dlist
		
					
	def findFounders(self):
		"""
			return a list of founders
		"""
		Fids=[]
		all_IDs = self.getIDs()
		for ID in all_IDs:
			# are the parents in the list of IDs as well?
			PARENTS=self.par_table[ID]
			if PARENTS==[0,0]:
				Fids.append(ID)
		return Fids
		
	def asList(self):
		""" export as a list """
		output=[]
		for ID in self.ids:
			entry = [ID] ; entry.extend(self.par_table[ID])
			output.append(entry)
		return output

# -----------------------------------
class pedFixer_random_sex(pedFixer):
	"""
		this class manages the sex column as well
	"""
	def __init__(self,ped_name=None):
		#super(pedFixer_random_sex,self).__init__(ped_name)
		pedFixer.__init__(self,ped_name)
		self.sex_table=dict() # sex table for each ID (1 male, 2 female, 0 unknown)
		
	def addSex(self,ID,SEX=0,FORCE=False):
		"""
			add a sex indicator (1 male, 2 female, 0 unknown)
			algorithm
				if FORCE:
					add sex
				if ID has sex already
					return the existing sex
				else
					add sex
					
		here are the possible states (assuming FORCE is False)
		SEX		sex_table[ID]	action
		0		0				return
		0		1/2				return
		1/2		0				sex_table[ID]=SEX
		1/2		1/2, 2/1		return	(we do not want to change what was written)
		
		0/1/2	None			sex_table[ID]=SEX
		
		"""
		if FORCE:
			self.sex_table[ID]=SEX
			return SEX
			
		# if there's an entry for this ID 
		if ID in self.sex_table:
			if SEX==0: # there is nothing to update in this case
				return self.sex_table[ID]
			if self.sex_table[ID]==0: # an update might be needed
				self.sex_table[ID]=SEX
				return SEX
			else: # sex_table[ID]!=0 -> return it
				return self.sex_table[ID]
		# there isn't an entry for this ID
		self.sex_table[ID]=SEX
		return SEX
		
	def GetSex(self,ID):
		if ID in self.sex_table:
			return self.sex_table[ID]
		return 0
	
	def checkUnionSex(self,P1,P2):
		""" check that a pair sex is consistent with male, female and 
			return the corresponding sexes  
			return an indicator suggesting that they are consistent with the sex table
		"""
		if P1==0 and P2==0:
			return 1,2,True
			
		consistent = True
		S1=0; S2=0 # the sexes are unknown at the moment
		S1 = self.GetSex(P1); S2 = self.GetSex(P2)
		if S1!=0 and S1==S2:
			consistent=False
			print '\nWarning: <pedInterface::checkUnionSex> inconsistent sexes for parents:'
			print 'parent:sex'
			print P1,':',S1
			print P2,':',S2,'\n'
		return S1, S2, consistent
		
	def determineUnionSex(self,P1,P2):
		"""
			F,M = determineUnionSex(P1,P2)
			
			determine sexes (arbitrary if needed) and return
			Ids for Father and Mother
			
			these are the options
			S1	S2	action
			---------------------
			0	0	arbitrary
			0	X	S1=flip(X)
			X	0	S2=flip(X)
			X	X	S1=flip(X) 	not consistent
			X	Y	---			else do nothing ...
			
		"""
			
		# helper
		def flipSex(Si):
			""" if Si is male return female and vice-versa """
			if Si==1:
				return 2
			return 1
			
		# code	
		S1,S2,consistent = self.checkUnionSex(P1,P2)
		if S1==0 and S2==0:
			#both unknown -> determine arbitrary
			S1=1; S2=2
		elif S1==0 and S2>0:
			# one unknown - flip the other
			S1=flipSex(S2)
			
		elif S2==0 and S1>0:
			S2=flipSex(S1)
		elif not(consistent):
			S1=flipSex(S2)
		# update the sex table
		self.addSex(P1,SEX=S1,FORCE=True)
		self.addSex(P2,SEX=S2,FORCE=True)			

		# return Father,Mother	
		if S1==1:
			return P1,P2
		return P2,P1
			
			
	def addEntry(self,ID,P1,P2):
		"""	
			decide arbitrary who is the male and who is the female
			add a dictionary of sex status which is consistent with the order of the parents
			(i.e. child father mother)
		
		"""
		self.addSex(ID) # add as unknown or leave current value unchanged
		F,M = self.determineUnionSex(P1,P2) # also adds the sex entries
		#super(pedFixer_random_sex,self).addEntry(ID,F,M)
		pedFixer.addEntry(self,ID,F,M)

	def asList(self):
		""" export as a list """
		output=[]
		for ID in self.ids:
			entry = [ID] ; entry.extend(self.par_table[ID])
			entry.append(self.sex_table[ID])
			#entry.append('1') # affected status
			output.append(entry)
		return output		

	def exportPed(self):
		"""
		TBD fix this module 
		need to change IDs from strings to integers!!!
		
		
			ped_lines, id_map = exportPed()
		
			--- added sex field
			export pedigree as a list(list(integers))
			ped_lines[i]=[id father mothers]
			id_map = mapping IDs (strings) to integers using self.mapIDs
		
		"""
		
		id_map = self.mapIDs()
		# get the pedigree as a list
		ped_lines_aslist = self.asList()
		ped_lines_int=[]		
		for line in ped_lines_aslist:
			ID=line[0]; father=line[1]; mother=line[2]; sex=line[3] 
			entry = [id_map.orig2new[ID],id_map.orig2new[father],id_map.orig2new[mother],int(sex)]
			ped_lines_int.append(entry)
		
		return ped_lines_int, id_map
		

