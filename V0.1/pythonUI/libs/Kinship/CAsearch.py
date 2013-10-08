"""
this is a general module for conducting searches in pedigrees
it requires an implementation of an interface class that walks the pedigree by retrieving children or parents for each entry (ID)
the interface class is defined in pedInterface.py

"""

from numpy import *

DEBUG=False
DEBUG1=False

# helpers
class Ancestrors():
	"""  
		a class to hold a list of ancestors and their distance (generations) from a given individual 
		ID should be INT
	"""
	def __init__(self,ID):
		self.ID=ID # of the individual
		self.generations=[]
		self.ancestor_ids=[]
	
	def __str__(self):
		S='Ancestry for ID:%d\n'%(self.ID,)
		for g in range(self.getMaxGeneration(),-1,-1):
			S+='%d\t%s\n'%(g,str(self.getUniqAncestors(g)[0]))
		return S	
			
	def getMaxGeneration(self):
		return max(self.generations)
		
	def addAncestors(self,list_of_ancestors,generation):
		""" given the generation (a num) and given a list of ancestors add them to the object data 
			note that if there are loops in the pedigree then some ancestors will appear twice
		"""
		self.ancestor_ids.extend(list_of_ancestors)
		self.generations.extend([generation]*len(list_of_ancestors))
	
	def getAncestors(self,generation=None):
		""" return a list of Ancestors for a given generation, 
			if generation is not given then return all 
			note that if there are loops in the pedigree then some ancestors will appear twice
			Do a unique() on the returned value for unique IDs
			"""
		if generation == None:
			return array(self.ancestor_ids)
		# else	
		g_arr=array(self.generations)
		a_arr=array(self.ancestor_ids)
		return a_arr[g_arr==generation]
	
	def getUniqAncestors(self,generation=None):
		"""
		same as getAncestors() with a unique set of ancestors IDs
		"""
		non_unq_anc = self.getAncestors(generation)
		(unq_anc,loc)=unique(non_unq_anc,return_index=True)
		gen = array(self.generations)[loc]
		# sort by generation
		I=argsort(gen)				
		return unq_anc[I], gen[I]
	
	def getAtorGeneration(self,AID):
		"""  return the generation found for the given Ancestor ID """
		I = self.ancestor_ids.index(AID)
		return self.generations[I]
			
# main classes
class AncestrySearch():
	def __init__(self,pedI_obj):
		self.ped=pedI_obj # the interface class for transversing the pedigree with the data loaded already
	
	def getParents(self,ID):
		return self.ped.getParents(ID)
		
		
	def getAncestors(self,ID,g_th):
		"""
			return all ancestors of ID up to generation < g_th
		"""
		anc = Ancestrors(ID)
		anc.addAncestors([ID],0) # add ID as generation 0
		for g in range(1,g_th+1):
			for PID in anc.getAncestors(g-1):
				try:
					list_of_ancestors=self.ped.getParents(PID)
				except:
					print ' '
					print 'ERROR: in ped.getParents(%d)'%PID
					print 'This id was not found! you might need to send string instead of numbers ...'
					print ' '
					raise ValueError
				anc.addAncestors(list_of_ancestors,g)
		return anc
	
	def getAncIntersect(self,ID1,G1,ID2,G2):
		"""
			return the intersection of ancestors for ID1 and ID2
			limiting the search by G1 and G2 generations
		"""
		anc_obj_1 = self.getAncestors(ID1,G1)
		anc_obj_2 = self.getAncestors(ID2,G2)
		# intersect the unique anc ids
		uniq1,gen1 = anc_obj_1.getUniqAncestors()
		uniq2,gen2 = anc_obj_2.getUniqAncestors()
		common_ids = intersect1d(uniq1,uniq2)
		# get the corresponding generation for each founder
		gen1 = [anc_obj_1.getAtorGeneration(AID) for AID in common_ids]
		gen2 = [anc_obj_2.getAtorGeneration(AID) for AID in common_ids]
		return common_ids,gen1,gen2
	
	# ----------------------
	# main interface method
	# ----------------------	
	def exportAncPedigree(self,ID1,G1,ID2,G2):
		"""	
			common_ids, ped_fixer = exportAncPedigree(self,ID1,G1,ID2,G2)
			
			get the ancestry of IDi up to Gi 
			calc the ancestry intersection 
		
			build a ped format array
			ID FATHER MOTHER	
			
			the output is sorted so that an individual parents appear before him
			the output should be a valid pedigree for computing IBD (sorted with no missing parents)
			
			return:
				common_ids = 
					None if no common ancestry was found
					the original ids for common ancestors 
				the pedigree as a pedInterface::pedFixer Object 
				
		"""
		# get unique ancestors for each ID
		anc_obj_1 = self.getAncestors(ID1,G1)
		anc_obj_2 = self.getAncestors(ID2,G2)
		# intersect the unique anc ids
		uniq1,gen1 = anc_obj_1.getUniqAncestors()
		uniq2,gen2 = anc_obj_2.getUniqAncestors()			
		common_ids = intersect1d(uniq1,uniq2)
		
		if DEBUG:
			#print '<CAsearch::exportAncPedigree> common_ids = ',common_ids
			# get the corresponding generation for each founder
			gen1 = [anc_obj_1.getAtorGeneration(AID) for AID in common_ids]
			gen2 = [anc_obj_2.getAtorGeneration(AID) for AID in common_ids]
			tot_g=array(gen1)+array(gen2)
			I=argsort(tot_g)
			print '-----------------------------------------'
			print '<CAsearch::exportAncPedigree> common_ids:'	
			print 'AID         \tGEN1\tGEN2\tTOTG'
			for i in I:	
				print common_ids[i],'\t',gen1[i],'\t',gen2[i],'\t',gen2[i]+gen1[i]
			print '-----------------------------------------'
			
		# if there are no common ancestors then return null
		if len(common_ids)==0:
			return None, None
		# build the pedigree
		# make a unique set of IDs
		ped_ids = unique(append(uniq1,uniq2))

		# debug
		if DEBUG1:
			print 'Debug info <CAsearch::exportAncPedigree>'
			print '<CAsearch::exportAncPedigree> ped_ids=',array(ped_ids)
			#print 'uniq1=',uniq1
			#print 'uniq2=',uniq2
			#print 'common_ids=',common_ids
			
			 
		# export to linkage format
		ped_fixer = self.ped.exportLinkage(ped_ids) # TBD add filename
		return common_ids, ped_fixer
					

