"""
	pair IBD calculator
	methods for IBD calculations
	given a pair of IDs (database IDs) calculate the IBD coeffs	
	
"""
import CAsearch
import gDBpedUI
import idcoefInterface

from numpy import *
import os
import sys
import traceback


	


class pairIBDcalcUI():
	"""
		obj = pairIBDcalcUI(ID1,ID2,tmp_dir)
	"""
	def __init__(self,ID1,ID2,tmp_dir,db_name=None, user_name=None, password=None):
		# database parameters
		self.db_name=db_name
		self.user_name=user_name
		self.password=password
		
		self.ID1=ID1 # original IDs from the DB
		self.ID2=ID2
		# tmp dir for idcoef runs
		self.tmp_dir=tmp_dir
		# results from ancestry search
		self.common_ids=[]
		self.ped_fixer=None
		
	def extract_ancestry(self,G1=7,G2=None):
		"""
		extract_ancestry(G)
			G = num of generations to go back
	
		extract_ancestry(G1,G2)
			uses different generations for each indv
		
		return bool = is there common ancestry
		"""
		if G2==None:
			G2=G1
		# init a connection to the DB and return a pedigree object
		db_ped = gDBpedUI.init_ped(self.user_name, self.password, self.db_name)	
		# init an ancestry-search object
		ac_obj=CAsearch.AncestrySearch(db_ped)
		# export a pedigree 
		ID1=self.ID1; ID2=self.ID2
		common_ids, ped_fixer = ac_obj.exportAncPedigree(ID1,G1,ID2,G2)
		self.common_ids=common_ids
		self.ped_fixer=ped_fixer
		# is there common anc?
		return self.isThereCommonAnc()
		
	def isThereCommonAnc(self):
		if self.common_ids==None:
			return False
		return len(self.common_ids)>0

	def get_idcoeff_input_filenames(self):
		tmp_dir=self.tmp_dir
		ped_filename =  os.path.join(tmp_dir,'test.%d.%d.ped'%(self.ID1,self.ID2))
		study_filename =  os.path.join(tmp_dir,'test.%d.%d.study'%(self.ID1,self.ID2))
		output_filename = os.path.join(tmp_dir,'test.%d.%d.output'%(self.ID1,self.ID2))
		return ped_filename, study_filename, output_filename	
	
	def prepare_idcoef_input_files(self):
		"""
			given the output from extract_ancestry(ID1,ID2,G), prepare the files for idcoef 
		"""
		common_ids=self.common_ids
		ped_fixer =self.ped_fixer
		ID1=self.ID1; ID2=self.ID2
						
		if len(common_ids)==0:
			# this should not happen
			print 'Error: no common qcestry in the lineage!'
			print '		this should not happen ...'
			raise ValueError
		
		ped_filename, study_filename, output_filename = self.get_idcoeff_input_filenames()	
			
		# write ped and study files
		# -------------------------
		# map the original IDs to integers and export the pedigree
		ped_lines, id_map = ped_fixer.exportPed() 
		# id_map is a hash table: 
		#	id_map.new2orig[int] -> original_id
		#	id_map.orig2new[original_id] -> int
		# write ped file
		fid=open(ped_filename,'w')
		for _line in ped_lines:
			fid.write('\t'.join([str(l) for l in _line])+'\n')
		fid.close()		
	
		# save study file with the new IDs 
		# i.e. save this list to file: [id_map.orig2new[ID1], id_map.orig2new[ID2]]		
		savetxt(study_filename,[id_map.orig2new[ID1], id_map.orig2new[ID2]],fmt='%d',newline='\t')	
		return ped_filename, study_filename, output_filename 
 
	def run_idcoef(self,ram):
		launcher = idcoefInterface.progUI()
		# run the prog
		pedfile,studyfile,outfile = self.get_idcoeff_input_filenames()
		call_result = launcher.run(pedfile,studyfile,outfile,ram)
		if call_result!=0:
			print '\nError running idcoefs program!\n';
			raise RuntimeError
		# read the output
		idcoeff_table, uniq_id_pairs = launcher.readOutput(outfile)
		return idcoeff_table[uniq_id_pairs[0]]

	def main(self,G=10,RAM=1000):
		"""
		main pipeline
		G = num of generation back
		RAM = ram to use
		"""
		# search for ancestry
		anc_found = self.extract_ancestry(G1=G)
		
		# is there is no ancestry return zero IBD
		if not(anc_found):
			print '<pairIBDcalcUI> no anc for pair: ', (self.ID1,self.ID2)
			return idcoefInterface.condensedIBDcoeffs(None)
		
		# prepare input files for idcoef
		self.prepare_idcoef_input_files()
		
		# run the prog and collect the output
		return self.run_idcoef(RAM)		
				
		
	def __del__(self):
		"""
			delete tmp files
		"""
		ped_filename, study_filename, output_filename = self.get_idcoeff_input_filenames()
		try:
			[os.remove(item) for item in [ped_filename, study_filename, output_filename]]
		except:
			#print 'Warning: some files were not deleted!'
			pass

# ----------------------------------------------------------------------------------------
# interface functions
# ----------------------------------------------------------------------------------------

def UIcaller(ID1,ID2,G=7,RAM=2000,TMP_DIR='.',db_name=None, user_name=None,  password=None):
	"""
		A caller that wraps the class to a simple function call (for map function calls)
	"""
	try:
		obj = pairIBDcalcUI(ID1, ID2, TMP_DIR, db_name=db_name, user_name=user_name, password=password)
		ibd = obj.main(G,RAM)				
	except: # some specific exception you want to wrap
		exp_value = sys.exc_info()[1]
		exp_type = sys.exc_info()[0]
		print 'Exception_type:',exp_type
		print 'Exception_value:',exp_value
		traceback.print_tb(sys.exc_info()[2],file=sys.stdout)		
		raise
	return ibd
	
def parameterWrapper(params):
	"""
		unpacking the parameters and calling the calculation
		this is a good interface for multiprocessing map function
	"""
	ID1,ID2,G,RAM,TMP_DIR,DATABASE,USER,PASSWORD=params
	return UIcaller(ID1,ID2,G,RAM,TMP_DIR,DATABASE,USER,PASSWORD)
	
	

