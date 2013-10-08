"""
pwIBDscript: This is a test script that uses the FamiLinx python API for the calculation of pairwise kinship scores.

Usage:
python pwIBDscript.py -i id1, id2
	This will calculate the pairwise Identity coefficients for the pair (id1, id2)
	For example:
	>> python pwIBDscript.py -i 1,2

	trying to run IBD for the pair: 1 ,  2
	pedigree file: ./tmp/test.1.2.ped
	sample file: ./tmp/test.1.2.study
	output file: ./tmp/test.1.2.output
	ramsize = 2000
	npeop = 4
	Successful completion.
	successfully calculated the kinship scores for this pair:
	kinship =  0.5
	the 9 identity coefficients are:  [ 0.    0.    0.    0.    0.    0.    0.25  0.5   0.25]

If this procedure fails it will initiate a test procedure that will document each step of a similar calculation.
If something goes wrong you will be able to see at what stage of the calculation the error occurred.

You can invoke this test procedure with this command:
	>> python pwIBDscript.py -t

		-----------------------------
		 Performing a test procedure 
		-----------------------------

	Trying to establish a connection to the database ...
	Successfully connected to the database

	Trying to perform this query: 'select Parent_id from relationship where Child_id=1' 
	This should return [4, 5]

	Performing the query ...
	The returned value:  [4, 5]

	Trying to run idcoefs program on the example files
	The example files reside here: pythonUI/idcoef/Idcoefs2.1.1/Example
	pedigree file: idcoef/Idcoefs2.1.1/Example/ex.pedigree
	sample file: idcoef/Idcoefs2.1.1/Example/ex.study
	output file: idcoef/Idcoefs2.1.1/Example/output.new.txt
	ramsize = 1000
	npeop = 168
	Successful completion.

	The idcoef program exited normally!

		------------------------------------------------------------------
		 The test procedure is now complete.
		 Any error or warnings should be documented for further assistance
		------------------------------------------------------------------


	
	
	
	
"""

import sys
import os
import getopt
import traceback

import familinxUI



# main helpers
# ------------
def do_test():
	"""
		Performing a test run by
			connecting to the database
			trying to run idcoef program
	"""
	print '\n\t-----------------------------'
	print '\t Performing a test procedure '
	print '\t-----------------------------'

	print '\nTrying to establish a connection to the database ...'
	try:
		conn=familinxUI.DBconnection.getDBconnection()
	except:
		print 'Could not connect to the database!'
		print 'Check that the database is configured correctly. See the database INSTALL file and the pythonUI/README file.'	
		sys.exit(2)
	if conn==None:
		print 'Could not connect to the database!'	
		print 'Make sure your username, password or database name are correctly configured'
		print 'For more help on this topic see the README and INSTALL files'
		sys.exit(2)
	print 'Successfully connected to the database'	
	print "\nTrying to perform this query: 'select Parent_id from relationship where Child_id=1' "
	print "This should return [4, 5]"
	try:
		cur=familinxUI.DBconnection.makeQuery('select Parent_id	from relationship where Child_id=1')
		res=[int(item[0]) for item in  cur.fetchall()]
		print '\nPerforming the query ...'
		print 'The returned value: ',res
	except:
		print '\n ... An error occurred!'
		print 'There seem to be a problem with the tables. \nSee the Installation file for more details on how to check database integrity.\n'
		traceback.print_tb(sys.exc_info()[2],file=sys.stdout)
		sys.exit(2)
	print '\nTrying to run idcoefs program on the example files'
	print 'The example files reside here: pythonUI/idcoef/Idcoefs2.1.1/Example'
	pedfile=os.path.join('idcoef/Idcoefs2.1.1/Example','ex.pedigree')
	studyfile=os.path.join('idcoef/Idcoefs2.1.1/Example','ex.study')
	outputfile=os.path.join('idcoef/Idcoefs2.1.1/Example','output.new.txt')
	res=familinxUI.pwIBDctor.idcoefInterface.globalrunner([pedfile,studyfile,outputfile])
	if res!=0:
		print 'There was a problem running the idcoef program!'
		print 'Some suggestions:'
		print '1) Make sure you have read/write permission to the idcoef directory'
		print '2) Try running the program from the command line'
		print '3) Make sure the relative path to the program is pythonUI/idcoef/Idcoefs2.1.1'
		exit(2)
	
	print '\nThe idcoef program exited normally!'
	print '\n\t------------------------------------------------------------------'
	print '\t The test procedure is now complete.' 
	print '\t Any error or warnings should be documented for further assistance' 
	print '\t------------------------------------------------------------------'
	

	
def run_ibd(id1,id2):
	"""
		running ibd calculation for (id1,id2)
	"""
	print 'trying to run IBD for the pair:',id1,', ', id2
	try:
		K, CID = familinxUI.IBDctorUI.calcPW(id1,id2)
		print 'successfully calculated the kinship scores for this pair:'
		print 'kinship = ',K
		print 'the 9 identity coefficients are: ',CID
		
	except: # some specific exception you want to wrap
		print '\n\n\nSome errors occurred during the program execution!'
		#print 'Here is the traceback ... \n' 
		exp_value = sys.exc_info()[1]
		exp_type = sys.exc_info()[0]
		#print 'Excpetion_type:',exp_type
		#print 'Excpetion_value:',exp_value
		#traceback.print_tb(sys.exc_info()[2],file=sys.stdout)
		print 'A test will now run that will try to pin-point the problem ...\n'
		do_test()
		sys.exit(2)
	
def parse_ids(arg):
	id1=id2=-1
	try: 
		id1,id2 = [int(item) for item in arg.split(',')]
	except:
		print 'Error parsing the input: ', arg
		print 'Expecting to get 2 integers.'
		print 'Try the -h option for help ...'
	return id1,id2
	
	
# -------------------------------------------------
# this is main
# -------------------------------------------------
	
def main(argv):
	"""
		parsing input parameters and calling a run function
	"""
	
	try:
		opts, args = getopt.getopt(argv,"hti:")
	except getopt.GetoptError:
		print __doc__
		sys.exit(2)
		
	if len(opts)==0:
		print __doc__
		sys.exit()
		
	# options
	for opt, arg in opts:
		if opt in ('-h',"--help"):
			print __doc__
			sys.exit()
			
		elif opt in ('-t', "--test"):
			do_test()
			sys.exit()
			
		elif opt in ('-i', "--ids"):
			id1,id2=parse_ids(arg)
			run_ibd(id1,id2)
		else:
			print 'else ...'
			print __doc__
			
if __name__ == "__main__":
	#	print __doc__
  	main(sys.argv[1:])			
