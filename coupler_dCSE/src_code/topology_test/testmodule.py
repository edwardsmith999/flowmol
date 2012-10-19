#! /usr/bin/env python
"""
	TOPOLOGY TEST MODULE

	Contents:
		
		prevjobID (integer)   -   Parameter incremented every time an object of
		                          type RunClass is created, then saved to
		                          that object's self.jobID
		
		RunClass  (class)     -   Everything required for a job object.

"""
prevjobID = 0

class RunClass:

	# Constructor
	def __init__(self, npx_md,  npy_md,  npz_md, 
	                   npx_cfd, npy_cfd, npz_cfd):
		from re import sub 
		from os import mkdir

		global prevjobID

		# Create new ID for every job initialised
		self.jobID = prevjobID + 1
		prevjobID  = self.jobID

		# Read and store input file in "lines"
		with open('input','r') as source:
			lines = source.readlines()

		# Clean input file and rewrite with new procs
		with open('input','w') as output:
			# Store procesors in an array for loop below
			procs = [npx_md, npy_md, npz_md, npx_cfd, npy_cfd, npz_cfd]
			for line in range(len(lines)):
				if (line < 6):
					output.write(sub('\d{1,}', str(procs[line]), 
					                    lines[line]))
				else:
					output.write(lines[line])

		# Store nprocs
		self.nproc = ( npx_md  * npy_md  * npz_md 
		             + npx_cfd * npy_cfd * npz_cfd )

	# Print statement
	def __str__(self):

		string =  'jobID: ' + str(self.jobID) + '\n'
		string += '-----------------------------------------\n'

		# Add input file contents
		with open( 'input', 'r' ) as source:
			lines = source.readlines()
		for line in lines:
			string += line

		string += '-----------------------------------------'

		return( string )

	
	def execute(self):
		import os	
	
		if (os.path.exists('./logs')):

			logfile = './logs/' + str(self.jobID) + '_log'
			errfile = './logs/' + str(self.jobID) + '_errlog'

			cmd = ('mpiexec -n ' + str(self.nproc) + ' ./a.out > ' + logfile 
				   + ' 2> ' + errfile)

			print(cmd)
			os.system(cmd)

		else:

			raise


	def concatenate(self):
		from os import system 


		cmd = 'cat fort.1* > info_realms && rm fort.1*'
		system(cmd)

		cmd = 'cat fort.2* > info_MD_recv && rm fort.2*'
		system(cmd)

		cmd = 'cat fort.3* > info_graph && rm fort.3*'
		system(cmd)

		cmd = 'cat fort.4* > info_MD_send && rm fort.4*'
		system(cmd)

		cmd = 'cat fort.5* > info_CFD_recv && rm fort.5*'
		system(cmd)

		cmd = 'cat fort.6* > info_maps && rm fort.6*'
		system(cmd)

		cmd = 'cat fort.7* > info_scatter_md && rm fort.7*'
		system(cmd)

		cmd = 'cat fort.8* > info_gather_cfd && rm fort.8*'
		system(cmd)

		cmd = 'cat fort.9* > info_CFD_send && rm fort.9*'
		system(cmd)

		
	def checkvalues(self):
		import check

		all_success = True


		success = check.vals('info_scatter_md')
		if success is not 0:
			print('Error in gather values')
			all_success = False
		success = check.vals('info_gather_cfd')
		if success is not 0:
			print('Error in scatter values')
			all_success = False
		success = check.vals('info_CFD_send')
		if success is not 0:
			print('Error in CFD send values')
			all_success = False
		success = check.vals('info_MD_send')
		if success is not 0:
			print('Error in MD send values')
			all_success = False
		success = check.vals('info_MD_recv')
		if success is not 0:
			print('Error in MD recv values')
			all_success = False
		success = check.vals('info_CFD_recv')
		if success is not 0:
			print('Error in CFD recv values')
			all_success = False

	
		if all_success:
			print('All value checks successful')
			print('=========================================')
