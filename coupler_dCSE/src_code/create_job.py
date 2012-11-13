#!/usr/bin/python
import sys
import subprocess

# Store number of command line arguments
nargs = len(sys.argv)

# default values
job_name   = 'default_jobname'
nproc_md   = 0
nproc_cfd  = 0
walltime   = '00:00:30'
queue      = ''
icib       = 'false'

# If wrong number of input arguments
if ( nargs < 3 ): 

	# If user asked for help
	if ( nargs == 2 and sys.argv[1] == 'help'):
		message  = 'create_job.py performs two operations: \n'
		message += '	1. Check simulation input files agree, and \n'
		message += '	2. Submit the simulation via PBS. \n'
		message += '\n'
		message += 'create_job.py requires the arguments nproc_md and \n'
		message += 'nproc_cfd. Other arguments are optional, but all \n'
		message += 'have the same syntax: \n'
		message += '	./create_job.py argname=argvalue argname=argvalue...\n'
		message += '\n'
		message += 'Available arguments: \n'
		message += '	"nproc_md"   or "m",\n'
		message += '	"nproc_cfd"  or "c",\n'
		message += '	"job_name"   or "j",\n'
		message += '	"walltime"   or "w",\n'
		message += '	"queue"      or "q",\n'
		message += '	"icib"       or "i".\n'
		message += '\n'
		message += 'Example:\n'
		message += '	./create_job.py m=4 c=1 j=debug_job w=00:01:00 q=pqtzaki i=false\n\n'
		print(message)
		sys.exit(1)

	else:
		message = ('create_job.py requires at least two input arguments \n' +
                   'that should be nproc_md and nproc_cfd. Type \n'         +
		           '"./create_job.py help" for more information on \n'      +
		           'available arguments and their syntax.')
		print(message)
		sys.exit(1)

# Run through and store arguments
for arg in sys.argv[1:nargs]:

	argname = arg.split('=')[0].strip()

	if   ( argname == 'nproc_md'  or argname == 'm' ):
		nproc_md = int(arg.split('=')[1])

	elif ( argname == 'nproc_cfd' or argname == 'c' ):
		nproc_cfd = int(arg.split('=')[1])

	elif ( argname == 'job_name'  or argname == 'j' ):
		job_name = arg.split('=')[1].strip()

	elif ( argname == 'icib'      or argname == 'i' ):
		icib = arg.split('=')[1].strip()

	elif ( argname == 'walltime'  or argname == 'w' ):
		walltime = arg.split('=')[1].strip()

	elif ( argname == 'queue'     or argname == 'q' ):
		queue = arg.split('=')[1].strip()
	
	else:
		message = ('Unrecognised argument "' + argname + '" to \n' +
		           'create_job.py. Aborting.')
		print(message)
		sys.exit(1)


# Check simulation is going to be OK
ierr = subprocess.call( ["./checksims.py", str(nproc_md), str(nproc_cfd)] )
if ( ierr == 0 ):
	print('Simulation check successful. Initialising job creation...\n')
else:
	print('Simulation check failed. Aborting.\n')
	sys.exit(1)

# Create header information
nproc = nproc_md + nproc_cfd
header  = '#!/bin/bash\n' 
header += '#PBS -l select='   + str(nproc) + ':icib=' + icib + '\n'
header += '#PBS -l walltime=' + walltime   + '\n'
header += '#PBS -N '          + job_name   + '\n'
if ( queue != '' ):
	header += '#PBS -q '      + queue      + '\n'
header += 'module load intel-suite\n'
header += 'module load mpi\n'

# Create script
script  = 'outfile=$PBC_JOBNAME_$PBS_JOBID\n\n'
script += 'cd $PBS_O_WORKDIR\n\n'
script += '#Make results directory if not present\n'
script += 'mkdir -p ./couette_data/results\n\n'
script += '#Clean it up if it is present\n'
script += 'rm ./couette_data/data ./couette_data/ucvcwc.dble.* ./couette_data/uuvvww.dble.* ./couette_data/conold.dble.* ./couette_data/pressure_ph.* ./couette_data/pres_p* ./couette_data/archive* ./couette_data/report ./couette_data/SubDom_dble*\n\n'
script += 'echo "  0" > data\n'
script += 'mv ./data ./couette_data/\n'
script += 'cp ./../../CFD_dCSE/src_code/main_code/input ./couette_data/\n'
script += 'cp ./../../CFD_dCSE/src_code/main_code/parallel_couette.exe ./couette_data/\n'
script += 'cp ./../../CFD_dCSE/src_code/setup/ucvcwc.dble.000000 ./couette_data/\n'
script += 'cp ./../../CFD_dCSE/src_code/setup/uuvvww.dble.000000 ./couette_data/\n'
script += 'cp ./../../CFD_dCSE/src_code/setup/conold.dble.000000 ./couette_data/\n'
script += 'cp ./../../CFD_dCSE/src_code/setup/pressure_ph.000000 ./couette_data/\n'
script += 'cp ./../../CFD_dCSE/src_code/setup/grid.data ./couette_data/\n'
script += 'cp ./../../CFD_dCSE/src_code/setup/archive ./couette_data/\n'
script += 'cp ./../../CFD_dCSE/src_code/setup/report ./couette_data/\n'
script += 'cp ./couette_data/archive ./couette_data/archive.000000\n'
script += 'cp ./couette_data/report ./couette_data/report.000000\n'
script += 'date\n\n'
script += 'mpiexec heterostart ' + str(nproc_cfd) + ' ./couette_data/parallel_couette.exe ' + str(nproc_md) + ' ./../../MD_dCSE/src_code/md.exe -i ./md_data/MD.in >> $outfile\n\n'
script += 'date\n'

# Ask user if they wish to proceed with submission
print('PBS submission header as follows:\n')
print(header)
proceed = raw_input('\nIs this the submission header you wish to proceed with? (y/n): ')

if ( proceed == 'y' or proceed == 'Y' ):
	# Submit job
	job_file = 'job_submission'
	fobj = open(job_file,'w')
	fobj.write(header)
	fobj.write(script)	
	fobj.close()
	ierr = subprocess.call( ["qsub", job_file] )
	sys.exit(ierr)
else:
	# Abort
	print('Aborting job submission')
	sys.exit(1)
