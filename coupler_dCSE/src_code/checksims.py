#!/usr/bin/python
import sys

# Define custom grep wrapper
def grep(args): 
	import subprocess
	cmd = ["grep"] + args
	obj = subprocess.Popen( cmd, stdout=subprocess.PIPE )
	out = obj.communicate()[0]
	return out 

def diff(args):
	import subprocess
	cmd = ["diff"] + args
	obj = subprocess.Popen( cmd, stdout=subprocess.PIPE )
	out = obj.communicate()[0]
	return out 

md_input_file = "./md_data/MD.in"
cfd_dir = "../../CFD_dCSE/src_code/"
cfd_param_file = cfd_dir + "main_code/param.inc"

# Info about processors
grep_out = grep(["-A3","PROCESSORS",md_input_file])
npx_md   = int( grep_out.split()[1] )
npy_md   = int( grep_out.split()[2] )
npz_md   = int( grep_out.split()[3] )
nproc_md = npx_md * npy_md * npz_md

message  = ('===================================================\n' +
            'MD processors:\n'  + 
            'npx: ' + str(npx_md) + '\n' +
            'npy: ' + str(npy_md) + '\n' +
            'npz: ' + str(npz_md) + '\n' )

grep_out = grep( ["parameter (npx", cfd_param_file ] )
npx_cfd  = int( grep_out.split('npx=')[1].split(',')[0] )
npy_cfd  = int( grep_out.split('npy=')[1].split(',')[0] )
npz_cfd  = int( grep_out.split('npz=')[1].split(',')[0] )
nproc_cfd = npx_cfd * npy_cfd * npz_cfd

message += ('==================================================\n' +
            'CFD processors:\n'  + 
            'npx: ' + str(npx_cfd) + '\n' +
            'npy: ' + str(npy_cfd) + '\n' +
            'npz: ' + str(npz_cfd) + '\n' + 
            '==================================================')

print(message)

# If not one argument to ./checksims.py then exit with error
# (sys.argv counts "./checksims.py" as an argument, hence the 2)
if ( len(sys.argv) != 2 ): 

	message = ('checksims.py requires a single input argument that \n' +
	           'should be the total number of MD processors.')
	print(message)
	sys.exit(1)


# Check user has specified correct number of MD procs
nproc_md_IN = int(sys.argv[1])
if ( nproc_md_IN != nproc_md ):
	message = ('User specified nproc_md differs to that in the \n' +
	           'input file. Aborting.')
	print(message)
	sys.exit(1)	

# Grid information
grid_infile = cfd_dir + "grid_generation/input.file"
message  = 'CFD grid info:\n'
message += '==================================================\n'

floatparams = ["Lx","Ly"]
for param in floatparams:
	grep_out = grep([param,grid_infile])
	vars()[param] = float(grep_out.split()[0])
	message += param + ': ' + str(vars()[param]) + '\n'

intparams = ["ngx","ngy"]
for param in intparams:
	grep_out = grep([param,grid_infile])
	vars()[param] = int(grep_out.split()[0])
	message += param + ': ' + str(vars()[param]) + '\n'

message += '=================================================='
print(message)

# Setup information
grid_infile = cfd_dir + "setup/input.setup"
message  = 'CFD setup info:\n'
message += '==================================================\n'

floatparams = ["xL","yL"]
for param in floatparams:
	grep_out = grep([param,grid_infile])
	vars()[param] = float(grep_out.split()[0])
	message += param + ': ' + str(vars()[param]) + '\n'

intparams = ["nix","niy","niz"]
for param in intparams:
	grep_out = grep([param,grid_infile])
	vars()[param] = int(grep_out.split()[0])
	message += param + ': ' + str(vars()[param]) + '\n'

message += '=================================================='
print(message)

# Check setup and grid info are the same 
if ( xL != Lx or yL != Ly or nix != ngx or niy != ngy ):
	print('Grid and setup input files do not match. Aborting.')
	sys.exit(1)
else:
	print('Grid and setup input files match. Proceeding...')

# Diff the data files
gengrid   = cfd_dir + '/grid_generation/grid.data'
setupgrid = cfd_dir + '/setup/grid.data'
if (diff([gengrid,setupgrid]) != ""):
	print('grid.data files differ. Aborting.')
	sys.exit(1)
else:
	print('grid.data files are identical. Proceeding...')
	sys.exit(1)

sys.exit(0)
