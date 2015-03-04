#!/usr/bin/python
import sys

# Define custom grep wrapper, "out" is a string
def grep(args): 
	import subprocess
	cmd = ["grep"] + args
	obj = subprocess.Popen( cmd, stdout=subprocess.PIPE )
	out = obj.communicate()[0]
	return out 

# Define custom diff wrapper, "out" is a string
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
if ( len(sys.argv) != 3 ): 

	message = ('checksims.py requires two input arguments that \n' +
	           'should be (in order) the total number of MD and \n' +
	           'CFD processors.')
	print(message)
	sys.exit(1)


# Check user has specified correct number of MD and CFD procs
nproc_md_IN = int(sys.argv[1])
nproc_cfd_IN = int(sys.argv[2])
if ( nproc_md_IN != nproc_md ):
	message = ('User specified nproc_md differs to that in the \n' +
	           'input file. Aborting.')
	print(message)
	sys.exit(1)	
if ( nproc_cfd_IN != nproc_cfd ):
	message = ('User specified nproc_cfd differs to that in the \n' +
	           'input file. Aborting.')
	print(message)
	sys.exit(1)	

# Grid information
grid_infile = cfd_dir + "grid_generation/input.file"
xL_grid  = float( grep( ['Lx',  grid_infile] ).split()[0] )
yL_grid  = float( grep( ['Ly',  grid_infile] ).split()[0] )
ngx_grid =   int( grep( ['ngx', grid_infile] ).split()[0] )
ngy_grid =   int( grep( ['ngy', grid_infile] ).split()[0] )

message  = 'CFD grid info:\n'
message += '==================================================\n'
message += ('xL: '  + str(xL_grid)  + '\n' +
            'yL: '  + str(yL_grid)  + '\n' +
            'ngx: ' + str(ngx_grid) + '\n' +
            'ngy: ' + str(ngy_grid) + '\n' )

message += '=================================================='
print(message)

# Setup information
setup_infile = cfd_dir + "setup/input.setup"
xL_setup  = float( grep( ['xL',  setup_infile] ).split()[0] )
yL_setup  = float( grep( ['yL',  setup_infile] ).split()[0] )
ngx_setup =   int( grep( ['nix', setup_infile] ).split()[0] )
ngy_setup =   int( grep( ['niy', setup_infile] ).split()[0] )
ngz_setup =   int( grep( ['niz', setup_infile] ).split()[0] )

message   = 'CFD setup info:\n'
message  += '==================================================\n'
message  += ('xL: '  + str(xL_setup)  + '\n' +
             'yL: '  + str(yL_setup)  + '\n' +
             'ngx: ' + str(ngx_setup) + '\n' +
             'ngy: ' + str(ngy_setup) + '\n' +
             'ngz: ' + str(ngz_setup) + '\n' )
message += '=================================================='
print(message)

# Check setup and grid info are the same 
if ( xL_setup  != xL_grid  or yL_setup  != yL_grid or 
     ngx_setup != ngx_grid or ngy_setup != ngy_grid ):
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

sys.exit(0)
