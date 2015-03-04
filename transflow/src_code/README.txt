!TJ: 20/03/2012

DNS_grid_generation_Couette:
	- contains grid generation routines

DNS_main_code_Couette:
	- contains source code for Couette DNS
	- compile on cx1 using: "make -f makefile.planes_fftz_fftx all"

DNS_setup_Couette:
	- contains setup routines for initial field
	- compile on cx1 using: "make simple"

DNS_TJ_data_Couette:
	- contains all the required files to submit a job on cx1
	- submit using the suitably ammended run_script.sh

!---------------------------------------------------------------
