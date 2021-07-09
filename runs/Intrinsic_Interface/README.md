These are the run files to reproduce the results from the paper:

Edward Smith (edwardsmith.co.uk, ORCID 0000-0002-7434-5912)
2021
The Importance of Reference Frame for Pressure at the Liquid-Vapour Interface
Accepted for publication in Molecular Simulation
https://arxiv.org/abs/2107.00499

and more generally, to study a full pressure profile for the intrinsic interface.

The latest version of the Flowmol code will likely work, but the paper outputs were obtained from git commit
c4a52d434053d676c0281449b0fce7112116fd54

To run the code, download a copy (on linux)

`git clone https://github.com/edwardsmith999/flowmol`

navigate to the source directory

`cd flowmol/src`

You need to have gcc complilers with fortran support (gfortran) or intel compliers (ifort) and MPI (tested with mpich)

`make PLATFORM=gfortran p`

or

`make PLATFORM=intel p`

If this works, you should get an executable "parallel_md.exe"
Move this to a runs/Intrinsic_interface directory and to run the case call,

`mpiexec -n 1 ./parallel_md.exe -i ./MD_2_phase.in -r ./initial_state`

This will take a long time to get good statistics and generate large files.
Once finished, you can analyse results by first calling

`python Analyse_data_and_pickle.py`

Not you need numpy and matplotlib. The default behaviour of this script is to read
the summary.p file (which is included from the data used in the paper above). If
you want to analyse your own data, delete or rename this summary.p folder and a new
one will be created from your data.
To generate the plots from the paper, you can then run

`python figures_for_paper.py`

Any questions, please contact the author on edward.smith@brunel.ac.uk or edwardsmith999@hotmail.com
