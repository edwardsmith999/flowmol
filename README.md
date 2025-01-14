<img src="out-1.webp" alt="Flowmol Logo" width="200"/>

A molecular dynamics solver for molecular fluid dynamics simulation.

You only need a Fortran compiler and MPI (tested with [MPICH](https://www.mpich.org/)). Assuming you have gfortran installed (from gcc), build it by going to the src directory and calling

    make PLATFORM=gfortran p
    
The Intel compliers should also work, and the code has been tested on a range of specialised supercomputer architectures, check the platform file to see some examples of build flags (and adding new ones is fairly straightforward). There is an assembly language version of the Heaviside function used for efficiency, which might cause problems on some platform. This can be disabled using

    make NO_ASSMBLY_HEAVISIDES=1 p
    
Simulations can then be run using,

    mpiexec -n 1 ./parallel_md.exe -i default.in
    
where default.in is an example input file to run a simple NVE case starting from an FCC lattice. 
You can change simulation parameters by changing the input file, which uses a keyword lookup system with a word in capitals describing the input variable to change and the following numbers being the values you want to set. For example, to make the simulation bigger you would change the following part,

    INITIALNUNITS
    8
    8
    8
    
where the default is 8 FCC units by 8 by 8. For documentation on the different inputs, please see `setup_read_inputs.f90` which includes notes on the form of inputs and what they do. This input file is then parsed to create documentation in flowmol_input, an experimental user interface provided in python,

    python3 flowmol_input.py

Which allows the user to tweak input files, run multiple cases (started in subprocess but can exist beyond the lifetime of flowmol_input)  and analyse results using an inbuilt version of PyDataView.

![image](https://user-images.githubusercontent.com/13366737/133908501-c3cbe270-0df5-4396-8bd0-ade8b5fa389c.png)

It requires Python3, wxPython, matplotlib, Numpy as well as an inclusion of visualisation software [PyDataView](https://github.com/edwardsmith999/pyDataView) and thread running software [SimWrapLib](https://github.com/edwardsmith999/SimWrapPy). You will also have to have a working version of Flowmol compiled to run cases.

For more information on the Flowmol code, including verification/validation and features, please see
http://www.edwardsmith.co.uk/content/flowmol/MD_coding.html
