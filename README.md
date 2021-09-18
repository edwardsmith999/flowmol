# flowmol
A molecular dynamics solver for molecular fluid dynamics simulation.

You only need a Fortran compiler. Assuming you have gfortran installed (from gcc), build it but going to the src directory and calling

    make PLATFORM=gfortran p
    
Intel or other specialised supercomputer architectures are supported, check the platform file (and adding new ones is fairly straightforward). 

An experimental user interface is provided in python,

    python3 flowmol_input.py

Which allows the user to tweak input files, run multiple cases (started in subprocess but can exist beyond the lifetime of flowmol_input)  and analyse results using an inbuilt version of PyDataView.

![image](https://user-images.githubusercontent.com/13366737/133908501-c3cbe270-0df5-4396-8bd0-ade8b5fa389c.png)

I will write up install instructions soon, but it requires Python3, wxPython, matplotlib, Numpy as well as inclusion of visualisation software [PyDataView](https://github.com/edwardsmith999/pyDataView) and thread running software [SimWrapLib](https://github.com/edwardsmith999/SimWrapPy).  
You will also have to have a working version of Flowmol compiled.

Please see my website for more information on the Flowmol code, including verification/validation and features

http://www.edwardsmith.co.uk/content/flowmol/MD_coding.html
