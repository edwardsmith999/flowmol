cp ../../Couette_serial/input_continuum ./couette_serial
cp ../../Couette_serial/continuum.exe ./couette_serial
mpiexec -n 1 ./../../MD_dCSE/src_code/parallel_md.exe -i ./md_data/MD.in : -n 1 ./couette_serial/continuum.exe
