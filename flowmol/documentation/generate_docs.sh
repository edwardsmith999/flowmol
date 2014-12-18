# Routine to generate documentation for the coupler code using f90doc
# Requires f90doc to be in user's paths
f90doc ./../coupler_dCSE/src_code/coupler*.f90
mv coupler*.html ./f90doc_output/
