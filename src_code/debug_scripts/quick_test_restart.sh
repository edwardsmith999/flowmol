#! /usr/bin/env bash

EXPECTED_ARGS=2
if [ $# -ne $EXPECTED_ARGS ]
then
	echo ""
	echo " ./test_restart.sh requires 2 arguments: "
	echo ""
	echo "    1) INFILE: input file name (string)"
	echo "    2) N1:     number of iterations to test (integer)"
	echo ""
	echo "The script will perform the following simulations based on the input file INFILE:"
	echo ""
	echo "    1) Serial initialisation simulation for \$N1 iterations"
	echo "           * ./results/final_state -> ./\$N1.state"
	echo "    2) Serial restart from \$N1.state for a further \$N1 iterations"
	echo "           * ./results/final_state -> ./\$N1-(2*\$N1).state"
	echo "    3) Continuous serial simulation for (2*\$N1) iterations"
	echo "           * ./results/final_state -> ./(2*\$N1).state"
	echo "    4) Parallel restart from \$N1.state for a further \$N1 iterations"
	echo "           * ./results/final_state -> ./\$N1-(2*\$N1)_parallel.state"
	exit
fi

INFILE=$1
N1=$2
N2=`echo "${N1}*2" | bc`

# Clean
make clean

# Make sure VMD output is on
sed -i '/VMD_OUTFLAG/{n; s/.*/1/}' $INFILE
# Make sure macro outflag is on 2
sed -i '/MACRO_OUTFLAG/{n; s/.*/2/}' $INFILE

# Set NSTEPS to N1 and run simulation
sed -i '/NSTEPS/{n; s/.*/'$N1'/}' $INFILE
make debug_s
./md.exe -i $INFILE
# Collect final_state, macro and vmd files
mv results/final_state ./$N1.state
mv results/macroscopic_properties ./$N1.macro
mv results/vmd_out.dcd ./${N1}.dcd

# Restart sim
./md.exe -i $INFILE -r $N1.state
# Collect final_state, macro and vmd files
mv results/final_state ./$N1-$N2.state
mv results/macroscopic_properties ./$N1-$N2.macro
mv results/vmd_out.dcd ./$N1-$N2.dcd

# Run complete simulation without restart
sed -i '/NSTEPS/{n; s/.*/'$N2'/}' $INFILE
./md.exe -i $INFILE
# Collect final_state, macro and vmd files
mv results/final_state ./$N2.state
mv results/macroscopic_properties ./$N2.macro
mv results/vmd_out.dcd ./$N2.dcd

# Clean and test parallel
make clean
make debug_p

# Set NSTEPS to N1 again
sed -i '/NSTEPS/{n; s/.*/'$N1'/}' $INFILE
# Set processors to 4, 2 in x and 2 in y
sed -i '/PROCESSORS/{n; s/.*/2/; n; s/.*/2/; n; s/.*/1/}' $INFILE
mpiexec -n 4 ./md.exe -i $INFILE -r $N1.state
mv results/final_state ./${N1}-${N2}_parallel.state
mv results/macroscopic_properties ./${N1}-${N2}_parallel.macro
mv results/vmd_out.dcd ./${N1}-${N2}_parallel.dcd
