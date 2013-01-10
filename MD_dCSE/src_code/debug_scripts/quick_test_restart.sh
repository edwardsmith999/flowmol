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
	exit
fi

INFILE=$1
N1=$2
N2=`echo "${N1}*2" | bc`

# Make sure VMD output is on
sed -i '/VMD_OUTFLAG/{n; s/.*/1/}' $INFILE
# Make sure macro outflag is on 2
sed -i '/MACRO_OUTFLAG/{n; s/.*/2/}' $INFILE

################################################################
#    S E R I A L   O N L Y     R E S T A R T 

make clean
make debug_s

# First run
sed -i '/NSTEPS/{n; s/.*/'$N1'/}' $INFILE
./md.exe -i $INFILE
mv results/final_state ./${N1}_s.state
mv results/macroscopic_properties ./${N1}_s.macro
mv results/vmd_out.dcd ./${N1}_s.dcd

# Restart from first run 
./md.exe -i $INFILE -r ${N1}_s.state
mv results/final_state ./${N1}-${N2}_s2s.state
mv results/macroscopic_properties ./${N1}-${N2}_s2s.macro
mv results/vmd_out.dcd ./${N1}-${N2}_s2s.dcd

# Run complete simulation without restart
sed -i '/NSTEPS/{n; s/.*/'$N2'/}' $INFILE
./md.exe -i $INFILE
mv results/final_state ./${N2}_s.state
mv results/macroscopic_properties ./${N2}_s.macro
mv results/vmd_out.dcd ./${N2}_s.dcd

################################################################
#    P A R A L L E L   O N L Y    R E S T A R T 

make clean
make debug_p

# Set processors to 4, 2 in x and 2 in y
sed -i '/PROCESSORS/{n; s/.*/2/; n; s/.*/2/; n; s/.*/1/}' $INFILE
# Set NSTEPS to N1 again
sed -i '/NSTEPS/{n; s/.*/'$N1'/}' $INFILE

# First run 
mpiexec -n 4 ./parallel_md.exe -i $INFILE
mv results/final_state ./${N1}_p.state
mv results/macroscopic_properties ./${N1}_p.macro
mv results/vmd_out.dcd ./${N1}_p.dcd

# Restart from first run 
mpiexec -n 4 ./parallel_md.exe -i $INFILE -r ${N1}_p.state
mv results/final_state ./${N1}-${N2}_p2p.state
mv results/macroscopic_properties ./${N1}-${N2}_p2p.macro
mv results/vmd_out.dcd ./${N1}-${N2}_p2p.dcd

# Run complete simulation without restart
sed -i '/NSTEPS/{n; s/.*/'$N2'/}' $INFILE
mpiexec -n 4 ./parallel_md.exe -i $INFILE
mv results/final_state ./${N2}_p.state
mv results/macroscopic_properties ./${N2}_p.macro
mv results/vmd_out.dcd ./${N2}_p.dcd


################################################################
#    S E R I A L / P A R A L L E L    R E S T A R T S

sed -i '/NSTEPS/{n; s/.*/'$N1'/}' $INFILE

mpiexec -n 4 ./parallel_md.exe -i $INFILE -r ${N1}_s.state
mv results/final_state ./${N1}-${N2}_s2p.state
mv results/macroscopic_properties ./${N1}-${N2}_s2p.macro
mv results/vmd_out.dcd ./${N1}-${N2}_s2p.dcd

make clean
make debug_s

./md.exe -i $INFILE -r ${N1}_p.state
mv results/final_state ./${N1}-${N2}_p2s.state
mv results/macroscopic_properties ./${N1}-${N2}_p2s.macro
mv results/vmd_out.dcd ./${N1}-${N2}_p2s.dcd

#################################################################
## P L O T T I N G
# Find out which columns KE and PE are in
if grep -q "FENE" *.macro
then
	tcol=2
	KEcol=6
	PEcol=9
else
	tcol=2
	KEcol=6
	PEcol=7
fi

plotfile=./debug_scripts/quick_test_restart.gp
sed -i "1{s/[0-9].*/${tcol}/};2{s/[0-9].*/${KEcol}/};3{s/[0-9].*/${PEcol}/}" $plotfile
cat $plotfile | gnuplot --persist
