NPROC=$1

mpiexec -n $NPROC ./parallel_md.exe -i concentric.in > concentric.out
mv results/vmd_out.dcd ./concentric.dcd
mv results/final_state ./concentric.state
mv results/cylinders ./concentric.cyl

mpiexec -n $NPROC ./parallel_md.exe -i fill.in -c concentric.cyl > fill.out
mv results/vmd_out.dcd ./fill.dcd
mv results/final_state ./fill.state

mpiexec -n $NPROC ./parallel_md.exe -i rotate.in -r fill.state -c concentric.cyl > rotate.out
mv results/vmd_out.dcd ./rotate.dcd
mv results/mbins ./rotate.mbins
mv results/vbins ./rotate.vbins
mv results/macroscopic_properties ./rotate.macro
mv results/final_state ./rotate.state

mpiexec -n $NPROC ./parallel_md.exe -i calcs.in -r rotate.state -c concentric.cyl > calcs.out
#cp results/mbins ./calcs.mbins
#cp results/vbins ./calcs.vbins
mv results/macroscopic_properties ./calcs.macro
