#!/bin/sh
#
mpiexe="mpiexec"
left=" -n 1 xterm -geometry "
right=" -hold -e $2 ./md.exe $3 $4 $5 $6 :"
end=" -hold -e $2 ./md.exe $3 $4 $5 $6"
xwidth=150
ywidth=10
xpos=1000
ypos=0

for (( i=1; i<$1; i++ ))
do
	size_loc="${xwidth}x${ywidth}+${xpos}+${ypos}"
	fullline="$fullline$left$size_loc$right"
	ypos=$((150 * $i))
done

size_loc="${xwidth}x${ywidth}+${xpos}+${ypos}"
fullline="$mpiexe$fullline$left$size_loc$end"

$fullline

#mpiexec -n 1 xterm -geometry 150x20+1000+0 -hold -e ./md.exe :  -n 1 xterm -geometry 150x20+100+220 -hold -e ./md.exe :  -n 1 xterm -geometry 150x20+100+240 -hold -e ./md.exe :  -n 1 xterm -geometry 150x20+100+260 -hold -e ./md.exe
#mpiexec -n $1 xterm -geometry 150x20 -hold -e ./md.exe
