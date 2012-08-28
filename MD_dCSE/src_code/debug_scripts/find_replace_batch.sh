#!/bin/bash

#ARRAYS ($A)

A[1]="r"
A[2]="rtrue"
A[3]="rinitial"
A[4]="rijsum"
A[5]="v"
A[6]="vtrue"
A[7]="vmagnitude"
A[8]="a"
A[9]="aold"
A[10]="theta"
A[11]="aD"
A[12]="aR"
A[13]="tag"
A[14]="fix"
A[15]="slidev"
A[16]="thermostat"

#MOLECULAR COLUMN ($M)
M[1]='molno'
M[2]='molnoi'
M[3]='molnoj'
M[4]='molnoX'
M[5]='n'
M[6]='m'
M[7]='i'
M[8]='j'
M[9]='k'
M[10]='nl'
M[11]='np+n'
M[12]='np+m'
M[13]='np+i'
M[14]='np+j'
M[15]='np+k'
M[16]=':'
M[17]='np+1:np+halo_np'
M[18]='1:np'
M[19]='molnopush'
M[20]='ip'
M[21]='np+new_np'

#DIMENSIONAL COLUMN ($D)
D[1]=':'
D[2]='1'
D[3]='2'
D[4]='3'
D[5]='ixyz'
D[6]='jxyz'
D[7]='kxyz'
D[8]='le_sp'
D[9]='le_rp'
D[10]='np'
D[11]='copyplane'

for k in {1..16}
do
	for j in {1..11}
	do
		for i in {1..21}
		do
			x="${A[k]}(${M[i]},${D[j]})"
			y="${A[k]}(${D[j]},${M[i]})"
			s="s/$x/$y/g"
	    	#echo $s
	    	#echo "=========================================="
	    	#echo "$x => $y"
			grep $x *.f90
	    	#echo "=========================================="
			#grep -rl $x *.f90 | xargs sed -i -e $s

		done
	done
done

