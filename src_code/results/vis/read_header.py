#! /usr/bin/env python
f = open('simulation_header','r')
headerdata = f.readlines()
for i in range(len(headerdata)):
	headerdata[i] = headerdata[i].split(';')
	for j in range(len(headerdata[i])):
		headerdata[i][j] = headerdata[i][j].strip().replace('(','').replace(')','')
	vars()[headerdata[i][1]] = headerdata[i][2]
f.close()

f = open('simulation_progress','r')
Nsteps = f.read().strip()
f.close()
