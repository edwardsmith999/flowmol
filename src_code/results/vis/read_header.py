#! /usr/bin/env python
#########################################################################################
# read_header.py
# author: David Trevelyan
#
# This script is designed for importing into other python routines that require
# information from the file 'simulation_header'
#
# ---------------------------------------------------------------------------------------

f = open('simulation_header','r')						# Create ascii file object
headerdata = f.readlines()								# Read all lines into list
for i in range(len(headerdata)):						# Loop through all elements in list (i.e. lines in simulation_header) 
	headerdata[i] = headerdata[i].split(';')			# Split list elements into sub-lists according to ; character in simulation_header
	for j in range(len(headerdata[i])):					# Loop through all elements in sub-list
		headerdata[i][j] = headerdata[i][j].strip().replace('(','').replace(')','') # Strip white space and remove brackets (can't create arrays using vars)
	vars()[headerdata[i][1]] = headerdata[i][2]			# Create variables with name headerdata[i][1] and value headerdata[i][2]
#	print(str(headerdata[i][1]),str(headerdata[i][2]))
f.close()												# Close file object

f = open('simulation_progress','r')						# Creat ascii file object
Nsteps = f.read().strip()								# Read number of timesteps 
f.close()												# Close file object
