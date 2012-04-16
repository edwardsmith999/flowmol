#! /usr/bin/env python

s = open('./debug_scripts/serial.macro',   'r')
p = open('./debug_scripts/parallel.macro', 'r')
d = open('./debug_scripts/diff.macro',     'w')

header = s.readline()
header = p.readline()
d.write(header)

while 1:
	serial = s.readline()
	parallel = p.readline()
	if not serial or not parallel:
		break
	serial   = serial.strip().split(';')
	parallel = parallel.strip().split(';')

	diff = [0.0]*len(serial)
	diff[0] = int(serial[0])

	outstring = str(diff[0]).rjust(12)
	for i in range(1,len(serial)):
		diff[i] = "%22.15f" % (float(serial[i]) - float(parallel[i]))
		outstring += diff[i]
	d.write(outstring+'\n')

d.write('\n')
d.close()
s.close()
p.close()
