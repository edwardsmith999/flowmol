#!/usr/bin/env python

from sim import *
import math

sim1 = Simulation()

sim1.changeLine('ENSEMBLE',1 )
sim1.changeLine('SLRC_FLAG',1 )

pressurelist = []
#f = open('./output', "w+b")
for i in range(0, 10):
	density = 0.1 + i*0.1
	sim1.changeLine('DENSITY',density)
	for j in range(0, 15):
		temperature = 0.2 + j * 0.1
		sim1.changeLine('MACRO_OUTFLAG',1 )
		sim1.changeLine('TEMPERATURE',temperature )	
		sim1.changeLine('NSTEPS',1000)

		#Run to generate restart file
		sim1.restart = False
		sim1.run()

		#sim1.changeLine('ENSEMBLE',0 )
		sim1.changeLine('MACRO_OUTFLAG',2 )
		sim1.changeLine('NSTEPS',10000)
		sim1.restart = True
		sim1.run()

		Pressure = sim1.macrolist('Pressure')
		print Pressure
		#f.writelines(["%s\n" % Pressure])
		#f = open('./output', "w+b")
		#pickle.dump(Pressure, f)
		#f.close()
		a = (density,temperature,mean(Pressure),std(Pressure))
		pressurelist.append(a)
		open('./try_sim_output', 'w').writelines(["%s\n" % pressurelist])
		#print mean(Pressure)
		#print std(Pressure)

#f.close()
#Write to ouput file


#sim1.run()
#sim1.writeSummary()

#sim1.changeLine('DENSITY',0.8)
#sim1.changeLine('TEMPERATURE',1.2 )
#sim1.run()
#sim1.writeSummary()


#plot=True
#analytic=True
