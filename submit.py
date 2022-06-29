#!/usr/bin/env python
import os
ifnprods = [0,1,3,6,10]

isstochastic=[["homogenous",0,280]]
#initial_virus = [[1,1,1],[20,20,10],[20,20,100],[200,500,10],[200,500,100]]
#initial_virus = [[1,1,1]]
initial_virus = [[20,20,10]]
for name,Is,totalstep in isstochastic:
	for s1,s2,vnumb in initial_virus:
		for iifnprods in ifnprods:
			path = "./initial_virus/{}/{}_{}_{}/{}/".format(name,s1,s2,vnumb,iifnprods)
			os.chdir(path)
			os.system("python COVID.py {} {} {} {} {} {}&".format(iifnprods,s1,s2,vnumb,Is,totalstep))
			os.chdir("../../../../")
