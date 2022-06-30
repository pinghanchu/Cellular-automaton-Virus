#!/usr/bin/env python
import os
ifnprods = [0,1,3,6,10]
isstochastic=["homogenous","stochastic"]
initial_virus = [[1,1,1],[20,20,10],[20,20,100],[200,500,10],[200,500,100]]
os.system("mkdir initial_virus")
for Is in isstochastic:
	os.system("mkdir initial_virus/{}".format(Is))
	for s1,s2,vnumb in initial_virus:
		os.system("mkdir initial_virus/{}/{}_{}_{}/".format(Is,s1,s2,vnumb))
		for iifnprods in ifnprods:
			os.system("mkdir ./initial_virus/{}/{}_{}_{}/{}/".format(Is,s1,s2,vnumb,iifnprods))
		        os.system("mkdir ./initial_virus/{}/{}_{}_{}/{}/data/".format(Is,s1,s2,vnumb,iifnprods))	
			path = "./initial_virus/{}/{}_{}_{}/{}/".format(Is,s1,s2,vnumb,iifnprods)
			os.system("cp COVIDi.py {}".format(path))
			os.system("cp COVID.py {}".format(path))
