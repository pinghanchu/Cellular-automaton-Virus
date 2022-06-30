#!/usr/bin/env python
import os
ifnprods = [0,1,3,6,10]

isstochastic=[["homogenous",0,280]]
initial_virus = [[1,1,1],[20,20,10]]
#initial_virus = [[1,1,1]]
#initial_virus = [[200,500,10]]
for name,Is,totalstep in isstochastic:
	for s1,s2,vnumb in initial_virus:
		for iifnprods in ifnprods:
			pathin = "./Archive/{}_{}_{}/{}/data".format(s1,s2,vnumb,iifnprods)
			pathout = "./initial_virus/{}/{}_{}_{}/{}/data/".format(name,s1,s2,vnumb,iifnprods)
			#os.system("ls {}/*.csv".format(pathin))
			os.system("mv {}/*.csv {}".format(pathin,pathout))
