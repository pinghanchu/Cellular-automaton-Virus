#!/usr/bin/env python
import os

ifnprods = [0,1,3,6,10]
initial_virus = [[1,1,1],[20,20,10],[20,20,100],[200,500,10],[200,500,100]]
isstochastic=[["stochastic",1,1000],["homogenous",0,280]]
isimage= 0
totalrun=100

os.system("mkdir data_py figure_py")
for name,Is,totalstep in isstochastic:
    for s1,s2,vnumb in initial_virus:
        for iifnprods in ifnprods:
            os.system("python COVID.py {} {} {} {} {} {} {} {}&".format(iifnprods,s1,s2,vnumb,Is,isimage,totalstep,totalrun))
