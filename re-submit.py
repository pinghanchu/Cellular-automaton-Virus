#!/usr/bin/env python                                                                                                 
import os
import numpy as np
ifnprods = [0,1,3,6,10]
isstochastic=[["homogenous",0,280]]
initial_virus = [[1,1,1]]
for name,Is,totalstep in isstochastic:
        for s1,s2,vnumb in initial_virus:
                for iifnprods in ifnprods:
                        os.system("ls ./initial_virus/{}/{}_{}_{}/{}/data/*.csv > list.txt".format(Is,s1,s2,vnumb,iifnprods))
                            
                        f = open("list.txt", "r")
                        Lines = f.readlines()
                        if len(Lines)>0:
                            run = []
                            for line in Lines:
                                    a1 = line.split("_")
                                    a2 = a1[len(a1)-1].split("run")
                                    a3 = a2[1].split(".csv")
                                    run.append(int(a3[0]))
                            runarr = np.asarray(run)
                            runarr = np.sort(runarr)
                            lastrun = runarr[len(runarr)-1]
                            path ="./initial_virus/{}/{}_{}_{}/{}".format(name,s1,s2,vnumb,iifnprods)
                            filename = "data_virusprods2_virusdiff1_ifnprods{}_ifndiff5_ifnprob0.1_virusreduct0.5_isstochastic{}_initialvirus{}_{}_{}_run{}.csv" .format(iifnprods,Is,s1,s2,vnumb,lastrun)
                            os.system("rm {}/data/{}".format(path,filename))
                            os.chdir(path)
                            for irun in range(lastrun,100):
                                os.system("python COVIDi.py {} {} {} {} {} {} {}&".format(iifnprods,s1,s2,vnumb,Is,totalstep,irun))
                            os.chdir("../../../../")
                        f.close()
                        os.system("rm list.txt")
