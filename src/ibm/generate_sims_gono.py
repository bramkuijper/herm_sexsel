#!/usr/bin/env python

import datetime
import numpy as np

avals = [ 0.3 ]
rvals = list(np.arange(0,2,0.05))
dvals = list(np.arange(0,2,0.05))

exe = "./gonochorist_mutual_direct_benefits.exe"

ctr = 0

date = datetime.datetime.now()

base_name = "sim_gono_" +\
        f"{date:%d}_{date:%m}_{date:%Y}_{date:%H}{date:%M}{date:%S}"

for a in avals:
    for r in rvals:
        for d in dvals:
            print("echo " + str(ctr))

            print(exe + " " + str(a) + " " + str(d) + " 0.001 0.2 0.001 " + str(r) + " 0.0 0.01 0.01 0.01 0.05 0.05 0.05 1.0 1.0 " + base_name + "_" + str(ctr))

            ctr += 1


