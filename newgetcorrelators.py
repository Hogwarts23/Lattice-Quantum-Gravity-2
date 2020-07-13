#!/home/mdai07/application/miniconda3/bin/python
from configclass import LQGConfigs
import os
import numpy as np
from pathlib import Path
import sys
import math
from multiprocessing.dummy import Pool as ThreadPool
# import time

maxlen = 50
numberofjobs = 50
nsources = 20
f = "/home/mdai07/EDT_configs/4b0/"
fs = '/home/mdai07/run/4b0/sources'
# ff = "./4b0"

#the program will do the i-th job of all the "njobs" submitted
def jobdistributor(njobs,i,folder):#i goes from 0 to njobs -1 
	dirs = sorted(os.listdir(Path(folder)))
	# print(type(dirs))
	l = len(dirs)
	nconfigs = math.floor(l/njobs)
	if i != (njobs-1):
		return dirs[int(i*nconfigs):int((i+1)*nconfigs)]
	else:
		return dirs[int(i*nconfigs):l]

def getcorrelators(fname):
	m = LQGConfigs(f+fname,np.arange(0.001,0.051,0.001),nsources,fs)
	totalcor,totalcor2 = m.correlators(True)
	np.save(Path("/home/mdai07/run/4b0/correlatordata/%s_correlator"%fname),totalcor)
	np.save(Path("/home/mdai07/run/4b0/correlatordata/%s_twoparcorr"%fname),totalcor2)
	# np.save(Path("/home/mdai07/run/4b0/sources/%s"%fname),m.sources)

	# np.save(Path("/home/mdai07/run/32b0/correlatordata/%s_correlator"%fname),totalcor)
	# np.save(Path("/home/mdai07/run/32b0/correlatordata/%s_twoparcorr"%fname),totalcor2)
	# np.save(Path("/home/mdai07/run/32b0/sources/%s"%fname),m.sources)

# files = jobdistributor(50,1,ff)
files = jobdistributor(50,int(sys.argv[1]),f)
# pool = ThreadPool(2)
# pool.map(getcorrelators,files)
# pool.close()
# pool.join()

for file in files:
	getcorrelators(file)



