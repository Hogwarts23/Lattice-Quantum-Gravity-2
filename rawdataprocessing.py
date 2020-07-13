import numpy as np
import os
from pathlib import Path
nsources = 20

dirs = os.listdir(Path("./rawcorrelatordata"))
l = len(dirs)
m = np.arange(0.001,0.051,0.001)
for i in range(50):
	data1 = np.zeros((50,int(l*nsources/2)))
	data2 = np.zeros((50,int(l*nsources/2)))
	j1=0
	j2=0
	for f in dirs:
		if "_correlator" in f:
			x1 = np.load(Path("./rawcorrelatordata")/f)
			data1[:,int(j1*nsources):int(j1*nsources+nsources)] = x1[:,int(i*nsources):int(i*nsources+nsources)]
			j1 = j1 + 1
		if "_two" in f:
			x2 = np.load(Path("./rawcorrelatordata")/f)
			data2[:,int(j2*nsources):int(j2*nsources+nsources)] = x2[:,int(i*nsources):int(i*nsources+nsources)]
			j2 = j2 + 1
	np.save(Path('./correlatordata/allcorrelators_m0=%6f'%m[i]),data1)
	np.save(Path('./correlatordata/alltwoparticlecorrelators_m0=%6f'%m[i]),data2)
