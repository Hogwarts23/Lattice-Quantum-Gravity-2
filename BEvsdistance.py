import numpy as np
import matplotlib.pyplot as plt
import correlatefitter as cf
import time
from inspect import signature
import numpy.linalg as lin
import random
from pathlib import Path

maxlen = 50
corperconfig = 20
binsize = 20
start=1
end=5

def fittedfunc(x,a,m,b):
	return a - m * x - b*np.log(x)

def fittedfunc2(x,a,m,b):
	return a*np.exp(-m*x)/(x**b)

def bindingenergy(x,a1,m,b1,a2,M,b2):
	return -(2*a1-a2-2*m*x+M*x+(2*b1-b2)*np.log(x))

sig = signature(fittedfunc2)
numpara = len(sig.parameters) - 1

m=0.001
totalcor1,allpara1,parastderr1,allchisq1 = cf.fits(fittedfunc,Path('./correlatordata/allcorrelators_m0=%f.npy'%m),maxlen,corperconfig,binsize,m,start,end,numpara,'log',0)
L = len(allpara1[1,:])


for bmass in np.arange(0.001,0.051,0.001):
# totalcor1,finalpara1,parastderr,finalchisq = cf.fits(fittedfunc,Path('testt.npy'),maxlen,corperconfig,binsize,bmass,start,end,numpara,'log')
	totalcor1,finalpara1,parastderr1,allchisq1 = cf.fits(fittedfunc,Path('./correlatordata/allcorrelators_m0=%f.npy'%bmass),maxlen,corperconfig,binsize,bmass,start,end,numpara,'log',1)
	p1 = np.mean(totalcor1,axis=1)
	std1 = np.std(totalcor1,axis=1)/np.sqrt(totalcor1.shape[1]-1)
	totalcor2,finalpara2,parastderr2,allchisq2 = cf.fits(fittedfunc,Path('./correlatordata/alltwoparticlecorrelators_m0=%f.npy'%bmass),maxlen,corperconfig,binsize,bmass,start,end,numpara,'log',1)
	p2 = np.mean(totalcor2,axis=1)
	std2 = np.std(totalcor2,axis=1)/np.sqrt(totalcor2.shape[1]-1)
	x = np.arange(1,15,0.1)
	y = bindingenergy(x,finalpara1[0],finalpara1[1],finalpara1[2],finalpara2[0],finalpara2[1],finalpara2[2])
	x1 = np.arange(1,30)
	y1 = np.zeros(29)
	for i in range(1,30):
		y1[i-1] = np.log(p2[i]/p1[i]/p1[i])
	fig1 = plt.figure()
	plt.plot(x,y)
	plt.errorbar(x1,y1,marker='*',yerr=np.sqrt(std1[1:30]**2+std2[1:30]**2))
	plt.xlabel('Geodesic Distance')
	plt.ylabel('Log of Binding Energy')
	plt.title('start %d, end %d'%(start,end))
	plt.savefig('./befigure15/%6f.png'%bmass)
	plt.close()
	del fig1,x,y,x1,y1,totalcor1,finalpara1,parastderr1,allchisq1,totalcor2,finalpara2,parastderr2,allchisq2