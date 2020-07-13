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
rpm = False
start=1
end=5

def fittedfunc(x,a,m,b):
	return a - m * x - b*np.log(x)

def fittedfunc2(x,a,m,b):
	return a*np.exp(-m*x)/(x**b)

def bindingenergy(x,a1,m,b1,a2,M,b2):
	return (2*a1-a2-2*m*x+M*x+(2*b1-b2)*np.log(x))/x

sig = signature(fittedfunc2)
numpara = len(sig.parameters) - 1

m=0.001
totalcor1,allpara1,parastderr1,allchisq1 = cf.fits(fittedfunc,Path('./correlatordata/allcorrelators_m0=%f.npy'%m),maxlen,corperconfig,binsize,m,start,end,numpara,'log',0,rpm)
L = len(allpara1[1,:])

# E= np.zeros(L)
# for bmass in np.arange(0.001,0.051,0.001):
# 	m,M,E[int(bmass*1000-1)] = smearefit(fittedfunc,Path('./correlatordata/allcorrelators_m0=%f.npy'%bmass),maxlen,corperconfig,binsize,bmass,start,end,numpara,'log',1)



onemass = np.zeros((50,L))
twomass = np.zeros((50,L))
for bmass in np.arange(0.001,0.051,0.001):
# totalcor1,finalpara1,parastderr,finalchisq = cf.fits(fittedfunc,Path('testt.npy'),maxlen,corperconfig,binsize,bmass,start,end,numpara,'log')
	totalcor1,allpara1,parastderr1,allchisq1 = cf.fits(fittedfunc,Path('./correlatordata/allcorrelators_m0=%f.npy'%bmass),maxlen,corperconfig,binsize,bmass,start,end,numpara,'log',0,rpm)
	onemass[int(1000*bmass-1),:]=allpara1[1,:]
	totalcor2,allpara2,parastderr2,allchisq2 = cf.fits(fittedfunc,Path('./correlatordata/alltwoparticlecorrelators_m0=%f.npy'%bmass),maxlen,corperconfig,binsize,bmass,start,end,numpara,'log',0,rpm)
	twomass[int(1000*bmass-1),:]=allpara2[1,:]

np.save(Path('./mass/onemass_4b0_%d-%d'%(start,end)),onemass)
np.save(Path('./mass/twomass_4b0_%d-%d'%(start,end)),twomass)





# print('chi:',finalchisq1)
# totalcor,finalpara,parastderr,finalchisq = cf.fits(fittedfunc2,Path('./correlatordata/allcorrelators_m0=%f.npy'%bmass),maxlen,corperconfig,binsize,bmass,start,end,numpara,'exp')
# m=finalpara1[1] 
# totalcor,finalpara2,parastderr,finalchisq = cf.fits(fittedfunc2,Path('./correlatordata/alltwoparticlecorrelators_m0=%f.npy'%bmass),maxlen,corperconfig,binsize,bmass,start,end,numpara,'exp')
# M=finalpara2[1]


exit()
nonzero1 = 10
fig1 = plt.figure()
#Plot the correlator with error bars
x = np.arange(1,nonzero1)
#y = fittedfunc2(x,finalpara[0],finalpara[1],finalpara[2])
y = np.exp(fittedfunc(x,finalpara1[0],finalpara1[1],finalpara1[2]))
y1 = np.mean(totalcor1[1:nonzero1],axis = 1)
err = np.std(totalcor1[1:nonzero1],axis = 1)/np.sqrt(len(totalcor1[0,:])-1)
#err1 = err[0:nonzero]
plt.errorbar(x,y1,yerr = err, ecolor = 'red',linewidth = 1, capsize = 1)
#plt.yscale('log')
plt.plot(x,y,'red')
plt.grid(True)
plt.xlabel('Geodesic Distance')
plt.ylabel('Correlator')
plt.title('start %d, end %d'%(start,end))
plt.show()


exit()
plotlenth = 33
plotdata(totalcortmp1,plotlenth)
plotfit(finalpara,form,plotlenth)