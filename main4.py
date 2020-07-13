import numpy as np
import matplotlib.pyplot as plt
import correlatefitter as cf
from inspect import signature
from pathlib import Path
maxlen = 50
corperconfig = 20
binsize = 20
start=2
end=9

def fittedfunc(x,a,m,b):
	return a - m * x - b*np.log(x)

def fittedfunc2(x,a,m,b):
	return a*np.exp(-m*x)/(x**b)

def bindingenergy(x,a1,m,b1,a2,M,b2):
	return (2*a1-a2-2*m*x+M*x+(2*b1-b2)*np.log(x))/x

sig = signature(fittedfunc2)
numpara = len(sig.parameters) - 1

onemass = np.load('./mass/onemass_4b0_2-9.npy')
twomass = np.load('./mass/twomass_4b0_2-9.npy')
N = len(onemass[1,:])
print(N)

E = np.zeros(50)
Estd = np.zeros(50)
for bmass in np.arange(0.001,0.051,0.001):
	m1,M1,E[int(bmass*1000-1)],Estd[int(bmass*1000-1)] = cf.smearefit(fittedfunc,Path('./correlatordata/allcorrelators_m0=%f.npy'%bmass),maxlen,corperconfig,binsize,bmass,start,end,numpara,'log',0)
#E[6] = (E[5]+E[7])/2
plt.figure()
x = np.arange(0.001,0.051,0.001)
plt.errorbar(x,E,yerr=Estd*np.sqrt(N-1),label = 'Smeared source',marker='o',capsize = 0.5)
# plt.show()


m=0.001
totalcor1,allpara1,parastderr1,allchisq1 = cf.fits(fittedfunc,Path('./correlatordata/allcorrelators_m0=%f.npy'%m),maxlen,corperconfig,binsize,m,start,end,numpara,'log',0)
L = len(allpara1[1,:])

# onemass = np.zeros((50,L))
# twomass = np.zeros((50,L))
# for bmass in np.arange(0.001,0.051,0.001):
# 	totalcor1,allpara1,parastderr1,allchisq1 = cf.fits(fittedfunc,Path('./correlatordata/allcorrelators_m0=%f.npy'%bmass),maxlen,corperconfig,binsize,bmass,start,end,numpara,'log',0)
# 	onemass[int(1000*bmass-1),:]=allpara1[1,:]
# 	totalcor2,allpara2,parastderr2,allchisq2 = cf.fits(fittedfunc,Path('./correlatordata/alltwoparticlecorrelators_m0=%f.npy'%bmass),maxlen,corperconfig,binsize,bmass,start,end,numpara,'log',0)
# 	twomass[int(1000*bmass-1),:]=allpara2[1,:]

# np.save(Path('./mass/onemass_0.5b0-l'),onemass)
# np.save(Path('./mass/twomass_0.5b0-l'),twomass)

m = np.mean(onemass,axis=1)
mstd = np.std(onemass,axis=1)
M = np.mean(twomass,axis=1)
Mstd = np.std(twomass,axis=1)
#m[6] = (m[5]+m[7])/2
#M[6] = (M[5]+M[7])/2
binding = 2*m - M
bstd = np.sqrt(4*(mstd**2)+Mstd**2)
# tot = np.zeros((2,50))
# tot[0,:]=binding
# tot[1,:]=E
#np.save('2-9 binding energy',tot)
plt.errorbar(x,binding,yerr=bstd*np.sqrt(N-1),label='Binding Energy',capsize = 0.5,marker='*')
plt.grid()
plt.legend()
plt.xlabel('bare mass m_0')
plt.ylabel('Binding Energy')
plt.title('Binding energy using smeaered sources vs binding energy')
plt.show()



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
plt.errorbar(x,y1,yerr = err,ecolor = 'red',linewidth = 1, capsize = 1)
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