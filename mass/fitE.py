import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import correlatefitter as cf
import ODR as cf1
start=1
end=5
def func(x,a,b,c):
	return a*x**2+b*x+c

onemass = np.load('onemass_4b0_%d-%d.npy'%(start,end))
twomass = np.load('twomass_4b0_%d-%d.npy'%(start,end))
# print(onemass)
l = onemass.shape[1]
w = onemass.shape[0]

m = np.mean(onemass,axis=1)
mstd = np.std(onemass,axis=1)/np.sqrt(l-1)
be = 2*onemass - twomass
bindingE = np.mean(be,axis=1)
bstd = np.std(be,axis=1)/np.sqrt(l-1)

resamplemean = np.zeros((w,l))
resamplemass = np.zeros((w,l))
for i in range(w):
	resamplemean[i,:] = cf.jackkniferesamplemean(be[i,:])
	resamplemass[i,:] = cf.jackkniferesamplemean(onemass[i,:])
fitstart = 9
fitend = 13
part = resamplemean[fitstart:fitend,:]
mpart = resamplemass[fitstart:fitend,:]
mmcov=np.cov(mpart,bias=True)*(l-1)
mcov = np.cov(part,bias=True)*(l-1)
# print(mcov)
# print(mcov.shape)


# fittedpara,allchisq = cf.fit(func,part,mcov,l,fitstart,fitend,3,[1,1,1],m[fitstart:fitend])
fittedpara = np.zeros((3,l))
allchisq = np.zeros(l)
for i in range(l):
	fitter = cf1.Fitter()
	fitter.Model(func)
	fitter.Beta0([1,1,1])
	fitter.Data(m[fitstart:fitend],part[:,i],covx=mmcov,covy=mcov)
	fitter.Run()
	fittedpara[:,i] = fitter.out.x[0:3]
	allchisq[i] = fitter.chi2
# print(fittedpara)

finalpara = np.mean(fittedpara,axis=1)
print(allchisq)
x1 = np.arange(0,0.2,0.01)
y1 = func(x1,finalpara[0],finalpara[1],finalpara[2])

fig = plt.figure()
plt.errorbar(m,bindingE,xerr=mstd,yerr=bstd,marker="*",capsize=0.5)
plt.plot(x1,y1,label='fitted curve')
plt.grid()
plt.xlabel('Bare Mass')
plt.ylabel('Binding Energy')
plt.show()