import numpy as np
import matplotlib.pyplot as plt
import correlatefitter as cf
from pathlib import Path
maxlen = 50
corperconfig = 20
binsize = 20
start=2
end=9
numpara = 3
def fittedfunc(x,a,m,b):
	return a - m * x - b*np.log(x)
for m in np.arange(0.001,0.051,0.001):
	# totalcor,allpara1,parastderr1,allchisq1 = cf.fits(fittedfunc,Path('./correlatordata/alltwoparticlecorrelators_m0=%f.npy'%m),maxlen,corperconfig,binsize,m,start,end,numpara,'log',0)
	# numfile = len(totalcor[1,:])
	plt.figure()
	# x = np.arange(0,50,1)
	# for i in range(numfile):
	# 	y = totalcor[:,i]
	# 	plt.plot(x,y)
	# plt.xlabel('Distance')
	# plt.ylabel('One-Particle-Correlator')
	# plt.savefig('correlators4b0.png',dpi=130)
	# plt.clf()
	# x = np.arange(1,10,1)
	# for i in range(numfile):
	# 	y = np.log(totalcor[1:10,i])
	# 	plt.plot(x,y)
	# plt.xlabel('Distance')
	# plt.ylabel('Log of Correlator')
	# plt.savefig('correlators4b0log1-10.png',dpi=130)
	# plt.clf()
	# x = np.arange(1,25,1)
	# for i in range(numfile):
	# 	y = np.log(totalcor[1:25,i])
	# 	plt.plot(x,y)
	# plt.xlabel('Distance')
	# plt.ylabel('Log of Correlator')
	# plt.savefig('correlators4b0log1-25.png',dpi=130)
	# plt.clf()

	# totalcor,finalpara1,parastderr1,finalchisq1 = cf.fits(fittedfunc,Path('./correlatordata/allcorrelators_m0=%f.npy'%m),maxlen,corperconfig,binsize,m,start,end,numpara,'log',1)

	# nonzero1 = 10
	# x = np.arange(1,nonzero1)
	# #y = fittedfunc2(x,finalpara[0],finalpara[1],finalpara[2])
	# y = np.exp(fittedfunc(x,finalpara1[0],finalpara1[1],finalpara1[2]))
	# y1 = np.mean(totalcor[1:nonzero1],axis = 1)
	# err = np.std(totalcor[1:nonzero1],axis = 1)/np.sqrt(len(totalcor[0,:])-1)
	# #err1 = err[0:nonzero]
	# plt.errorbar(x,y1,yerr = err,linewidth = 1, capsize = 1)
	# #plt.yscale('log')
	# plt.plot(x,y,'red')
	# plt.grid(True)
	# plt.xlabel('Geodesic Distance')
	# plt.ylabel('Correlator')
	# # plt.title('Distance %d to %d'%(start,end))
	# plt.title('m_0=0.001, Ren. mass = %6f, distance %d to %d'%(finalpara1[1],start,end))
	# plt.savefig('cor2-9.png',dpi=130)
	# plt.clf()
	start = 1
	end = 5
	totalcor,finalpara1,parastderr1,finalchisq1 = cf.fits(fittedfunc,Path('./correlatordata/alltwoparticlecorrelators_m0=%f.npy'%m),maxlen,corperconfig,binsize,m,start,end,numpara,'log',1)
	nonzero1 = 10
	nonzero2 = 30
	x = np.arange(1,nonzero1)
	x1 = np.arange(1,nonzero2)
	#y = fittedfunc2(x,finalpara[0],finalpara[1],finalpara[2])
	y = np.exp(fittedfunc(x,finalpara1[0],finalpara1[1],finalpara1[2]))
	y1 = np.mean(totalcor[1:nonzero2],axis = 1)
	err = np.std(totalcor[1:nonzero2],axis = 1)/np.sqrt(len(totalcor[0,:])-1)
	
	#err1 = err[0:nonzero]
	plt.errorbar(x1,y1,yerr = err,linewidth = 1, capsize = 1)
	#plt.yscale('log')
	plt.plot(x,y,'red')
	plt.grid(True)
	plt.xlabel('Geodesic Distance')
	plt.ylabel('Correlator')
	# plt.title('Distance %d to %d'%(start,end))
	plt.title('m_0=%6f, Ren. mass = %6f, distance %d to %d'%(m,finalpara1[1],start,end))
	plt.savefig('befigure15/two%6f.png'%m)
	plt.close()
	# plt.savefig('cor3-8.png',dpi=130)




