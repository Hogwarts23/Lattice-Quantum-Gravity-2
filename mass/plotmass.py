import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
start = 1
end = 5

onemass = np.load('onemass_4b0_%d-%d.npy'%(start,end))
twomass = np.load('twomass_4b0_%d-%d.npy'%(start,end))
alljm = np.load('renorm_masses_sqcorr_bin20_4kb0_n11500_m00p001-mf0p05_be0-f2-7_rm0-f2-7.npy')
m0 = np.arange(0.001,0.051,0.001)
m = np.mean(onemass,axis=1)
print(len(onemass[1,:]))

M = np.mean(twomass,axis=1)
binding = 2*m - M
fig = plt.figure()
l1=plt.plot(m0,m,'blue',label='1-particle-mass')
l2=plt.plot(m0,M,'green',label='2-particle-mass')
l3=plt.plot(m0,binding,'red',label='binding energy')
plt.grid(ls='--')
plt.legend()
plt.xlabel('bare mass m_0')
plt.ylabel('mass parameters')
plt.title('E vs m_0 (ensemble 4b0, distance 3 to 8)')
plt.savefig('mass4b0_3-8.png')


# jm=np.mean(alljm,axis = 1)
# fig = plt.figure()
# plt.plot(m0,m,label='My mass')
# plt.plot(m0,jm,label='Judah\'s mass')
# plt.grid(ls='--')
# plt.legend()
# plt.xlabel('bare mass m_0')
# plt.ylabel('mass parameters')
# plt.title('Comparison with Judah\'s results(4b0)')
# plt.savefig('compareresults.png')
