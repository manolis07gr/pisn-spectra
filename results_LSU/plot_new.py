import numpy as np
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

#t1,L1 = loadtxt('full_lightcurve1D.dat',usecols=(0,1),unpack=True,skiprows=0)
#t2,L2 = loadtxt('full_lightcurve_2D1DAA.dat',usecols=(0,1),unpack=True,skiprows=0)
#t3,L3 = loadtxt('full_lightcurve_2D1DAASNEC.dat',usecols=(0,1),unpack=True,skiprows=0)
#t4,L4 = loadtxt('full_lightcurve_SNEC.dat',usecols=(0,1),unpack=True,skiprows=0)

l1,F1 = loadtxt('spectrum_158d_1D.dat',usecols=(0,1),unpack=True,skiprows=0)
l2,F2 = loadtxt('spectrum_169d_2D1DAA.dat',usecols=(0,1),unpack=True,skiprows=0)
l3,F3 = loadtxt('spectrum_167d_2D1DAASNEC.dat',usecols=(0,1),unpack=True,skiprows=0)

fig, ax = plt.subplots()
#plt.plot(t1,L1,label='1D spherical',linewidth=2.0,color='k')
#plt.plot(t2,L2,label='2D angle-averaged',linewidth=2.0,color='r')
#plt.plot(t3,L3,label='2D angle-averaged (SNEC evol)',linewidth=2.0,color='b')
#plt.plot(t4,L4,label='SNEC',linewidth=2.0,color='g')

plt.plot(l1,F1,label='1D spherical',linewidth=2.0,color='k')
plt.plot(l2,F2,label='2D angle-averaged',linewidth=2.0,color='r')
plt.plot(l3,F3,label='2D angle-averaged (SNEC evol)',linewidth=2.0,color='b')

#plt.xlim([0,500])
#plt.ylim([0.001e+44,2e+44])

plt.xlim([3000,12000])
#plt.ylim([0.001e+44,2e+44])

plt.xlabel("t [days]",fontsize=14)
plt.ylabel("L [erg/s] " +  r'$\left(\mathregular{10^{44}}\right)$',fontsize=14)
ax.yaxis.offsetText.set_visible(False)
plt.minorticks_on()
legend = plt.legend(loc='upper right', shadow=False)
plt.show()




