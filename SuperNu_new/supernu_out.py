from numpy import *
import numpy as np
import math
import pandas as pd
import sys
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pylab import *

#----------
#Users Guide:
#First input: time in days (float) or 'max'
#to plot spectrum during maximum light
#Second input: type 'animate' to produce
#movie of spectral evolution. Any other
#string will not generate a movie

#Control spectral evolution animation y-scale:
#If set true animation is scaled in y-axis
#based on peak spectrum
ScaleToPeakLum = False

#----------
#Constants
cm2A = 1.e+08 # cm to angstroms
day  = 86400. # days to seconds

#----------
#Use pandas to read output.flx_grid file as CSV
data = pd.read_csv('output.flx_grid')
#Extract data for wavelength, mu(=cos(theta)), polar angle phi
#and timestep
[wl,mu,phi,ts] = [data.iloc[:,0].values,data.iloc[:,0].values,data.iloc[:,0].values,data.iloc[:,0].values]

#Truncate previous lines
ts = ts[3:]
for i in range(0,len(ts)):
	ts[i] = eval(ts[i])/day
#Split long strings into a lists of strings
wl = wl[0].split()
mu = mu[1].split()
phi = phi[2].split()

#Convert strings to floats with eval
#and convert from cm to Angstroms
for i in range(0,len(wl)):
	wl[i] = eval(wl[i])*cm2A
	
#Convert wavelength array to cell-centered values
wl2=[wl[0]/2.0,]
for i in range(1,len(wl)-1):
	wl2.append(i)
	wl2[i] = wl[i-1]+(wl[i]-wl[i-1])/2.0
	
for i in range(0,len(mu)):
	mu[i] = eval(mu[i])	
	
#Convert mu(=cos(theta)) array to cell-centered values
mu2=[mu[0]/2.0,]
for i in range(1,len(mu)-1):
	mu2.append(i)
	mu2[i] = mu[i-1]+(mu[i]-mu[i-1])/2.0	
	
for i in range(0,len(phi)):
	phi[i] = eval(phi[i])	
	
#Convert polar angle (phi) array to cell-centered values
phi2=[phi[0]/2.0,]
for i in range(1,len(phi)-1):
	phi2.append(i)
	phi2[i] = phi[i-1]+(phi[i]-phi[i-1])/2.0
	
#Convert time array to cell-centered values
ts2=[ts[0]/2.0,]
for i in range(1,len(ts)-1):
	ts2.append(i)
	ts2[i] = ts[i-1]+(ts[i]-ts[i-1])/2.0		
	
#-----------
#Use pandas to read output.grd_grid file
data2 = pd.read_csv('output.flx_luminos',header=None)

#Save spectral data in array
#each row corresponds to array of the spectrum of that
#particular timestep (same index as timesteps index)
sp=[]
for i in range(0,len(data2)):
	sp.append(i)
	sp[i] = data2.iloc[i,:].values[0].split()
	
#Append spectral data in 2D array where the rows
#correspond to different times and the columns
#to different luminosities and calculate total
#luminosity in each time-step (LC)
lum=[]
sp2 = [[0 for y in range(len(wl2))] for x in range(len(ts2))]
for i in range(0,len(ts2)):
	lum.append(i)
	for j in range(0,len(wl2)):
		sp2[i][j] = eval(sp[i][j])
		
	lum[i] = sum(sp2[i][:])
	
#Figure out if user wants peak spectrum or a spectrum
#at a specific day
diff2 = [(abs(max(lum) - x),idx) for (idx,x) in enumerate(lum)]
diff2.sort()
maxL_ind = diff2[0][1]
tmax = ts2[maxL_ind]

if sys.argv[1] == 'max':
	time = tmax
else:
	time = eval(sys.argv[1])	
	
#Plot selected spectrum
#First find index in ts2 array
#that corresponds to input time (in days)
diff = [(abs(time - x),idx) for (idx,x) in enumerate(ts2)]
diff.sort()
index = diff[0][1]

#Plot spectrum for input time
plt.title('Spectrum at '+str(ts2[index])+' days')
plt.xlabel('Wavelength [A]')
plt.ylabel('Luminosity [erg/s]')
plt.plot(wl2,sp2[index][:],'r',linewidth=2.0,label='Spectrum')
legend = plt.legend(loc='upper right', shadow=True)
plt.xlim(0.,10000.)
plt.minorticks_on()
plt.show()

#Plot lightcurve
plt.title('Full lightcurve')
plt.xlabel('Time since explosion [d]')
plt.ylabel('Luminosity [erg/s]')
plt.plot(ts2,lum,'k',marker='o',markersize=5.0,linestyle='none',label='Model')
legend = plt.legend(loc='upper right', shadow=True)
plt.minorticks_on()
plt.show()

#Make animation of spectral evolution if 
#desired by the user
anim = sys.argv[2]

if anim != 'animate': sys.exit()
	
elif anim == 'animate':
	fig,ax = plt.subplots()
	line, = ax.plot(wl2,sp2[i][:])

	xlabel('Wavelength [A]')
	ylabel('Luminosity [erg/s]')
	xlim(1000.,10000.)
	if ScaleToPeakLum:
		ylim(0.,max(sp2[maxL_ind][:]))
	minorticks_on()

	def animate(i):
		line.set_xdata(wl2)
		line.set_ydata(sp2[i][:])
		title('Spectral Evolution, t = '+str(ts2[i])+' days')
		if not ScaleToPeakLum:
			ylim(0.,max(sp2[i][:]))
		return line
	
	ani = animation.FuncAnimation(fig, animate, interval=1,repeat = False,blit=False)
	ani.save('SpecEvol.mp4',fps=7)	
	



	

