import h5py
import numpy as np
import matplotlib.pyplot as plt
import time
import math
import sys
from decimal import Decimal
from numpy import *


vsize = 120

nx = vsize
ny = 1
nz = 1
nabund = 31 # 32
ncol = nabund + 3

CV = 1.247223E+08
eV2K = 11604.5221

#Input file creation
#is adjusted to correspond
#to the elements present in FLASH'
#Aprox19 network

print '# spherical'
print '#       nr      nabund'
print '#       '+str(nx)+'          '+str(ny)+' '+str(nz)+' '+str(ncol)+' '+str(nabund)
print '# rightvel        mass           temp           h          he          li          be           b           c           n           o           f          ne          na          mg          al          si           p           s          cl          ar           k          ca          sc          ti           v          cr          mn          fe          co          ni        cu       zn       ni56'
#print '# rightvel        mass           temp          h          he          c          n           o           ne           mg           si           s          ar          ca          ti          cr          fe           ni          ni56'
#print '# rightvel        mass            temp           h          he          c          n           o           ne           mg           si           s          ar          ca          ti          cr          fe           ni'

'''
file = h5py.File('pisn_341d_1D_lref_5_exp_0033_hdf5_chk_0025','r') 

density = file['dens']
velocity = file['velx']
temperature = file['temp']
coord = file['bounding box']
node = file['node type']
ref = file['refine level']
blksize=file['block size']
#Elements
prot_i = file['prot']
neut_i = file['neut']
h1_i = file['h1  ']
he3_i = file['he3 ']
he4_i = file['he4 ']
c12_i = file['c12 ']
n14_i = file['n14 ']
o16_i = file['o16 ']
ne20_i = file['ne20']
mg24_i = file['mg24']
si28_i = file['si28']
s32_i  = file['s32 ']
ar36_i = file['ar36']
ca40_i = file['ca40']
ti44_i = file['ti44']
cr48_i = file['cr48']
fe52_i = file['fe52']
fe54_i = file['fe54']
ni56_i = file['ni56']

nblocks = len(coord)
nzones = density.shape[3]

#print nzones,nblocks

l = -1
dx = []
rad = []
rho = []
temp = []
blk = []
vel = []
#Element arrays
prot = []
neut = []
h1 = []
he3 = []
he4 = []
c12 = []
n14 = []
o16 = []
ne20 = []
mg24 = []
si28 = []
s32  = []
ar36 = []
ca40 = []
ti44 = []
cr48 = []
fe52 = []
fe54 = []
ni56 = []
#First determine zone center coordinates using block boundary coordinates
for i in range(0,nblocks):
    blk.append(i)
    blk[i] = blksize[i][0]
    dd=(blk[i]/16)
    dd = dd/2
    if node[i] == 1:
        for j in range(0,nzones):
            l = l + 1
            rad.append(l)
            rho.append(l)
            vel.append(l)
	    temp.append(l)
            
            prot.append(l)
            neut.append(l)
            h1.append(l)
            he3.append(l)
            he4.append(l)
            c12.append(l)
            n14.append(l)
            o16.append(l)
            ne20.append(l)
            mg24.append(l)
            si28.append(l)
            s32.append(l)
            ar36.append(l)
            ca40.append(l)
            ti44.append(l)
            cr48.append(l)
            fe52.append(l)
            fe54.append(l)
            ni56.append(l)            
            
            rad[l] = coord[i][0][0]+(2*j+1)*dd
            rho[l] = density.value[i][0][0][j]
            vel[l] = velocity.value[i][0][0][j]
	    temp[l] = temperature.value[i][0][0][j]
            prot[l] = prot_i.value[i][0][0][j]
            neut[l] = neut_i.value[i][0][0][j]
            h1[l] = h1_i.value[i][0][0][j]
            he3[l] = he3_i.value[i][0][0][j]
            he4[l] = he4_i.value[i][0][0][j]
            c12[l] = c12_i.value[i][0][0][j]
            n14[l] = n14_i.value[i][0][0][j]
            o16[l] = o16_i.value[i][0][0][j]
            ne20[l] = ne20_i.value[i][0][0][j]
            mg24[l] = mg24_i.value[i][0][0][j]
            si28[l] = si28_i.value[i][0][0][j]
            s32[l] = s32_i.value[i][0][0][j]
            ar36[l] = ar36_i.value[i][0][0][j]
            ca40[l] = ca40_i.value[i][0][0][j]
            ti44[l] = ti44_i.value[i][0][0][j]
            cr48[l] = cr48_i.value[i][0][0][j]
            fe52[l] = fe52_i.value[i][0][0][j]
            fe54[l] = fe54_i.value[i][0][0][j]
            ni56[l] = ni56_i.value[i][0][0][j]
	    

h_tot  = []
he_tot = []
fe_tot = []
for i in range(0,len(rad)):
    h_tot.append(i)
    he_tot.append(i)
    fe_tot.append(i)
    h_tot[i] = neut[i]+prot[i]+h1[i]
    he_tot[i] = he3[i]+he4[i]
    fe_tot[i] = fe52[i]+fe54[i]
    #print rad[i],rho[i],temp[i],vel[i],h_tot[i],he_tot[i],c12[i],n14[i],o16[i],ne20[i],mg24[i],si28[i],s32[i],ar36[i],ca40[i],ti44[i],cr48[i],fe_tot[i],fe52[i],ni56[i]
'''
rad,rho,temp,vel,h_tot,he_tot,c12,n14,o16,ne20,mg24,si28,s32,ar36,ca40,ti44,cr48,fe_tot,fe52,ni56=loadtxt('P250_RAGE_8hrs_trunc.dat',usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19),unpack=True,skiprows=1)
#temp = temp * eV2K
eint = CV * temp


#Calculate shell volumes
vol = [(4./3.)*math.pi*rad[0]**3,]
mass = [vol[0]*rho[0],]
for i in range(1,len(rad)):
    vol.append(i)
    mass.append(i)
    vol[i] = (4./3.)*math.pi*(rad[i]**3-rad[i-1]**3)
    mass[i] = rho[i]*vol[i]

#Define element mass list
#m_prot = []
#m_neut = []
#m_h1 = []
m_hyd = []
#m_he3 = []
#m_he4 = []
m_hel = []
m_c12 = []
m_n14 = []
m_o16 = []
m_ne20 = []
m_mg24 = []
m_si28 = []
m_s32 = []
m_ar36 = []
m_ca40 = []
m_ti44 = []
m_cr48 = []
m_fe52 = []
#m_fe54 = []
m_fe = []
m_ni = []
m_ni56 = []
for i in range(0,len(rad)):
   # m_prot.append(i)
   # m_neut.append(i)
   # m_h1.append(i)
    m_hyd.append(i)
   # m_he3.append(i)
   # m_he4.append(i)
    m_hel.append(i)
    m_c12.append(i)
    m_n14.append(i)
    m_o16.append(i)
    m_ne20.append(i)
    m_mg24.append(i)
    m_si28.append(i)
    m_s32.append(i)
    m_ar36.append(i)
    m_ca40.append(i)
    m_ti44.append(i)
    m_cr48.append(i)
    m_fe52.append(i)
  #  m_fe54.append(i)
    m_fe.append(i)
    m_ni.append(i)
    m_ni56.append(i)
    
   # m_prot[i] = mass[i]*prot[i]
   # m_neut[i] = mass[i]*neut[i]
   # m_h1[i] = mass[i]*h1[i]
    m_hyd[i] = mass[i]*h_tot[i]
   # m_he3[i] = mass[i]*he3[i]
   # m_he4[i] = mass[i]*he4[i]
    m_hel[i] = mass[i]*he_tot[i]
    m_c12[i] = mass[i]*c12[i]
    m_n14[i] = mass[i]*n14[i]
    m_o16[i] = mass[i]*o16[i]
    m_ne20[i] = mass[i]*ne20[i]
    m_mg24[i] = mass[i]*mg24[i]
    m_si28[i] = mass[i]*si28[i]
    m_s32[i] = mass[i]*s32[i]
    m_ar36[i] = mass[i]*ar36[i]
    m_ca40[i] = mass[i]*ca40[i]
    m_ti44[i] = mass[i]*ti44[i]
    m_cr48[i] = mass[i]*cr48[i]
    m_fe52[i] = mass[i]*fe52[i]
  #  m_fe54[i] = mass[i]*fe54[i]
    m_fe[i] = mass[i]*fe_tot[i]
    m_ni56[i] = mass[i]*ni56[i]
    m_ni[i] = m_ni56[i]

#Generate homologous velocity grid
min_vel = 0
max_vel = max(vel)
vel_step = (max_vel-min_vel)/vsize
vel_grid = [vel_step,]
for i in range(1,vsize):
    vel_grid.append(i)
    vel_grid[i] = vel_grid[i-1]+vel_step
    
    
mass2 = []
m_h_2 = []
m_he_2 = []
m_c_2 = []
m_n_2 = []
m_o_2 = []
m_ne_2 = []
m_mg_2 = []
m_si_2 = []
m_s_2 = []
m_ar_2 = []
m_ca_2 = []
m_ti_2 = []
m_cr_2 = []
m_fe_2 = []
m_ni_2 = []
m_fe52_2 = []
m_ni56_2 = []
X_h = []
X_he = []
X_c = []
X_n = []
X_o = []
X_ne = []
X_mg = []
X_si = []
X_s = []
X_ar = []
X_ca = []
X_ti = []
X_cr = []
X_fe = []
X_ni = []
X_fe52 = []
X_ni56 = []
n_bin = []

mass2 = []
eint2 = []
temp2 = []

for i in range(0,len(vel_grid)):
    mass2.append(i)
    eint2.append(i)
    temp2.append(i)
    m_h_2.append(i)
    m_he_2.append(i)
    m_c_2.append(i)
    m_n_2.append(i)
    m_o_2.append(i)
    m_ne_2.append(i)
    m_mg_2.append(i)
    m_si_2.append(i)
    m_s_2.append(i)
    m_ar_2.append(i)
    m_ca_2.append(i)
    m_ti_2.append(i)
    m_cr_2.append(i)
    m_fe_2.append(i)
    m_ni_2.append(i)
    m_fe52_2.append(i)
    m_ni56_2.append(i) 

    mass2[i] = 0.
    m_h_2[i] = 0.
    m_he_2[i] = 0.
    m_c_2[i] = 0.
    m_n_2[i] = 0.
    m_o_2[i] = 0.
    m_ne_2[i] = 0.
    m_mg_2[i] = 0.
    m_si_2[i] = 0.
    m_s_2[i] = 0.
    m_ar_2[i] = 0.
    m_ca_2[i] = 0.
    m_ti_2[i] = 0.
    m_cr_2[i] = 0.
    m_fe_2[i] = 0.
    m_ni_2[i] = 0.
    m_fe52_2[i] = 0.
    m_ni56_2[i] = 0.

    X_h.append(i)
    X_he.append(i)
    X_c.append(i)
    X_n.append(i)
    X_o.append(i)
    X_ne.append(i)
    X_mg.append(i)
    X_si.append(i)
    X_s.append(i)
    X_ar.append(i)
    X_ca.append(i)
    X_ti.append(i)
    X_cr.append(i)
    X_fe.append(i)
    X_ni.append(i)
    X_fe52.append(i)
    X_ni56.append(i)
    n_bin.append(i)

for i in range(0,len(vel)):
    for j in range(1,len(vel_grid)):
        if vel_grid[j-1] < vel[i] <= vel_grid[j]:
            mass2[j] = mass2[j]+mass[i]
	    eint2[j] = eint2[j]+eint[i]
	    temp2[j] = eint2[j]/CV
            m_h_2[j] = m_h_2[j]+m_hyd[i]
            m_he_2[j] = m_he_2[j] + m_hel[i]
            m_c_2[j] = m_c_2[j] + m_c12[i]
            m_n_2[j] = m_n_2[j] + m_n14[i]
            m_o_2[j] = m_o_2[j] + m_o16[i]
            m_ne_2[j] = m_ne_2[j] + m_ne20[i]
            m_mg_2[j] = m_mg_2[j] + m_mg24[i]
            m_si_2[j] = m_si_2[j] + m_si28[i]
            m_s_2[j] = m_s_2[j] + m_s32[i]
            m_ar_2[j] = m_ar_2[j] + m_ar36[i]
            m_ca_2[j] = m_ca_2[j] + m_ca40[i]
            m_ti_2[j] = m_ti_2[j] + m_ti44[i]
            m_cr_2[j] = m_cr_2[j] + m_cr48[i]
            m_fe_2[j] = m_fe_2[j] + m_fe[i]
            m_ni_2[j] = m_ni_2[j] + m_ni[i]
            m_fe52_2[j] = m_fe52_2[j] + m_fe52[i]
            m_ni56_2[j] = m_ni56_2[j] + m_ni56[i]

    if vel[i] < vel_grid[0]:
        mass2[0] = mass2[0] + mass[i]
	eint2[0] = eint2[0] + eint[i]
	temp2[0] = eint2[0]/CV
        m_h_2[0] = m_h_2[0]+m_hyd[i]
        m_he_2[0] = m_he_2[0] + m_hel[i]
        m_c_2[0] = m_c_2[0] + m_c12[i]
        m_n_2[0] = m_n_2[0] + m_n14[i]
        m_o_2[0] = m_o_2[0] + m_o16[i]
        m_ne_2[0] = m_ne_2[0] + m_ne20[i]
        m_mg_2[0] = m_mg_2[0] + m_mg24[i]
        m_si_2[0] = m_si_2[0] + m_si28[i]
        m_s_2[0] = m_s_2[0] + m_s32[i]
        m_ar_2[0] = m_ar_2[0] + m_ar36[i]
        m_ca_2[0] = m_ca_2[0] + m_ca40[i]
        m_ti_2[0] = m_ti_2[0] + m_ti44[i]
        m_cr_2[0] = m_cr_2[0] + m_cr48[i]
        m_fe_2[0] = m_fe_2[0] + m_fe[i]
        m_ni_2[0] = m_ni_2[0] + m_ni[i]
        m_fe52_2[0] = m_fe52_2[0] + m_fe52[i]
        m_ni56_2[0] = m_ni56_2[0] + m_ni56[i]    


#Ni_mass = []
for i in range(0,len(vel_grid)):
    #Ni_mass.append(i)
    
    X_h[i] = '%.4E' % Decimal(m_h_2[i]/mass2[i])
    X_he[i] = '%.4E' % Decimal(m_he_2[i]/mass2[i])
    X_c[i] = '%.4E' % Decimal(m_c_2[i]/mass2[i])
    X_n[i] = '%.4E' % Decimal(m_n_2[i]/mass2[i])
    X_o[i] = '%.4E' % Decimal(m_o_2[i]/mass2[i])
    X_ne[i] = '%.4E' % Decimal(m_ne_2[i]/mass2[i])
    X_mg[i] = '%.4E' % Decimal(m_mg_2[i]/mass2[i])
    X_si[i] = '%.4E' % Decimal(m_si_2[i]/mass2[i])
    X_s[i]  = '%.4E' % Decimal(m_s_2[i]/mass2[i])
    X_ar[i] = '%.4E' % Decimal(m_ar_2[i]/mass2[i])
    X_ca[i] = '%.4E' % Decimal(m_ca_2[i]/mass2[i])
    X_ti[i] = '%.4E' % Decimal(m_ti_2[i]/mass2[i])
    X_cr[i] = '%.4E' % Decimal(m_cr_2[i]/mass2[i])
    X_fe[i] = '%.4E' % Decimal(m_fe_2[i]/mass2[i])
    X_ni[i] = '%.4E' % Decimal(m_ni_2[i]/mass2[i])
    X_fe52[i] = '%.4E' % Decimal(m_fe52_2[i]/mass2[i])
    X_ni56[i] = '%.4E' % Decimal(m_ni56_2[i]/mass2[i])

    #Ni_mass[i] = (X_ni[i]*mass2[i])/1.99e+33
        

#print sum(mass)/1.99e+33
#print sum(Ni_mass)

for i in range(0,len(vel_grid)):
	mass2[i]='%.4E' % Decimal(mass2[i])
	vel_grid[i] = '%.4E' % Decimal(vel_grid[i])
    
for i in range(0,len(vel_grid)):
    #print vel_grid[i],'',mass2[i],'',temp2[i],'',X_h[i],'',X_he[i],'',X_c[i],'',X_n[i],'',X_o[i],'',X_ne[i],'',X_mg[i],'',X_si[i],'',X_s[i],'',X_ar[i],'',X_ca[i],'',X_ti[i],'',X_cr[i],'',X_fe[i],'',X_ni[i],'',X_ni56[i]#,'',X_fe52[i],'',X_ni56[i]
    print vel_grid[i],mass2[i],temp2[i],X_h[i],X_he[i],0.0,0.0,0.0,X_c[i],X_n[i],X_o[i],0.0,X_ne[i],0.0,X_mg[i],0.0,X_si[i],0.0,X_s[i],0.0,X_ar[i],0.0,X_ca[i],0.0,X_ti[i],0.0,X_cr[i],0.0,X_fe[i],0.0,X_ni[i],0.0,0.0,X_ni56[i]

