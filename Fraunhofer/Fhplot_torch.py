import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('xtick', direction='in', top=True)
mpl.rc('ytick', direction='in', right=True)
mpl.rc('xtick.minor', visible=True)
mpl.rc('ytick.minor', visible=True)
plt.rcParams['font.size'] = 13
plt.rcParams['figure.figsize'] = [8,6]
import datetime as dt
import sys
import Fhcalc_torch
#import tensorflow as tf
import torch

if (len(sys.argv)!=4):
    print('Wrong inputs!')
    print('Usage: python FRBplot_torch.py output_filename N_screen N_theta_interp')
    sys.exit()

name = sys.argv[1] #file name
N_screen = np.int(sys.argv[2])
N_theta = np.int(sys.argv[3])

starttime = dt.datetime.now()
screen = {}
screen['N'] = N_screen
screen['D'] = 10
Fhcalc_torch.Initialize_tf(screen)
center = (screen['D']/2, screen['D']/2)
Fhcalc_torch.MultByGaussian_tf(screen, center, 1.0)
Fhcalc_torch.InCircle_tf(screen, center, 2.0)
Fhcalc_torch.ScreenFFT_tf(screen)

lam = 0.002 #mm wavelength
kphot = 2*torch.pi/lam
thetamaxdeg = 2.0
thetamax = 2*torch.pi/180. # 2 degrees in radians
theta_vec = torch.linspace(-thetamax,thetamax,N_theta) 
II0 = Fhcalc_torch.Project_I_on_thetagrid_tf(theta_vec,screen,lam) 
# fig = plt.figure(figsize=(8,6))
# plt.imshow(screen['E'],extent=(0,screen['D'],screen['D'],0))
# plt.colorbar()
endtime = dt.datetime.now()
t = (endtime-starttime).seconds
print('Time = {}s'.format(t))
#plt.savefig('/home/gemma/Fraunhofer/{}_E0.png'.format(name))


starttime = dt.datetime.now()
screen = {}
screen['N'] = N_screen
screen['D'] = 10
Fhcalc_torch.Initialize_tf(screen)
center = (screen['D']/2, screen['D']/2)
Fhcalc_torch.MultByGaussian_tf(screen, center, 1.0) #sigma 1.0
Fhcalc_torch.InCircle_tf(screen, center, 2.0) #radius 2.0
c2 = (center[0]-0.2, center[1]-0.1)
#c2 = (center[0], center[1])
Fhcalc_torch.CircleAtten_tf(screen,c2,0.1,1.0*np.exp(np.pi * 1j))
#c2 = (center[0]+0.1, center[1]+0.3)
#FRBcalc.CircleAtten(screen,c2,0.1,1/1.3)
#c2 = (center[0]+0.3, center[1]-0.1)
#FRBcalc.CircleAtten(screen,c2,0.1,1.3)
#c2 = (center[0]-0.3, center[1]+0.2)
#FRBcalc.CircleAtten(screen,c2,0.1,1/1.3)
Fhcalc_torch.ScreenFFT_tf(screen)

lam = 0.002 #mm wavelength
kphot = 2*torch.pi/lam
thetamaxdeg = 2.0
thetamax = 2*torch.pi/180. # 2 degrees in radians
theta_vec = torch.linspace(-thetamax,thetamax,N_theta) 
II1 = FRBcalc_torch.Project_I_on_thetagrid_tf(theta_vec,screen,lam) 
# fig = plt.figure(figsize=(8,6))
# plt.imshow(torch.abs(screen['E']),extent=(0,screen['D'],screen['D'],0))
# plt.colorbar()
endtime = dt.datetime.now()
t = (endtime-starttime).seconds
print('Time = {}s'.format(t))
#plt.savefig('/home/gemma/Fraunhofer/{}_E1.png'.format(name))


# Idiff = II1 - II0
# fig = plt.figure(figsize=(16,6))
# plt.subplot(1,3,1)
# plt.imshow(torch.log10(torch.abs(II0)),interpolation = None,extent=[-thetamaxdeg,thetamaxdeg,thetamaxdeg,-thetamaxdeg])
# plt.colorbar()
# plt.subplot(1,3,2)
# plt.imshow(torch.log10(torch.abs(II1)),interpolation = None,extent=[-thetamaxdeg,thetamaxdeg,thetamaxdeg,-thetamaxdeg])
# plt.colorbar()
# plt.subplot(1,3,3)
# plt.imshow(torch.log10(torch.abs(Idiff)),interpolation = None,extent=[-thetamaxdeg,thetamaxdeg,thetamaxdeg,-thetamaxdeg])
# plt.colorbar()
# plt.savefig('/home/gemma/{}_I.png'.format(name))

