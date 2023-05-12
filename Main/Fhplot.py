### uses Fhcalc.py 
### generates a perfectly Gaussian illuminated aperture, truncated by a circle, place among a square screen.
### interpolates the beam and generates an on-sky beam

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
import Fhcalc

if (len(sys.argv)!=4):
    print('Wrong inputs!')
    print('Usage: python Fhplot.py output_filename N_screen N_theta_interp')
    sys.exit()

# name of the output file
name = sys.argv[1]
# number of pixels in the E-field screen - N_screen by N_screen
N_screen = np.int(sys.argv[2])
# number of points to interpolate - N_theta by N_theta
N_theta = np.int(sys.argv[3])

# starttime = dt.datetime.now() # timeit
screen = {}
screen['N'] = N_screen
screen['D'] = 10
Fhcalc.Initialize(screen)
center = (screen['D']/2, screen['D']/2)
Fhcalc.MultByGaussian(screen, center, 1.0)
Fhcalc.InCircle(screen, center, 2.0)
Fhcalc.ScreenFFT(screen)

lam = 0.002 #mm wavelength
kphot = 2*np.pi/lam
thetamaxdeg = 2.0
thetamax = thetamaxdeg*np.pi/180. # 2 degrees in radians
theta_vec = np.linspace(-thetamax,thetamax,N_theta) 
II0 = Fhcalc.Project_I_on_thetagrid(theta_vec,screen,lam)
theta_vec = np.linspace(0, 2*thetamax, N_theta) 
fig = plt.figure(figsize=(12,5))
plt.subplot(1,2,1)
plt.imshow(screen['E'],extent=(0,screen['D'],screen['D'],0))
plt.colorbar()
plt.xlabel('X (m)')
plt.ylabel('Y (m)')
plt.subplot(1,2,2)
plt.imshow(II0, extent=(-thetamax*180/np.pi,thetamax*180/np.pi,thetamax*180/np.pi, -thetamax*180/np.pi))
plt.xlabel('degree')
plt.ylabel('degree')
plt.colorbar()
#plt.show()
# endtime = dt.datetime.now()
# t = (endtime-starttime).seconds
# print('Time = {:.2f}s'.format(t))
plt.savefig('./{}_E0_I0.png'.format(name))