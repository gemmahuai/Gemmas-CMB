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
    print('Usage: python FRBplot.py output_filename N_screen N_theta_interp')
    sys.exit()

name = sys.argv[1] #file name
N_screen = np.int(sys.argv[2])
N_theta = np.int(sys.argv[3])

starttime = dt.datetime.now()
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
thetamax = 2*np.pi/180. # 2 degrees in radians
theta_vec = np.linspace(-thetamax,thetamax,N_theta) 
II0 = Fhcalc.Project_I_on_thetagrid(theta_vec,screen,lam) 
#fig = plt.figure(figsize=(8,6))
#plt.imshow(screen['E'],extent=(0,screen['D'],screen['D'],0))
#plt.colorbar()
endtime = dt.datetime.now()
t = (endtime-starttime).seconds
print('Time = {:.2f}s'.format(t))
#plt.savefig('/home/gemma/{}_E0.png'.format(name))



# Idiff = II1 - II0
# fig = plt.figure(figsize=(16,6))
# plt.subplot(1,3,1)
# plt.imshow(np.log10(np.abs(II0)),interpolation = None,extent=[-thetamaxdeg,thetamaxdeg,thetamaxdeg,-thetamaxdeg])
# plt.colorbar()
# plt.subplot(1,3,2)
# plt.imshow(np.log10(np.abs(II1)),interpolation = None,extent=[-thetamaxdeg,thetamaxdeg,thetamaxdeg,-thetamaxdeg])
# plt.colorbar()
# plt.subplot(1,3,3)
# plt.imshow(np.log10(np.abs(Idiff)),interpolation = None,extent=[-thetamaxdeg,thetamaxdeg,thetamaxdeg,-thetamaxdeg])
# plt.colorbar()
# plt.savefig('/home/gemma/{}_I.png'.format(name))

