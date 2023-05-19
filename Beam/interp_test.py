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
    print('Usage: python interp_test.py output_filename N_screen N_theta_interp')
    sys.exit()

name = sys.argv[1] #file name
N_screen = np.int(sys.argv[2])
N_theta = np.int(sys.argv[3])

screen0 = {}
screen0['N'] = N_screen
screen0['D'] = 10
Fhcalc.Initialize(screen0)
center = (screen0['D']/2, screen0['D']/2)
Fhcalc.MultByGaussian(screen0, center, 1.0)
Fhcalc.InCircle(screen0, center, 2.0)
Fhcalc.ScreenFFT(screen0)

lam = 0.002 #mm wavelength
kphot = 2*np.pi/lam
thetamaxdeg = 2.0
thetamax = thetamaxdeg*np.pi/180. # 2 degrees in radians
theta_vec = np.linspace(-thetamax,thetamax,N_theta) 
II0 = Fhcalc.Project_I_on_thetagrid(theta_vec,screen0,lam) 

#fft to ell space
fft_iI = np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(II0))) # ifft of sky intensity 
fft_I = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(II0)))   # fft of sky intensity (1/rad)
fig = plt.figure(figsize=(8,6))
plt.subplot(2,1,1)
plt.imshow(np.abs(fft_iI))
plt.colorbar()
plt.subplot(2,1,2)
plt.imshow(np.abs(fft_I))
plt.colorbar()
plt.savefig('./{}.png'.format(name))
