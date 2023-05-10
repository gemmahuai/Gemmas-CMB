### 01/31/2023 updated
### this python script compares the power spectrum of the perturbed and unperturbed beams
### analytical and numerical solutions matched!!!

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import binned_statistic
import scipy.optimize as opt
mpl.rc('xtick', direction='in', top=True)
mpl.rc('ytick', direction='in', right=True)
mpl.rc('xtick.minor', visible=True)
mpl.rc('ytick.minor', visible=True)
plt.rcParams['font.size'] = 13
plt.rcParams['figure.figsize'] = [8,6]
import sys
import Fhcalc
import ErrMask

# input
if (len(sys.argv)!=9):
    print('Wrong inputs!')
    print('Usage: python BEAM_error_new.py N_screen N_theta_interp screenD sigma angle(deg) truncation(y/n) phase/amp amplitude')
    sys.exit()
# input parameters
N_screen = np.int(sys.argv[1])
N_theta = np.int(sys.argv[2])
D = float(sys.argv[3])
sigma = float(sys.argv[4])
maxdeg = float(sys.argv[5])
trunc = str(sys.argv[6])
option = str(sys.argv[7])
amp = float(sys.argv[8])

# measure sigma in angle (rad)
def gaussian(x, A, sigma, x0): 
  g = A*np.exp(-(x-x0)**2/(2*sigma**2))
  return(g)

def analytical(l, sigma): # blm
  fft = bin_mean0.max()* np.exp(-l*(l+1)*sigma**2/2)
  return(fft)
 

# perfect gaussian create screen [m]
screen = {}
screen['N'] = N_screen
screen['D'] = D
Fhcalc.Initialize(screen)
center = (screen['D']/2, screen['D']/2)
Fhcalc.MultByGaussian(screen, center, sigma)
if trunc=='y':
  Fhcalc.InCircle(screen, center, 2.0)
Fhcalc.ScreenFFT(screen)

# perturbed screen
screen1 = {}
screen1['N'] = N_screen
screen1['D'] = D
Fhcalc.Initialize(screen1)
Fhcalc.MultByGaussian(screen1, center, sigma)
if trunc=='y':
  Fhcalc.InCircle(screen1, center, 2.0)

if option=='phase':
  ErrMask.filter_annulus_phase(screen1, amp, 10, 12)
elif option=='amp':
  ErrMask.filter_annulus_amp(screen1, amp, 10, 12)
else: print('Choose phase or amplitude errors')

Fhcalc.ScreenFFT(screen1)

# interpolation in sky intensity [rad]
lam = 0.002 #mm wavelength
thetamaxdeg = maxdeg
thetamax = thetamaxdeg*np.pi/180. # in radians
theta_vec = np.linspace(-thetamax,thetamax,N_theta) 
II0 = Fhcalc.Project_I_on_thetagrid(theta_vec,screen,lam)   # unperturbed
II1 = Fhcalc.Project_I_on_thetagrid(theta_vec, screen1, lam) # perturbed
# shift the beam from being centered at theta=0 to theta=thetamax so that the beam spans from 0 deg to 2*thetamax deg
theta_vec = np.linspace(0,2*thetamax,N_theta) #rad

# FT of sky intensity
fft_I0 = np.abs(np.fft.fftshift(np.fft.fft2(np.fft.fftshift(II0)))) # in 1/rad space
fft_I1 = np.abs(np.fft.fftshift(np.fft.fft2(np.fft.fftshift(II1))))
# fig = plt.figure(figsize=(12,5))
# plt.subplot(1,3,1)
# plt.imshow(fft_I0,extent=(-thetamaxdeg, thetamaxdeg,-thetamaxdeg, thetamaxdeg)) # plot absolute value of fft_I
# plt.colorbar()
# plt.subplot(1,3,2)
# plt.imshow(fft_I1,extent=(-thetamaxdeg, thetamaxdeg,-thetamaxdeg, thetamaxdeg)) # plot absolute value of fft_I
# plt.colorbar()
# plt.subplot(1,3,3)
# plt.imshow(np.abs(fft_I0 - fft_I1),extent=(-thetamaxdeg, thetamaxdeg,-thetamaxdeg, thetamaxdeg)) # plot absolute value of fft_I
# plt.colorbar()
# plt.savefig('/home/gemma/Beam/FFT_{}.png'.format(name))

#calculate ell
n = theta_vec.shape[0]
dl = 2*np.pi/theta_vec.max() # dl in 1/rad space
l_vec = np.fft.fftshift(dl * np.fft.fftfreq(n)*n)
(l_x, l_y) = np.meshgrid(l_vec,l_vec) # 1/rad 
l = np.sqrt(l_x**2 + l_y**2)

# bin sky intensity radially
(thetax, thetay) = np.meshgrid(theta_vec, theta_vec) # rad 
theta_r = np.sqrt((thetax-thetax.max()/2)**2 + (thetay-thetay.max()/2)**2)
bins = np.linspace(0,theta_r.max(),int(len(theta_vec)/1.2))
theta_flatten = theta_r.flatten()
I_flatten = II0.flatten()
bin_I0, bin_edge, bin_num = binned_statistic(theta_flatten, I_flatten, statistic='mean', bins=bins) 

### fit 1D in sky (rad space)
### this uses the horizontal slice of data
# (fit, err) = opt.curve_fit(gaussian, theta_vec, II0[int(n/2)], p0=np.array([II0[int(n/2)].max(), 0.01, theta_vec.max()/2]), absolute_sigma=True)
# beam_sky = gaussian(theta_vec, fit[0], fit[1], fit[2])
### this uses the binned sky intensity
(fit, err) = opt.curve_fit(gaussian, bins[:-1], bin_I0, p0=np.array([II0[int(n/2)].max(), 0.01, 0.0]), absolute_sigma=True)
beam_sky = gaussian(theta_vec, fit[0], fit[1], theta_vec.max()/2)
print('beam sigma is {:.5f} rad'.format(fit[1]))

# average (FT of II0)^2 radially 
bin_edges = np.linspace(0,l.max(),int(len(l_vec)/2))
l_flatten = l.flatten()
fft_numerical0 = fft_I0.flatten()
fft_numerical1 = fft_I1.flatten()
bin_mean0, bin_edge, bin_num = binned_statistic(l_flatten, fft_numerical0, statistic='mean', bins=bin_edges) 
bin_mean1, bin_edge, bin_num = binned_statistic(l_flatten, fft_numerical1, statistic='mean', bins=bin_edges) # bin_mean is the binned numerical beam
l_vec = bin_edges[0:-1] # ell 1D vector


# plot the numerical beam and analytical beam
fig = plt.figure(figsize=(22,5))
plt.subplot(1,4,1)
plt.imshow(fft_I1 - fft_I0, interpolation=None) 
plt.colorbar()
plt.title('FT(I) difference')
plt.subplot(1,4,2)
plt.imshow(fft_I1 - fft_I0, interpolation=None) 
plt.colorbar()
plt.xlim(int(N_theta*0.43),int(N_theta*0.57)) 
plt.ylim(int(N_theta*0.43),int(N_theta*0.57))
plt.subplot()
plt.title('FT(I) difference')
plt.subplot(1,4,3)
plt.plot(l_vec,bin_mean0,color='blue', lw=2, label='binned unperturbed') # normalize the numerical solution
plt.plot(l_vec,bin_mean1,color='red', lw=2, label='binned perturbed')
plt.xlabel(r'$\ell$')
plt.plot(bin_edges, analytical(bin_edges, fit[1]), color='black', lw=3, ls=':', label='analytical blm')
plt.xlim(0,2e4)
plt.subplot(1,4,4)
plt.plot(l_vec,bin_mean0,color='blue', lw=2, label='binned unperturbed') # normalize the numerical solution
plt.plot(l_vec,bin_mean1,color='red', lw=2, label='binned perturbed')
plt.xlabel(r'$\ell$')
analytic = analytical(bin_edges, fit[1])
plt.plot(bin_edges, analytic, lw=3, color='black', ls=':', label='analytical blm')
plt.xlim(0,500)

if bin_mean0.max()<bin_mean1.max():
  plt.ylim(0.9*bin_mean0.max(), 1.1*bin_mean1.max())
else: 
  plt.ylim(0.7*bin_mean1.max(), 1.1*bin_mean0.max())

plt.legend()
plt.savefig('/home/gemma/Beam/Spectrum_N{}N{}_D{}_sig{}_ang{}_{}{}.png'.format(N_screen,N_theta,D,sigma,maxdeg,option,amp))



