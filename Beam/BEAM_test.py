# this python script tests if our numerical fourier transform of a interpolated perfect gaussian beam (in sky) 
# agrees with the analytical solution of a gaussian beam given by exp(-l(l+1)sigma^2/2), given size of the screen 
# and the interpolation size. It returns the fitted gaussian beam width in sky [rad] and outputs three plots - 
# screen&interpolation, abs(fft of gaussian beam in sky) [1/rad], and 1D power spectrum as a function of ell. 

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

# input
if (len(sys.argv)!=4):
    print('Wrong inputs!')
    print('Usage: python Beam_test.py output_filename N_screen N_theta_interp')
    sys.exit()

name = sys.argv[1] #file name
N_screen = np.int(sys.argv[2])
N_theta = np.int(sys.argv[3])

# define beam and gaussian functions
def b_lm(l, sigma):
  blm = np.exp(-l*(l+1)*sigma**2/2)
  return(blm)
def sigma_calc(r,sigma):
  gaussian = bin_mean.max()*np.exp(-r**2/(2*sigma**2))
  return(gaussian)


# create screen [m]
screen = {}
screen['N'] = N_screen
screen['D'] = 10
Fhcalc.Initialize(screen)
center = (screen['D']/2, screen['D']/2)
Fhcalc.MultByGaussian(screen, center, 1.0)
Fhcalc.InCircle(screen, center, 2.0)
Fhcalc.ScreenFFT(screen)

# interpolation in sky intensity [rad]
lam = 0.002 #mm wavelength
thetamaxdeg = 3.0
thetamax = thetamaxdeg*np.pi/180. # 2 degrees in radians
theta_vec = np.linspace(-thetamax,thetamax,N_theta) 
II0 = Fhcalc.Project_I_on_thetagrid(theta_vec,screen,lam) 
# plot the screen and interpolation 
fig = plt.figure(figsize=(14,6))
plt.subplot(1,2,1)
plt.imshow(screen['E'],extent=(0, screen['D'], screen['D'], 0))
plt.colorbar()
plt.subplot(1,2,2)
plt.imshow(np.log10(II0))
plt.colorbar()
#plt.xlim(460,564)
#plt.ylim(460,564)
plt.xlim(int(0.4*N_theta),int(0.6*N_theta))
plt.ylim(int(0.4*N_theta),int(0.6*N_theta))
plt.savefig('/home/gemma/Beam/ScreennInterpolation_{}.png'.format(name))

# FFT the sky intensity to ell space (ifft for comparison)
fft_iI = np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(II0))) # ifft of sky intensity 
fft_I = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(II0)))   # fft of sky intensity (1/rad) 
blm_numerical = np.abs(fft_I)**2 #numerically calculated beam in [1/rad] k space
fig = plt.figure(figsize=(7,5))
# plt.imshow(np.abs(fft_iI))
# plt.colorbar()
# plt.show()
plt.imshow(np.abs(fft_I),extent=(-thetamaxdeg, thetamaxdeg,-thetamaxdeg, thetamaxdeg)) # plot absolute value of fft_I
plt.colorbar()
plt.savefig('/home/gemma/Beam/{}_FFT_I.png'.format(name))

# calculate ell
n = theta_vec.shape[0]
dk = 1/theta_vec.max() # dk in 1/rad space
k_vec = dk*np.fft.fftshift(np.fft.fftfreq(n))*n
(k_x, k_y) = np.meshgrid(k_vec,k_vec) # 1/rad 
k_r = np.sqrt(k_x**2 + k_y**2)
l = k_r * 2 * np.pi
#plt.imshow(l)

# fit sigma in sky (rad)
(thetax, thetay) = np.meshgrid(theta_vec, theta_vec) # rad 
theta_r = np.sqrt(thetax**2 + thetay**2)
bin_edges = np.linspace(0,theta_r.max(),int(0.5*N_theta)) # bin radially
theta_flatten = theta_r.flatten()
I_flatten = II0.flatten()
bin_mean, bin_edge, bin_num = binned_statistic(theta_flatten, I_flatten, statistic='mean', bins=bin_edges) # calculate binned sky intensity
(fit, err) = opt.curve_fit(sigma_calc, bin_edges[0:-1], bin_mean, p0=1e-4, absolute_sigma=True) # fit the binned intensity for the beam width sigma
print('Fitted beam sigma in sky is {:.5f} radians'.format(fit[0]))
fig = plt.figure(figsize=(7,5))
plt.plot(bin_edges[0:-1],bin_mean, lw=3, color='blue',label='binned')
plt.xlabel('k')
plt.plot(bin_edges, sigma_calc(bin_edges, fit[0]), lw=1.5, color='orange',label='fit')
plt.legend()
plt.xlim(0,0.005)
plt.savefig('/home/gemma/Beam/{}_Fitted_skyII0'.format(name))

# average (FT of II0)^2 radially 
bin_edges = np.arange(0,l.max(),1000)
l_flatten = l.flatten()
beam_numerical_flatten = blm_numerical.flatten()
bin_mean, bin_edge, bin_num = binned_statistic(l_flatten, beam_numerical_flatten, statistic='mean', bins=bin_edges) # bin_mean is the binned numerical beam
l_vec = bin_edges[0:-1] # ell 1D vector
# plot the numerical beam and analytical beam
fig = plt.figure(figsize=(8,6))
plt.plot(l_vec,bin_mean/bin_mean.max(),label='binned') # normalize the numerical solution
plt.xlabel(r'$\ell$')
blm = b_lm(l_vec, fit[0]) # analytical solution
plt.plot(l_vec, blm,label='blm') # plot the analytical solution
plt.xlim(0,3e4)
plt.legend()
plt.savefig('/home/gemma/Beam/{}_Beam_NumAna'.format(name))



