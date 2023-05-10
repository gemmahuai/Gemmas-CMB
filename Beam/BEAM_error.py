# this python script compares the power spectrum of the perturbed and unperturbed beams


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
if (len(sys.argv)!=6):
    print('Wrong inputs!')
    print('Usage: python BEAM_error.py output_filename N_screen N_theta_interp phase/amp amplitude')
    sys.exit()

name = sys.argv[1] #file name
N_screen = np.int(sys.argv[2])
N_theta = np.int(sys.argv[3])
option = str(sys.argv[4])
amp = float(sys.argv[5])

# define beam and gaussian functions
def b_lm(l, sigma):
  #l = screen['kap']*np.pi*2
  #blm = (np.sqrt(2*l+1)/(4*np.pi)) * np.exp(-l*(l+1)*sigma**2/2)
  blm = np.exp(-l*(l+1)*sigma**2/2)
  return(blm)
def sigma_calc(r,sigma):
  gaussian = bin_mean.max()*np.exp(-r**2/(2*sigma**2))
  return(gaussian)


# perfect gaussian create screen [m]
screen = {}
screen['N'] = N_screen
screen['D'] = 10
Fhcalc.Initialize(screen)
center = (screen['D']/2, screen['D']/2)
Fhcalc.MultByGaussian(screen, center, 1.0)
Fhcalc.InCircle(screen, center, 2.0)
Fhcalc.ScreenFFT(screen)

# perturbed screen
screen1 = {}
screen1['N'] = N_screen
screen1['D'] = 10
Fhcalc.Initialize(screen1)
Fhcalc.MultByGaussian(screen1, center, 1.0)
Fhcalc.InCircle(screen1, center, 2.0)
if option=='phase':
  ErrMask.filter_annulus_phase(screen1, amp, 10, 12)
elif option=='amp':
  ErrMask.filter_annulus_amp(screen1, amp, 10, 12)
else: print('Choose phase or amplitude errors')

Fhcalc.ScreenFFT(screen1)

# interpolation in sky intensity [rad]
lam = 0.002 #mm wavelength
thetamaxdeg = 3.0
thetamax = thetamaxdeg*np.pi/180. # in radians
theta_vec = np.linspace(-thetamax,thetamax,N_theta) 
II0 = Fhcalc.Project_I_on_thetagrid(theta_vec,screen,lam)   # unperturbed
II1 = Fhcalc.Project_I_on_thetagrid(theta_vec, screen1, lam) # perturbed

# FFT the sky intensity to ell space (ifft for comparison)
fft_I0 = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(II0)))  # fft of sky intensity (1/rad) 
fft_I1 = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(II1)))
blm_numerical0 = np.abs(fft_I0)**2 #numerically calculated beam in [1/rad] k space
blm_numerical1 = np.abs(fft_I1)**2
fig = plt.figure(figsize=(12,5))
plt.subplot(1,3,1)
plt.imshow(np.abs(blm_numerical0),extent=(-thetamaxdeg, thetamaxdeg,-thetamaxdeg, thetamaxdeg)) # plot absolute value of fft_I
plt.colorbar()
plt.subplot(1,3,2)
plt.imshow(np.abs(blm_numerical1),extent=(-thetamaxdeg, thetamaxdeg,-thetamaxdeg, thetamaxdeg)) # plot absolute value of fft_I
plt.colorbar()
plt.subplot(1,3,3)
plt.imshow(np.abs(blm_numerical1 - blm_numerical0),extent=(-thetamaxdeg, thetamaxdeg,-thetamaxdeg, thetamaxdeg)) # plot absolute value of fft_I
plt.colorbar()
plt.savefig('/home/gemma/Beam/FFT_{}.png'.format(name))

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
I0_flatten = II0.flatten()
bin_mean, bin_edge, bin_num = binned_statistic(theta_flatten, I0_flatten, statistic='mean', bins=bin_edges) # calculate binned sky intensity
(fit, err) = opt.curve_fit(sigma_calc, bin_edges[0:-1], bin_mean, p0=1e-4, absolute_sigma=True) # fit the binned intensity for the beam width sigma
print('Fitted beam sigma in sky is {:.5f} radians'.format(fit[0]))
# fig = plt.figure(figsize=(7,5))
# plt.plot(bin_edges[0:-1],bin_mean, lw=3, color='blue',label='binned')
# plt.xlabel('k')
# plt.plot(bin_edges, sigma_calc(bin_edges, fit[0]), lw=1.5, color='orange',label='fit')
# plt.legend()
# plt.xlim(0,0.005)
# plt.savefig('/home/gemma/Beam/Fitted_skyII0_{}'.format(name))

# average (FT of II0)^2 radially 
bin_edges = np.arange(0,l.max(),1000)
l_flatten = l.flatten()
beam_numerical0_flatten = blm_numerical0.flatten()
beam_numerical1_flatten = blm_numerical1.flatten()
bin_mean0, bin_edge, bin_num = binned_statistic(l_flatten, beam_numerical0_flatten, statistic='mean', bins=bin_edges) # bin_mean is the binned numerical beam
bin_mean1, bin_edge, bin_num = binned_statistic(l_flatten, beam_numerical1_flatten, statistic='mean', bins=bin_edges) # bin_mean is the binned numerical beam
l_vec = bin_edges[0:-1] # ell 1D vector
# plot the numerical beam and analytical beam
fig = plt.figure(figsize=(10,6))
plt.subplot(1,2,1)
plt.plot(l_vec,bin_mean0/bin_mean0.max(),color='mediumblue', lw=5, label='binned unperturbed') # normalize the numerical solution
plt.plot(l_vec,bin_mean1/bin_mean1.max(),color='orangered', lw=2, label='binned perturbed')
plt.xlabel(r'$\ell$')
blm = b_lm(l_vec, fit[0]) # analytical solution
plt.plot(l_vec, blm, color='black', label='blm') # plot the analytical solution
plt.xlim(0,3e4)
plt.subplot(1,2,2)
plt.plot(l_vec,bin_mean0/bin_mean0.max(),color='mediumblue', lw=5, label='binned unperturbed') # normalize the numerical solution
plt.plot(l_vec,bin_mean1/bin_mean1.max(),color='orangered', lw=2, label='binned perturbed')
plt.xlabel(r'$\ell$')
blm = b_lm(l_vec, fit[0]) # analytical solution
plt.plot(l_vec, blm, color='black', label='blm') # plot the analytical solution
plt.xlim(0,1000)
plt.ylim(0.98,1)
plt.legend()
plt.savefig('/home/gemma/Beam/Spectrum_{}'.format(name))



