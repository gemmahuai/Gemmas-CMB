### 01/30/2023 updated
### this python script tests if our numerical fourier transform of a interpolated perfect gaussian beam (in sky) 
### agrees with the analytical solution of a gaussian beam given by exp(-l(l+1)sigma^2/2), 
### numerical & analytical results match in this version!

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
import csv

# input
if (len(sys.argv)!=7):
    print('Wrong inputs!')
    print('Usage: python BEAM_test_new.py N_screen N_theta_interp screenD screen_sigma truncate(y/n) angle(deg)')
    #for example, python BEAM_test_new.py 4096 1024 10 1.0 y 3.0
    sys.exit()

N_screen = np.int(sys.argv[1])
N_theta = np.int(sys.argv[2])
D = np.float(sys.argv[3])
sigma = np.float(sys.argv[4])
trunc = str(sys.argv[5])
maxdeg = np.float(sys.argv[6]) # in degrees

# measure sigma in angle (rad)
def gaussian(x, A, sigma, x0): 
  g = A*np.exp(-(x-x0)**2/(2*sigma**2))
  return(g)

def analytical(l, sigma): # blm
  fft = bin_mean0.max()* np.exp(-l*(l+1)*sigma**2/2)
  return(fft)

# create screen [m]
screen = {}
screen['N'] = N_screen
screen['D'] = D
Fhcalc.Initialize(screen)
center = (screen['D']/2, screen['D']/2)
Fhcalc.MultByGaussian(screen, center, sigma)
if trunc=='y':
  Fhcalc.InCircle(screen, center, 2.0)
Fhcalc.ScreenFFT(screen)

# interpolation in sky intensity [rad]
lam = 0.002 #mm wavelength
thetamax = maxdeg*np.pi/180. # degrees to radians
theta_vec = np.linspace(-thetamax,thetamax,N_theta)  #rad
II0 = Fhcalc.Project_I_on_thetagrid(theta_vec,screen,lam) 
# shift the beam from being centered at theta=0 to theta=thetamax so that the beam spans from 0 deg to 2*thetamax deg
theta_vec = np.linspace(0,2*thetamax,N_theta) #rad


# FT of sky intensity
fft_I0 = np.abs(np.fft.fftshift(np.fft.fft2(np.fft.fftshift(II0)))) # in 1/rad space

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
### this uses a horizontal slice of data
# (fit, err) = opt.curve_fit(gaussian, theta_vec, II0[int(n/2)], p0=np.array([II0[int(n/2)].max(), 0.01, theta_vec.max()/2]), absolute_sigma=True)
# beam_sky = gaussian(theta_vec, fit[0], fit[1], fit[2])
### this uses the binned sky intensity
(fit, err) = opt.curve_fit(gaussian, bins[:-1], bin_I0, p0=np.array([II0[int(n/2)].max(), 0.001, 0.0]), absolute_sigma=True)
beam_sky = gaussian(theta_vec, fit[0], fit[1], theta_vec.max()/2)
print('beam sigma is {:.5f} rad'.format(fit[1]))


# numerical: average (FT of II0)^2 radially 
bin_edges = np.linspace(0,l.max(),int(len(l_vec)/2))
l_flatten = l.flatten()
fft_numerical = fft_I0.flatten()
bin_mean0, bin_edge, bin_num = binned_statistic(l_flatten, fft_numerical, statistic='mean', bins=bin_edges) 


# plot screen, II0, and fft(II0)
fig = plt.figure(figsize=(17,5))
plt.subplot(1,3,1)
plt.imshow(screen['E'], extent=(0, screen['D'], screen['D'], 0))
plt.colorbar()
plt.title('screenE')
plt.subplot(1,3,2)
plt.imshow(np.log10(II0), extent=(0, 2*maxdeg, 2*maxdeg, 0))
plt.colorbar()
plt.title('Interpolated sky intensity')
plt.subplot(1,3,3)
plt.imshow(fft_I0, extent=(-l_vec.min(), l_vec.max(), l_vec.max(), -l_vec.min()))
plt.colorbar()
plt.title('FFT (II0)')
# if trunc=='y':
#   plt.savefig('/home/gemma/Beam/beam_imshow_N{}N{}_D{}_sig{}_trunc_ang{}.png'.format(N_screen,N_theta,D,sigma,maxdeg))
# else:
#   plt.savefig('/home/gemma/Beam/beam_imshow_N{}N{}_D{}_sig{}_ang{}.png'.format(N_screen,N_theta,D,sigma,maxdeg))


# plot fit
fig = plt.figure(figsize=(17,5))
plt.subplot(1,3,1)
plt.plot(theta_vec, II0[int(n/2)], label='horizontal slice')
plt.plot(theta_vec, II0[:,int(n/2)], label='vertical slice')
plt.plot(theta_vec, beam_sky, ls='--', lw=5, color='black', label='gaussian fit')
plt.title('beam sigma is {:.5f} rad'.format(fit[1]))
plt.xlabel('angle (radians)')
plt.xlim(0.45*theta_vec.max(), 0.55*theta_vec.max())
plt.legend()
#plot analytical & numerical solutions
plt.subplot(1,3,2)
plt.plot(bin_edges[:-1], bin_mean0/bin_mean0.max(), color='black', label='numerical')
plt.plot(bin_edges, analytical(bin_edges, fit[1])/bin_mean0.max(), color='red', label='analytical')
plt.xlabel(r'$\ell$')
plt.xlim(0,20000)
plt.legend()
plt.subplot(1,3,3)
plt.plot(bin_edges[:-1], bin_mean0/bin_mean0.max(), color='black', label='numerical')
plt.plot(bin_edges, analytical(bin_edges, fit[1])/bin_mean0.max(), color='red', label='analytical')
plt.xlabel(r'$\ell$')
plt.xlim(0,1000)
plt.ylim(0.95, 1.0)
plt.legend()
#plt.show()
if trunc=='y':
  plt.savefig('./beam_test_N{}N{}_D{}_sig{}_trunc_ang{}.png'.format(N_screen,N_theta,D,sigma,maxdeg))
else:
  plt.savefig('./beam_test_N{}N{}_D{}_sig{}_ang{}.png'.format(N_screen,N_theta,D,sigma,maxdeg))

### write to a text file
# a = bin_mean0 #numerical
# b = analytical(bin_edges, fit[1]) #analytical
# zip(a,b)
# with open('output_angle{}.csv'.format(maxdeg), 'w') as f:
#    writer = csv.writer(f, delimiter='\t')
#    writer.writerows(zip(a,b))
