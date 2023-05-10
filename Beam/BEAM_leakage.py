### 04/04/2023 updated
### this python script compares the power spectrum of the perturbed and unperturbed beams;
### takes the difference in theta space (sky), and then FFT to ell space;
### generates a plot of leakage level (error^2) as a function of error amplitude for a fixed k_in and k_out and after noise normalization

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import binned_statistic
import scipy.optimize as opt
mpl.rc('xtick', direction='in', top=True)
mpl.rc('ytick', direction='in', right=True)
mpl.rc('xtick.minor', visible=True)
mpl.rc('ytick.minor', visible=True)
plt.rcParams['font.size'] = 15
import Fhcalc
import ErrMask
import sys
import csv

# input & Initialize parameters
if (len(sys.argv)!=8):
    print('Wrong inputs!')
    print('Usage: python BEAM_leakage.py N_screen N_theta_interp screenD sigma angle(deg) truncation(y/n) phase/amp')
    # ex. python BEAM_leakage.py 8192 1024 10.0 1.0 5.0 y amp
    sys.exit()
# input parameters
N_screen = np.int(sys.argv[1])
N_theta = np.int(sys.argv[2])
D = float(sys.argv[3])
sigma = float(sys.argv[4])
maxdeg = float(sys.argv[5])
trunc = str(sys.argv[6])
option = str(sys.argv[7])
error = np.array([])
noise_norm = np.array([])
#amp = np.array([0.5,1.0]) # test
RMS_real = np.array([])
amp = np.array([0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.3, 1.6, 2.0])
#amp = np.array([0.005, 0.3, 0.6]) #
kin = 15
kout = 25

def gaussian(x, A, sigma, x0): 
    g = A*np.exp(-(x-x0)**2/(2*sigma**2))
    return(g)

def analytical(l, sigma): # blm
    fft = bin_mean0.max()* np.exp(-l*(l+1)*sigma**2/2)
    return(fft)

def rms(x):
    return np.sqrt(np.mean(x**2))


# unperturbed perfect gaussian create E screen [m]
screen = {}
screen['N'] = N_screen
screen['D'] = D
Fhcalc.Initialize(screen)
center = (screen['D']/2, screen['D']/2)
Fhcalc.MultByGaussian(screen, center, sigma)
if trunc=='y':
    Fhcalc.InCircle(screen, center, 2.0)
Fhcalc.ScreenFFT(screen)

# interpolation in sky intensity [rad] (unperturbed beam)
lam = 0.002 #mm wavelength
thetamaxdeg = maxdeg
thetamax = thetamaxdeg*np.pi/180. # in radians
theta_vec = np.linspace(-thetamax,thetamax,N_theta) 
II0 = Fhcalc.Project_I_on_thetagrid(theta_vec,screen,lam)   # unperturbed
fft_I0 = np.abs(np.fft.fftshift(np.fft.fft2(np.fft.fftshift(II0)))) # in ell space
theta_vec = np.linspace(0,2*thetamax,N_theta) #rad

#calculate ell
n = theta_vec.shape[0]
dl = 2*np.pi/theta_vec.max() # dl in 1/rad space
l_vec = np.fft.fftshift(dl * np.fft.fftfreq(n)*n)
(l_x, l_y) = np.meshgrid(l_vec,l_vec) # 1/rad 
l = np.sqrt(l_x**2 + l_y**2)

for i in range(len(amp)):
    
    # perturbed E screen
    screen1 = {}
    screen1['N'] = N_screen
    screen1['D'] = D
    Fhcalc.Initialize(screen1)
    Fhcalc.MultByGaussian(screen1, center, sigma)
    if trunc=='y':
        Fhcalc.InCircle(screen1, center, 2.0)
    
    if option=='phase':
        emap_E = ErrMask.filter_annulus_phase(screen1, amp[i], kin, kout)
    elif option=='amp':
        emap_E = ErrMask.filter_annulus_amp(screen1, amp[i], kin, kout)
    else: print('Choose phase or amplitude errors')
    
    Fhcalc.ScreenFFT(screen1)
    
    # interpolate the perturbed beam 
    theta_vec = np.linspace(-thetamax,thetamax,N_theta) 
    II1 = Fhcalc.Project_I_on_thetagrid(theta_vec, screen1, lam) # perturbed
    II1 = II1 * (np.sum(II0)/np.sum(II1)) # normalize the perturbed beam
    # shift the beam from being centered at theta=0 to theta=thetamax so that the beam spans from 0 deg to 2*thetamax deg
    theta_vec = np.linspace(0,2*thetamax,N_theta) #rad
    
    # FT of sky intensity
    I_diff = II1 - II0 # take the difference in theta space (sky)
    fft_Idiff = np.abs(np.fft.fftshift(np.fft.fft2(np.fft.fftshift(I_diff)))) # FFT the difference to ell space
    fft_I1 = np.abs(np.fft.fftshift(np.fft.fft2(np.fft.fftshift(II1)))) # in ell space
    
    # average (FT of II0)^2 radially 
    bin_edges = np.linspace(0,l.max(),int(len(theta_vec)/2))
    l_flatten = l.flatten()
    fft_numerical0 = fft_I0.flatten()
    # fft_numerical1 = fft_I1.flatten()
    fft_numerical_diff = fft_Idiff.flatten()
    bin_mean0, bin_edge, bin_num = binned_statistic(l_flatten, fft_numerical0, statistic='mean', bins=bin_edges) 
    # bin_mean1, bin_edge, bin_num = binned_statistic(l_flatten, fft_numerical1, statistic='mean', bins=bin_edges) # bin_mean is the binned numerical beam
    bin_mean_diff, bin_edge, bin_num = binned_statistic(l_flatten, fft_numerical_diff, statistic='mean', bins=bin_edges) # bin the beam difference
    l_vec = bin_edges[0:-1] # ell 1D vector
    
    beam_diff_rela = bin_mean_diff/bin_mean0 # relative beam difference
    error = np.append(error, np.mean(beam_diff_rela[(l_vec>=49)&(l_vec<=251)])) # average the beam difference at l=85 and l=170 (2nd and 3rd elements)
    print('Error amplitude = {}, leakage (error^2) = {:.3e}'.format(amp[i], error[i]**2))

    # Calculates the noise level in fourier space (V^2/Hz)
    emap_fft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(emap_E))) # scaled error map in fourier space
    norm_fft = np.abs(emap_fft**2)/screen1['dk']**2/screen1['N']**4 # V^2/Hz^2
    idx = np.where((screen1['kap']<kout) & (screen1['kap']>kin)) # inside the filter
    avg_fft = np.mean(norm_fft[idx]) # V^2/Hz^2
    total_Vf = np.sum(norm_fft)*screen['dk']**2
    rms_sq = np.abs(rms(emap_E)**2)
    print('V^2/Hz^2 = ',total_Vf)
    print('RMS^2 = ', rms_sq)
    noise_norm = np.append(noise_norm, avg_fft)
    RMS_real = np.append(RMS_real, np.sqrt(rms_sq))


### write to a csv file

with open('/home/gemma/Beam/data/Replot_leakrms_{}.csv'.format(option), 'w') as f:
   writer = csv.writer(f, delimiter='\t')
   writer.writerows(zip(amp, error**2, noise_norm, RMS_real))

print('leakage = ', error**2)
print('noise k = ', noise_norm)
print('rms noise (real) = ', RMS_real)
fig = plt.figure(figsize=(10,8))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
ax1.semilogy(RMS_real, error**2, 'o--', ms=7, color='black')
ax1.set_xlabel('RMS noise in real space')
ax1.set_ylabel('Leakage level')
ax2.set_xlim(noise_norm.min(), noise_norm.max())
ax2.set_xlabel('Noise in k space')
#plt.show()
plt.grid()
plt.savefig('/home/gemma/Beam/leak_rms_{}.png'.format(option))
