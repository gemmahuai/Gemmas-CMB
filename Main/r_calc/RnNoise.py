### updated on 03/16/2023
### This script reads in actual CMB T and P power spectra, iterates over an array of error amplitude, 
### calculates the leakage power spectrum at each amplitude, calculates T to P leakage spectra
### measures r of the T to P spectra by comparing to the theoretical spectra with r=1 generated from the Goddard Space Flight Center web
### Returns: a plot of tensor to scalar ratio r as a function of normalized noise
###         + a text file containing error amp, normalized noise, ell vectors, r, leakage power spectra, and rms

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import binned_statistic
import scipy.optimize as opt
import scipy.interpolate as interp
mpl.rc('xtick', direction='in', top=True)
mpl.rc('ytick', direction='in', right=True)
mpl.rc('xtick.minor', visible=True)
mpl.rc('ytick.minor', visible=True)
plt.rcParams['font.size'] = 13
import sys
sys.path.insert(0, '/home/gemma/Beam/')
import Fhcalc
import ErrMask
import csv


# input & Initialize parameters
if (len(sys.argv)!=8):
    print('Wrong inputs!')
    print('Usage: python RnNoise.py N_screen N_theta_interp screenD sigma angle(deg) truncation(y/n) phase/amp')
    # ex. python RnNoise.py 4096 1024 10.0 1.0 5.0 y amp
    sys.exit()

# input parameters
N_screen = np.int(sys.argv[1])
N_theta = np.int(sys.argv[2])
D = float(sys.argv[3])
sigma = float(sys.argv[4])
maxdeg = float(sys.argv[5])
trunc = str(sys.argv[6])
option = str(sys.argv[7])

leak = []
#amp = np.arange(0.1,1.21,0.05) # phase errors
amp = np.arange(0.01, 0.3, 0.02)
#amp = np.array([0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.24, 0.3])
#amp = np.arange(0.05, 2.0, 0.1) # amplitude errors

#amp = np.array([0.1,1.0])
r = np.array([])
kin = 15
kout = 25
noise_norm = np.zeros_like(amp)
txt_leak = np.array([]) # 1D array containing all leakage spectrum, will be printed to a text file.
RMS = np.array([])

# CMB power spectra from https://lambda.gsfc.nasa.gov/toolbox/camb_online.html
(ll, T, E, B) = np.loadtxt('camb_96896687_totcls.dat.txt', unpack=True, usecols=(0,1,2,3))

def rms(error_map):
    return np.sqrt(np.mean(error_map**2))

# perfect gaussian create E screen [m]
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

# fig = plt.figure(figsize=(13,6))
# plt.subplot(1,2,1)

for i in range(len(amp)):
    print('amp =', amp[i])
    
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
    II1 = II1 * (np.sum(II0)/np.sum(II1)) ## normalize the perturbed beam
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
    fft_numerical1 = fft_I1.flatten()
    fft_numerical_diff = fft_Idiff.flatten()
    bin_mean0, bin_edge, bin_num = binned_statistic(l_flatten, fft_numerical0, statistic='mean', bins=bin_edges) 
    bin_mean1, bin_edge, bin_num = binned_statistic(l_flatten, fft_numerical1, statistic='mean', bins=bin_edges) # bin_mean is the binned numerical beam
    bin_mean_diff, bin_edge, bin_num = binned_statistic(l_flatten, fft_numerical_diff, statistic='mean', bins=bin_edges) # bin the beam difference
    l_vec = bin_edges[0:-1] # ell 1D vector
    
    beam_diff_rela = bin_mean_diff/bin_mean0 # relative beam difference
    print(l_vec[1:3]) # making sure we average the leakage at ell=85 and ell=170
    leak = np.append(leak, np.mean(beam_diff_rela[(l_vec>=49)&(l_vec<=251)])) # average the beam difference at l=85 and l=170 (2nd and 3rd elements)
    txt_leak = np.append(txt_leak, beam_diff_rela**2) # leakage power spectrum stored in the 1D array leak_spec
    print('Error amplitude is {}, the corresponding leakage^2 is {:.5f}'.format(amp[i], leak[i]**2))
    
    
    ### plot the leakage spectrum and T->B power spectrum
    spline = interp.InterpolatedUnivariateSpline(l_vec, beam_diff_rela**2)
    leakage = spline(ll)
    # plt.loglog(ll,T*leakage, label='leakage {}'.format(amp[i])) 
    
    ### calculate r
    # B0 = B[np.where(ll==100)] # BB power at ell=100 with r=1
    # r = np.append(r, 1 * (T[np.where(ll==100)]*spline(100) / B0)) # r of the leakage at ell=100
    # print('r=',r[i])
    want = np.where((ll<250)&(ll>50))
    B0 = np.mean(B[want]) # BB power at ell=100 with r=1
    B_leak = np.mean(T[want]*spline(ll[want]))
    r = np.append(r, 1 * (B_leak / B0)) # r of the leakage at ell=100
    print('r=',r[i])
    
    ### calculate normalized errors  
    emap_fft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(emap_E))) # scaled error map in fourier space
    norm_fft = np.abs(emap_fft**2)/screen1['dk']**2/screen1['N']**4 # V^2/Hz^2
    idx = np.where((screen1['kap']<kout) & (screen1['kap']>kin)) # inside the filter
    avg_fft = np.mean(norm_fft[idx]) # V^2/Hz^2
    total_Vf = np.sum(norm_fft)*screen1['dk']**2
    rms_sq = np.abs(rms(emap_E)**2)
    noise_norm[i] = avg_fft
    print('average error is ', avg_fft)
    RMS = np.append(RMS, np.sqrt(rms_sq))
    
# plt.loglog(ll, T, label='TT')
# plt.loglog(ll, B, label='BB')
# plt.xlabel(r'$\ell$')
# plt.ylabel(r'$\mu K^2$')
# plt.legend()
# plt.xlim(2,2200)
# plt.legend()

# plt.subplot(1,2,2)
# plt.semilogy(noise_norm, r, 'o--', ms=3)
# plt.xlabel('Noise level')
# plt.ylabel('r')
# plt.savefig('/home/gemma/Beam/r/fig_LeakSpec_r.png')

### write to a csv file
txt_amp = np.repeat(amp, len(beam_diff_rela)) # 0th column - noise amplitude
txt_noise = np.repeat(noise_norm, len(beam_diff_rela)) # 1st column - normalized noise level in k space
txt_ell = np.tile(l_vec, len(amp)) # 2nd column - ell vector
txt_r = np.repeat(r, len(beam_diff_rela)) # 3rd column - tensor to scalar ratio r
txt_rms = np.repeat(RMS, len(beam_diff_rela)) # 5th column - rms error in real space 
# 4th column - txt_leakage stored in the loop, already a 1D array
with open('./R_noise_{}_r.csv'.format(option), 'w') as file:
    writer = csv.writer(file, delimiter='\t')
    writer.writerows(zip(txt_amp, txt_noise, txt_ell, txt_r, txt_leak, txt_rms))