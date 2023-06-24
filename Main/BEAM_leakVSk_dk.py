### 04/23/2023 updated
### this python script fixes the noise level in k space
### and plots leakage level as a function of dk, fixing k_in
###

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
import Fhcalc
import ErrMask
import sys
import csv

# input & Initialize parameters
if (len(sys.argv)!=9):
    print('Wrong inputs!')
    print('Usage: python BEAM_leakVSk_dk.py N_screen N_theta_interp screenD sigma angle(deg) truncation(y/n) "phase/amp" amplitude')
    # ex. python BEAM_leakVSk_dk.py 4096 1024 10.0 1.0 5.0 y amp 0.1
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
kin = 0
#k_in = np.array([1, 10]) # test
dk = np.arange(1,80,4)
error = np.zeros(len(dk)) #leakage = error^2
noise_norm = np.zeros_like(error)

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

for i in range(len(dk)):
    kout = kin + dk[i]

    print('error amp =', amp, 'k_in =', kin, 'k_out =', kout)
    
    # perturbed E screen
    screen1 = {}
    screen1['N'] = N_screen
    screen1['D'] = D
    Fhcalc.Initialize(screen1)
    Fhcalc.MultByGaussian(screen1, center, sigma)
    if trunc=='y':
        Fhcalc.InCircle(screen1, center, 2.0)

    if option=='phase':
        emap_E = ErrMask.filter_annulus_phase(screen1, amp, kin, kout)
    elif option=='amp':
        emap_E = ErrMask.filter_annulus_amp(screen1, amp, kin, kout)
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
    #fft_numerical1 = fft_I1.flatten()
    fft_numerical_diff = fft_Idiff.flatten()
    bin_mean0, bin_edge, bin_num = binned_statistic(l_flatten, fft_numerical0, statistic='mean', bins=bin_edges) 
    #bin_mean1, bin_edge, bin_num = binned_statistic(l_flatten, fft_numerical1, statistic='mean', bins=bin_edges) # bin_mean is the binned numerical beam
    bin_mean_diff, bin_edge, bin_num = binned_statistic(l_flatten, fft_numerical_diff, statistic='mean', bins=bin_edges) # bin the beam difference
    l_vec = bin_edges[0:-1] # ell 1D vector

    beam_diff_rela = bin_mean_diff/bin_mean0 # relative beam difference
    #error[i] = np.mean(beam_diff_rela[1:3]) ## maxdeg=3.0
    #error[i] = np.mean(beam_diff_rela[4]) ## maxdeg=10.0
    error[i] = np.mean(beam_diff_rela[(l_vec>=49)&(l_vec<=251)]) # average the beam difference over 0<l<500
    print('leakage (error^2) = {}'.format(error[i]**2))

    # Calculates the noise level in fourier space (V^2/Hz)
    emap_fft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(emap_E))) # scaled error map in fourier space
    norm_fft = np.abs(emap_fft**2)/screen1['dk']**2/screen1['N']**4 # V^2/Hz^2
    idx = np.where((screen1['kap']<kout) & (screen1['kap']>kin)) # inside the filter
    avg_fft = np.mean(norm_fft[idx]) # V^2/Hz^2
    total_Vf = np.sum(norm_fft)*screen1['dk']**2
    rms_sq = np.abs(rms(emap_E)**2)
    print('V^2/Hz^2 = ',total_Vf)
    print('RMS^2 = ', rms_sq)
    noise_norm[i] = avg_fft


### write to a csv file
with open('data/LeakvsK_{}{}_dk.csv'.format(option,amp), 'w') as f:
   writer = csv.writer(f, delimiter='\t')
   writer.writerows(zip(dk, error**2))

# fig = plt.figure(figsize=(10,8))
# for i in range(len(k_in)):
#     plt.plot(k_in, error**2, 'o--', label='amp={}'.format(amp))
# plt.xlabel('k_in') ### not the actual noise level in V^2/Hz
# plt.ylabel('leakage level')
# plt.legend()
# plt.savefig('/home/gemma/Beam/leak_vs_k.png')
