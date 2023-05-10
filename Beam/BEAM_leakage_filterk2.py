### 02/15/2023 updated
### this python script compares the power spectrum of the perturbed and unperturbed beams with varying k_in and k_out
### generates a plot of leakage level (error^2) as a function of error amplitude 
### with varying k_in and k_out of the noise annulus 
### calculates the normalized noise level in the fourier space as the error amplitude (x-axis) --> noise_norm

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
    print('Usage: python BEAM_leakage_filterk.py N_screen N_theta_interp screenD sigma angle(deg) truncation(y/n) phase/amp, filter_radius')
    # ex. python BEAM_leakage_filterk2.py 4096 1024 10.0 1.0 3.0 y amp 2.0
    sys.exit()
# input parameters
N_screen = np.int(sys.argv[1])
N_theta = np.int(sys.argv[2])
D = float(sys.argv[3])
sigma = float(sys.argv[4])
maxdeg = float(sys.argv[5])
trunc = str(sys.argv[6])
option = str(sys.argv[7])
#radius = float(sys.argv[8])
amp = np.array([0.1, 0.3, 0.6, 1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0])
#amp = np.array([1.0, 2.0]) # test
#k_in = np.array([1, 10]) # test
k_in = np.arange(1,52,5)
A = 150 # area of the noise filter in k space
error = np.zeros((len(k_in),len(amp))) #leakage = error^2
noise_norm = np.zeros_like(error)

def gaussian(x, A, sigma, x0): 
    g = A*np.exp(-(x-x0)**2/(2*sigma**2))
    return(g)

def analytical(l, sigma): # blm
    fft = bin_mean0.max()* np.exp(-l*(l+1)*sigma**2/2)
    return(fft)

def rms(x):
    return np.sqrt(np.mean(x**2))

# fix area of filter at A, given an inner radius, find the corresponding outer radius
def r_out(r_in, A):
    rout = np.sqrt((A/np.pi) + r_in**2)
    return(rout)


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

for j in range(len(k_in)):
    kin = k_in[j]
    kout = r_out(kin, A)
    
    for i in range(len(amp)):
        print('amp =', amp[i], 'k_in =', kin, 'k_out =', kout)
        
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
        #fft_numerical1 = fft_I1.flatten()
        fft_numerical_diff = fft_Idiff.flatten()
        bin_mean0, bin_edge, bin_num = binned_statistic(l_flatten, fft_numerical0, statistic='mean', bins=bin_edges) 
        #bin_mean1, bin_edge, bin_num = binned_statistic(l_flatten, fft_numerical1, statistic='mean', bins=bin_edges) # bin_mean is the binned numerical beam
        bin_mean_diff, bin_edge, bin_num = binned_statistic(l_flatten, fft_numerical_diff, statistic='mean', bins=bin_edges) # bin the beam difference
        l_vec = bin_edges[0:-1] # ell 1D vector

        beam_diff_rela = bin_mean_diff/bin_mean0 # relative beam difference
        #error[j,i] = np.mean(beam_diff_rela[1:3])
        error[j,i] = np.mean(beam_diff_rela[l_vec<500]) # average the beam difference for 0<l<500
        print('leakage (error^2) = {}'.format(error[j,i]**2))

        # Calculates the noise level in fourier space (V^2/Hz)
        emap_fft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(emap_E))) # scaled error map in fourier space
        norm_fft = np.abs(emap_fft**2)/screen1['dk']**2/screen1['N']**4 # V^2/Hz^2
        idx = np.where((screen1['kap']<kout) & (screen1['kap']>kin)) # inside the filter
        avg_fft = np.mean(norm_fft[idx]) # V^2/Hz^2
        total_Vf = np.sum(norm_fft)*screen['dk']**2
        rms_sq = np.abs(rms(emap_E)**2)
        print('V^2/Hz^2 = ',total_Vf)
        print('RMS^2 = ', rms_sq)
        noise_norm[j,i] = avg_fft


### write to a csv file
k_in_output = np.repeat(k_in, len(amp))
#radius_output = np.repeat(radius, len(amp))
amp_output = np.tile(amp, len(k_in))
error_output = (error**2).flatten()
noise_output = noise_norm.flatten()
with open('data_k_area150.csv', 'w') as f:
   writer = csv.writer(f, delimiter='\t')
   writer.writerows(zip(k_in_output, amp_output, noise_output, error_output))
