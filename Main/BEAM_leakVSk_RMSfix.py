### 04/22/2023 updated
### this python script fixes the RMS error in real space by varying the error scaling factor 
### and plots leakage level as a function of k_in

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
import Fhcalc
import ErrMask
import sys
import csv

# input & Initialize parameters
if (len(sys.argv)!=9):
    print('Wrong inputs!')
    print('Usage: python BEAM_leakVSk_RMSfix.py N_screen N_theta_interp screenD sigma angle(deg) truncation(y/n) "phase/amp" RMS_want')
    # ex. python BEAM_leakVSk_RMSfix.py 4096 1024 10.0 1.0 5.0 y amp 0.01
    # RMS_want: specify an RMS level which will be fixed in the calculation
    sys.exit()
# input parameters & initialize
N_screen = np.int(sys.argv[1])
N_theta = np.int(sys.argv[2])
D = float(sys.argv[3])
sigma = float(sys.argv[4])
maxdeg = float(sys.argv[5])
trunc = str(sys.argv[6])
option = str(sys.argv[7])
RMS_want = float(sys.argv[8])
#k_in = np.arange(5,30,5)
k_in = np.arange(1,72,2)
amp = 0.1 # the scaling factor C_1 or C_2. this is used for interpolation and generating Fig.18 in my thesis. May vary depending on the choice of RMS_want since we want C_1 or C_2 to roughly match the scaling factor corresponding to RMS_want
dk = 2 # width of the annular filter 
error = np.zeros(len(k_in))  # initialize the relative beam window difference
noise_norm = np.zeros_like(error) # initialize leakage level
#k_in = np.arange(0,screen1['kx'].max()-2,2)
k_out = k_in+dk
RMS = np.array([]) # initialize RMS of error maps used for generating Fig.18 
RMS_sq = np.array([])
RMS_test = np.array([]) # initialize RMS of error maps used for actual calculation of leakage (they should be fixed...)

### define functions
def gaussian(x, A, sigma, x0): 
    g = A*np.exp(-(x-x0)**2/(2*sigma**2))
    return(g)
def analytical(l, sigma): # blm
    fft = bin_mean0.max()* np.exp(-l*(l+1)*sigma**2/2)
    return(fft)
def linear(x, k, b):
    return(k*x+b)
def sqrt(x, k, b, c):
    """square root fit to the RMS vs. k_in"""
    return(k*np.sqrt(x+c)+b)

"""
line 65 to 101: generate error masks with different k_in and measure RMS of each. Then do a square root fit of RMS vs. k_in, which will be used later to fix RMS 
"""
### create a screen and measure RMS of each k_in
screen1 = {}
screen1['N'] = N_screen
screen1['D'] = D
Fhcalc.Initialize(screen1)
center = (screen1['D']/2, screen1['D']/2)
Fhcalc.MultByGaussian(screen1, center, sigma)
if trunc=='y':
    Fhcalc.InCircle(screen1, center, 2.0)
for i in range(len(k_in)):
    if option=='phase':
        emap_E = ErrMask.filter_annulus_phase(screen1, amp, k_in[i], k_out[i]) ### scaled error map in real space and unscaled error in k space
    elif option=='amp':
        emap_E = ErrMask.filter_annulus_amp(screen1, amp, k_in[i], k_out[i])
    rms_sq = np.abs(ErrMask.rms(emap_E)**2) #rms square
    RMS_sq = np.append(RMS_sq, rms_sq)
    rms_nosq = ErrMask.rms(emap_E) # rms
    RMS = np.append(RMS, rms_nosq)

### fit RMS vs. k_in
(fit1, err1) = opt.curve_fit(linear, k_in, RMS_sq, absolute_sigma=True) #rms square: fit a line
(fit2, err2) = opt.curve_fit(sqrt, k_in, RMS,absolute_sigma=True) # rms: fit sqrt
scaling = RMS_want * (amp/sqrt(k_in, fit2[0], fit2[1], fit2[2]))
sp = interp.InterpolatedUnivariateSpline(k_in, sqrt(k_in, fit2[0], fit2[1], fit2[2])) # spline for the RMS error
sp_scale = interp.InterpolatedUnivariateSpline(k_in, scaling) # spline for the scaling factor
### plot RMS, RMS^2 vs. k_in and fits
# fig = plt.figure(figsize=(14,5))
# plt.subplot(1,2,1)
# plt.plot(k_in, RMS, 'o', ms=5, color='black', label='data')
# plt.plot(k_in, sqrt(k_in, fit2[0], fit2[1], fit2[2]), label='fit')
# plt.xlabel(r'$k_{in}$')
# plt.ylabel('RMS in real space')
# plt.subplot(1,2,2)
# plt.plot(k_in, RMS_sq, 'o', ms=5, color='black', label='data')
# plt.plot(k_in, linear(k_in, fit1[0], fit1[1]), label='fit')
# plt.xlabel(r'$k_{in}$')
# plt.ylabel(r'$RMS^2$')



"""
line 109 to 198: 
    generate an unperturbed screen and a perturbed one. 
    The RMS level of the error mask used to perturbed the screen remains constant by using the sqrt fit derived above and scaling the factors C_1 or C_2.
    Then the leakage level is found as a function of k_in.
"""
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

for i in range(len(k_in)):
    kin = k_in[i]
    #kin = 5
    kout = kin + dk

    print('error amp =', sp_scale(kin), 'k_in =', kin, 'k_out =', kout)
    
    # perturbed E screen
    screen1 = {}
    screen1['N'] = N_screen
    screen1['D'] = D
    Fhcalc.Initialize(screen1)
    Fhcalc.MultByGaussian(screen1, center, sigma)
    if trunc=='y':
        Fhcalc.InCircle(screen1, center, 2.0)

    if option=='phase':
        emap_E = ErrMask.filter_annulus_phase(screen1, sp_scale(kin), kin, kout)
    elif option=='amp':
        emap_E = ErrMask.filter_annulus_amp(screen1, sp_scale(kin), kin, kout)
    else: print('Choose phase or amplitude errors')
    rms_test = ErrMask.rms(emap_E)
    RMS_test = np.append(RMS_test, rms_test) # check if RMS is kept constant
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
    error[i] = np.mean(beam_diff_rela[(l_vec>=50)&(l_vec<=250)]) # average the beam difference over 0<l<500
    print('leakage (error^2) = {}'.format(error[i]**2))

    # Calculates the noise level in fourier space (V^2/Hz)
    emap_fft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(emap_E))) # scaled error map in fourier space
    norm_fft = np.abs(emap_fft**2)/screen1['dk']**2/screen1['N']**4 # V^2/Hz^2
    idx = np.where((screen1['kap']<kout) & (screen1['kap']>kin)) # inside the filter
    avg_fft = np.mean(norm_fft[idx]) # V^2/Hz^2
    total_Vf = np.sum(norm_fft)*screen1['dk']**2
    rms_sq = np.abs(ErrMask.rms(emap_E)**2)
    print('V^2/Hz^2 = ',total_Vf)
    print('RMS^2 = ', rms_sq)
    noise_norm[i] = avg_fft


### write to a csv file
with open('./RMS{}_leak_k_{}.csv'.format(RMS_want, option), 'w') as f:
#with open('data/test.csv', 'w') as f:
   writer = csv.writer(f, delimiter='\t')
   writer.writerows(zip(k_in, error**2, scaling, RMS_test))
   # 1st col: inner radius k_in (fixing dk)
   # 2nd col: leakage level
   # 3rd col: scaling factor of error maps (C_1 or C_2) for a fixed given RMS level
   # RMS of each filtered error map (should be the same with <1% variation)

# fig = plt.figure(figsize=(10,8))
# for i in range(len(k_in)):
#     plt.plot(k_in, error**2, 'o--', label='amp={}'.format(amp))
# plt.xlabel('k_in') ### not the actual noise level in V^2/Hz
# plt.ylabel('leakage level')
# plt.legend()
# plt.savefig('/home/gemma/Beam/leak_vs_k_RMSfix.png')
