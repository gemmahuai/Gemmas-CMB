### updated on 03/17
### This script takes in the output csv file from RnNoise.py (data_r.csv)
### and makes a plot of the tensor-to-scalar ratio r as a function of normalized error
### r (4th column (idx 3rd)) vs. Error (2nd column (idx 1st))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.interpolate as interp
mpl.rc('xtick', direction='in', top=True)
mpl.rc('ytick', direction='in', right=True)
mpl.rc('xtick.minor', visible=True)
mpl.rc('ytick.minor', visible=True)
plt.rcParams['font.size'] = 15
import csv
import sys

# input csv file
if (len(sys.argv)!=3):
    print('Wrong inputs!')
    print('Usage: python plot_r.py csv_filename option(amp/phase')
    # ex. python plot_r.py R_noise_amp_r.csv amp
    sys.exit()
filename = str(sys.argv[1])
option = str(sys.argv[2])
(amp, noise_norm, r, RMS_real) = np.loadtxt('./{}'.format(filename), unpack=True, usecols=(0,1,3,5))

# remove repeaters --- merge 3 for loops
r_1d = np.array([])
for r_i in r:
    if r_i not in r_1d:
        r_1d = np.append(r_1d, r_i)
noise_1d = np.array([])
for noise_i in noise_norm:
    if noise_i not in noise_1d:
        noise_1d = np.append(noise_1d, noise_i)
amp_1d = np.array([])
for amp_i in amp:
    if amp_i not in amp_1d:
        amp_1d = np.append(amp_1d, amp_i)
rms_1d = np.array([])
for rms_i in RMS_real:
    if rms_i not in rms_1d:
        rms_1d = np.append(rms_1d, rms_i)

fig = plt.figure(figsize=(7,6))
ax1 = fig.add_subplot(111)
ax1.semilogy(rms_1d, r_1d, 'o--', ms=9, color='black', label='Simulation')
plt.axvline(x=0.0127, ls=':', lw=3, color='black')
plt.axhline(0.001, ls=':', lw=3, label='CMB-S4: r=0.001', color='black')

ax1.set_xlabel('RMS noise in real space')
ax1.set_ylabel('Bias in r')
#plt.ylim(1e-8, 1e-5)
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.legend()
ax2 = ax1.twiny()
ax2.set_xlim(noise_1d.min(), noise_1d.max())
ax2.set_xlabel('Noise in k space')
#plt.show()
plt.savefig('./fig_{}_r.png'.format(option))

# # plot power spectrum
# (ll, T, E, B) = np.loadtxt('camb_96896687_totcls.dat.txt', unpack=True, usecols=(0,1,2,3))
# fig = plt.figure(figsize=(10,7))
