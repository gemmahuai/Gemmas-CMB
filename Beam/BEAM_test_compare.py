### 01/31/2023 updated
### this python script reads several text files containing numerical solution (1st column) and analytical solution (2nd column)
### and plot a single power spectrum including all files.


import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.optimize as opt
mpl.rc('xtick', direction='in', top=True)
mpl.rc('ytick', direction='in', right=True)
mpl.rc('xtick.minor', visible=True)
mpl.rc('ytick.minor', visible=True)
plt.rcParams['font.size'] = 10
plt.rcParams['figure.figsize'] = [8,6]
import sys
import Fhcalc

# input
if (len(sys.argv)!=3):
    print('Wrong inputs!')
    print('Usage: python BEAM_test_compare.py N_theta_interp angle(deg)')
    sys.exit()

N_theta = np.int(sys.argv[1])
maxdeg = np.float(sys.argv[2]) # in degrees
#screen = [5.0, 10.0, 15.0, 20.0]
#sigma = [0.1, 0.5, 1.0, 2.0]
angle = [1.0, 2.0, 3.0, 4.0]
colors = ['red', 'orange', 'green', 'blue']


#generate coordiantes
lam = 0.002 #mm wavelength
thetamax = maxdeg*np.pi/180. # degrees to radians
# shift the beam from being centered at theta=0 to theta=thetamax so that the beam spans from 0 deg to 2*thetamax deg
theta_vec = np.linspace(0,2*thetamax,N_theta) #rad
#calculate ell
dl = 2*np.pi/theta_vec.max() # dl in 1/rad space
l_vec = np.fft.fftshift(dl * np.fft.fftfreq(N_theta)*N_theta)
(l_x, l_y) = np.meshgrid(l_vec,l_vec) # 1/rad 
l = np.sqrt(l_x**2 + l_y**2)
bin_edges = np.linspace(0,l.max(),int(len(l_vec)/2))

fig = plt.figure(figsize=(10,4))
plt.subplot(1,2,1)
for i in range(len(angle)):
    (numerical, analytical) = np.loadtxt('/home/gemma/Beam/output_angle{}.csv'.format(angle[i]), unpack=True)
    # plot power spectrum -- analytical & numerical solutions
    thetamax = angle[i]*np.pi/180. 
    theta_vec = np.linspace(0,2*thetamax,N_theta) #rad
    dl = 2*np.pi/theta_vec.max() # dl in 1/rad space
    l_vec = np.fft.fftshift(dl * np.fft.fftfreq(N_theta)*N_theta)
    (l_x, l_y) = np.meshgrid(l_vec,l_vec) # 1/rad 
    l = np.sqrt(l_x**2 + l_y**2)
    bin_edges = np.linspace(0,l.max(),int(len(l_vec)/2))
    plt.plot(bin_edges[:-1], numerical/numerical.max(), color=colors[i], label='num_sig{}'.format(angle[i]))
    plt.plot(bin_edges[:-1], analytical/numerical.max(), ls=':', lw=2, color=colors[i], label='ana_sig{}'.format(angle[i]))
    plt.xlabel(r'$\ell$')
    plt.xlim(0,13000)
    plt.legend()
plt.subplot(1,2,2)
for i in range(len(angle)):
    (numerical, analytical) = np.loadtxt('/home/gemma/Beam/output_angle{}.csv'.format(angle[i]), unpack=True)
    thetamax = angle[i]*np.pi/180. 
    theta_vec = np.linspace(0,2*thetamax,N_theta) #rad
    dl = 2*np.pi/theta_vec.max() # dl in 1/rad space
    l_vec = np.fft.fftshift(dl * np.fft.fftfreq(N_theta)*N_theta)
    (l_x, l_y) = np.meshgrid(l_vec,l_vec) # 1/rad 
    l = np.sqrt(l_x**2 + l_y**2)
    bin_edges = np.linspace(0,l.max(),int(len(l_vec)/2))
    plt.plot(bin_edges[:-1], numerical/numerical.max(), color=colors[i], label='num_sig{}'.format(angle[i]))
    plt.plot(bin_edges[:-1], analytical/numerical.max(), ls=':', lw=2, color=colors[i], label='ana_sig{}'.format(angle[i]))
    plt.xlabel(r'$\ell$')
    plt.xlim(0,1000)
    plt.ylim(0.96, 1.0)
    plt.legend()
plt.savefig('/home/gemma/Beam/beam_compare_angle.png')

