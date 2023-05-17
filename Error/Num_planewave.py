"""
usage: python Amp_planewave.py <number of screens> 
ex. python Amp_planewave.py 3
    3 means that three screens will be generated. For the first one, there will be a single wave added, and a figure is created plotting FT of the screen
    For the second one, there will be two waves added, and another figure created.
    For the third one, there will be three waves added, and another figure created. 
    Figure name is wave_num_plot_{}.png'.format(i) plotting the FT and wave_num_im_{}.png'.format(i) showing 2D image of FT
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('xtick', direction='in', top=True)
mpl.rc('ytick', direction='in', right=True)
mpl.rc('xtick.minor', visible=True)
mpl.rc('ytick.minor', visible=True)
plt.rcParams['font.size'] = 13
plt.rcParams['figure.figsize'] = [8,6]
import datetime as dt
import sys
import Fhcalc


def Plane_wave(k,size,A,theta):
  """A = amplitude
  theta = angle in radians"""
  x = np.linspace(0, 10, size)
  y = np.linspace(0, 10, size)
  [xx,yy] = np.meshgrid(x,y)
  kx = k * np.cos(theta)
  ky = k * np.sin(theta)
  plane_wave = A*np.sin(kx*xx + ky*yy)
  #plane_wave += plane_wave.max() # offset
  return(plane_wave)

B = 10
# input number of plane waves
n_pw = int(sys.argv[1])
N_pw = np.arange(1,n_pw) # 1, 2, 3, plane waves
for i in N_pw:
  pw = []
  screeni = {}
  screeni['N'] = 4096
  screeni['D'] = 10
  Fhcalc.Initialize(screeni)
  center = (screeni['D']/2, screeni['D']/2)
  Fhcalc.MultByGaussian(screeni, center, 1.0)
  Fhcalc.InCircle(screeni, center, 2.0)
  for j in range(i):
    A_j = np.abs(np.random.normal(0,1)) # generate a random amplitude of a plane wave centered at 0, std=1
    k_j = np.abs(np.random.normal(11,1)) # generate a random k of a plane wave centered at 11, std=1
    print('k value of {} waves is {}'.format(i,k_j))
    print('Amplitude of {} waves is {}'.format(i,A_j))
    pw_j = Plane_wave(k_j,screeni['N'],A_j,0)
    screeni['E'] = screeni['E'] * np.exp(1j*B*pw_j)
    if j == i:
      break
  fig = plt.figure(figsize=(8,6))
  #plt.imshow(np.real(screen['E']))
  #plt.colorbar()
  #plt.show()
  #plt.rcParams['figure.figsize'] = [15, 10]
  ffti = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(screeni['E'])))
  plt.plot(screeni['X'][0], np.real(ffti[2048]),label='real')
  plt.plot(screeni['X'][0], np.imag(ffti[2048]),label='imag')
  plt.plot(screeni['X'][0], np.abs(ffti[2048]),label='abs')
  plt.xlim(3.5,6.5)
  #plt.ylim(-8e3,8e3)
  plt.legend()
  plt.savefig('./wave_num_plot_{}.png'.format(i))

  fig = plt.figure(figsize=(18,6))
  #plt.rcParams['figure.figsize'] = [25, 7]
  plt.subplot(1,3,1)
  plt.imshow(np.real(ffti))
  plt.xlim(1900,2200)
  plt.ylim(1900,2200)
  plt.title('real')
  plt.colorbar()
  plt.subplot(1,3,2)
  plt.imshow(np.imag(ffti))
  plt.xlim(1900,2200)
  plt.ylim(1900,2200)
  plt.title('imag')
  plt.colorbar()
  plt.subplot(1,3,3)
  plt.imshow(np.abs(ffti))
  plt.xlim(1900,2200)
  plt.ylim(1900,2200)
  plt.title('abs')
  plt.colorbar()  
  plt.savefig('./wave_num_im_{}.png'.format(i))
