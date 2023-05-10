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

n = int(sys.argv[1]) # number of amplitudes from 0.1 to 100
B = 10**np.linspace(-1,2,n)
# input number of plane waves
N_pw = 5
for b_i in B:
  pw = []
  screeni = {}
  screeni['N'] = 1024
  screeni['D'] = 10
  Fhcalc.Initialize(screeni)
  center = (screeni['D']/2, screeni['D']/2)
  Fhcalc.MultByGaussian(screeni, center, 1.0)
  Fhcalc.InCircle(screeni, center, 2.0)
  for j in range(N_pw):
    A_j = np.abs(np.random.normal(0,1)) # generate a random amplitude of a plane wave centered at 0, std=1
    k_j = np.abs(np.random.normal(11,1)) # generate a random k of a plane wave centered at 11, std=1
    pw_j = Plane_wave(k_j,screeni['N'],A_j,0)
    screeni['E'] = screeni['E'] * np.exp(1j*b_i*pw_j)
    if j == N_pw:
      break
  fig = plt.figure(figsize=(8,6))
  #plt.imshow(np.real(screen['E']))
  #plt.colorbar()
  #plt.show()
  #plt.rcParams['figure.figsize'] = [15, 10]
  ffti = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(screeni['E'])))
  plt.plot(screeni['X'][0], np.real(ffti[512]),label='real')
  plt.plot(screeni['X'][0], np.imag(ffti[512]),label='imag')
  plt.plot(screeni['X'][0], np.abs(ffti[512]),label='abs')
  plt.xlim(3.5,6.5)
  #plt.ylim(-8e3,8e3)
  plt.legend();
  plt.savefig('/home/gemma/Error/wave_amp_{}.png'.format(b_i))
