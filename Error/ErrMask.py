# this python script generates white noise in real space, fourier transforms to k space, applies an annulus filter, 
# and iFT to real space, returning our error mask. Then modify the amplitude of the original gaussian beam in function
# filter_annulus_amp, and the phase of the original gaussian beam in function filter_annulus_phase. Or, returns the error mask in
# function filter_WN

import numpy as np

def filter_annulus_amp(screen, A, k_in, k_out): 
    """
    k_in, k_out: filter radius in the fourier space
    A = amplitude of the amplitude error mask
    Modify screen['E']
    Return filtered WN in real space
    """
    N = screen['N']
    k = screen['kap'] # radius in k space
    # white noise
    WN = np.random.normal(0,1,(N,N))
    WN_FT = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(WN)))
    #filter
    cut1 = np.where(k<k_in, 0, 1)
    cut2 = np.where(k>k_out, 0, 1)
    WN_FT_fil = WN_FT*cut1*cut2 # filtered WN in fourier space
    WN_fil = np.real(np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(WN_FT_fil)))) # filtered WN in real space
    screen['E'] = (1 + A*WN_fil) * screen['E']
    return(A*WN_fil) #return amplitude
    

def filter_annulus_phase(screen, B, k_in, k_out):
    """
    k_in, k_out: filter radius in the fourier space
    B = amplitude of the phase error mask
    Modify screen['E']
    Return filtered WN in real space
    """
    N = screen['N']
    k = screen['kap']
    WN = np.random.normal(0,1,(N,N))
    WN_FT = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(WN)))
    cut1 = np.where(k<k_in, 0, 1)
    cut2 = np.where(k>k_out, 0, 1)
    WN_FT_fil = WN_FT*cut1*cut2 # filtered WN in fourier space
    WN_fil = np.real(np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(WN_FT_fil)))) # filtered WN in real space
    E_complex = screen['E']*np.exp(1j*B*WN_fil) 
    #screen['E'] = np.sqrt(np.real(E_complex**2) + np.imag(E_complex)**2)
    screen['E'] = E_complex
    return(B*WN_fil) 
    

def filter_WN(screen, A, k_in, k_out):
    """
    k_in, k_out: filter radius in the fourier space
    A = amplitude of error mask
    Returns the filtered WN in real space
    """
    N = screen['N']
    k = screen['kap'] # radius in k space
    # white noise
    WN = np.random.normal(0,1,(N,N))
    WN_FT = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(WN)))
    #filter
    cut1 = np.where(k<k_in, 0, 1)
    cut2 = np.where(k>k_out, 0, 1)
    WN_FT_fil = WN_FT*cut1*cut2 # filtered WN in fourier space
    WN_fil = np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(WN_FT_fil))) # filtered WN in real space
    return(WN_fil)
    
def rms(error_map):
    """
    Calculates the RMS of an error mask.
    """
    return np.sqrt(np.mean(error_map**2))