import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import scipy.interpolate as spi
c = 3e8
import sys
np.set_printoptions(threshold=sys.maxsize)
import scipy.interpolate as interp

def Initialize(screen):
    # must run this first
    # It modifies the screen dictionary in place.
    # Before calling this, you have to first define 
    #   screen['N'] == number of pixels across the screen, and
    #   screen['D'] == the physical length across the screen.
    # These two values are used to populate various new elements in the dictionary 
    # associated with the screen pixels and k-values of the FFT.
    #
    n = screen['N']
    xvec = np.linspace(0,screen['D'],n)
    screen['X'] = np.tile(xvec,(len(xvec),1))
    screen['Y'] = screen['X'].T #np.flipud(screen['X'].T)
    screen['E'] = np.ones((n,n))
    dk = 1/screen['D'] #the separation in k-space corresponding to the separation in physical space # 2pi/D
    screen['dk']=dk
    kvec = dk*np.fft.fftfreq(n)*n #largest k corresponds to smallest physical resolution on the aperture
    screen['kvec'] = kvec
    screen['kx']= np.fft.fftshift(np.tile(kvec,(kvec.size,1)))
    screen['ky']= screen['kx'].T #np.flipud(screen['kx'].T)
    screen['kap']= np.sqrt(screen['kx']**2 + screen['ky']**2)  # k vector in aperture plane (#radial distance at each point)

def MultByGaussian(screen, center, sigma):
    # center must be a tuple, (xcenter,ycenter)
    x0 = center[0]
    y0 = center[1]
    R = np.sqrt((screen['X']-x0)**2 + (screen['Y']-y0)**2)
    screen['E']=screen['E']*np.exp(-(R**2)/(2*sigma**2))
    
def InCircle(screen,center,radius):
    # Must have called makeXY first.
    # center must be a tuple, (xcenter,ycenter) in meters
    # radius in meters
    x0 = center[0]
    y0 = center[1]
    R = np.sqrt((screen['X']-x0)**2 + (screen['Y']-y0)**2)
    cut_ap = np.where(R<radius,1,0) # multiplication factor: in-circle: 1; out-circle: 0
    screen['E']=screen['E']*cut_ap

def CircleAtten(screen,center,radius,factor):
    x0 = center[0]
    y0 = center[1]
    R = np.sqrt((screen['X']-x0)**2 + (screen['Y']-y0)**2)
    atten_ap = np.where(R<radius,factor,1) #factor: multiplication factor of in-circle.
    screen['E']=screen['E']*atten_ap
    
def ScreenFFT(screen):
    screen['FFT_E'] = np.fft.fftshift(np.fft.fft2(screen['E']))
    screen['I'] = np.abs(screen['FFT_E'])**2

def singleslit(screen,w,h,factor):
    x1 = screen['D']/2-w/2
    x2 = screen['D']/2+w/2
    y1 = screen['D']/2-h/2
    y2 = screen['D']/2+h/2
    slit = np.where((screen['X']>x1)&(screen['X']<x2)&(screen['Y']>y1)&(screen['Y']<y2), factor, 1)
    screen['E']=screen['E']*slit

def doubleslit(screen, w, h ,d, factor):
    x1l = screen['D']/2-d/2-w/2
    x1r = screen['D']/2-d/2+w/2
    x2l = screen['D']/2+d/2-w/2
    x2r = screen['D']/2+d/2+w/2
    y1 = screen['D']/2-h/2
    y2 = screen['D']/2+h/2
    xposition = ((screen['X']>x1l)&(screen['X']<x1r)) | ((screen['X']>x2l)&(screen['X']<x2r))
    yposition = ((screen['Y']>y1)&(screen['Y']<y2))
    slits = np.where((xposition&yposition), factor, 1)
    screen['E']=screen['E']*slits

def CDpanelgap_mask(screen):
    # SO secondary (from image) has 9 panels in "x" across about a 6m aperture, 
    # so they are 0.667 meters in size.
    # We build a 5m mirror with similar panel sizes here, not trying to be exact.
    #
    slit_centers1 = np.arange(0.0,2.5,0.667)
    slit_centers2 = -slit_centers1[1:]
    slit_centers = np.concatenate((slit_centers2[::-1],slit_centers1))
    print(slit_centers)
    slit_width = 0.0007 #mm
    #
    pix_size = screen['D']/screen['N']
    xvec = np.linspace(-screen['D']/2,screen['D']/2,screen['N'])  # centered at zero
    onerow = np.ones(len(xvec))
    
    if pix_size >= slit_width:
        print("Pixel size >= slit width:  Removing one column/row for each slit.")
        for ii in range(len(slit_centers)):
            slitlocator = np.abs(xvec - slit_centers[ii])
            slitloc = np.argmin(slitlocator)
            onerow[slitloc] = 0
    else:
        print("Pixel size <= slit width:  Removing multiple columns/rows for each slit.")
        for ii in range(len(slit_centers)):
            slitlocator = np.abs(xvec - slit_centers[ii])
            slitloc = np.where(slitlocator < (pix_size/1.9))
            onerow[slitloc] = 0
    
    plt.plot(onerow)
    gapmask_x = np.tile(onerow,(len(onerow),1)) 
    gapmask_y = gapmask_x.T
    gapmask = gapmask_x*gapmask_y
    #screen['E']=screen['E']*gapmask
    return gapmask

def Project_I_on_thetagrid(theta_vec,screen,lam):
    # Run after you've calculated the intensity as a function of kx, ky.
    # thetavec = 1D array of positions to be used for thetax, and for thetay
    #
    # Make 2D arrays of thetax, thetay coordinate for our map
    # We use the "_grid" suffix to indicate things that are 2D arrays associated
    # with the (new) thetax_grid, thetay_grid created next.
    N_thetagrid = len(theta_vec)
    thetax_grid, thetay_grid = np.meshgrid(theta_vec,theta_vec)
    #thetay_grid = np.flipud(thetay_grid)
    #
    kphot = 1/lam    # photon k vector. # 2pi/lambda
    #
    # Find kx and ky at each (thetax,thetay) grid spot
    #
    #kx_grid = kphot*thetax_grid/(np.pi/2) # fix this
    #ky_grid = kphot*thetay_grid/(np.pi/2)
    kx_grid = kphot*np.sin(thetax_grid) 
    ky_grid = kphot*np.sin(thetay_grid)
    #
    # Now look at the original kx,ky map of our FT screen.
    #
    I_grid = np.zeros((N_thetagrid, N_thetagrid))
    
    dk = screen['kx'][0,1]-screen['kx'][0,0]
    for xx in range(N_thetagrid):
        for yy in range(N_thetagrid):
            kx = kx_grid[yy,xx]
            ky = ky_grid[yy,xx]
            #
            nx = kx/dk + screen['N']/2  # what element is this in the screen's kx,ky
            ny = ky/dk + screen['N']/2
            #
            # find the four values of kx,ky that surround
            nx1 = int(np.floor(nx))
            nx2 = int(np.ceil(nx))
            ny1 = int(np.floor(ny))
            ny2 = int(np.ceil(ny))
            # Go through four points surrounding our grid point's position and take weighted average.
            numsum = 0
            denomsum = 0
            for pt in [(nx1,ny1), (nx1,ny2), (nx2,ny1), (nx2,ny2)]:
                rr2 = (kx - screen['kx'][pt[0],pt[1]])**2 + (ky - screen['ky'][pt[0],pt[1]])**2
                w = 1/rr2
                numsum += w*screen['I'][pt[0],pt[1]]
                denomsum += w
                #
                #numsum += screen['I'][pt[0],pt[1]]
                #denomsum += 1
            I_grid[xx,yy] = numsum/denomsum
            #I_grid[xx,yy] = screen['I'][nx2,ny1]
        
    return I_grid

