# Effects of CMB Telescope Beam Systematics on Measurements of B-Mode Polarization

## Main - main calculation


### Fhcalc.py: 
Functions used for generating a screen, including initialization, Gaussian illumination, circular truncation, FFT, interpolation, etc.

### ErrMask.py: 
Functions used for generating errors, including amplitude and phase errors, and rms calculation

### ErrorMask.ipynb: 
Plots the annular noise filter and filtered noise, and compares unperturbed and perturbed screens.

### BEAM_test_new.py: 
Tests if numerical & analytical beam window functions are consistent (no errors yet)

### BEAM_error_difference.py: 
Compares the beam window difference between a perturbed and unperturbed one.

### BEAM_leakage.py: 
Generates a plot of leakage level (error^2) as a function of error amplitude for the fixed $k_{in}=15$ and $k_{out}=25$ and after noise normalization

### r_calc: 
A directory containing calculation of bias in r.

camb_96896687_totcls.dat.txt: CMB power spectrum downloaded from the CAMB webpage

RnNoise.py: Reads in the CMB text file 'camb_96896687_totcls.dat.txt', calculates the bias in r as a function of normalized noise, and outputs a text file including error amplitude ($C_1$ or $C_2$), normalized noise in k space, ell vectors, $\Delta r$, leakage power spectra, and RMS in real space. 

plot_r.py: Reads in the output .csv file from RnNoise.py and plots $\Delta r$ as a function of noise

run_amp.sh: shell script that runs RnNoise.py and plot_r.py altogether for amplitude perturbation.

run_phase.sh: shell script that runs RnNoise.py and plot_r.py altogether for phase perturbation.

### BEAM_leakage_filterk1.py:
Iterates through different k's


## From CMB computer:

Fraunhofer

Beam

Error

## My Jupyter notebooks:

Notebooks

## Google collab notebooks:

Collab
