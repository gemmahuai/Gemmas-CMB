# Gemma's BS&T - Effects of CMB Telescope Beam Systematics on Measurements of B-Mode Polarization

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
Calculates the leakage level as a function of normalized noise in $k$ space, by scanning the inner radius of the annular noise filter $k_{in}$ from 1 to ~70 and fixing width $dk=2$. Generates a csv file "data_k_r2.csv", which is read and plotted in the notebook ErrorMap_norm.ipynb and Beam_norm0403.ipynb.

### BEAM_leakage_filterk3.py:
Calculates the leakage level as a function of normalized noise in $k$ space, by fixing $k_{in}$ and increasing radius of the annular noise filter $dk$ from 1 to ~50. Generates a csv file "data_k_kin{}.csv", which is read and plotted in the notebook ErrorMap_norm.ipynb and Beam_norm0403.ipynb.

### Beam_norm0403.ipynb:
Normalize unperturbed & perturbed beams by matching their sums. 

Also generates some thesis plots - amplitude/phase perturbed screen, beam window difference, leakage vs. noise (fixing noise in k space) for different $k_{in}$ values, leakage vs. $k_{in}$, leakage vs. $dk$.

### ErrorMap_norm.ipynb:
Containing some tests normalizing error maps based on the Parseval theorem. Also reads BEAM_leakage_filterk1.py and BEAM_leakage_filterk3.py output csv files and plots leakage level vs. $k_{in}$. 



## From CMB computer:

Fraunhofer

Beam

Error

## My Jupyter notebooks:

Notebooks

## Google collab notebooks:

Collab
