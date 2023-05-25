# Gemma's BS&T - Effects of CMB Telescope Beam Systematics on Measurements of B-Mode Polarization

## Main - main calculation


- ### Fhcalc.py: 
    Functions used for generating a screen, including initialization, Gaussian illumination, circular truncation, FFT, interpolation, etc.

- ### ErrMask.py: 
    Functions used for generating errors, including amplitude and phase errors, and rms calculation

- ### ErrorMask.ipynb: 
    Plots the annular noise filter and filtered noise, and compares unperturbed and perturbed screens.

- ### BEAM_test_new.py: 
    Tests if numerical & analytical beam window functions are consistent (no errors yet)

- ### BEAM_error_difference.py: 
    Compares the beam window difference between a perturbed and unperturbed one.

- ### BEAM_leakage.py: 
    Generates a plot of leakage level (error^2) as a function of error amplitude for the fixed $k_{in}=15$ and $k_{out}=25$ and after noise normalization

- ### BEAM_LeakVSk.py:
    Fixes the error amplitude (noise level in k space), and plots leakage level as a function of $k_{in}$. Also generates a csv file LeakvsK_{}{}.csv'.format(option,amp). 

- ### r_calc: 
    A directory containing calculation of bias in r.

    - camb_96896687_totcls.dat.txt: CMB power spectrum downloaded from the CAMB webpage

    - RnNoise.py: Reads in the CMB text file 'camb_96896687_totcls.dat.txt', calculates the bias in r as a function of normalized noise, and outputs a text file including error amplitude ($C_1$ or $C_2$), normalized noise in k space, ell vectors, $\Delta r$, leakage power spectra, and RMS in real space. 

    - plot_r.py: Reads in the output .csv file from RnNoise.py and plots $\Delta r$ as a function of noise

    - run_amp.sh: shell script that runs RnNoise.py and plot_r.py altogether for amplitude perturbation.

    - run_phase.sh: shell script that runs RnNoise.py and plot_r.py altogether for phase perturbation.

- ### BEAM_leakage_filterk1.py:
    Calculates the leakage level as a function of normalized noise in $k$ space, by scanning the inner radius of the annular noise filter $k_{in}$ from 1 to ~70 and fixing width $dk=2$. Generates a csv file "data_k_r2.csv", which is read and plotted in the notebook ErrorMap_norm.ipynb and Beam_norm0403.ipynb.

- ### BEAM_leakage_filterk3.py:
    Calculates the leakage level as a function of normalized noise in $k$ space, by fixing $k_{in}$ and increasing radius of the annular noise filter $dk$ from 1 to ~50. Generates a csv file "data_k_kin{}.csv", which is read and plotted in the notebook ErrorMap_norm.ipynb and Beam_norm0403.ipynb.

- ### Beam_norm0403.ipynb:
    Normalize unperturbed & perturbed beams by matching their sums. 

    Also generates some thesis plots - amplitude/phase perturbed screen, beam window difference, leakage vs. noise (fixing noise in k space) for different $k_{in}$ values, leakage vs. $k_{in}$, leakage vs. $dk$.

- ### ErrorMap_norm.ipynb:
    Contains some tests normalizing error maps based on the Parseval theorem. Also reads BEAM_leakage_filterk1.py and BEAM_leakage_filterk3.py output csv files and plots leakage level vs. $k_{in}$. 

- ### BEAM_leakVSk_RMSfix.py:
    By fixing RMS level in real space, this script calculates the leakage level as a function of $k_{in}$, and outputs to RMS{}_leak_k_{}.csv, including $k_{in}$, leakage level, scaling factors used for normalization, and RMS, which is read and plotted in Leak_k_fixRMS.ipynb.

- ### Leak_k_fixRMS.ipynb:
    Includes some initial calculation of normalizing RMS errors. Also generates figures of leakage vs. $k_{in}$ (Fig.12) and Fig.18 in my thesis.

- ### PowerSpectrum.ipynb:
    Generates leakage power spectra with measured $\Delta r$ labeled. 

- ### thesis_plot.ipynb:
    Generates some figures used in my thesis: Fig.13, Fig.17, Fig.16, Fig.14


## From CMB computer:

- ### Fraunhofer
    A directory containing basic Fraunhofer calculation and visualization.

    - Fhcalc.py: copied to Main

    - Fhplot.py: copied to Main

    - Fhcalc_torch.py: pytorch version of Fraunhofer calculation, may be faster...

    - Fhplot_torch.py: pytorch version of Fraunhofer visualization, may be faster...

    - test_torch_E0.png / test_torch_E1.png: generated using the pytorch version.

- ### Error
    A directory including some error mask calculation and tests.

    - ErrMask.py: copied to Main

    - Fhcalc.py: same as before

    - Amp_planewave.py: To understand why the on-sky beam of a phase perturbed screen has shell structures by adding sine waves with different frequency to a screen. This script takes a slice and generates a plot of the FT of a screen with multiple sine waves of different k added. output images: 'wave_amp_{}.png'.format(b_i).

    - Num_planewave.py: Similar to Amp_planewave.py but generates multiple screens with added sine waves of different k. Output images: wave_num_plot_{}.png'.format(i) plotting the FT and wave_num_im_{}.png'.format(i) showing 2D image of FT.

    - Num_planewave_single.py: Same as Num_planewave.py except that only one screen is generated with specified number of waves added. Output images: singlewave_num_plot_{}.png'.format(n_pw) plotting the FT and singlewave_num_im_{}.png showing 2D image of the FT

 - ### Beam
    A directory including beam calculations.

    - data: a directory including all csv data files. 

    - r: copied to Main (a directory of $\Delta r$ calculation).

    - Fhcalc.py: Fraunhofer calculation

    - ErrMask.py: Error mask calculation

    - BEAM_error_difference.py: copied to Main

    - BEAM_error_new.py: plots the beam window difference in $\ell$ space by taking FT of the perturbed and unperturbed on-sky beams and differencing in $\ell$ space, as a test version. For my thesis, I took the on-sky beam difference and then FT the difference. 

    - BEAM_error.py: generates a perturbed and unperturbed beam, and compares them with analytical solution of a Gaussian beam window. No matched in this version...

    - BEAM_test.py: tests if the numerical and analytical solution of a Gaussian beam agrees - not really in this version.

    - BEAM_test_compare.py: this python script reads several text files containing numerical solution (1st column) using different interpolation angles and analytical solution (2nd column), and plots all beam windows in the same figure for comparison.

    - BEAM_test_new.py: copied to Main

    - BEAM_leakage.py: copied to Main

    - BEAM_LeakVSk.py: copied to Main

    - BEAM_leakage_filterk1.py: copied to Main

    - BEAM_leakage_filterk2.py: studies the leakage level as a function of error scaling factor $C_1$ or $C_2$ with different filter radii. $k_{in}$ increases while fixing the total area of the annulus and finding the corresponding $k_{out}$. Generates a csv file, including $k_{out}$, $C_1$ or $C_2$, leakage level, and normalized noise level. 

    - BEAM_leakage_filterk3.py: copied to Main

    - BEAM_leakVSk_dk.py: fixes the noise level in k space, and plots the leakage level as a function of the annular width $dk$, fixing $k_{in}$. Generates a csv file including width and leakage.

    - BEAM_leakVSk_RMSfix.py: copied to Main

    - BEAM_leakVSk_RMSfix_dk.py: fixes the noise level in k space, and plots leakage level as a function of the width $dk$, fixing $k_{in}$. Generates a csv file. Not very helpful...

    - interp_test.py: tests the difference between FT and iFT of an on-sky beam...

    - run.sh: runs BEAM_error_difference.py with different error amplitudes, given phase or amp option.

## Notebooks - Jupyter notebooks on my local machine

- Data: A directory containing output data files and CAMB power spectra. These are plotted in some notebooks.

- ErrorMask.ipynb: copied to Main

- ErrorMap_norm.ipynb: copied to Main

- Beam_norm0403.ipynb: copied to Main

- Leak_k_fixRMS.ipynb: copied to Main

- PowerSpectrum.ipynb: copied to Main

- thesis_plot.ipynb: copied to Main

- ApertureDiffraction.ipynb: the initial Fraunhofer calculation notebook written by Prof. Ruhl.

- Fraunhofer.py: same as Fhcalc.py, copied to Main

- ErrMask.py: copied to Main

- BEAM_error_difference.ipynb: compares the beam window difference between a perturbed and unperturbed screen. Also generates a plot of leakage level as a function of error amplitude for a fixed k_in and k_out (before noise normalization)... 

- CMB_SummerSchool.ipynb: some codes from https://github.com/CMB-S4/CMBAnalysis_SummerSchool

- Fraunhofer.ipynb: includes a Fraunhofer calculation I wrote for fun when taking optics in the spring of 2022, and some tests of a screen perturbed by a dot representing local boost/reduction.

- Interpolation_test.ipynb: tests the interpolation algorithm by identifying the 1st, 2nd, ... minima and compare those to the expected values - $1.22\frac{\lambda}{D}$, $2.23\frac{\lambda}{D}$, $3.24\frac{\lambda}{D}$. Matched.

- Pytorch.ipynb: tests FT using pytorch for optimization...

- scratch.ipynb: play with Fraunhofer calculation...



## Collab - Google collab notebooks:

- the initial Fraunhofer calculation notebook written by Prof. Ruhl, with some tests.

- beam difference.ipynb: checks the interpolation centering, beam verification; compares the FT difference...

- beam_test.ipynb: try to get the numerical solution agree with the analytical solution. But they don't match in this notebook. 

- beam_error.ipynb: more tests to verify the beam calculation. no success in this notebook either.

- blm_error.ipynb: more tests. Finally matched at the very end!

- CPU_benchmark.ipynb: compares performace of CPU/GPU Fraunhofer calculation. Contains pytorch version of all Fraunhofer functions. 

- error_test.ipynb: adding plane waves to amplitude/phase, and study the FT behavior. 

- FraunhoferCalc.ipynb: some tests regarding our Fraunhofer calculation.

- Interpolation_test.ipynb: tests the interpolation algorithm by identifying minima. Also tests beam verification (failed...) with different screen parameters.

- PhaseShift.ipynb: tests a screen with plane wave added. 

- test_gpu.ipynb: test gpu implementation with cupy.

- Untitled0.ipynb: Fraunhofer function I wrote for the optics class...



