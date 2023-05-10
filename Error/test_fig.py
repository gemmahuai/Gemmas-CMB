

from skimage.measure import block_reduce
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.linalg as la
import scipy.interpolate as interp
import scipy.special as sf
import matplotlib as mpl

fig = plt.figure(figsize=(20, 20), constrained_layout=True)
subfigs = fig.subfigures(1, 3, width_ratios=[1, 1, 1], wspace=.15)