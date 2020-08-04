#!/usr/bin/env python

# Author: Dogacan S. Ozturk

# Import default Python libraries.
import os
import sys
from glob import glob
import tables
import numpy as np
import datetime as dt
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.dates as mdates
import matplotlib.ticker as mtick
myFmt = mdates.DateFormatter('%H:%M')

# Import custom Python libraries.
sys.path.insert(0, '../Code/')
from spacepy.pybats import gitm
import hime_helper_functions
from downsample_data import downsample_pfisr_data
from merge_potentials import merge_pfisr_with_gitm_potentials

# Enter filename for the PFISR 2D VEF estimates.
filename = '../Examples/Files/PFISR_Data/20191026.002_lp_1min-fitcal_2dVEF_001001-geo600km.h5'

# Enter desired grid resolution.
gridRes = 0.75

# Downsample the grid and calculate the potential differences.
PhiX, PhiY, Ex_downsampled, Ey_downsampled, Ex_calculated, Ey_calculated, XnewGrids, YnewGrids, experimentTimes = downsample_pfisr_data(filename, gridRes)

# Define the path to global potential values.
weimerSimulationList = glob('../Examples/Files/Simulations/3D*.bin')

#  Define the merge parameter.
mergeParameter = 0.6

# Set plot potentials to True for saving plots.
plotPotentials = True

# Set save potentials to True for saving output ASCII files.
savePotentials = True

# Merge the local and global potentials together.
phiXhime, phiYhime, himeEx, himeEy, xHimeMesh, yHimeMesh, himeTimes = merge_pfisr_with_gitm_potentials(PhiX, PhiY, XnewGrids, YnewGrids, experimentTimes, weimerSimulationList, gridRes, mergeParameter, plotPotentials, savePotentials)
