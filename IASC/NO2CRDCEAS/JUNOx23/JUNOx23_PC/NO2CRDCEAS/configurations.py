#########################################################################################
### This is the USER CONFIGURATION file
### Variables will be stored here to be called by main BBCEAS routines (blank and sample)
### Comment header has a short description of variables

import numpy as np

### Need to import the codes for the Andor Camera
from pyAndorSDK2 import atmcd_codes as codes

### Both Blank and Sample store data to a subdirectory, windows and linux system use 
### different symbols, "\\" for windows, "/" for linux
folder_symbol = "\\"
savepath = "C:\\JUNOx23_PC\\Data"


### The following are instrument configurations (used in both Blank and Sample)
# Acquisition mode (SINGLE_SCAN, ACCUMULATE ... check codes for more)
acqMode = codes.Acquisition_Mode.ACCUMULATE

# Read Mode
readMode = codes.Read_Mode.FULL_VERTICAL_BINNING

# Trigger Mode
trigMode = codes.Trigger_Mode.INTERNAL

# Wavelength calibration factors (we calculate wavelengths from pixels ourselves)
calfactors=(3.10322588e+02-3.3, 1.72905279e-01,
            -1.48816626e-05,1.95619706e-09)

# Number of accumulations 
accums = 3

# Accumulation cycle delay (Exp + Delay = Cycle time)
delay = 0.5

# Exposure time. 
exptime_sample = 5
exptime_blank = 1

# Camera temperature
temp = -30

# Measurement to start average
start_avg = 2

# Number of backgrounds to take
bckg_shots = 10

### The following are signal analysis parameters (not used in Blank, only Sample)
# Distance, the optical length of the sample (this is a physcal cavity parameter)
distance = 44.17

# The resonance window, lower (start) and upper (end) wavelengths
lower_wavelength = 395
upper_wavelength = 420

# Reference and background filenames to load, they should be located in the local dir
back_filename = "background.npy"
no2_refname = "NO2_JUNOx23.npy"
chocho_refname = "CHOCHO_JUNOx23.npy"
zero_filename = "zero.npy"

# Effective Reflectivity, check lab notebook on more information, can be a vector 
# to call using np.load, or it can be a constant
Reff = np.load('Reff.npy')
scale_index = 58

# Dilution factor, units are sccm
tflow = 5200
n2flow = 0 


