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
savepath = "C:\\CRAC\\Data\\Instruments\\HONO_IBBCEAS"


### The following are instrument configurations (used in both Blank and Sample)
# Acquisition mode (SINGLE_SCAN, ACCUMULATE ... check codes for more)
acqMode = codes.Acquisition_Mode.ACCUMULATE

# Read Mode
readMode = codes.Read_Mode.FULL_VERTICAL_BINNING

# Trigger Mode
trigMode = codes.Trigger_Mode.INTERNAL

# Wavelength calibration factors (we calculate wavelengths from pixels ourselves)
#FOR NO2 at 405 (JUNOx23)
#calfactors=(3.10322588e+02-3.3, 1.72905279e-01,
#            -1.48816626e-05,1.95619706e-09)
#FOR NO3 11/09/2023
#calfactors=(6.9664720862e+02, -1.2308600988e-01,
#            -8.6814165183e-06,8.0507623387e-10)
#FOR HONO/NO2 29/09/2023
# note: cannot use calibrator for this, needs simple offset
offset=0.35


# Number of accumulations 
accums = 100

# Accumulation cycle delay (Exp + Delay = Cycle time)
# if your exposure time is short, this can add a lot of overhead so choose a number that is small (5-10 miliseconds maybe)
delay = 0.005

# Exposure time. 
exptime_sample = .1
exptime_blank = .1

# Camera temperature
temp = -30

# Measurement to start average
start_avg = 2

# Number of backgrounds to take
bckg_shots = 10

### The following are signal analysis parameters (not used in Blank, only Sample)
# Distance, the optical length of the sample (this is a physcal cavity parameter)
distance = 500

# The resonance window, lower (start) and upper (end) wavelengths
lower_wavelength = 364
upper_wavelength = 390

# Reference and background filenames to load, they should be located in the local dir
back_filename = "background.npy"
no2_refname = "NO2_IASC_1.npy"
no3_refname = "NO3_IASC.npy"
hono_refname = "HONO_IASC.npy"
zero_filename = "zero.npy"

# Effective Reflectivity, check lab notebook on more information, can be a vector 
# to call using np.load, or it can be a constant
Reff = np.load('Reff.npy')
#Reff = 0.998182
scale_index = 58

# Dilution factor, units are sccm
tflow = 5200
n2flow = 0 


