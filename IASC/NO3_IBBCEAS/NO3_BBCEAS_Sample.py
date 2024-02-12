#########################################################################################
# This script configures the DV401 CCD plugged via the CCI-001 PCI Card                 #
# to take FVB single scans or accumulations for real time trace gas                     #
# concentration analysis using SVD.                                                     #
#                                                                                       #
# IMPORTANT. PLEASE MAKE SURE YOU KNOW ABOUT THE WARNING REGARDING COOLING BELOW -20C.  #
# YOU MUST LET IT COOL TO >-20C (PREFERABLY 0C) USING COOLEROFF() BEFORE USING SHUTDOWN #
# AND KEEP THIS IN MIND IF THERE IS AN EXCEPTION OR EXIT TO DESKTOP EVENT               #
#                                                                                       #
# Created by Mixtli Campos on 22/2/2023                                                 #
# mcampos@ucc.ie                                                                        #
#########################################################################################

# Python packages
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import datetime as dt
import numpy as np
from time import sleep
from matplotlib.dates import DateFormatter

# Local
import configurations as conf
import andorfunctions as andor
import CESfunctionsJUNOx23 as cf

# The pyAndorSDK2 is a proprietary package from the ANDOR SDK
from pyAndorSDK2 import atmcd
from pyAndorSDK2 import atmcd_codes as codes
from pyAndorSDK2 import atmcd_errors as errors

#########################################################################################
#####                       PARAMETER CONFIGURATION                                 #####
# These parameters can be changed manually or with a configuration file

### Instrument 
temp = conf.temp                                    # Camera temperature
exptime = conf.exptime_sample                       # Exposure time in seconds
bckg_shots = conf.bckg_shots                        # Number of background shots
                                                    # (for averaging in analysis)
acqMode = conf.acqMode        
                                                    # Acquisition mode
                                                    # e.g. SINGLE_SCAN, ACCUMULATE
                                                    # check codes for more
accum_number = conf.accums                          # Number of accumulations (if needed)
accum_cycle = exptime + conf.delay                  # Exp + Delay = Cycle time
                                                    # (only for internal trigger)
readMode = conf.readMode                            # Read mode
trigMode = conf.trigMode                            # Trigger Mode

### Signal analysis
# Cavity parameters
distance = conf.distance                    # Sample optical length

# Resonance window 
lower_wavelength=conf.lower_wavelength      # Starting wavelength of resonance window
upper_wavelength=conf.upper_wavelength      # Ending wavelength of resonance window

# Reference and background files
back_filename = conf.back_filename
no3_refname = conf.no3_refname

# Start average from measure #
start_avg = conf.start_avg

# Reff : Either a number conf.Reff or a vector np.load(conf.Reff_matrix)
#Reff= conf.Reff
Reff = 0.99999
# Dilution factor --> SET TO 1 for IASC
dfactor = 1
#dfactor = 1-(conf.n2flow/conf.tflow)

### Path for saving data
savepath = conf.savepath

#########################################################################################
### Reference and Background file loading for analisis                                ###
no3reference = np.load(no3_refname)
background = np.load(back_filename)

#########################################################################################
##### Making a subdirectory for generated files %Y%m%d

directory = dt.datetime.now().strftime('%Y%m%d')
path = os.path.join(savepath,directory)   # path of the folder wherein to store data

# Tries to make a new folder with the current date, does nothing if folder already
# exists

try:
    os.mkdir(path)
except:
    pass

path_file = path + conf.folder_symbol   # full path to append filename when writing

#########################################################################################
#####               INSTRUMENT PREPARATION FOR SAMPLING                             #####
# Starts the instrument object and initializes it
# Stabilizes to the set temperature
# Configures acquisition parameters

sdk = atmcd()  # Load the atmcd library

ret = sdk.Initialize(r"c:\Program Files\Andor SDK\\")   # Initialize camera, path points
                                                        # to DETECTOR.INI
print("Function Initialize returned {}".format(ret))

if errors.Error_Codes.DRV_SUCCESS != ret:
    print("...Could not initialize camera with error {}, will exit".format(ret))
    sys.exit()

# Uncomment the following if you really want to see the serial number
#(ret, iSerialNumber) = sdk.GetCameraSerialNumber()
#print("Function GetCameraSerialNumber returned {} Serial No: {}".format(
#        ret, iSerialNumber))

# Configure the acquisition, lines outsourced to AndorFunctions.py
try:
    andor.prepare_temperature(sdk,temp)
except:
    sys.exit()

xpixels = andor.prepare_camera(sdk,acqMode,readMode,trigMode,
        accum_number,accum_cycle,exptime)


#########################################################################################
### Calculating the wavelengths with the calibration factors from configuration file  ###

wavelengths = cf.andor_calibrator(xpixels,*conf.calfactors)

#########################################################################################
#####                                   SAMPLING                                    #####
# Initializes an interactive plot
# Performs an acquisition loop
# Analyzes data
# Plots

### Initialize plot
fig = plt.figure()              # Figure initialization
ax1 = fig.add_subplot(211)      # Axes 1 : Signal
ax2 = fig.add_subplot(212)      # Axes 2 : Concentration 1
#ax3 = ax2.twinx()               # Axes 3 : Concentration 2

# Initialize empty plots
minwave,maxwave = cf.segment_indices(background[:,0:2],lower_wavelength,
            upper_wavelength)
ax1.set_ylim([0,500])
ax2.set_ylim([0,50])
#ax3.set_ylim([0,50])
xs = background[minwave:maxwave,0]
#print(xs.shape)
ys = [0] * xs
#print(ys.shape)
line, = ax1.plot(xs,ys,'-k')
line1, = ax1.plot(xs,ys,'-g')
line2, = ax2.plot(xs,ys,'-b')
#line3, = ax3.plot(xs,ys,'-r')
#t0 = dt.datetime.now() # testing for total elapsed time

### Initializing timestamp and concentration list, and measurement array
measurements = np.array(background[:,0]).reshape(len(background[:,0]),1)
meastime = []
meastime2 = []
ppbs = []

# Perform Acquisition loop as an animate function

def init_func():
    return line, line1, line2, 

def animate(i):
    global measurements,meastime,meastime2,ppbs
    # Perform Acquisition
    # Uncomment the print statements for verbosity
    print("Acquisition number",i)
    
    ret = sdk.StartAcquisition()
    #print("Function StartAcquisition returned {}".format(ret))
    
    #tt1=dt.datetime.now() # for testing wait delay 
    
    ret = sdk.WaitForAcquisition()
    #print("Function WaitForAcquisition returned {}".format(ret))
    
    #tt2=dt.datetime.now() # for testing wait delay
    #print((tt2-tt1).total_seconds())

    (ret, arr, validfirst, validlast) = sdk.GetImages16(1, 1, xpixels)
    #print("Function GetImages16 returned {} first pixel = {} size = {}".format(
    #    ret, arr[0], xpixels))
    #print(arr.shape)

    ### Calculating number density
    counts = np.copy(arr).reshape(len(arr),1)
    minwave,maxwave = cf.segment_indices(no3reference,lower_wavelength,
            upper_wavelength)
    bckg = np.copy(background[minwave:maxwave,:])
    no3ref = np.copy(no3reference[minwave:maxwave,:])
    I_sample = np.copy(counts[minwave:maxwave,:])
    I_0 = np.average(bckg[:,1:],axis=1).reshape(len(bckg),1)
    #print(bckg.shape,no3ref.shape,I_sample.shape,I_0.shape)
    pPa = 101335
    tK = 293.15
    
    ### This one does everything (see recursive_fit_2ref function in CESfunctions.py)
    try:
        alpha,fl,a,b,ndensity1 = cf.fit_alg_1A_it(I_sample, I_0, Reff, distance, 
        no3ref,pPa,tK,parameters=1)
    except Exception as e:
        print("fit_alg_1A failed with exception:")
        print(e)
        pass
    
    ### The timestamp for this measurement is now
    timenow = dt.datetime.now()
    stamp = timenow.strftime('%y%m%d%H%M%S')
    meastime2.append(timenow)

    ### Add sample to measurements array and save individual sample datafile
    measurements = np.concatenate((measurements,counts.reshape(len(counts),1)),axis=1)
    
    np.savetxt(path_file+'Im'+stamp+'.txt',measurements[:,[0,-1]],fmt='%s')

    ### Populate ppbs and meastime arrays with currents sample, make/overwrite datafile
    conc = (ndensity1*1e15*1.380649e-23*tK/pPa)
    ppbs.append(conc)
    meastime.append(timenow.strftime('%Y/%m/%d-%H:%M:%S'))
       
    np.savetxt(path_file+'Mtemp.txt',np.column_stack((meastime,ppbs)),fmt='%s')

    # Print calculated NO3 in ppb
    print('NO3 ppb: ', ppbs[-1])

    ### Plotting
    # Plot 1 : Axes 1
    # Plot 1 : Axes 1
    #print('Getting first plot')
    ax1.set_ylim([min(alpha),max(alpha)])
    line.set_ydata(alpha)
    line1.set_ydata(a+b*fl+no3ref[:,1]*ndensity1)
    
    # Plot 2 : Axes 2
    #print('Getting second plot')
    if i!=0:
        ax2.set_xlim([min(meastime2),max(meastime2)])
        ax2.set_ylim([min(ppbs),max(ppbs)])
        #ax3.set_ylim([min(ppbs2),max(ppbs2)])
    line2, = ax2.plot(meastime2,ppbs,'-b')
    #line3, = ax3.plot(meastime2,ppbs2,'-r')
    ax2.xaxis.set_major_formatter(DateFormatter('%H:%M'))
    #line2.set_ydata(ppbs)
    #line3.set_ydata(ppbs2)    
    return line, line1, line2, 

# call animation
ani = animation.FuncAnimation(fig,animate, init_func=init_func,
                              interval=1,blit=False,cache_frame_data=False)
plt.show()

t1 = dt.datetime.now()                  # End time
#print("Seconds elapsed: ",(t1-t0).total_seconds())

# We save all measurements in a numpy file
np.save(path_file + "Imeas" + t1.strftime("%y%m%d%H%M"), measurements)

# We save all concentrations in a datafile
np.savetxt(path_file + "M" + t1.strftime('%y%m%d%H%M') + '.txt',
        np.column_stack((meastime,ppbs)), fmt='%s')

#########################################################################################

#########################################################################################
#####                               SHUTTING DOWN                                   #####
# This part will prepare the camera for shutting down
# "Clears" the status
# Let's the camera reach around 0C
# Shuts down camera object to free the resource
# Lines outsourced to AndorFunctions.py

andor.shutdown_camera(sdk)

#########################################################################################
print('Shape of measurements array:',measurements.shape)
print('bye')
