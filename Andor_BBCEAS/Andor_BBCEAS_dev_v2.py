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

# Local
import configurations as conf
import AndorFunctions as andor
sys.path.append('..')
import py.CESfunctions_dev as cf

# The pyAndorSDK2 is a proprietary package from the ANDOR SDK
from pyAndorSDK2 import atmcd
from pyAndorSDK2 import atmcd_codes as codes
from pyAndorSDK2 import atmcd_errors as errors

#########################################################################################
#####                       PARAMETER CONFIGURATION                                 #####
# These parameters can be changed manually or with a configuration file

temp = conf.temp                                    # Camera temperature
exptime = conf.exptime                              # Exposure time in seconds
#shots = 5                                           # Number of acquisitions (deprecated)
acqMode = conf.acqMode        
                                                    # Acquisition mode
                                                    # e.g. SINGLE_SCAN, ACCUMULATE
                                                    # check codes for more
accum_number = conf.accums                          # Number of accumulations (if needed)
accum_cycle = exptime + conf.delay                  # Exp + Delay = Cycle time
                                                    # (only for internal trigger)
readMode = conf.readMode                            # Read mode
trigMode = conf.trigMode                            # Trigger Mode

#########################################################################################

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

# Configure the acquisition, lines outsourced to the AndorFunctions.py
try:
    andor.prepare_temperature(sdk,temp)
except:
    sys.exit()

xpixels = andor.prepare_camera(sdk,acqMode,readMode,trigmode,
        accum_number,accum_cycle,exptime)


#########################################################################################

#########################################################################################
#####                                   SAMPLING                                    #####
# Initializes an interactive plot
# Performs an acquisition loop
# Analyzes data
# Plots

### initialize plot
fig = plt.figure()              # Figure initialization
ax1 = fig.add_subplot(111)      # Axes 1 : Signal
ax1.set_ylim([0,500])
x_len = xpixels
xs = list(range(0,xpixels))
ys = [0] * x_len
line, = ax1.plot(xs,ys,'-k')
t0 = dt.datetime.now()

# Perform Acquisition loop as an animate function
# animate function

def animate(i):
    # Perform Acquisition
    # Uncomment the print statements for verbosity
    print("Acquisition number",i)
    
    ret = sdk.StartAcquisition()
    #print("Function StartAcquisition returned {}".format(ret))
    
    #tt1=dt.datetime.now() # for testing wait delay 
    
    ret = sdk.WaitForAcquisition()
    #print("Function WaitForAcquisition returned {}".format(ret))
    
    ### BEGIN status ready with loop
    #cont = 1
    #(ret, status) = sdk.GetStatus()
    #while status != errors.Error_Codes.DRV_IDLE:
    #    (ret, status) = sdk.GetStatus()
    #    cont += 1
    ### END status ready with loop
    #print('ran {} times'.format(cont))

    #tt2=dt.datetime.now() # for testing wait delay
    #print((tt2-tt1).total_seconds())

    (ret, arr, validfirst, validlast) = sdk.GetImages16(1, 1, xpixels)
    #print("Function GetImages16 returned {} first pixel = {} size = {}".format(
    #    ret, arr[0], xpixels))
    #print(arr.shape)

    ### Plotting every 10 samples
    #ttt1 = dt.datetime.now() # for testing plotting time
    if (i%10) == 0:
        ax1.set_ylim([min(arr)-10,max(arr)+10])
        line.set_ydata(arr)
    #ttt2 = dt.datetime.now() # for testing plotting time
    #print((ttt2-ttt1).total_seconds())

    return line,

# call animation
ani = animation.FuncAnimation(fig,animate,
                              interval=1,blit=True,cache_frame_data=False)
plt.show()

t1 = dt.datetime.now()
print("Seconds elapsed: ",(t1-t0).total_seconds())

#########################################################################################

#########################################################################################
#####                               SHUTTING DOWN                                   #####
# This part will prepare the camera for shutting down
# "Clears" the status
# Let's the camera reach around 0C
# Shuts down camera object to free the resource

(ret, status) = sdk.GetStatus()
print("Function GetStatus returned {} with status {}".format(ret,status))

while status == errors.Error_Codes.DRV_ACQUIRING:
    print("...Device still acquiring, will wait a bit")
    sleep(5)
    (ret, status) = sdk.GetStatus()
    print("Function GetStatus returned {} with status {}".format(ret,status))

#### BEGIN Returning to OC for shutdown so that it doesn't go beyond specs 
ret = sdk.CoolerOFF()
print("Function CoolerOFF returned {}".format(ret))

s = True
print('Waiting for temperature to reach 0 C')
while s:
    sleep(15)
    (ret, temperature) = sdk.GetTemperature()
    if temperature >=0:
        s = False
#### END Returning to 0C for shutdown

ret = sdk.ShutDown()
print("Function ShutDown returned {}".format(ret))

#########################################################################################

print('bye')
