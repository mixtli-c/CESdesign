##########################################################################################
# This is the companion script for the NO2 CRD/CEAS instrument. The "service" consists on
# listening to rs232 calls from the main sscript in the MIC1816 and running the CEAS in
# different measurement schemes.
#
# Created by Mixtli Campos on 02.06.2023
#
##########################################################################################

########## General packages
import serial, os, sys
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import datetime as dt
import numpy as np
from time import sleep
from matplotlib.dates import DateFormatter

########## Local
import configurations as conf
import andorfunctions as andor
import CESfunctionsJUNOx23 as cf

########## pyAndorSDK2 is a proprietary package from the ANDOR SDK
from pyAndorSDK2 import atmcd
from pyAndorSDK2 import atmcd_codes as codes
from pyAndorSDK2 import atmcd_errors as errors

############################## Measurement functions #####################################
def ceas_blank(meastype,accums,exposure,shots):
    '''
    This function executes the blank measurements.
    meastype = ['b','z']
    '''
    exptime = float(exposure)                      # Exposure time in seconds
    blnk_shots = int(shots)                      # Number of background shots
                                            # (for averaging in analysis)
    acqMode = conf.acqMode        
                                            # Acquisition mode
                                            # e.g. SINGLE_SCAN, ACCUMULATE
                                            # check codes for more
    accum_number = int(accums)                   # Number of accumulations (if needed)
    accum_cycle = exptime + conf.delay      # Exp + Delay = Cycle time
                                            # (only for internal trigger)
    readMode = conf.readMode                # Read mode
    trigMode = conf.trigMode                # Trigger Mode
    ### Path for saving data
    savepath = conf.savepath
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
    ### Prepare the camera
    xpixels = andor.prepare_camera(sdk,acqMode,readMode,trigMode,
        accum_number,accum_cycle,exptime)

    ### Calculating the wavelengths with the calibration factors from configuration file
    wavelengths = cf.andor_calibrator(xpixels,*conf.calfactors)
    
    ### Initialize plot
    plt.ion()                       # Interactive plot
    fig = plt.figure()              # Figure initialization
    ax1 = fig.add_subplot(111)      # Axes 1 : Signal
    xs = list(range(0,xpixels))     # x axis
    ys = [0] * xpixels              # y axis
    ax1.plot(xs,ys,'-k')            

    t0 = dt.datetime.now() # testing for total elapsed time

    ### Initialize measurement array
    measurements = np.copy(wavelengths).reshape(len(wavelengths),1)

    # Perform Acquisition loop as an interactive pyplot figure
    for i in range(blnk_shots):
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

        ### Plotting
        ax1.cla()
        ax1.plot(wavelengths,arr,'-k')
        fig.canvas.draw()
        fig.canvas.flush_events()

        ### Making arrays
        counts = np.copy(arr).reshape(len(arr),1)
        measurements = np.concatenate((measurements,counts),axis=1)

    t1 = dt.datetime.now()                  # End time
    plt.close('all')

    print("Seconds elapsed: ",(t1-t0).total_seconds()) #testing for total elapsed time
    # we generate a name to save the background
    blank_archive = "I" + meastype + t1.strftime("%y%m%d%H%M") +".txt"
    
    if meastype == 'b':
        measname = 'background'
    else:
        measname = 'zero'

    np.save(measname, measurements)     # for use by other functions
    
    np.savetxt(path_file + blank_archive, measurements)    # for further analysis

def ceas_measure(meastype,accums,exposure,shots):
    '''
    This function executes the sample air measurements
    NEEDS WORK!
    '''
    ### Instrument 
    exptime = float(exposure)                      # Exposure time in seconds
    meas_shots = int(shots)                      # Number of measurement shots
                                            # (for averaging in analysis)
    acqMode = conf.acqMode        
                                            # Acquisition mode
                                            # e.g. SINGLE_SCAN, ACCUMULATE
                                            # check codes for more
    accum_number = int(accums)                   # Number of accumulations (if needed)
    accum_cycle = exptime + conf.delay      # Exp + Delay = Cycle time
                                            # (only for internal trigger)
    readMode = conf.readMode                # Read mode
    trigMode = conf.trigMode                # Trigger Mode

    ### Signal analysis
    # Cavity parameters
    distance = conf.distance                    # Sample optical length

    # Resonance window  
    lower_wavelength=conf.lower_wavelength      # Starting wavelength of resonance window
    upper_wavelength=conf.upper_wavelength      # Ending wavelength of resonance window

    # Reference and background files
    back_filename = conf.back_filename
    no2_refname = conf.no2_refname
    zero_filename = conf.zero_filename

    # Start average from measure #
    start_avg = conf.start_avg
    
    # Reff : Either a number conf.Reff or a vector np.load(conf.Reff_matrix)
    Reff= conf.Reff

    # Dilution factor --> SET TO 1 for IASC
    dfactor = 1
    #dfactor = 1-(conf.n2flow/conf.tflow)

    ### Path for saving data
    savepath = conf.savepath

    ### Reference and Background file loading for analisis             
    no2reference = np.load(no2_refname)
    background = np.load(back_filename)
    zeroair = np.load(zero_filename)
    
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
    ## Prepare the camera and wavelength vector
    xpixels = andor.prepare_camera(sdk,acqMode,readMode,trigMode,
            accum_number,accum_cycle,exptime)
    wavelengths = cf.andor_calibrator(xpixels,*conf.calfactors)

    ### Initialize plot
    plt.ion()
    fig = plt.figure()              # Figure initialization
    ax1 = fig.add_subplot(211)      # Axes 1 : Signal
    ax2 = fig.add_subplot(212)      # Axes 2 : Concentration 1

    # Initialize empty plots
    minwave,maxwave = cf.segment_indices(no2reference,lower_wavelength,
                 upper_wavelength)
    xs = background[minwave:maxwave,0]
    #print(xs.shape)
    ys = [0] * xs
    #print(ys.shape)
    ax1.plot(xs,ys,'-k')
    ax1.plot(xs,ys,'-g')
    ax2.plot(xs,ys,'-b')

    t0 = dt.datetime.now() # testing for total elapsed time

    ### Initializing timestamp and concentration list, and measurement array
    measurements = np.array(no2reference[:,0]).reshape(len(no2reference[:,0]),1)
    meastime = []
    meastime2 = []
    ppbs = []

    # Perform Acquisition loop as interactive pyplot
    for i in range(meas_shots):
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
        minwave,maxwave = cf.segment_indices(no2reference,lower_wavelength,
                upper_wavelength)
        bckg = np.copy(background[minwave:maxwave,:])
        zero = np.copy(zeroair[minwave:maxwave,:])
        no2ref = np.copy(no2reference[minwave:maxwave,:])
        I_s = np.copy(counts[minwave:maxwave,:])
        I_b = np.average(bckg[:,1:],axis=1).reshape(len(bckg),1)
        I_not = np.average(zero[:,start_avg:],axis=1).reshape(len(zero),1)
        I_0 = np.subtract(I_not,I_b)
        I_sample = np.subtract(I_s,I_b)
    
        ### This one does everything (see CESfunctionsJUNOx23.py)
        try:
            alpha,fl,a,b,ndensity1 = cf.fit_alg_1A_it(I_sample, I_0, Reff, distance, 
            no2ref,parameters=1)
        except Exception as e:
            print("fit_alg_1A failed with exception:")
            print(e)
            continue
    
        ### The timestamp for this measurement is now
        timenow = dt.datetime.now()
        stamp = timenow.strftime('%y%m%d%H%M%S')
        meastime2.append(timenow)

        ### Add sample to measurements array and save individual sample datafile
        measurements = np.concatenate((measurements,counts.reshape(len(counts),1)),axis=1)
    
        np.savetxt(path_file+'Im'+stamp+'.txt',measurements[:,[0,-1]],fmt='%s')

        ### Populate ppbs and meastime arrays with currents sample, make/overwrite datafile
        ppbs.append((ndensity1/2.504e10)/dfactor)
        meastime.append(timenow.strftime('%Y/%m/%d-%H:%M:%S'))
       
        np.savetxt(path_file+'Mtemp.txt',np.column_stack((meastime,ppbs)),fmt='%s')

        # Print calculated NO2 in ppb
        print('NO2 ppb: ', ppbs[-1])

        ### Plotting
        ax1.cla()
        ax2.cla()
        # Plot 1 : Axes 1
        #print('Getting first plot')
        #ax1.plot(bckg[:,0],(I_0/I_sample)-1,'-k')
        ax1.plot(bckg[:,0],alpha,'-k')
        ax1.plot(bckg[:,0],a+b*fl+no2ref[:,1]*ndensity1,'-g')
    
        # Plot 2 : Axes 2
        #print('Getting second plot')
        ax2.plot(meastime2,ppbs,'-b')
        ax2.xaxis.set_major_formatter(DateFormatter('%H:%M'))
        fig.canvas.draw()
        fig.canvas.flush_events()
        

    plt.close('all')
    t1 = dt.datetime.now()                  # End time
    print("Seconds elapsed: ",(t1-t0).total_seconds())
    
    
    # We save all measurements in a numpy file
    np.save(path_file + "Im" + t1.strftime("%y%m%d%H%M"), measurements)

    # We save all concentrations in a datafile
    np.savetxt(path_file + "M" + t1.strftime('%y%m%d%H%M') + '.txt',
            np.column_stack((meastime,ppbs)), fmt='%s')


############################## End of Measurement functions ##############################

############################## CEAS and COM port initialization ##########################

### Starting the communication port
ceascom = serial.Serial('COM1',19200)

### Instrument parameters
temp = conf.temp                                    # Camera temperature

### CCD initialization

sdk = atmcd()  # Load the atmcd library
ret = sdk.Initialize(r"c:\Program Files\Andor SDK\\")   # Initialize camera, path points
print("Function Initialize returned {}".format(ret))

if errors.Error_Codes.DRV_SUCCESS != ret:
    print("...Could not initialize camera with error {}, will exit".format(ret))
    sys.exit()

# Configure the acquisition, lines outsourced to AndorFunctions.py
try:
    andor.prepare_temperature(sdk,temp)
except Exception as e:
    print('Will exit due to following error:',e)
    sys.exit()


############################## Listener and operation ####################################
while True:
    try:
        sleep(2)
        msg_in = ceascom.read(ceascom.in_waiting).decode('utf-8')
        while msg_in == '':
            sleep(0.5)
            msg_in = ceascom.read(ceascom.in_waiting).decode('utf-8')
    
        ceascom.write(b'k')
    
        params = msg_in.split(',')

        if params[0] == 'm':
            #ceas_measure(*tuple(params))
            try:
                ceas_measure(*tuple(params))
            except Exception as e:
                print('Error running measure')
                print(e)
                ceascom.write(b'e')
                continue
        else:
            #ceas_blank(*tuple(params))
            try:
                ceas_blank(*tuple(params))
            except Exception as e:
                print('Error running blank')
                print(e)
                ceascom.write(b'e')
                continue

        ceascom.write(b'd')
    except KeyboardInterrupt:
        print("Loop stopped by keyboard.")
        choice = input('R - Reload conf file. Otherwise will end')
        if choice.lower() == 'r':
            import configurations as conf
            continue
        else:
            break
    except Exception as e:
        print('Loop ended by exception:')
        print(e)

############################## SHUTDOWN ##################################################
andor.shutdown_camera(sdk)
ceascom.close()
