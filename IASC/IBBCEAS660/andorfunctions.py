# These are some functions for the Andor CCI-001 and cameras that use it
# Mostly to simplify the main script.

from pyAndorSDK2 import atmcd
from pyAndorSDK2 import atmcd_codes as codes
from pyAndorSDK2 import atmcd_errors as errors
from time import sleep

def take_measurement(sdk,xpixels):
    '''
    Takes sdk object and xpixels, starts acquisition, waits until done, gets array
    and returns as arr
    '''
    ret = sdk.StartAcquisition()
    
    ret = sdk.WaitForAcquisition()

    (ret,arr,validfirst,validlast) = sdk.GetImages16(1,1,xpixels)

    return arr

def prepare_temperature(sdk,temp):
    '''
    Takes the sdk object and the temperature, sets the temperature and 
    '''
    ret = sdk.SetTemperature(temp)
    print("Function SetTemperature returned {} Set Temperature: {}".format(ret, temp))

    ret = sdk.CoolerON()
    print("Function CoolerON returned {}".format(ret))

    #### BEGIN Temperature stabilization
    print('Waiting for temperature to stabilize:')
    s = True
    p = 1
    while s:
        sleep(15)
        (ret, temperature) = sdk.GetTemperature()
        print("...Current temperature is {}C".format(temperature))
        if ret == errors.Error_Codes.DRV_TEMP_STABILIZED:
            print('...Temperature stabilized')
            s = False
        if p > 30:
            print('...Instrument did not stabilize')
            if (temperature <= (temp+5)) and (temperature >= (temp-5)):
                print("......Temperature ({}C) is close to set temp, will continue".format(
                    temperature))
                s = False
            else:
                print("......Temperature ({}C) is outside Set Temp +-5C, will exit".format(
                    temperature))
                print("......Will run CoolerOFF for a bit for safety")
                ret = sdk.CoolerOFF()
                print("......Function CoolerOFF returned {}".format(ret))
                sleep(120)
                ret = sdk.ShutDown()
                print("......Function ShutDown returned {}".format(ret))
                sys.exit()
        p += 1
    #### END Temperature stabilization

def prepare_camera(sdk,acqMode,readMode,trigMode,accum_number,accum_cycle,exptime):
    ret = sdk.SetAcquisitionMode(acqMode)
    print("Function SetAcquisitionMode returned {}".format(ret))

    #### BEGIN Set Number and Cycles for ACCUMULATION MODE
    if acqMode == codes.Acquisition_Mode.ACCUMULATE:
        ret = sdk.SetNumberAccumulations(accum_number)
        print("Function SetNumberAccumulations returned {}".format(ret))
    
        ret = sdk.SetAccumulationCycleTime(accum_cycle)
        print("Function SetAccumulationCycleTime returned {}".format(ret))
    #### END Set Number and Cycles for ACCUMULATION MODE

    ret = sdk.SetReadMode(readMode)
    print("Function SetReadMode returned {}".format(ret))

    ret = sdk.SetTriggerMode(trigMode)
    print("Function SetTriggerMode returned {}".format(ret))

    (ret, xpixels, ypixels) = sdk.GetDetector()
    print("Function GetDetector returned {} xpixels = {} ypixels = {}".format(
        ret, xpixels, ypixels))

    ret = sdk.SetExposureTime(exptime)
    print("Function SetExposureTime returned {} time = {}s".format(ret,exptime))

    (ret, fminExposure, fAccumulate, fKinetic) = sdk.GetAcquisitionTimings()
    print("Function GetAcquisitionTimings returned",
            "{} exposure = {} accumulate = {} kinetic = {}".format(ret, fminExposure, 
                                                             fAccumulate, fKinetic))

    ret = sdk.PrepareAcquisition()
    print("Function PrepareAcquisition returned {}".format(ret))
    return xpixels

def shutdown_camera(sdk):
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

