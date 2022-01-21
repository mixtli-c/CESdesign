import numpy as np
import datetime as dt
import time
import matplotlib.pyplot as plt
import scipy.signal as scs
from matplotlib.dates import DateFormatter
import configurations as conf
import CESfunctions as cf
import avaspec as av

#########################################################################################
### USER CONFIGURATION (from configuration file)                                      ###

### Instrument
accums = conf.accums                    # Accumulations
checkrate = conf.checkrate              # Query freq. to check if ready (calls per itime)
samples = conf.samples                  # Number of accumulated measurements
itime = conf.itime                      # Integration time in ms
averages = conf.averages                # Number of integrations to average 

### Signal analysis
# Cavity parameters
distance = conf.distance                # Sample optical length

# Resonance window 
lower_wavelength=conf.lower_wavelength  # Starting wavelength of resonance window
upper_wavelength=conf.upper_wavelength  # Ending wavelength of resonance window

# Reference and background files
back_filename = conf.back_filename
no2_refname = conf.no2_refname
chocho_refname = conf.chocho_refname

# Reff : Either number or np.load(Reff_matrix)
Reff= conf.Reff

#########################################################################################
### Reference and Background file loading for analisis                                ###
no2reference = np.load(no2_refname)
#chochoref = np.load(chocho_refname)
background = np.load(back_filename)


#########################################################################################
### Instrument initialization                                                         ###

av.AVS_Init(0) 

ret = av.AVS_GetNrOfDevices()
req = 75*ret
a_pList = av.AvsIdentityType * ret 
retn=av.AVS_GetList(req,req,a_pList)
handle = av.AVS_Activate(retn[1])

config = av.DeviceConfigType
ret = av.AVS_GetParameter(handle, 63484, 63484, config)
pixels = ret[1].m_Detector_m_NrPixels
wavelengths = []
ret = av.AVS_GetLambda(handle,wavelengths)
for pixel in range(pixels):
    wavelengths.append(ret[pixel])
    
#########################################################################################
### Configurations                                                                    ###

measconfig = av.MeasConfigType 
measconfig.m_StartPixel = 0 
measconfig.m_StopPixel = pixels - 1 
measconfig.m_IntegrationTime = itime 
measconfig.m_IntegrationDelay = 0 
measconfig.m_NrAverages = averages 
measconfig.m_CorDynDark_m_Enable = 0 
measconfig.m_CorDynDark_m_ForgetPercentage = 100 
measconfig.m_Smoothing_m_SmoothPix = 0 
measconfig.m_Smoothing_m_SmoothModel = 0 
measconfig.m_SaturationDetection = 0 
measconfig.m_Trigger_m_Mode = 0 
measconfig.m_Trigger_m_Source = 0 
measconfig.m_Trigger_m_SourceType = 0 
measconfig.m_Control_m_StrobeControl = 0 
measconfig.m_Control_m_LaserDelay = 0 
measconfig.m_Control_m_LaserWidth = 0 
measconfig.m_Control_m_LaserWaveLength = 0.0 
measconfig.m_Control_m_StoreToRam = 0

#########################################################################################
### Instrument Preparation for Measurements                                           ###

av.AVS_PrepareMeasure(handle, measconfig)
scans = 1                       # Number of scans (only last is recorded)

#########################################################################################
### LOOPS for accumulations and sample logging                                        ###

plt.ion()                       # Interactive plot
fig = fig=plt.figure()          # Figure initialization

t0=dt.datetime.now()            # Start time

# Initializing timestamp and concentration list, and measurement array
measurements = np.array(wavelengths).reshape(len(wavelengths),1)
timestamp = []
ppbs = []

for n in range(samples):        # Sample loop
    print("Measurement number: ", n+1)
    counts = np.zeros(len(wavelengths))

    for i in range(accums):     # Accumulations loop
        #print(i+1)
        av.AVS_Measure(handle, 0, scans)

        ### Query instrument ready
        dataready = False
        while not dataready:
            time.sleep(measconfig.m_IntegrationTime/(checkrate*1000))
            dataready = av.AVS_PollScan(handle)

        ### Get data and accumulate
        timestamp = 0
        spectraldata = []
        ret = av.AVS_GetScopeData(handle, timestamp, spectraldata)
        timestamp = ret[0]
        spectra = np.array(ret[1][0:len(wavelengths)])
        counts = np.add(counts,spectra) 
    
    ### Plotting
    fig.clear()
    plt.plot(wavelengths,counts)
    plt.show()
    plt.pause(0.001)

    ### Measurements array
    measurements = np.concatenate((measurements,counts.reshape(len(counts),1)),axis=1)
    
    ### Calculating number density
    sample = np.copy(counts.reshape(len(counts),1))
    minwave,maxwave = cf.segment_indices(measurements[:,0:2],lower_wavelength,upper_wavelength)
    bckg = np.copy(background[minwave:maxwave,:])
    no2ref = np.copy(no2reference[minwave:maxwave,:])
    I_sample = np.copy(sample[minwave:maxwave,:])
    I_0 = np.average(bckg[:,1:],axis=1).reshape(len(bckg),1)
    #print(sample.shape,I_sample.shape,I_0.shape,no2ref.shape)

    # This one does everything (see recursive_fit function for CESfunctions.py)
    ndensity = cf.recursive_fit(I_sample, I_0, Reff, distance, no2ref)
    ppb = ndensity/2.5e10
    print('NO2 ppb: ', ppb)
    
    # Add timestamp and ppb to lists
    ppbs.append(ppb)
    timestamp.append(dt.datetime.now)
    

t1=dt.datetime.now()            # End time

# we generate a name to save the measurements
sample_archive = "Isample" + t1.strftime("%y%m%d%H%M")

np.save(sample_archive, measurements)   # for archiving (further analysis)

# we generate a name to save the timestamped concentrations
conc_filename = "measurements" + t1.strftime("%y%m%d%H%M")

conc = np.concatenate((np.asarray(timestamp).reshape(len(timestamp),1),
    np.asarray(ppbs).reshape(len(ppbs),1)),axis=1)
np.save(conc_filename, conc)


print("Seconds elapsed: ",(t1-t0).total_seconds())
print("Shape of measurements array: ",measurements.shape)
print("Shape of concentrations array: ",conc.shape)
print("EOF")
