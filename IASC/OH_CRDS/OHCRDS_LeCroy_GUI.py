import tkinter as tk
import sys
import time
from pyvisa import constants,ResourceManager
from time import sleep
from scipy.optimize import least_squares
import numpy as np
import datetime as dt

class App:
    """Define the application class."""
    def __init__(self,root=tk.Tk()):
        ### Initialising GUI
        self.root = root
        self.root.resizable(False, False)
        ### The following grid_columnconfigure statemets make all the columns to belong to the same group, will have the
        ### size of the largest widget
        #self.root.grid_columnconfigure(0, weight=1, uniform="fred")
        #self.root.grid_columnconfigure(1, weight=1, uniform="fred")
        #self.root.grid_columnconfigure(2, weight=1, uniform="fred")
        #self.root.grid_columnconfigure(3, weight=1, uniform="fred")
        #self.root.grid_columnconfigure(4, weight=1, uniform="fred")
        #self.root.grid_columnconfigure(5, weight=1, uniform="fred")
        #self.root.grid_columnconfigure(6, weight=1, uniform="fred")
        self.root.title('OH CRDS LeCroy')


        ########################## PARAMETERS ###################################################

        # Instument parameters
        self.size = 4                        # the size of the HEX chunks to process the waveform
                                        # depends on the precision, 2 for BYTE, 4 for WORD
                                        # check COMM_FORMAT command for LeCroy instruments

        self.scanmatevi = 'ASRL1::INSTR'     # pyvisa name for INSTR at serial COM1
        self.lecroyvi = 'GPIB5::1::INSTR'    # pyvisa name for INSTR at GPIB channel 1
                                        # you can check them by the list.resources() method
                                        # check pyvisa DOCS for further info

        self.cycles = 10000                  # number of meas cycles before switching wavelength
        self.wavemeas = 307.921              # wavelength for tau_sample
        self.waveblank = 308                 # wavelength for tau_0

        # Fitting parameters
        self.len_offset = 300                # Length of baseline offset for logarithm fitting
        self.start_fit = 600
        self.end_fit = 5000
        self.x0 = np.array([0,2,-50000])     # x0 for least_squares()

        # Cycle parameters
        self.sweeps = 30
        self.runtime = 1.1 * (sweeps/10)

        ######################### END OF PARAMETERS #############################################
        ### Initialising dictionary for variables
        self.keyvars = {}

        ### Frame 1: Radiobuttons
        self.frame1 = tk.Frame(self.root)
        #self.radio = Radiobutton(self.frame1,'ID;Set A0;Set A1;Set PWM;Ground Outputs;Read Inputs;Diff. Measurement', 'print(self.mssg)')
        # We add a button to test our setup
        self.test_button = tk.Button(self.frame1, text="Start", command=self.run_mode)
        self.test_button.grid(row=1,column=0)
        self.frame1.grid(row=1,column=0,columnspan=15,sticky='we')

        # Add windows where we are going to write the std output.
        self.console_text = tk.Text(self.root, state='disabled', height=10)
        self.console_text.grid(row=2,column=0,rowspan=20,columnspan=15,sticky='we')

        # We redirect sys.stdout -> TextRedirector
        self.redirect_sysstd()

        #### Resource initialization and opening

        self.rm = ResourceManager()
        #print(rm.list_resources()) ### uncomment to look at the names
        self.lecroy = rm.open_resource(lecroyvi)
        #self.scanmate = rm.open_resource(scanmatevi,
        #                            stop_bits = constants.StopBits.two,
        #                            read_termination = '\r',
        #                            write_termination = '\r')

        #### LeCroy setup
        #lecroy.write('DISP OFF')   # uncomment to stop display
        self.lecroy.write('COMM_HEADER OFF')
        self.lecroy.write('COMM_FORMAT OFF,WORD,HEX')
        self.lecroy.write('TA:DEF EQN,\'AVGS(C1)\',MAXPTS,100000,SWEEPS,%i' %sweeps)
        #lecroy.write('CLM M1')     # clear memory M1

        #### Good to go queries
        print('LECROY:',lecroy.query('ALST?')) # ALL STATUS just to make sure it works
        #print('SCANMATE:',scanmate.query('S?')) # STATUS to make sure it gets R

        #### Waveform info
        self.counts,self.vertgain,self.vertoff,self.horint,self.horoff = getWaveformParams(self.lecroy)
        self.tdivs = np.arange(self.counts)*self.horint+self.horoff

    def redirect_sysstd(self):
        # We specify that sys.stdout point to TextRedirector
        sys.stdout = TextRedirector(self.console_text, "stdout")
        sys.stderr = TextRedirector(self.console_text, "stderr")

    def getFloat(res):
        '''
        Gets a splitted string response to a query and converts what it can to a float
        '''
        for ele in res:
            try:
                number=float(ele)
            except:
                pass
        return number

    def getInt(res):
        '''
        Gets a splitted string response to a query and converts what it can to an integer
        '''
        for ele in res:
            try:
                number=int(ele)
            except:
                pass
        return number

    def getChunks(hexes,size,vertgain,vertoff):
        '''
        Converts the waveform HEX DAT1 string to values by chunking the string,
        converting each chunk into HEX, then into big endian signed integers, and
        then using the equation from the manual
        V = Vertical Gain * INT + Vertical offset
        '''
        chunks = [hexes[i:i+size] for i in range(0,len(hexes),size)]
        #print(chunks)
        vals = []
        for ele in chunks:
            hhex = bytes.fromhex(ele)
            hint = int.from_bytes(hhex,byteorder='big',signed=True)
            #print(hint)
            vals.append(vertgain*hint+vertoff)
        return vals

    def getWaveformParams(self):
        '''
        Gets needed waveform parameters from a series of queries that are converted
        to integer or float depending on the parameters
        '''
        counts = getInt(self.lecroy.query('TA:INSPECT? \"WAVE_ARRAY_COUNT\"').split())
        print('Counts:',counts)

        vertgain = getFloat(self.lecroy.query('TA:INSPECT? \"VERTICAL_GAIN\"').split())
        print('Vertical gain:',vertgain)

        vertoff = getFloat(self.lecroy.query('TA:INSPECT? \"VERTICAL_OFFSET\"').split())
        print('Vertical offset:',vertoff)

        horint = getFloat(self.lecroy.query('TA:INSPECT? \"HORIZ_INTERVAL\"').split())
        print('Horizontal interval:',horint)

        horoff = getFloat(self.lecroy.query('TA:INSPECT? \"HORIZ_OFFSET\"').split())
        print('Horizontal offset:',horoff)

        return counts,vertgain,vertoff,horint,horoff

    def gendata(t,a,b,c):
        '''
        Generates the exponential fit
        '''
        return a+b*np.exp(t*c)

    def funct(x,t,y):
        '''
        Function to fit with optimize.least_squares
        '''
        return x[0] + x[1]*np.exp(x[2]*t)-y

    def run_mode(self):
        #### Begin measurements
        ## Moves to blank wavelength and sets check variable
        #scanmate.write('X=%.3f' %waveblank)
        #isblank = True
        # Initializes lists
        pid = subprocess.Popen(["python","OHCRDS_plotter.py"]).pid
        taus = []
        tstamp = []
        i=0
        waveform = np.arange(self.counts-1)
        while True:
            try:
                with open('endme','r') as f:
                    read=f.read()
                os.remove("endme")
                t1 = dt.datetime.now()                  # End time
                #print("Seconds elapsed: ",(t1-t0).total_seconds())
                # Generates arrays to export to npy or txt
                taumat = np.array(taus)
                tstampmat = np.array(tstamp,dtype='U19')

                # Saves NPY files
                #np.save('taus.npy', taumat) ### The taus as npy
                #np.save('waveform.npy',waveform) ### The last waveform as npy
                #np.save('timestamps.npy',tstampmat) ### The timestamps as npy

                # Saves TXT file: [timestamp tau]
                np.savetxt("M" + t1.strftime('%y%m%d%H%M') + '.txt',np.column_stack((tstampmat,taumat)),fmt='%s')

                ### Some rough statistical data in case you want to know
                print('Average TAU: %.2f' %np.average(taumat))
                print('Stdev: %.2f' %np.std(taumat))
                break

            except:
                pass

            ## Resets the trace, waits for ready
            ta = dt.datetime.now()
            self.lecroy.write('TA:FRST')
            tb = dt.datetime.now()
            sleep(self.runtime - (tb-ta).total_seconds())
            #waitReady(lecroy,0.01,timeout=3,verbose=0) ### Deprecated, problems with WaitReady()

            ## Main loop, stores to M1, resets trace, queries data, gets values,
            ## fits paraments, gets residuals, waits for ready
            ## sends wavelength change to scanmate after n cycles
            ## is presented as a function for the animation class

            print('RUN:',i+1)
            t1 = dt.datetime.now()
            self.lecroy.write('STO TA,M1')
            self.lecroy.write('TA:FRST')
            #lecroy.write('CLM M1')
            data = self.lecroy.query('M1:WF? DAT1').split()
            t2 = dt.datetime.now()
            print((t2-t1).total_seconds())
            values = getChunks(data[-1],self.size,self.vertgain,self.vertoff)
            waveform = np.array(values[:-1])

            # Prepare data for fitting
            ts=self.tdivs[self.start_fit:self.end_fit]
            ys=waveform[self.start_fit:self.end_fit]

            #### LEAST SQUARES

            # Levenberg - Marquardt
            #res_lsq = least_squares(fun,x0, method = 'lm', ftol=1e-12,xtol=1e-12,gtol=1e-12,
            #                       args=(ts,ys))

            # TRF
            res_log = least_squares(funct,self.x0, ftol=1e-12,xtol=1e-12,gtol=1e-12,
                                    loss = 'cauchy', f_scale=0.1, args=(ts,ys))

            # The results
            xs = res_log.x
            #### END OF LEAST SQUARES

            #### NATURAL LOG
            # Fit the baseline and get a function
            #coef = np.polyfit(tdivs[:len_offset+1],values[:len_offset+1],1)
            #poly1d_fn = np.poly1d(coef)
            # Remove the baseline
            #offset_values = values[:-1]-poly1d_fn(tdivs[:-1])
            # Take the log
            #logs = np.log(offset_values[start_fit:end_fit])
            # Fit the log
            #coef2 = np.polyfit(tdivs[start_fit:end_fit],logs,1)
            # The results
            #xs = (coef[1],np.exp(coef2[1]),coef2[0])
            #### END OF NATURAL LOG

            #print(xs)              # in case you want to see the fit parameters

            # Shows and saves tau to a list, timestamps to list as well
            tau = -1e6/xs[2]
            taus.append(tau)
            timenow = dt.datetime.now()
            stamp = timenow.strftime('%Y/%m/%d-%H:%M:%S')
            tstamp.append(stamp)
            print('TAU is: %.2f' %tau)

            # Builds the fit = non fitted waveform + generated data from params
            fita = waveform[:self.start_fit]
            fitb = gendata(self.tdivs[self.start_fit:-1],*xs)
            fit = np.concatenate((fita,fitb))

            np.save("ax1data", np.concatenate((self.tdivs,waveform),axis=1))
            np.save("ax1adata", np.concatenate((self.tdivs,fit),axis=1))
            np.save("ax2data", np.concatenate((self.tdivs,waveform-fit),axis=1))


            # wavelength change
            #if (i+1)% cycles == 0:
            #    if isblank:
            #        scanmate.write('X=%.3f' %wavemeas)
            #        isblank = False
            #    else:
            #        scanmate.write('X=%.3f' %waveblank)
            #        isblank = True
            # end

            # Total computing time in cycle, rest will be sleeping
            t3 = dt.datetime.now()
            print((t3-t1).total_seconds())
            sleep(runtime-(t3-t1).total_seconds())
            #waitReady(lecroy,0.01) ### Deprecated due to problems with WaitReady()

            # Total cycle time
            t4 = dt.datetime.now()
            print((t4-t1).total_seconds())


class LabelEntry:
    def __init__(self,place,keyvars={'Entry':'0'}):
        sef.val = tk.StringVar()
        for k,v in keyvars.times():
            r = tk.Radiobutton(place, text=item, variable=self.val,
                                value=i, command=self.cb, **kwargs)
            r.grid(row=1,column=i)#,sticky='w')
            #r.pack(side=tk.LEFT)

class Radiobutton:
    """Create a list-based Radiobutton object."""

    def __init__(self, place, items='Radiobutton', cmd='', val=0, **kwargs):
        self.items = items.split(';')
        self.cmd = cmd
        self.val = tk.IntVar()
        self.val.set(val)
        self.item = self.items[self.val.get()]
        for i, item in enumerate(self.items):
            r = tk.Radiobutton(place, text=item, variable=self.val,
                                value=i, command=self.cb, **kwargs)
            r.grid(row=1,column=i)#,sticky='w')
            #r.pack(side=tk.LEFT)
    def cb(self):
        """Evaluate the cmd string in the Radiobutton context."""
        self.item = self.items[self.val.get()]
        self.mssg = 'Mode '+self.item+' selected'
        exec(self.cmd)

class TextRedirector(object):
    def __init__(self, widget, tag):
        self.widget = widget
        self.tag = tag

    def write(self, text):
        self.widget.configure(state='normal') # Edit mode
        self.widget.insert(tk.END, text, (self.tag,)) # insert new text at the end of the widget
        self.widget.configure(state='disabled') # Static mode
        self.widget.see(tk.END) # Scroll down
        self.widget.update_idletasks() # Update the console

    def flush(self):
        pass

if __name__ == "__main__":
    app = App()
    app.root.mainloop()
