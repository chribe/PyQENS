import numpy as np
import h5py
import datetime
import matplotlib.pyplot as plt
import mantid.simpleapi as MTD
import os
from scipy.optimize import curve_fit
class Data:
    """this is a container for the data collected in multiple or single run and information regarding that measurement
        Attributes
        ----------
        sample: str
            The name given to the sample during the data acquisition (default is empty string)
        type: str
            If the measurement is a qens spectra or a fws (default is empty string)
        run_num: array
            the run numbers of the files (default is empty string)
        temperature: float
            this is the temperature recorded during the acquisition (default is Not a Number)
        time: date value
           The time at which the measurement started (default the current date and time, meaning the date and time at which the initialisation happens )
        end_time: date value
           The time at which the measurement ended (default the current date and time, meaning the date and time at which the initialisation happens )
        time_delta: float
            this is the time in second since the time at which the set temperature was changed (default is Not A Number)
        time_initial: date value
            The time at which the set temperature was changed (default the current date and time, meaning the date and time at which the initialisation happens )
        temperature_initial: float
            The temperature at which the set temperature was changed (default is Not a Number)
        temperature_final: float
            The temperature of the set temperature(default is Not a Number)
        x:?? depending on what analysis is performed the value of x and y can change, it is safer to not use it
        y:??depending on what analysis is performed the value of x and y can change, it is safer to not use it
        sqw: array
            The dynamic structure factor (default is Not a Number)
        dsqw: array
            The error of dynamic structure factor (default is Not a Number)
        q: array
            the value of the scattering vector at which the data were collected
        hw: the energy transfer for this measurement

        Fit:list of fits

        Methods
        -------
        -load data
        -calculate_msd
        -clean_nonsense
        -get_values_qens
        """

    # initialising the class. written so you can also put values in it if you have already loaded the data in another variable by writing:

    # Measurement(sample:'Ig with whatever',type:'fws', temperature:300)
    # with this an instance of the class Measurement is created and has default values for every attribute except  sample,type and temperature

    def __init__(self,
                 namesample='',
                 sample=[],
                 type=[],
                 temperature=[],
                 time=[],
                 end_time=[],
                 time_delta=[],
                 time_initial=[],
                 temperature_initial=[],
                 temperature_final=[],
                 set_temperature=np.NaN,
                 x=None,
                 y=None,
                 sqw=None,
                 dsqw=None,
                 q=None,
                 hw=None,
                 run_num=[],
                 MSD=[],
                 dMSD=[]):
        # initialising all the variables
        self.namesample = namesample
        self.sample = sample
        self.type = type
        self.temperature = temperature
        self.set_temperature = set_temperature
        self.time = time
        self.end_time = end_time
        self.time_initial = time_initial  # not in load data
        self.time_delta = time_delta  # not in load data
        self.temperature_initial = temperature_initial  # not in load data
        self.temperature_final = temperature_final
        self.x = [] if x is None else x
        self.y = [] if y is None else y
        self.sqw = [] if sqw is None else sqw
        self.dsqw = [] if dsqw is None else dsqw
        self.q = [] if q is None else q
        self.hw = [] if hw is None else hw
        self.run_num = run_num
        self.MSD=MSD
        self.dMSD = dMSD

    def find_Vanadium(pathraw):
        'find the runnumbers of samples called vanadium or Vanadium in a path pathraw'
        files = os.listdir(pathraw)
        Vanadium = [];
        for file in files:
            if '.nxs' in file:
                f_raw = h5py.File(pathraw + file + '.nxs', 'r')
                name = [(f_raw['/entry0/subtitle'][0]).decode()]
                run_num = [(f_raw['/entry0/run_number'][0])]
                if name == 'vanadium' or name == 'Vanadium':
                    Vanadium += [run_num]

        return Vanadium

    def find_EmptyCan(pathraw):
        'find the runnumbers of samples called EmptyCan or EC in a path pathraw'
        files = os.listdir(pathraw)
        EmptyCan = []
        for file in files:
           if '.nxs' in file:
                f_raw = h5py.File(pathraw + file + '.nxs', 'r')
                name = [(f_raw['/entry0/subtitle'][0]).decode()]
                run_num = [(f_raw['/entry0/run_number'][0])]
                if name == 'EC' or name == 'EmptyCan':
                    EmptyCan += [run_num]
        return EmptyCan

    def load_data(self, pathraw, filepath_analysed, Vanadium, EC):
        """Load data in the class Measurement
        pathraw is the file path where the raw data are
        filepath_analysed is the file path where the mantid reduced data are

        an instance of the class Data has to be already present, the minimal information to pass to make the code to work is
        to provide:
        -the run number (or run numbers) as a list of int (in self.run_num)
        -the type of measurement ('qens' or 'fws') in self.type ... this will be changed so that you don't have to say the type

        - you can give a name to the selected runs, in samplename.(the name that was written during the measurement is
         in anycase saved under sample ) it is not required
         all the other informations will be stored automatically
                """
        if os.path.exists(filepath_analysed) == False:
            os.makedirs(filepath_analysed)
            print("Attention! The folder was not existing. Created new folder" + filepath_analysed)

        file_list = os.listdir(filepath_analysed)
        # -----------------------------------------------------------------------------------------------------------------------
        # this is the definition of the name of the saved reduced file:
        # it is one number if the file is only one and it is first run number_last run number if there are more files
        name_save_reduced_file = str(self.run_num[0])
        if len(self.run_num) > 1:
            name_save_reduced_file = str(self.run_num[0]) + '_' + str(
                self.run_num[len(self.run_num) - 1])
        #runs is a string containing all the files, which is an input of the function to reduce the data   MTD.IndirectILLReductionFWS
        runs = ''
        for r in self.run_num:
            runs = runs+pathraw+str(r)+'.nxs,'
        runs = runs[:-1]
        import mantid.simpleapi as MTD
        # -----------------------------------------------------------------------------------------------------------------------
        # if the file was not reduced yet it will be reduced
        if name_save_reduced_file not in file_list:
            if 'fws' in self.type:
                print('Reducing the FWS runs')
                ws = MTD.IndirectILLReductionFWS(Run=runs, MapFile=filepath_analysed+'IN16B_grouping.xml',
                                                 Observable='start_time')  # BackgroundRun='265627:265638'
                ws = MTD.ConvertSpectrumAxis(ws, Target='theta', Version=1)
                MTD.SaveNexus(ws, name_save_reduced_file)

            if 'qens' in self.type:
                print('Reducing the QENS runs')
                ws = MTD.IndirectILLReductionQENS(Run=runs, AlignmentRun=Vanadium, SumRuns=True,
                                              UnmirrorOption=7, Version=1,
                                              MapFile=filepath_analysed+'IN16B_grouping.xml',
                                              BackgroundRun=EC)
                ws = MTD.ConvertSpectrumAxis(ws, Target="theta", Version=1)
                MTD.SaveNexus(ws, filepath_analysed+name_save_reduced_file)
    # -----------------------------------------------------------------------------------------------------------------------
        # load some data from the mantid analysed
        f = h5py.File(filepath_analysed + name_save_reduced_file, 'r')  # open hdf file for reading
        k = 1  # stored Mantid workspace index (one for each energy)
        h5path = '/mantid_workspace_' + str(k) + '/'  # internal hdf node name for that workspace
        while h5path + 'workspace/axis1/' in f.keys():  # loop over all energy offsets
            a=np.squeeze(np.transpose(np.array(f[h5path + 'workspace/axis1/'])))
            self.x    += [a]
            self.y    += [np.squeeze(np.transpose(np.array(f[h5path + 'workspace/axis2/'])))]
            self.sqw  += [np.squeeze(np.transpose(np.array(f[h5path + 'workspace/values/'])))]
            self.dsqw += [np.squeeze(np.transpose(np.array(f[h5path + 'workspace/errors/'])))]
            self.q    += [np.squeeze(4 * np.pi / 6.271 * np.sin(np.transpose(np.array(f[h5path + 'workspace/axis2/'])) * np.pi / 360.0))]
            self.hw   += [a[:-1] * 1e3]
            k = k + 1
            h5path = '/mantid_workspace_' + str(k) + '/'
        self.x = np.squeeze(self.x)
        self.y = np.squeeze(self.y)
        self.sqw = np.squeeze(self.sqw)
        self.dsqw = np.squeeze(self.dsqw)
        self.q = np.squeeze(self.q)
        self.hw = np.squeeze(self.hw)

        print('data raw loading')
        self.sample = []
        self.temperature = []
        self.temperature_final = []
        self.time = []
        self.time_delta = []
        self.end_time = []
        for run_nm in self.run_num:
            # load some data from the raw data
            f_raw = []

            f_raw = h5py.File(pathraw + str(run_nm) + '.nxs', 'r')
            print('loading raw data from run '+str(run_nm))
            self.sample += [(f_raw['/entry0/subtitle'][0]).decode()]
            self.temperature += [np.squeeze([(f_raw['/entry0/sample/temperature'])])]
            self.time += [datetime.datetime.strptime(str(f_raw['/entry0/start_time/'][0])[2:-1], "%d-%b-%y %H:%M:%S")]
            self.temperature_final += [np.squeeze( [(f_raw['/entry0/sample/setpoint_temperature'])])]
            self.end_time += [datetime.datetime.strptime(str(f_raw['/entry0/end_time/'][0])[2:-1], "%d-%b-%y %H:%M:%S")]
            # self.run_num += [(f_raw['/entry0/run_number'][0])]
            #maxEnergytrans = np.squeeze([np.array(f_raw['/entry0/instrument/Doppler/maximum_delta_energy/'])])

            #if maxEnergytrans < 20:
            #    self.type += ['fws']
            #else:
            #    self.type += ['qens']

    def calculate_MSD(self):
        if 'qens' in self.type:
            def line(m,q,x):
                return q-m*x

            self.MSD=[]
            self.dMSD = []

            plt.figure()
            popt, pcov = curve_fit(line, self.q[2:] ** 2, np.squeeze(self.sqw[self.hw ==0 ,2: ]),
                                 bounds=(0, [.1, 1]), ftol=1e-14, xtol=1e-14 )

            plt.errorbar( self.q[2:]**2, y=np.squeeze(self.sqw[[int(np.round(len(self.hw)/2)-1),int(np.round(len(self.hw)/2)),int(np.round(len(self.hw)/2)+1)],2:].mean(axis=0)),
                          yerr=np.squeeze(self.dsqw[self.hw ==0,2:]), marker='o', linestyle='--')
            plt.plot(self.q**2, popt[0] - popt[1] * self.q**2)
            plt.xlabel(r'q [\AA^{2}]')
            plt.ylabel(r'S(q,hw=0)')
            self.MSD += [popt[1] * 3]
            self.dMSD += [np.sqrt(np.diag(pcov))[1] * 3]
            print('the MSD is'+str([popt[1] * 3])+'\pm'+str([np.sqrt(np.diag(pcov))[1] * 3]) )

        if self.type=='fws':
            for i in range(len(self.hw)):
                if self.hw[i]!=0:
                    print('The fws scans used are not elastic. Choose an elastic scan to calculate the MSD')

                elif self.hw[i] == 0:

                    popt, pcov = curve_fit(line, self.q**2, ydata=np.squeeze(self.sqw[:, i]),sigma=np.squeeze(self.dsqw[:, i] **(-2)),
                                           bounds=(0, [.1, 1]), ftol=1e-14, xtol=1e-14 )
                    self.MSD +=[-popt[1]*3]
                    self.dMSD += [np.sqrt(np.diag(pcov))[1]*3]
                    #plt.errorbar( self.q**2, ydata=self.sqw[:, i],sigma=self.dsqw[:, i] **(-2), marker='o', linestyle='--')
                    #plt.plot(self.q**2, popt[0] + popt[1] * self.q**2)
                    #plt.xlabel(r'q [\AA^{2}]')
                    #plt.ylabel(r'S(q,hw=0)')

    def clean_nonsense(self):
        self.sqw[np.isnan(self.sqw)] = 0
        self.sqw[np.isinf(self.sqw)] = 0
        self.sqw[self.sqw < 5e-7 * self.sqw.max()] = 0  # reasonable signal-to-noise range
        self.dsqw[np.isnan(self.sqw)] = np.inf
        self.dsqw[np.isinf(self.sqw)] = np.inf
        self.dsqw[np.isnan(self.dsqw)] = np.inf
        self.dsqw[self.sqw == 0] = np.inf  # no signal at all should have infinite error
        self.dsqw[self.dsqw == 0] = np.inf  # "a measurement without error is nonsense"

    def get_values_qens(self):
        'returns the values of q,hw,sqw,dsqw,sqw_copy,dsqw_copy the last two are simply a copy of sqw,dsqw'
        q = self.q.copy()
        hw   = self.hw.copy()
        sqw  = self.sqw.copy()
        dsqw = self.dsqw.copy()
        sqw_copy = self.sqw.copy()
        dsqw_copy= self.dsqw.copy()
        return q,hw,sqw,dsqw,sqw_copy,dsqw_copy

    def subtract_sqw(self,data):
        'function to subtract the S(q,hw) of a measurement data from a dataset (self)'
        self.sqw = self.sqw-data.sqw
        self.dsqw =(self.dsqw**2+data.dsqw**2)**(1/2)

    def subtract_sqw_with_factor(self,data,factor):
        """function to subtract the S(q,hw) of a measurement data(instance of the class Data) from a dataset (self)

        factor is an array with 18 float or int (one each q) for a q dependent weight in the subtraction

        """
        for i in range(len(self.q)):
            self.sqw[:,i] = self.sqw[:,i] - factor[i]*data.sqw[:,i]
            self.dsqw[:,i] = (self.dsqw[:,i] ** 2 + (factor[i]*data.dsqw[:,i]) ** 2) ** (1 / 2)

