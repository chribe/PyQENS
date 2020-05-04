import numpy as np
import h5py
import datetime
import mantid.simpleapi as MTD
import os

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
        x:?? i don't know that it is (default set to Not a Number)
        y:?? i don't know that it is (default set to Not a Number)
        sqw: array
            The dynamic structure factor (default is Not a Number)
        dsqw: array
            The error of dynamic structure factor (default is Not a Number)
        q: array
            the value of the scattering vector at which the data were collected
        hw: the energy transfer for this measurement

        Methods
        -------
        calculate_MSD to be implemented
        """

    # initialising the class. written so you can also put values in it if you have already loaded the data in another variable by writing:

    # Measurement(sample:'Ig with whatever',type:'fws', temperature:300)
    # with this an instance of the class Measurement is created and has default values for every attribute except  sample,type and temperature

    def __init__(self, namesample='', sample=[], type=[], temperature=[], time=[],
                 end_time=[], time_delta=[], time_initial=[],
                 temperature_initial=[], temperature_final=[], set_temperature=np.NaN, x=None, y=None, sqw=None,
                 dsqw=None, q=None, hw=None, run_num=[]):
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

#    def find_Vanadium(pathraw):
#        files = os.listdir(pathraw)
#        Vanadium = [];
#        for file in files:
#            if '.nxs' in file:
#                f_raw = h5py.File(pathraw + file + '.nxs', 'r')
#                name = [(f_raw['/entry0/subtitle'][0]).decode()]
#                run_num = [(f_raw['/entry0/run_number'][0])]
#                if name == 'vanadium' or name == 'Vanadium':
#                    Vanadium += [run_num]
#
#       return Vanadium

#    def find_EmptyCan(pathraw):
#        files = os.listdir(pathraw)
#        EmptyCan = []
#        for file in files:
#            if '.nxs' in file:
#                f_raw = h5py.File(pathraw + file + '.nxs', 'r')
#                name = [(f_raw['/entry0/subtitle'][0]).decode()]
#                run_num = [(f_raw['/entry0/run_number'][0])]
#                if name == 'EC' or name == 'EmptyCan':
#                    EmptyCan += [run_num]
#        return EmptyCan

    def load_data(self,pathraw, filepath_analysed, Vanadium, EC):

        """Load data in the class Measurement
        pathraw is the file path where the raw data are
        filepath_analysed is the file path where the mantid reduced data are
        sample identifier is a dictionary with:
        - Sample Name
        - Set Temperature
        - run_num is an array of measurements' numbers
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
        #print(name_save_reduced_file)
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

        for run_nm in self.run_num:
            # load some data from the raw data
            f_raw = h5py.File(pathraw + str(run_nm) + '.nxs', 'r')
            print('loading raw data from run '+str(run_nm))
            self.sample += [(f_raw['/entry0/subtitle'][0]).decode()]
            self.temperature += np.squeeze([(f_raw['/entry0/sample/temperature'][0])])
            self.time += [datetime.datetime.strptime(str(f_raw['/entry0/start_time/'][0])[2:-1], "%d-%b-%y %H:%M:%S")]
            self.temperature_final += np.squeeze( [(f_raw['/entry0/sample/setpoint_temperature'][0])])
            self.end_time += [datetime.datetime.strptime(str(f_raw['/entry0/end_time/'][0])[2:-1], "%d-%b-%y %H:%M:%S")]
            # self.run_num += [(f_raw['/entry0/run_number'][0])]
            #maxEnergytrans = np.squeeze([np.array(f_raw['/entry0/instrument/Doppler/maximum_delta_energy/'])])

            #if maxEnergytrans < 20:
            #    self.type += ['fws']
            #else:
            #    self.type += ['qens']
