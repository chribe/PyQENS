import numpy as np
import h5py
import copy
import datetime
import lmfit
from shutil import copyfile
import pickle
import subprocess, os # for summary
from scipy.special import wofz
import matplotlib.pyplot as plt
# from lmfit.models import ExpressionModel
import matplotlib.style
import matplotlib
from Measurement import Measurement
import fws_utilities
import pylab
import os

class Experiment:

    def __init__(self,Data=Data(),Information='',Fit=[]):
        self.Data = Data
        self.Information = Information
        self.Fit = Fit




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

    def __init__(self, sample = '', type = '', temperature = np.NaN, time = datetime.datetime.now(),
                 end_time=  datetime.datetime.now(), time_delta = np.NaN, time_initial = datetime.datetime.now(),
                 temperature_initial = np.NaN, temperature_final = np.NaN, x = None, y = None, sqw = None,
                 dsqw = None, q = None, hw = None, run_num= np.NaN):
     # initialising all the variables
        self.sample = sample
        self.type = type
        self.temperature = temperature
        self.time = time
        self.end_time = end_time
        self.time_initial = time_initial# not in load data
        self.time_delta = time_delta    # not in load data
        self.temperature_initial = temperature_initial # not in load data
        self.temperature_final = temperature_final
        self.x = [] if x is None else x
        self.y = [] if y is None else y
        self.sqw = [] if sqw is None else sqw
        self.dsqw = [] if dsqw is None else dsqw
        self.q = [] if q is None else q
        self.hw = [] if hw is None else hw
        self.run_num = run_num

    @staticmethod
    def load_data_onebyoneIN16b(pathraw,filepath_analysed,sample_identifier):
        """Load data in the class Measurement
        pathraw is the file path where the raw data are
        filepath_analysed is the file path where the mantid reduced data are
        sample identifier is a dictionary with:
        - Sample Name
        - Set Temperature
        - Set Pressure
        ....
        - run_num is an array of measurements' numbers
        """

        # initialise an object, called new_measurement, with the class measurement
        new_measurement = Measurement()

        for run_num in sample_identifier['run_num']:
            # load some data from the raw data
            f_raw = h5py.File(pathraw + str(run_num)+'.nxs', 'r')

            new_measurement.sample            = sample_identifier['SampleName']
            new_measurement.temperature       = (f_raw['/entry0/sample/temperature'][0])
            new_measurement.time              = datetime.datetime.strptime(str(f_raw['/entry0/start_time/'][0])[2:-1],"%d-%b-%y %H:%M:%S")
            new_measurement.temperature_final = sample_identifier['SetTemperature']
            new_measurement.end_time          = datetime.datetime.strptime(str(f_raw['/entry0/end_time/'][0])[2:-1], "%d-%b-%y %H:%M:%S")

            maxEnergytrans = np.array(f_raw['/entry0/instrument/Doppler/maximum_delta_energy/'])

            if maxEnergytrans < 20:
                new_measurement.type = 'fws'
            else:
                new_measurement.type = 'qens'

            # load some data from the mantid analysed
            f = h5py.File(filepath_analysed + str(run_num), 'r')  # open hdf file for reading
            k = 1  # stored Mantid workspace index (one for each energy)
            h5path = '/mantid_workspace_' + str(k) + '/'  # internal hdf node name for that workspace
            while h5path + 'workspace/axis1/' in f.keys():  # loop over all energy offsets
                new_measurement.x.append( np.transpose(np.array(f[h5path + 'workspace/axis1/'])))
                new_measurement.y.append(np.transpose(np.array(f[h5path + 'workspace/axis2/'])))
                new_measurement.sqw.append(np.transpose(np.array(f[h5path + 'workspace/values/'])))
                new_measurement.dsqw.append(np.transpose(np.array(f[h5path + 'workspace/errors/'])))
                new_measurement.q.append(4 * np.pi / 6.271 * np.sin(np.transpose(np.array(f[h5path + 'workspace/axis2/'])) * np.pi / 360.0))
                new_measurement.hw.append(np.array(f[h5path + 'logs/Doppler.maximum_delta_energy/value']))
                k = k + 1
                h5path = '/mantid_workspace_' + str(k) + '/'

        return new_measurement
    @staticmethod
    def load_data_binning(pathraw, filepath_analysed, run_nums):
        new_measurement=Measurement()
# put here the mantid function to
        return new_measurement


class Fit:
    """ this is a class which contains the fits, this class is normally part of the class data,
     so that the corresponding data are associated to it"""
    def __init__(self, FitName='', FitModel = '', Parameters = [], Comments = [],Subfit=[],Background=''):
        # initialising all the variables
        self.FitName = FitName
        self.FitModel = FitModel #string
        self.Parameters = Parameters
        self.Comments = Comments
        self.Subfit=Subfit
    @staticmethod
    def fit_QENS(data,FitName):
        if FitName=='whatever':







