#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 21:15:23 2020

@author: beck
"""
import h5py
import datetime as date
import numpy as np
import sys
import os as os
import copy
from shutil import copyfile
from tqdm import tqdm
import Data as DT


class Sample:
    '''
    This class contains all the data of one sample.
    It contains the following information:
        Name
        Sample composition
        volume fraction of sample
        total run files
        start times of the runs
        end times of the runs
        duration of the runs
        type of runs
        runnumbers
        binned runnumbers
        list of class Data for all Data based on binned runnumbers
    
    function contained in this class:
        CreateSample
        binning
        
        
        TODO: volumefraction
        TODO: Samplecomposition
    '''
    #%% initialize the class
    def __init__(self,
                 Name='',
                 volumefraction=[],
                 Data=[]):
        # initialising all the variables
        self.Name = Name
        self.volumefraction=[]
        self.Samplecomposition=''
        self.runnumbers=[]
        self.type = []
        self.instrument=''
        self.bin=''
        self.filename=[]
        self.Data=[]
        self.ProtEntry={}
    #%% define functions 
    def CreateSample(self,ProtEntry):
        '''
        This function creates the Sample based on the Protocol Entry

        Parameters
        ----------
        ProtEntry : Dictionary

        Returns
        -------
        None.
        
        Uses
        ----
        binning_loading(self)

        '''
        self.Name=ProtEntry['Sample']
        self.volumefraction=np.NaN
        self.SampleComposition=''
        self.runnumbers=ProtEntry['Runfiles']
        self.starttime=[]
        self.endtime=[]
        self.duration=[]
        self.runnum_org=[]
        self.bin=ProtEntry['binning']
        self.instrument=ProtEntry['instrument']
        self.ProtEntry=ProtEntry
        self.binning_loading()
        return self
    def binning_loading(self):
        '''
        This function creates the list of samples based on the binning type 
        specified. 
        In case of in16b, this function is using the mantid routines
        for the desired binning. They are only executed, if the corresponding 
        sample name does not yet exist.
        For all other instruments, import from ascii is provided as single file 
        input

        Returns
        -------
        None.

        '''
        #% create samplenames
        if self.instrument!='in16b':
            self.bin='floatfile1'
            Warning('binning is only supported for IN16b data up to now!')# TODO
        # define which files will be binned together
        if 'floatfile' in self.bin:
            numberoffiles=int(self.bin[9:])
            for hi in range(1,len(self.runnumbers)-numberoffiles):
                self.runnum_org.append(self.runnumbers[hi:hi+numberoffiles])
        if self.bin=='all':
            self.runnum_org.append(self.runnumbers)
        if 'floathour' in self.bin:
            floathours=date.timedelta(hours=float(self.bin[9:]))
            with tqdm(total=len(self.runnumbers), file=sys.stdout,desc='get starting times') as pbar:
                for hi,num in enumerate(self.runnumbers):
                    starttime,endtime,duration=self.gettimes(num)
                    self.starttime.append(starttime)
                    self.endtime.append(endtime)
                    self.duration.append(duration)
                    pbar.update(1)
            with tqdm(total=len(self.runnumbers), file=sys.stdout,desc='binning run runnumbers') as pbar:
                for startnum,_ in enumerate(self.runnumbers):
                    for endnum in range(startnum,len(self.runnumbers)):
                        duration=self.endtime[endnum]-self.starttime[startnum]
                        if duration>floathours:
                            self.runnum_org.append(self.runnumbers[startnum:endnum])
                            break
                        else:
                            continue
                    pbar.update(1)
                    if floathours>self.endtime[-1]-self.starttime[startnum]:
                        self.runnum_org.append(self.runnumbers[startnum:])
                        pbar.update(len(self.runnumbers)-startnum-1)
                        break
        # load data for different instruments
        with tqdm(total=len(self.runnum_org), file=sys.stdout,desc='load data') as pbar:
            for _,runs in enumerate(self.runnum_org):
                ProtEntry=self.ProtEntry
                self.Data.append(DT.Data(ProtEntry=ProtEntry))
                if len(runs)==1:
                    filename=str(runs)
                else:
                    filename=str(runs[0]) + '_' + str(runs[-1])
                if self.instrument=='in16b':
                    # reduce via Mantid, if the reduced files do not exist
                    if os.path.exists(self.ProtEntry['filepath_analysed']) == False:
                        os.makedirs(self.ProtEntry['filepath_analysed'])
                        print("Attention! The folder was not existing. Created new folder" + self.ProtEntry['filepath_analysed'])
                    file_list = os.listdir(self.ProtEntry['filepath_analysed'])
                    if filename + 'QENS' not in file_list and filename + 'FWS' not in file_list:
                        self.Mantidreduction(runs,filename)
                    self.Data[-1].load_mantid(filename)
                else:
                    self.Data[-1].load_ascii(filename)
                pbar.update(1)
    #%% sub routines
    def gettimes(self,file):
        '''
        This function returns the start and end times read out from the nxs file 

        Parameters
        ----------
        file : runnumber of the file. 

        Returns
        -------
        starttime : datetime
        endtime : datetime
        duration : datetime

        '''
        copyfile( self.ProtEntry['filepath_raw'] + str(file) + '.nxs', "currentfile.nxs")
        f=h5py.File('currentfile.nxs', 'r')
        starttime_str=str(f['/entry0/start_time'][()])
        endtime_str=str(f['/entry0/end_time'][()])
        os.remove("currentfile.nxs")
        val0=starttime_str.find('0')
        val1=starttime_str.find('1')
        val2=starttime_str.find('2')
        startparam=[val0,val1,val2]
        startparam2=[x for x in startparam if x >= 1 ];
        index=min(startparam2)
        starttime_red=starttime_str[index:index+18]
        vale0=endtime_str.find('0')
        vale1=endtime_str.find('1')
        vale2=endtime_str.find('2')
        endparam=[vale0,vale1,vale2]
        endparam2=[x for x in endparam if x >= 1 ];
        indexe=min(endparam2)
        endtime_red=endtime_str[indexe:indexe+18]
        starttime=date.datetime.strptime(starttime_red,"%d-%b-%y %H:%M:%S")
        endtime=date.datetime.strptime(endtime_red,"%d-%b-%y %H:%M:%S")
        duration=endtime-starttime
        return starttime,endtime,duration
    
    def Mantidreduction(self,runnumbers,filename):
        '''
        This function is reducing the QENS data based on mantid. 
        Mantid has to be installed!

        Parameters
        ----------
        runnumbers : list of runnumbers. It can contain FWS and QENS
        filename : string of the file name

        Raises
        ------
        NameError
            if Mantid is not installed in the specified path

        Returns
        -------
        None.

        '''
        Vanadiumruns=self.ProtEntry['VanadiumRuns']
        if os.path.exists(flags['mantidpath'])==False:
            raise NameError('Is Mantid installed on your PC?! If so adapt path!')
        sys.path.append( flags['mantidpath'] + 'lib/')
        sys.path.append( flags['mantidpath'] + 'bin/')
        # then import Mantid function namespace:
        import mantid.simpleapi as MTD
        # add data directory to search path:
        MTD.config.appendDataSearchDir( self.ProtEntry['filepath_raw'] )
        #separate the FWS from the QENS
        FWS_runs=''
        QENS_runs=''
        for _,runs in enumerate(runnumbers):
            copyfile( self.ProtEntry['filepath_raw'] + str(runs) + '.nxs', "currentfile.nxs")
            f=h5py.File('currentfile.nxs', 'r')
            maxEnergytrans=np.array(f['/entry0/instrument/Doppler/maximum_delta_energy/'])
            #TODO write decision if FWS or QENS
            if maxEnergytrans>20:
               QENS_runs+=',' + str(runs)
            if maxEnergytrans<20:
                FWS_runs+=',' + str(runs)
# =============================================================================
            os.remove("currentfile.nxs")
        # Run reduction for FWS and QENS
        if len(FWS_runs)!=0:
            ws = MTD.IndirectILLReductionFWS(Run=FWS_runs[1:], MapFile=self.ProtEntry['groupingfile'],Observable='start_time')
            ws = MTD.ConvertSpectrumAxis(ws, Target='theta', Version=1)
            MTD.SaveNexus(ws,self.ProtEntry['filepath_analysed'] + filename + 'FWS')
            del ws
        if len(QENS_runs)!=0:
            ws = MTD.IndirectILLReductionQENS(Run=QENS_runs[1:], AlignmentRun=str(Vanadiumruns[0]), SumRuns=True,
                                          UnmirrorOption=7, Version=1,
                                          MapFile=self.ProtEntry['groupingfile'])
            ws = MTD.ConvertSpectrumAxis(ws, Target="theta", Version=1)
            MTD.SaveNexus(ws,self.ProtEntry['filepath_analysed'] + filename + 'QENS')
            
            
    def Fit(self,fitname,Vanadium,hwlimits=[],qlimits=[],hwindex=[],qindex=[]):
        '''
        This function is limiting the energy transfers and momentum transfers 
        and handing the arguments to the fit function
        '''
        hwlimits=np.array(hwlimits)
        qlimits=np.array(qlimits)
        with tqdm(total=len(self.Data), file=sys.stdout,desc='Fit Sample',position=0) as pbar:
            for _,Dat in enumerate(self.Data):
                # check for consistency
                if (not len(hwlimits)==0) and (not len(hwindex)==0):
                    Warning('Entries for "hwlimits" and "hwindex" found! Taking the Values of "hwindex"!')
                    hwindex=[]
                if (not len(qlimits)==0) and (not len(qindex)==0):
                    Warning('Entries for "qlimits" and "qindex" found! Taking the Values of "qindex"!')
                    qindex=[]
                for objecttotest in [hwindex,qindex,hwlimits,qlimits]:
                    if len(objecttotest)==1:
                        Warning('only one boundary in ' + str(objecttotest) + '! considering it as lower limit!')
                        objecttotest[1]=np.inf
                    if len(objecttotest)>2:
                        Warning('More then 2 Entries in ' + str(objecttotest) + '! Ignoring additional entries!')
                # copy data and reduce its size
                tempdat=copy.deepcopy(Dat)
                if len(hwlimits)>0 and len(hwindex)==0:
                    # convert limits to indices
                    if len(hwlimits)>0:
                        if tempdat.hw[0]>tempdat.hw[-1]:
                            raise Warning('Mirrowing data to have increasing hw!')
                            tempdat.sqw=tempdat.sqw[-1:-1:0, :]
                            tempdat.dsqw=tempdat.dsqw[-1:-1:0, :]
                            tempdat.hw=tempdat.hw[-1:-1:0]
                        if tempdat.hw[0]>hwlimits.min():
                            hwindex.append(0)
                        if tempdat.hw[-1]<hwlimits.min():
                            raise NameError('lower limit is higher than highest energy transfer!')
                        else:
                            for hiw1, _ in enumerate(tempdat.hw[0:-2]):
                                if (tempdat.hw[hiw1]<hwlimits.min()) and (tempdat.hw[hiw1+1]>hwlimits.min()):
                                    hwindex.append(hiw1)
                                    break
                        if tempdat.hw[-1]<hwlimits.max():
                            hwindex.append(len(tempdat.hw))
                        if tempdat.hw[0] > hwlimits.max():
                            raise NameError('upper limit is lower than lowest energy transfer!')
                        else:
                            for hiw2, _ in enumerate(tempdat.hw[0:-2]):
                                if (tempdat.hw[hiw2]<hwlimits.max()) and (tempdat.hw[hiw2+1]>hwlimits.max()):
                                    hwindex.append(hiw2+1)
                                    break
                if len(qlimits) > 0 and len(qindex) == 0:
                    if len(qlimits) > 0:
                        if tempdat.q[0] > tempdat.q[-1]:
                            raise Warning('Mirrowing data to have increasing q!')
                            tempdat.sqw = tempdat.sqw[:, -1:-1:0]
                            tempdat.dsqw = tempdat.dsqw[:, -1:-1:0]
                            tempdat.q = tempdat.q[-1:-1:0]
                        if tempdat.q[0] > qlimits.min():
                            qindex.append(0)
                        if tempdat.q[-1] < qlimits.min():
                            raise NameError('lower limit is higher than highest momentum transfer!')
                        else:
                            for hiq1, _ in enumerate(tempdat.q[0:-2]):
                                if (tempdat.q[hiq1] < qlimits.min()) and (tempdat.q[hiq1+1] > qlimits.min()):
                                    qindex.append(hiq1)
                                    break
                        if tempdat.q[-1] < qlimits.max():
                            qindex.append(len(tempdat.q))
                        if tempdat.q[0] > qlimits.max():
                            raise NameError('upper limit is lower than lowest momentum transfer!')
                        else:
                            for hiq2, _ in enumerate(tempdat.q[0:-2]):
                                if (tempdat.q[hiq2] < qlimits.max()) and (tempdat.q[hiw2+1]>qlimits.max()):
                                    qindex.append(hiq2+1)
                                    break
                # first reshape hw
                if not len(hwindex) == 0:
                    tempdat.sqw = Dat.sqw[hwindex[0]:hwindex[1], :]
                    tempdat.dsqw = Dat.dsqw[hwindex[0]:hwindex[1], :]
                    tempdat.hw = Dat.hw[hwindex[0]:hwindex[1]]
                # now reshape q
                if not len(qindex) == 0:
                    tempdat.sqw = tempdat.sqw[:, qindex[0]:qindex[1]]
                    tempdat.dsqw = tempdat.dsqw[:, qindex[0]:qindex[1]]
                    tempdat.q = tempdat.q[qindex[0]:qindex[1]]
                tempdat.Fit(fitname, Vanadium)
                Dat.Fits = copy.deepcopy(tempdat.Fits)
                pbar.update(1)

    def Subtract(self, subdata, scaling):
        for _, Dat in enumerate(self.Data):
            Dat.Subtract(subdata.Data, scaling)
