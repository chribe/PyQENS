#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 10:46:45 2020
This is the PyQENSLib Library which will take over all tasks concerning the 
analysis of the QENS data excluding the fitting procedure


File written by Christian Beck
beck@ill.eu
ORCID: 0000-0001-7214-3447

"""
import os
import time
import numpy as np
import h5py
import pickle
#%%
def open_pickle(input):
    '''
    Open pickle files (.p)
    '''
    with open(input, 'rb') as infile:
        output=pickle.load(infile)
    return output

def save_pickle(input,name_out):
    '''
    Convert a variable in pickle (.p) file with name read in name_out.
    '''
    with open(name_out, 'wb') as outfile:
        pickle.dump(input, outfile, protocol=pickle.HIGHEST_PROTOCOL)
    
#%% some short functions
def SampleProperties(hfile,h5path):
    data={'ExperimentTitle':str(hfile[ h5path+'logs/title/value/'][0])[2:-1],
          'SampleTitle':str(hfile[ h5path+'logs/subtitle/value/'][0])[2:-1],
          'StartTime':str(hfile[ h5path+'logs/start_time/value/'][0])[2:-1],
          'EndTime':str(hfile[ h5path+'logs/end_time/value/'][0])[2:-1],
          'SetTemperature':np.array( hfile[ h5path+'logs/sample.setpoint_temperature/value/']),
          'SampleTemperatureValue':np.array( hfile[ h5path+'logs/sample.temperature/value']),
          'SampleRegulationTemperature':np.array( hfile[ h5path+'logs/sample.regulation_temperature/value']),
          'SamplePressureValue':np.array( hfile[ h5path+'logs/sample.pressure/value']),
          'SampleSetPointPressure':np.array( hfile[ h5path+'logs/sample.setpoint_pressure/value']),
          'RunNumber':np.array( hfile[ h5path+'logs/run_number/value']),
          'ExperimentNumber':str(hfile[ h5path+'logs/experiment_identifier/value/'][0])[2:-1]}
    try:
        Timdep={'SampleTemperatureTime':np.array( hfile[ h5path+'logs/sample.temperature/time']),'SamplePressureTime':np.array( hfile[ h5path+'logs/sample.pressure/time'])}
        print(Timedep)
    except:
        Timedep={'SampleTemperatureTime':[],'SamplePressureTime':[]}
    data.update(Timedep)
    return data

def loadFWS(folder,name):
    h5file=folder + name
    hfile = h5py.File( h5file, 'r' )
    data=[]
    for k in range(1,len(hfile.keys())+1):
        h5path = '/mantid_workspace_'+str(k)+'/' 
        dat={'Type':'FWS',
              'Name':name,
              'ModTime':time.ctime(os.path.getmtime(folder + name))}
        SP=SampleProperties(hfile,h5path)
        dat.update(SP)
        data.append(dat)
    hfile.close()
    #TODO add specific data loading
    return data

def MergeData(Data,mergeQENSFWS):
    if mergeQENSFWS==1:
        for hiD,Dat in enumerate(Data):
            for _,DAT in enumerate(Dat):
                print(DAT)
    return Data

def loadQENS(folder,name):
    h5file=folder + name
    hfile = h5py.File( h5file, 'r' )
    data=[]
    for k in range(1,len(hfile.keys())+1):
        h5path = '/mantid_workspace_'+str(k)+'/' 
        dat={'Type':'QENS',
              'Name':name,
              'ModTime':time.ctime(os.path.getmtime(folder + name)),
              'sqw':np.matrix(hfile[ h5path+'workspace/values/'] ),
              'dsqw':np.matrix(hfile[ h5path+'workspace/errors/' ] ),
              'q':np.array( hfile[ h5path+'workspace/axis2/']),
              'hw':np.array( hfile[ h5path+'workspace/axis1/'][0:-1])
              }
        #dat['hw']=np.delete(dat['hw'],-1)
        SP=SampleProperties(hfile,h5path)
        dat.update(SP)
        data.append(dat)
    hfile.close()
    #TODO add specific data loading
    return data

def load_data(folder):
    FilesInFolder = [f for f in os.listdir(folder) if os.path.isfile(folder + f)]
    Data=[]
    for _,name in enumerate(FilesInFolder):
        if '_FWS_q.nxs' in name:
            data=loadFWS(folder,name)
        if  '_QENS_q.nxs' in name:
                    data=loadQENS(folder,name)
        Data.append(data)
    return Data

def FindEntry(Data,Samplename):
    Entry=[]
    for _,entry in enumerate(Data):
        if entry[0]['Name']==Samplename:
            Entry.append(entry)
    return Entry
#%% Background treatment
def BackgroundTreatment(Data,exclude=[]):
    Bkg=Background[0][0]
    #%%% check scaling of background: first scalar is treated
    if type(backgroundscaling)==str and backgroundscaling=='Paalman-Pings':
        #TODO write Paalman Pings correction TO BE DISCUSSED -- should be in the mantid routine!
        raise NameError('Paalman-Pings is not yet implemented!')
    elif type(backgroundscaling)==np.ndarray and size(backgroundscaling)==size(Bkg['q']):
        correction=backgroundscaling
    elif type(backgroundscaling)==int or type(backgroundscaling)==float:
        
        correction=np.array(np.ones(len(Bkg['q'])))*backgroundscaling
    else:
        raise NameError('The Entry backgroundscaling with the value' + str(backgroundscaling) + ' does not fit the requirements!')
    #%%% scale background data
    for hiq,_ in enumerate(Bkg['q']):
        Bkg['sqw'][hiq]=Bkg['sqw'][hiq]*correction[hiq]
        Bkg['dsqw'][hiq]=Bkg['dsqw'][hiq]*correction[hiq]
    #%%% now subtracting it from the different data
    if backgroundtreatment=='subtract' and backgrounduseage=='data':
        for _,data in enumerate(Data):
            if data[0]['Name'] not in exclude and data[0]['Type']=='QENS':
                data[0]['sqw']=data[0]['sqw']-Bkg['sqw']
                data[0]['dsqw']=np.power(np.power(data[0]['dsqw'],2)+np.power(Bkg['dsqw'],2),0.5)
    # TODO add option: add to fit
    # TODO            use model
    return Data
    
def SubtractSolvent(dat,datsub,scaling):
    dat[0]['sqw']+=-1*scaling*datsub[0]['sqw']
    dat[0]['dsqw']=np.power(np.power(dat[0]['dsqw'],2) +scaling*np.power(datsub[0]['dsqw'],2),0.5)
    return dat

def VolumeFraction(Protein='',cpreal=[],cpnom=[],T=[]):
    '''
    %==========================================================================
    %
    % For given preparative protein concentration cp [mg Protein / ml solvent]
    % and temperature T [K] calculate resulting volume fraction
    %
    % phi_hydr:     volume fraction of protein + hydration shell
    % dphi_hydr:    error of phi_hydr
    % phi_dry:      volume fraction of protein
    % dphi_dry:     error of phi_dry
    % n:            number density of proteins in [#/m^3]
    % dn:           error of n
    %
    %==========================================================================
    '''
    with open("ProteinProperties.txt") as f:
        data=f.readlines()
    for line in data:
        entries=line.split(';')
        if Protein==entries[0]:
            v0=float(entries[1])
            T0=float(entries[2])
            Mw=float(entries[3])
    '''
    % Hydration level of BSA, eta = 0.4 g H2O / g protein
    % Kuntz, I.D.; Kauzmann, W. Adu. Protein Chem. 1974,28, 239. 
    % In the softmatter lab I measured 0.44 g/g protein with D2O, The most
    % conservative value in the community is however 0.3 g/g protein (Doster & Co.)
    % Chen measured 0.34 g/g protein, J. Phys. Chem. lQ03,87, 1473-1477
    '''
    eta = 0.4 # [g D2O(H2O)/g protein]
    '''
      % ratio of densities of hyrdration shell to bulk, chi =
    % rho_hydr/rho_bulk, most scientisct in the field of biophysics believe
    % that the density of water on the protein surface is 3-22% higher than
    % in bulk water 
    '''
    p_min = 0.03
    p_max = 0.22
    chi = 1 + ( p_max + p_min ) / 2
    dchi = ( p_max - p_min ) / 2
    v = v0 - eta / DensityH2O( T0 ) * ( 1 / chi - 1 )
    dv = eta / DensityH2O( T0 ) * dchi / chi**2
    v1  = v  + eta / DensityD2O( T  ) * ( 1 / chi - 1 )
    dv1 = dv + eta / DensityD2O( T  ) * dchi / chi**2
    if cpreal!=[]:
        CPREAL=cpreal
        CPNOM=CPREAL/ ( 1 -CPREAL/1000* v1 )/1000
        dCPREAL=0
    elif cpnom!=[]:
        CPNOM=cpnom/1000
        CPREAL=CPNOM / ( 1 +CPNOM* v1 )*1000
        dCPREAL = CPNOM**2/ ( 1 + CPNOM* v1 )**2 * dv1 * 1000
    else:
        raise NameError('No concentration given!')
    NA = 6.022142E23
    n  = CPREAL  / Mw * NA * 1E3
    dn = dCPREAL / Mw * NA * 1E3
    phi_dry  = CPNOM * v / ( 1 + CPNOM * v1 )
    dphi_dry = CPNOM / ( 1 + CPNOM * v1 )**2 * ( CPNOM**2 * dv1**2 * v**2 + ( dv + CPNOM *dv*v1 )**2 )**0.5
    phi_hydr  = ( eta / DensityD2O( T ) + v1 ) * CPNOM / ( 1 + v1 * CPNOM )
    dphi_hydr = dv1 * CPNOM / ( 1 + CPNOM * v1 ) * ( 1 - phi_hydr )
    return phi_dry,dphi_dry,phi_hydr,dphi_hydr,n,dn


def ViscosityD2O( T, dT=0 ):
    '''
    %==========================================================================
    %
    % Calculates the viscosity of D2O in units of [Pa*s] for given temperature in
    % units of [K]
    % Dependence of viscosity on temperature was determined by Cho et al. (1999)
    % "Thermal offset viscosities of liquid H2O, D2O and T2O", J.Phys. Chem. B
    % 103(11):1991-1994
    %
    % eta(T) is valid from 280k up to 400K
    %
    % eta = ViscosityD2O( T ) or [ eta deta ] = ViscosityD2O( T, dT )
    %
    %==========================================================================
    '''
  
    C  = 885.60402;
    a  = 2.799 * 10**(-3)
    b  = -1.6342 * 10**(-5)
    c  = 2.9067 * 10**(-8)
    g  = 1.55255
    T0 = 231.832
    t = T - T0
    
    eta  = 1E-3 * C * ( t + a * t**2 + b * t**3 + c * t**4 )**(-g)
    detadT = -1E-3 * C * g * ( 1 + 2 * a * t + 3 * b * t**2 + 4 * c * t**3 ) *( t + a * t**2 + b * t**3 + c * t**4 ) ** ( - 1 - g )  
    if dT==[]:
        deta = detadT * dT
    else:
        deta=[]
    return eta, deta, detadT 

def DensityH2O(T,dT=0):
    t = T - 273.15
    rho = 1E-3 / (1 + 16.87985E-3 *  t ) * ( 999.83952 + 16.945176 * t - 7.9870401E-3 * t**2 - 46.17046E-6 * t**3 + 105.56302E-9 * t**4 - 280.54253E-12 * t**5)
    return rho

def DensityD2O( T ):
    '''
        
    %==========================================================================
    %
    % Calculates the density of D2O in units of [g/cm^3] for given temperature  
    % T in [K]
    % Handbook of Chemistry and Physics, 73rd Edition, 
    % Lide, D.R. Ed.; CRC Press: Boca Raton 1992; Chapter 6, pg. 13
    % http://physchem.uni-graz.ac.at/sm/Service/Water/D2Odens.htm
    %
    %==========================================================================

    '''
    T_table   = 273.15 + np.array([3.82,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105])
    rho_table = np.array([1.1053, 1.1055,1.1057, 1.1056, 1.1050, 1.1044, 1.1034, 1.1019,1.1001, 1.0979, 1.0957, 1.0931, 1.0905, 1.0875, 1.0847, 1.0815, 1.0783, 1.0748, 1.0712, 1.0673, 1.0635, 1.0598 ])
    rho = np.interp( T, T_table, rho_table)
    return rho

def IntegrateSQW(Data,suffix=''):
    for _,Dat in enumerate(Data):
        for _,DAT in enumerate(Dat):
            Sq=[]
            if DAT['Type']=='QENS':
                for hiq,_ in enumerate(DAT['q']):
                    print(DAT['sqw'][hiq])
                    Sq.append(np.sum(DAT['sqw'][hiq]))
                DAT.update({'integrSq'+suffix:Sq})
    return Data


def trim_axs(axs, N):
    """little helper to massage the axs list to have correct length..."""
    axs = axs.flat
    for ax in axs[N:]:
        ax.remove()
    return axs[:N]