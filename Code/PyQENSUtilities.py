#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 19:02:45 2020

@author: beck
"""
import h5py
import numpy as np
import pickle
import os as os
from shutil import copyfile
#%% functions for the protocol
def fillupProtocol(Protocol, beamtimeproperties):
    ''' 
    This function fills up the protocol with some general information, so that 
    the one which is written by hand is as short as possible, but the code 
    later has no problems with missing entries.
    '''
    #% check if paths are written in a good way:
    #TODO
    # complete Protocol
    VanadiumRuns=[]
    for _, Protentry in enumerate(Protocol):
        if 'beamtime' not in Protentry:
            Protentry.update({'beamtime':beamtimeproperties['beamtime']})
        if 'binning' not in Protentry:
            Protentry.update({'binning':'all'})
        if 'container' not in Protentry:
            Protentry.update({'container':beamtimeproperties['container']})
        if 'instrument' not in Protentry:
            Protentry.update({'instrument':beamtimeproperties['instrument']})
        if 'type' not in Protentry and Protentry['Sample']=='Vanadium':
            Protentry.update({'type':'Calibration'})
        elif 'type' not in Protentry and Protentry['Sample']=='EmptyCan':
            Protentry.update({'type':'Calibration'})
        if'type' not in Protentry:
            Protentry.update({'type':'Sample'})
        Protentry.update({'filepath_raw':beamtimeproperties['filepath_raw']})
        if 'groupingfile' in beamtimeproperties and 'groupingfile' not in Protentry:
            Protentry.update({'groupingfile':beamtimeproperties['groupingfile']})
        if 'filepath_analysed' in beamtimeproperties and 'filepath_analysed' not in Protentry:
            Protentry.update({'filepath_analysed':beamtimeproperties['filepath_analysed']})
        if 'ECtreatment' not in Protentry:
            Protentry.update({'ECtreatment':beamtimeproperties['ECtreatment']})
        if Protentry['Sample']=='Vanadium':
            for _,num in enumerate(Protentry['Runfiles']):
                copyfile(Protentry['filepath_raw'] + str(num) + '.nxs', "currentfile.nxs")
                f=h5py.File('currentfile.nxs', 'r')
                maxenergytransfer=np.array(f['/entry0/instrument/Doppler/maximum_delta_energy/'])
                if maxenergytransfer>20: # TODO
                    VanadiumRuns.append(num)
                os.remove("currentfile.nxs")
    for _, Protentry in enumerate(Protocol):
        Protentry.update({'VanadiumRuns':VanadiumRuns})
        Protentry.update({'Solvent':beamtimeproperties['Solvent']})
        for _,prop in enumerate(['hwlimits','hwindex','qindex','qlimits']):
            if prop in beamtimeproperties:
                Protentry.update({prop:eval('beamtimeproperties["' + prop + '"]')})
    return Protocol

def CheckProtocol(Protocol):
    Vana=[]
    sol=[]
    ec=[]
    for _, Protentry in enumerate(Protocol):
        if Protentry['type']=='Solvent':
            sol.append(Protentry)
        if Protentry['Sample']=='Vanadium':
            Vana.append(Protentry)
        if Protentry['Sample']=='EmptyCan':
            ec.append(Protentry)
    if len(Vana)==0:
        raise NameError('No Vanadium found!')
    if len(sol)==0:
        raise Warning('No Solvent data found!')
    if len(ec)==0:
        raise Warning('No empty Can found! Add an elastic contribution to the fit...')
        flag['ECtreatment']='model'

def CheckPaths(beamtimeproperties):
    for _ ,path in enumerate([flags['savepath'], flags['mantidpath'],
                              beamtimeproperties['filepath_raw'],
                              beamtimeproperties['filepath_analysed']]):
        while path[-1]==' ':
            path=path[0:-2]
        if path[-1]!='/':
            path+='/'
        while path[0]==' ':
            path=path[1:]
        if path[0:1]=='./':
            path=os.getcwd() + path[1:]
        if path[0]!='/':
            path=os.getcwd() + '/' + path
    if not os.path.isdir(flags['savepath']):
        try:
            os.mkdir(flags['savepath'])
        except OSError:
            print ("Creation of the directory %s failed" % flags['savepath'])
        else:
            print ("Successfully created the directory %s " % flags['savepath'])

#%% functions for writing summary
def ConvertModel2LaTeX(string):
    from collections import defaultdict
    import sympy
    class GenerateSymbols(defaultdict):
        def __missing__(self, key):
            return sympy.Symbol(key)
    return sympy.latex(eval(string,GenerateSymbols()))
#%% load and save
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
