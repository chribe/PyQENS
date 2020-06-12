#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 11:13:01 2020

@author: beck
"""

import Sample as SP


class Experiment:
    '''
    This class contains the information of an experiment.
    It contains the following information:
        beamtime
        instrument
        Container
        list of class Samples

    function contained in this class:
        Create

    '''
    #%% import utilities
    #%% initialize the class

    def __init__(self,
                 instrument='',
                 beamtime='',
                 container='',
                 Sample=[]):
        # initialising all the variables
        self.instrument = instrument
        self.container = container
        self.Vanadium = []
        self.EmptyCan = []
        self.Solvent = []
        self.Sample = []
        self.beamtime = ''
#%% defining the functions

    def CreatefromProtocol(self, Protocol):
        '''
        This function creates the experiment class specified in the protocol

        '''
        solvent = 0
        Vanadium = 0
        for _, Protentry in enumerate(Protocol):
            if Protentry['Sample'] == 'Vanadium':
                self.Vanadium.append(SP.Sample())
                self.Vanadium[-1].CreateSample(Protentry)
                Vanadium += 1
        if Vanadium == 0:
            raise NameError('No Vanadium Entry!')
        for hi, Protentry in enumerate(Protocol):
            if Protentry['Sample'] == 'Vanadium':
                1 == 1
            elif Protentry['type'] == 'Solvent':
                self.Solvent.append(SP.Sample())
                self.Solvent[-1].CreateSample(Protentry)
                solvent += 1
            elif Protentry['Sample'] == 'EmptyCan':
                self.EmptyCan.append(SP.Sample())
                self.EmptyCan[-1].CreateSample(Protentry)
            else:
                self.Sample.append(SP.Sample())
                self.Sample[-1].CreateSample(Protentry)
            # check of consistency
            if self.instrument == '':
                self.instrument = Protentry['instrument']
            elif self.instrument != Protentry['instrument']:
                Warning('Change of Instrument during one experiment!' +
                        ' Please verify! Change from "' + self.instrument +
                        '" to "' + Protentry['instrument'])
                self.instrument = Protentry['instrument']
            if self.beamtime == '':
                self.beamtime = Protentry['beamtime']
            elif self.beamtime != Protentry['beamtime']:
                Warning('Change of Beamtime during one experiment!' +
                        ' Please verify! Change from "' + self.beamtime +
                        '" to "' + Protentry['beamtime'])
                self.beamtime = Protentry['beamtime']
            if self.container == '':
                self.container = Protentry['container']
            elif self.container != Protentry['container']:
                Warning('Change of Container during one experiment!' +
                        ' Please verify! Change from "' + self.container +
                        '" to "' + Protentry['container'])
                self.container = Protentry['container']
#%% Fit

    def Fit(self, fitname):
        if type(fitname) == str:
            fitn = [fitname]
        else:
            fitn = fitname
        for _, FitName in enumerate(fitn):
            for _, prop in enumerate(['hwlimits', 'qlimits', 'hwindex', 'qindex']):
                if prop not in self.Vanadium[-1].ProtEntry:
                    self.Vanadium[-1].ProtEntry.update({prop: []})
            if FitName == 'Vanadium':
                self.Vanadium[-1].Fit(FitName, Vanadium=self.Vanadium,
                                      hwlimits=self.Vanadium[-1].ProtEntry['hwlimits'],
                                      qlimits=self.Vanadium[-1].ProtEntry['qlimits'],
                                      hwindex=self.Vanadium[-1].ProtEntry['hwindex'],
                                      qindex=self.Vanadium[-1].ProtEntry['qindex'])
            elif FitName == 'D2O':
                for _, Samp in enumerate(self.Solvent):
                    Samp.Fit(FitName, Vanadium=self.Vanadium,
                             hwlimits=self.Vanadium[-1].ProtEntry['hwlimits'],
                             qlimits=self.Vanadium[-1].ProtEntry['qlimits'],
                             hwindex=self.Vanadium[-1].ProtEntry['hwindex'],
                             qindex=self.Vanadium[-1].ProtEntry['qindex'])
            else:
                for _, Samp in enumerate(self.Sample):
                    Samp.Fit(FitName, Vanadium=self.Vanadium,
                             hwlimits=self.Vanadium[-1].ProtEntry['hwlimits'],
                             qlimits=self.Vanadium[-1].ProtEntry['qlimits'],
                             hwindex=self.Vanadium[-1].ProtEntry['hwindex'],
                             qindex=self.Vanadium[-1].ProtEntry['qindex'])
    #%% define subtraction
    def Subtract(self, name, Specification='', Scaling=[1],
                 method='simple', includealso=[]):
        '''
        This function subtracts a data set from another. as standard method,
        no scaling is considered. the scaling parameter should have either
        one entry or should correspond to the length of q. As method, there are
        simple and Paalman Pings

        Via "includealso" the subtraction is also performed for the calibration
        files mentioned as input.
        '''

        # Selection of the data set
        if name == 'Solvent':
            dataforsubtraction = self.Solvent
        elif name == 'EmptyCan':
            dataforsubtraction = self.EmptyCan
        elif name == 'Vanadium':
            dataforsubtraction = self.Vanadium
        else:
            dataforsubtraction = self.Sample
        if len(dataforsubtraction) > 1:
            found = 0
            for _, dat in enumerate(dataforsubtraction):
                if dat.Name == Specification:
                    found = 1
                    datasub = dat
                    break
            if found == 0:
                raise NameError("Data with " + Specification + " not found!")
        else:
            datasub = dataforsubtraction[0]
        #Paalman Pings 
        # TODO 

        # do subtraction
        for _, Sample in enumerate(self.Sample):
            Sample.Subtract(datasub, Scaling)
        for _, Samp in enumerate(includealso):
            Sampl = getattr(self, Samp)
            for _,Sample in enumerate(Sampl):
                Sample.Subtract(datasub, Scaling)
