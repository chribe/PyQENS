#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 17:40:32 2020
@author: beck
"""
import Experiment as EXP
import PyQENSUtilities as PQU
import builtins  # to set systemwide golbal parameters
import numpy as np
#%% set global flags
builtins.flags = {'savepath': './figures/',
                  'mantidpath': '/opt/mantidnightly/',
                  'plotresults': 'showsave',  # all, results, fits, silent, none, writeascii
                  'pictureextension': 'pdf',  # eps png pdf
                  'writesummary': 'none',  # all fits results none
                  'savename': 'BLG_CdCl2_Crystallization'
                  }
builtins.VanadiumFit = []
builtins.currq = np.NAN
builtins.Sampleidx = 'test'

#%% define protocol for one beamtime
beamtimeproperties = {'beamtime': '9-13-829',
                      'instrument': 'in16b',
                      'container': 'cylinder',
                      'filepath_raw': '/home/tofhr/beck/Documents/nonclassical_crystallization_BLG_CdCl2/exp_8-04-810_181_in16b/processed/rawdata_corrected_SANS_angles/',
                      'filepath_analysed': '/home/tofhr/beck/Documents/nonclassical_crystallization_BLG_CdCl2/exp_8-04-810_181_in16b/processed/reduceddataPyQENS_2/',
                      'ECtreatment': 'subtract',  # subtract model
                      'Solvent': {'treatment': 'subtract',  # subtract, model
                                  'scaling': 'fixed',
                                  'model': 'D2Ofixed'},
                      'groupingfile': '/home/tofhr/beck/Documents/nonclassical_crystallization_BLG_CdCl2/exp_8-04-810_181_in16b/processed/draft_grouping_IN16B_cycle201.xml',
                      'qindex': [3, 18]
                      }
Vanadium = {'bkg': 'flat',  # flat, slope, none
            'NumberOfGaussians': 2,
            'NumberOfLorentzians': 0,
            'NumberOfVoigts': 0}
Protocol = []
Protocol.append({'Sample': 'D2O', 'type': 'Solvent', 'Runfiles': list(range(216682, 216692))})
Protocol.append({'Sample': 'Vanadium', 'Runfiles': list(range(217518, 217560)), 'VanadiumProperties': Vanadium})
Protocol.append({'Sample': 'EmptyCan', 'Runfiles': list(range(216672, 216682))})
Protocol.append({'Sample': 'BLG ZnCl2', 'Runfiles': list(range(216692, 217517))})  # ,'binning':'floathour4'
Protocol = PQU.fillupProtocol(Protocol, beamtimeproperties)
PQU.CheckProtocol(Protocol)
PQU.CheckPaths(beamtimeproperties)
#%% create Experiment Class from Protocol
Experiment = []
Experiment.append(EXP.Experiment())
Experiment[-1].CreatefromProtocol(Protocol)
PQU.save_pickle(Experiment, flags['savename']+ '.p')
#%% perform subtraction
Experiment[-1].Subtract('EmptyCan', includealso=['Vanadium', 'Solvent'])
Experiment[-1].Subtract('D2O', Scaling=[0.7])
#%% perform fits
Experiment[-1].Fit(['Vanadium', 'twolorentz_free_dirac_q', 'immobile_jump_diffusion_internal_jump'])
PQU.save_pickle(Experiment, flags['savename'] + '.p')