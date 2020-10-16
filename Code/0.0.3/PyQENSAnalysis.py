#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 10:48:17 2020

This is PyQENSAnalysis. It takes over the fitting part of the QENS analysis.

File written by Christian Beck
beck@ill.eu
ORCID: 0000-0001-7214-3447

"""

from lmfit import Model,fit_report
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes as axes
import PyQENSModels as PQM
import builtins
from tqdm import tqdm
import os
#%% Analysis frameworks for FWS
def FWS_MSD(Data):
    for _,Dataentry in enumerate(Data):
       if Dataentry['Type']=='FWS':
           #TODO perform poly2 fit on FWS to extract MSD, Save in Structure
           print('test')
    return Data

def FWS_Ratio(Data):
    for _,Dataentry in enumerate(Data):
       if Dataentry['Type']=='FWS':
           print('test')
    return Data
#%% define main parts of fitfunction
def gaussian(x,bkg, amp, cen, wid):
    return bkg+ amp * np.exp(-(x-cen)**2 / wid)

#%% Analysis of QENS
def QENS_Fit(Data,Type='Sample',exclude=[],qrange=[]):
    if Type=='Calibration':
        for Numgaus in range(0,NumberOfGaussians):
            if Numgaus==0:
                ResolutionFitModel=Model(gaussian,prefix='g' + str(Numgaus) + '_')
            else:
                ResolutionFitModel+=Model(gaussian,prefix='g' + str(Numgaus) + '_')
        params = ResolutionFitModel.make_params()
        for _,par in enumerate(params.keys()):
            if 'wid' in par and Numgaus==0:
                params.add(par, value=4, min=1e-9)
            if 'wid' in par and not Numgaus==0:
                params.add(par, value=0.8, min=1e-9)
            if 'amp' in par:
                params.add(par, value=0.1, min=0)
            if 'cen' in par:
                params.add(par, value=0,vary=False)
            if 'bkg' in par:
                params.add(par, value=0,vary=False)
        params.add('g0_bkg', value=0.001,min=0,vary=True)
        for _,LE in enumerate(Calibration[0]):
            FitResults=[]
            if qrange!=[]:
                raise NameError("Resolution should be fitted for all 'q'!")
            for hiq,q in enumerate(LE['q']):
                XN,SQWN,DSQWN=extract_X_Sqw_Q(LE,hiq)
                FitResults.append(ResolutionFitModel.fit(SQWN, params, x=XN,weights=DSQWN**-2))
                plotFit(FitResults[-1],XN,'Vanadium/','q_'+ "%.2f"%q +'.pdf')
        Calibration[0][0].update({'Fitresults':[]})
        Calibration[0][0]['Fitresults'].append(FitResults)
    else:
        builtins.qrange=qrange
        params, FitModel,FitLevel=PQM.DefineFitModel(Type)
        if FitLevel=='q':
            with tqdm(total=len(Data)) as pbar:
                for _,Dat in enumerate(Data):
                    for _,DAT in enumerate(Dat):
                        if 'Fit' not in DAT.keys():
                            DAT.update({'Fit':[]})
                        FIT=[]
                        for hiq,q in enumerate(DAT['q']):
                            if qrange==[] or (q<qrange[1] and q>qrange[0]):
                                XN,SQWN,DSQWN=extract_X_Sqw_Q(DAT,hiq)
                                builtins.currqidx=hiq
                                FIT.append(FitModel.fit(SQWN, params, x=XN,weights=DSQWN**-2))
                                plotFit(FIT[-1],XN,DAT['Name'],'q_'+ "%.2f"%q +'.pdf')
                        plotfitresults(FIT,DAT['Name'],DAT['q'],qrange)
                        DAT['Fit'].append(FIT)
                    pbar.update(1)
        elif FitLevel=='S':
            for _,Dat in enumerate(Data):
                with tqdm(total=len(Data)) as pbarS:
                    for _,DAT in enumerate(Dat):
                        if 'Fit' not in DAT.keys():
                            DAT.update({'Fit':[]})
                        XN,SQWN,DSQWN=extract_X_Sqw(DAT)
                        DAT['Fit'].append(FitModel.fit(SQWN,params,x=XN,weights=DSQWN**-2))
                        plotFit2D(DAT,DAT['Name'],'globalfit_q_')
                        pbarS.update(1)
    return Data

def plotFit(FitResult,x,subfolder,name):
    if not os.path.isdir(figure_folder +subfolder + '/' ):
        os.mkdir(figure_folder +subfolder + '/' )
    FitResult.plot(xlabel='$\hbar\omega$',ylabel='S($q,\omega$) [arb.units]',yerr=FitResult.weights**-0.5,ax_fit_kws={'YScale':'log'},fit_kws={'zorder':5,'linewidth':5,'linestyle':'dashed'})
    dely = FitResult.eval_uncertainty(sigma=3)
    plt.fill_between(x, FitResult.best_fit-dely, FitResult.best_fit+dely, color="#ABABAB",label='3-$\sigma$ uncertainty band')
    comps = FitResult.eval_components(x=x)
    plt.autoscale(False)
    zval=10
    for _,val in enumerate(comps):
        zval=10*zval
        plt.plot(x,comps[val],zorder=zval)
    plt.savefig(figure_folder +subfolder + '/' + name)
    plt.close()
    
def plotfitresults(FIT,name,Q,qrange):
    Parameternames=FIT[0].params.keys()
    for _,parameter in enumerate(Parameternames):
        if FIT[0].params[parameter].expr == None:
            y=[]
            dy=[]
            qp=[]
            for hiq,q in enumerate(Q):
                if qrange==[] or (q>qrange[0] and qrange[1]>q):
                    y.append(FIT[hiq].params[parameter].value)
                    dy.append(FIT[hiq].params[parameter].stderr)
                    qp.append(q)
            qp=np.array(qp)
            plt.figure()
            plt.errorbar(qp, y, dy)
            plt.savefig(figure_folder + name + '/' + parameter + 'q.pdf')
            plt.close()
            plt.figure()
            plt.errorbar(qp**2, y, dy)
            plt.savefig(figure_folder + name + '/' + parameter + 'q2.pdf')
            plt.close()

def plotFit2D(DAT,subfolder,name):
    if not os.path.isdir(figure_folder +subfolder + '/' ):
        os.mkdir(figure_folder +subfolder + '/' )
    results=DAT['Fit'][-1].best_fit
    lr=int(len(results)/len(DAT['hw']))
    results.resize(lr,len(DAT['hw']))
    components=DAT['Fit'][-1].eval_components()
    for hi,comp in enumerate(components):
        components[comp].resize(lr,len(DAT['hw']))
    Q=np.array(Calibration[0][0]['q'])
    if qrange!=[]:
        q=Q[(Q>qrange[0])&(Q<qrange[1])]
    else:
        q=Q
    for hiq,q in enumerate(q):
        plt.figure()
        plt.errorbar(Calibration[0][0]['hw'],np.array(Calibration[0][0]['sqw'][q==Calibration[0][0]['q']])[0],np.array(Calibration[0][0]['dsqw'][q==Calibration[0][0]['q']])[0])
        plt.errorbar(DAT['hw'],np.array(DAT['sqw'][hiq])[0],np.array(DAT['dsqw'][hiq])[0])
        zval=10
        for _,val in enumerate(components):
            zval=10*zval
            plt.plot(DAT['hw'],components[val][hiq],zorder=zval)
        plt.plot(DAT['hw'],results[hiq],zorder=zval*10)
        plt.yscale('log')
        plt.ylim(bottom=1e-6)
        plt.savefig(figure_folder +subfolder + '/' + name+ "%.2f"%q +'.pdf')
        plt.close()
        
def extract_X_Sqw(LE):
    XN=np.array(LE['hw'])*1000
    x=XN
    Q=LE['q']
    sqw=np.array(LE['sqw'])
    dsqw=np.array(LE['dsqw'])
    sqwn=[]
    dsqwn=[]
    for hiq,q in enumerate(Q):
        if qrange==[] or (q<qrange[1] and q>qrange[0]):
            sqwn.append(sqw[hiq])
            dsqwn.append(dsqw[hiq])
    SQW=np.array(sqwn).flatten()
    DSQW=np.array(dsqwn).flatten()
    SQW[  np.isnan( SQW )  ] = 0
    SQW[  np.isinf( SQW )  ] = 0
    SQW[  SQW < 1e-6*SQW.max() ] = 0 # reasonable signal-to-noise range
    DSQW[ np.isnan( SQW )  ] = np.inf
    DSQW[ np.isinf( SQW )  ] = np.inf
    DSQW[ np.isnan( DSQW ) ] = np.inf
    DSQW[  SQW==0 ] = np.inf # no signal at all should have infinite error
    DSQW[ DSQW==0 ] = np.inf # "a measurement without error is nonsense"
    
    return x,SQW,DSQW


def extract_X_Sqw_Q(LE,hiq):
    x=np.array(LE['hw'])*1000
    sqw=np.array(LE['sqw'][hiq])
    dsqw=np.array(LE['dsqw'][hiq])
    SQW=sqw[0]
    DSQW=dsqw[0]
    SQW[  np.isnan( SQW )  ] = 0
    SQW[  np.isinf( SQW )  ] = 0
    SQW[  SQW < 1e-6*SQW.max() ] = 0 # reasonable signal-to-noise range
    DSQW[ np.isnan( SQW )  ] = np.inf
    DSQW[ np.isinf( SQW )  ] = np.inf
    DSQW[ np.isnan( DSQW ) ] = np.inf
    DSQW[  SQW==0 ] = np.inf # no signal at all should have infinite error
    DSQW[ DSQW==0 ] = np.inf # "a measurement without error is nonsense"
    SQWN=SQW[x!=0]
    DSQWN=DSQW[x!=0]
    XN=x[x!=0]
    return XN,SQWN,DSQWN