#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 11:32:03 2020

Code written by Christian Beck
christian.beck@uni-tuebingen.de
ORCID: 0000-0001-7214-3447
"""
import numpy as np
from lmfit import Model, Parameters
from scipy.special import wofz
import builtins

def gaussian(x,bkg, amp, cen, wid):
    return bkg+ amp * np.exp(-(x-cen)**2 / wid)

def V(x, alpha, gamma):
    """
    Return the Voigt line shape at x with Lorentzian component HWHM gamma
    and Gaussian component HWHM alpha.
    copied from https://scipython.com/book/chapter-8-scipy/examples/the-voigt-profile/

    """
    sigma = alpha / np.sqrt(2 * np.log(2))

    return np.real(wofz((x + 1j*gamma)/sigma/np.sqrt(2))) / sigma/np.sqrt(2*np.pi)

def Lorentz_convol(x,width,scaling):
    """
    x is the energy transfer
    width is the HWHM from the Lorentzian
    """
    y=np.zeros(len(x))
    for hi in range(0,NumberOfGaussians):
        Rparam=Calibration[0][0]['Fitresults'][0][currqidx].params
        cen=Rparam['g' + str(hi) + '_cen'].value
        amp=Rparam['g' + str(hi) + '_amp'].value
        wid=Rparam['g' + str(hi) + '_wid'].value
        y+=scaling*amp*V(x-cen,wid,width)
    return y
def Lorentz_convol_global(x,**parameterlist):
    yf=[]
    for hiq,qv in enumerate(q):
        builtins.currqidx=hiq
        scaling=parameterlist['scaling%i'%hiq]
        width=parameterlist['width%i'%hiq]
        yf.append(Lorentz_convol(x,width,scaling))
    yf=np.array(yf)
    return yf.flatten()


def Resolution(x):
    y=np.zeros(x.shape())
    for hi in range(0,NumberOfGaussians):
        Rparam=Calibration[0][0]['Fitresults'][currqidx][0]
        cen=Rparam['g' + str(hi) + '_cen'].value
        amp=Rparam['g' + str(hi) + '_amp'].value
        wid=Rparam['g' + str(hi) + '_wid'].value
        y+=gaussian(x,0,amp,cen,wid)
    return y

def DefineFitModel(name):
    #%%q dependent Models
    if name=='twofreeLorentz':
        Lorentz1=Model(Lorentz_convol,prefix='L1')
        Lorentz2=Model(Lorentz_convol,prefix='L2')
        FitModel=Lorentz1+Lorentz2
        FitLevel='q'
        params = FitModel.make_params()
        params.add('A0',min=0,value=0.5)
        params.add('beta',min=0,value=1)
        params.add('L1scaling',expr='beta*A0')
        params.add('L2scaling',expr='beta*(1-A0)')
        params.add('L1width',min=0,value=1)
        params.add('Gamma',min=0,value=10)
        params.add('L2width',expr='Gamma+L1width')
    elif name=='twofreeLorentz_R':
        Lorentz1=Model(Lorentz_convol,prefix='L1')
        Lorentz2=Model(Lorentz_convol,prefix='L2')
        Lorentz3=Model(Lorentz_convol,prefix='ec')
        FitModel=Lorentz1+Lorentz2+Lorentz3
        FitLevel='q'
        params = FitModel.make_params()
        params.add('A0',min=0,value=0.5)
        params.add('beta',min=0,value=1)
        params.add('L1scaling',expr='beta*A0')
        params.add('ecscaling',min=0,value=1)
        params.add('L2scaling',expr='beta*(1-A0)')
        params.add('L1width',min=0,value=1)
        params.add('ecwidth',value=0,vary=False)
        params.add('Gamma',min=0,value=10)
        params.add('L2width',expr='Gamma+L1width')
        #%% sample models
    elif name=='browndiff_jump':
        Q=np.array(Calibration[0][0]['q'])
        if qrange!=[]:
            builtins.q=Q[(Q>qrange[0])&(Q<qrange[1])]
        else:
            builtins.q=Q
        lambdastr='lambda x'
        paramstr=':Lorentz_convol_global(x'
        for hiq,_ in enumerate(q):
            lambdastr+=',scaling%i,width%i'%(hiq,hiq)
            paramstr+=',scaling%i=scaling%i,width%i=width%i'%(hiq,hiq,hiq,hiq)
        #LCG=eval(lambdastr+paramstr+')')
        FitLevel='S'
        params=Parameters()
        params.add('D',value=0.1,min=0)
        params.add('tauint',value=0.01,min=0)
        params.add('dint',value=1,min=0)
        for hiq,qv in enumerate(q):
            params.add('globwidth%i'%hiq,expr='D*%f**2'%qv)
            params.add('intwidth%i'%hiq,expr='D*%f**2+dint*%f**2/(1+tauint*dint*%f**2)'%(qv,qv,qv))
            params.add('EISF%i'%hiq,value=0.5,min=0,max=1)
            params.add('beta%i'%hiq,value=2)
            params.add('intscaling%i'%hiq,expr='beta%i*(1-EISF%i)'%(hiq,hiq))
            params.add('globscaling%i'%hiq,expr='beta%i*EISF%i'%(hiq,hiq))
        Lorentzglob=Model(eval(lambdastr+paramstr+')'),prefix='glob')
        Lorentzint=Model(eval(lambdastr+paramstr+')'),prefix='int')
        FitModel=Lorentzglob+Lorentzint
    elif name=='browndiff_free':
        Q=np.array(Calibration[0][0]['q'])
        if qrange!=[]:
            builtins.q=Q[(Q>qrange[0])&(Q<qrange[1])]
        else:
            builtins.q=Q
        lambdastr='lambda x'
        paramstr=':Lorentz_convol_global(x'
        for hiq,_ in enumerate(q):
            lambdastr+=',scaling%i,width%i'%(hiq,hiq)
            paramstr+=',scaling%i=scaling%i,width%i=width%i'%(hiq,hiq,hiq,hiq)
        LCG=eval(lambdastr+paramstr+')')
        FitLevel='S'
        params=Parameters()
        params.add('D',value=0.1,min=0)
        for hiq,qv in enumerate(q):
            params.add('globwidth%i'%hiq,expr='D*%f**2'%qv)
            params.add('intwidth%i'%hiq,value=10,min=0)
            params.add('EISF%i'%hiq,value=0.5,min=0,max=1)
            params.add('beta%i'%hiq,value=2)
            params.add('intscaling%i'%hiq,expr='beta%i*(1-EISF%i)'%(hiq,hiq))
            params.add('globscaling%i'%hiq,expr='beta%i*EISF%i'%(hiq,hiq))
        Lorentzglob=Model(LCG,prefix='glob')
        Lorentzint=Model(LCG,prefix='int')
        FitModel=Lorentzglob+Lorentzint
    elif name=='browndiff_free+R':
        Q=np.array(Calibration[0][0]['q'])
        if qrange!=[]:
            builtins.q=Q[(Q>qrange[0])&(Q<qrange[1])]
        else:
            builtins.q=Q
        lambdastr='lambda x'
        paramstr=':Lorentz_convol_global(x'
        for hiq,_ in enumerate(q):
            lambdastr+=',scaling%i,width%i'%(hiq,hiq)
            paramstr+=',scaling%i=scaling%i,width%i=width%i'%(hiq,hiq,hiq,hiq)
        LCG=eval(lambdastr+paramstr+')')
        FitLevel='S'
        params=Parameters()
        params.add('D',value=5,min=0)
        for hiq,qv in enumerate(q):
            params.add('globwidth%i'%hiq,expr='D*%f**2'%qv)
            params.add('intwidth%i'%hiq,value=10,min=0)
            params.add('ecwidth%i'%hiq,value=0,vary=False)
            params.add('ecscaling%i'%hiq,value=0.5,min=0)
            params.add('EISF%i'%hiq,value=1,min=0,max=1)
            params.add('beta%i'%hiq,value=2)
            params.add('intscaling%i'%hiq,expr='beta%i*(1-EISF%i)'%(hiq,hiq))
            params.add('globscaling%i'%hiq,expr='beta%i*EISF%i'%(hiq,hiq))
        Lorentzglob=Model(LCG,prefix='glob')
        Lorentzint=Model(LCG,prefix='int')
        LorentzEC=Model(LCG,prefix='ec')
        FitModel=Lorentzglob+Lorentzint+LorentzEC
    elif name=='a+browndiff_jump':
        Q=np.array(Calibration[0][0]['q'])
        if qrange!=[]:
            builtins.q=Q[(Q>qrange[0])&(Q<qrange[1])]
        else:
            builtins.q=Q
        lambdastr='lambda x'
        paramstr=':Lorentz_convol_global(x'
        for hiq,_ in enumerate(q):
            lambdastr+=',scaling%i,width%i'%(hiq,hiq)
            paramstr+=',scaling%i=scaling%i,width%i=width%i'%(hiq,hiq,hiq,hiq)
        LCG=eval(lambdastr+paramstr+')')
        FitLevel='S'
        params=Parameters()
        params.add('D',value=0.1,min=0)
        params.add('a',value=0.0)
        params.add('tauint',value=0.01,min=0)
        params.add('dint',value=1,min=0)
        for hiq,qv in enumerate(q):
            params.add('globwidth%i'%hiq,expr='a+D*%f**2'%qv)
            params.add('intwidth%i'%hiq,expr='a+D*%f**2+dint*%f**2/(1+tauint*dint*%f**2)'%(qv,qv,qv))
            params.add('EISF%i'%hiq,value=0.5,min=0,max=1)
            params.add('beta%i'%hiq,value=2)
            params.add('intscaling%i'%hiq,expr='beta%i*(1-EISF%i)'%(hiq,hiq))
            params.add('globscaling%i'%hiq,expr='beta%i*EISF%i'%(hiq,hiq))
        Lorentzglob=Model(LCG,prefix='glob')
        Lorentzint=Model(LCG,prefix='int')
        FitModel=Lorentzglob+Lorentzint
    elif name=='browndiff_jump_ecfree':
        Q=np.array(Calibration[0][0]['q'])
        if qrange!=[]:
            builtins.q=Q[(Q>qrange[0])&(Q<qrange[1])]
        else:
            builtins.q=Q
        lambdastr='lambda x'
        paramstr=':Lorentz_convol_global(x'
        for hiq,_ in enumerate(q):
            lambdastr+=',scaling%i,width%i'%(hiq,hiq)
            paramstr+=',scaling%i=scaling%i,width%i=width%i'%(hiq,hiq,hiq,hiq)
        LCG=eval(lambdastr+paramstr+')')
        FitLevel='S'
        params=Parameters()
        params.add('D',value=0.1,min=0)
        params.add('tauint',value=0.01,min=0)
        params.add('dint',value=1,min=0)
        for hiq,qv in enumerate(q):
            params.add('globwidth%i'%hiq,expr='D*%f**2'%qv)
            params.add('intwidth%i'%hiq,expr='D*%f**2+dint*%f**2/(1+tauint*dint*%f**2)'%(qv,qv,qv))
            params.add('EISF%i'%hiq,value=0.5,min=0,max=1)
            params.add('beta%i'%hiq,value=2)
            params.add('ecwidth%i'%hiq,value=0,vary=False)
            params.add('ecscaling%i'%hiq,value=0,min=0)
            params.add('intscaling%i'%hiq,expr='beta%i*(1-EISF%i)'%(hiq,hiq))
            params.add('globscaling%i'%hiq,expr='beta%i*EISF%i'%(hiq,hiq))
        Lorentzglob=Model(LCG,prefix='glob')
        Lorentzint=Model(LCG,prefix='int')
        EC=Model(LCG,prefix='ec')
        FitModel=Lorentzglob+Lorentzint+EC
    return params, FitModel,FitLevel
