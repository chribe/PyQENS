#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 11:40:29 2020

@author: beck
"""
import scipy as sci
import numpy as np


def L(gamma,x):
    #%%% convoluted Lorentz
    # write a Lorentzian in the model as a sum of several voigt functions and lorentian functions with parameters from Vanadium fit
    values = 0.0*x
    if 'NumberOfGaussians' in VanadiumFit.ProtEntry['VanadiumProperties']:
        for hiG in np.arange(1,VanadiumFit.ProtEntry['VanadiumProperties']['NumberOfGaussians']):
            CurScal=VanadiumFit.Data[-1].Fits[-1].Result[currq].params['Gscaling%i'%hiG]
            CurCenter=VanadiumFit.Data[-1].Fits[-1].Result[currq].params['Gcenter%i'%hiG]
            CurSigma=VanadiumFit.Data[-1].Fits[-1].Result[currq].params['Gsigma%i'%hiG]
            try:
                values += CurScal*np.real(sci.special.wofz((x -CurCenter + 1j*gamma) / CurSigma / np.sqrt(2.))) / CurSigma / np.sqrt( 2.*np.pi)
            except:
                print(x.shape)
                print("shouldn't be here....")
    if 'NumberOfLorentz' in VanadiumFit.ProtEntry['VanadiumProperties']:
        for hiL in np.arange(1,VanadiumFit.ProtEntry['VanadiumProperties']['NumberOfLorentz']): # TODO : check convolution (scaling parameters?)
            CurScal=VanadiumFit.Data[-1].Fits[-1].Result[currq].params['Lscaling%i'%hiL]
            CurCenter=VanadiumFit.Data[-1].Fits[-1].Result[currq].params['Lcenter%i'%hiL]
            HWHM=VanadiumFit.Data[-1].Fits[-1].Result[currq].params['LHWHM%i'%hiL]+gamma
            values += CurScal/np.pi*HWHM/((x-CurCenter)**2+(HWHM)**2)
    if 'NumberVoigt' in VanadiumFit.ProtEntry['VanadiumProperties']:
        for hiG in np.arange(1,VanadiumFit.ProtEntry['VanadiumProperties']['NumberVoigt']):# TODO : check convolution (scaling parameters?)
            CurScal=VanadiumFit.Data[-1].Fits[-1].Result[currq].params['Voigtscaling%i'%hiG]
            CurCenter=VanadiumFit.Data[-1].Fits[-1].Result[currq].params['Voigtcenter%i'%hiG]
            CurSigma=VanadiumFit.Data[-1].Fits[-1].Result[currq].params['VoigtSigma%i'%hiG]
            Curgamma=VanadiumFit.Data[-1].Fits[-1].Result[currq].params['Voigtgamma%i'%hiG]+gamma
            values += CurScal*np.real( sci.special.wofz( (x -CurCenter + 1j*Curgamma) / CurSigma / np.sqrt(2.) ) ) / CurSigma / np.sqrt( 2.*np.pi )
    return values


def G(sigma,x):
    #%%% convoluted Gauss
    values = 0.0*x
    if 'NumberOfGaussians' in VanadiumFit.ProtEntry['VanadiumProperties']:
        for hiG in np.arange(1,VanadiumFit.ProtEntry['VanadiumProperties']['NumberOfGaussians']):
            CurScal=VanadiumFit.Data[-1].Fits[-1].Result[currq].params['Gscaling%i'%hiG]
            CurCenter=VanadiumFit.Data[-1].Fits[-1].Result[currq].params['Gcenter%i'%hiG]
            CurSigma=VanadiumFit.Data[-1].Fits[-1].Result[currq].params['Gsigma%i'%hiG]+sigma
            values+=CurScal/CurSigma*(2*np.pi)**(-0.5)*np.exp(-0.5*(x-CurCenter)**2/CurSigma**2)
    if 'NumberOfLorentz' in VanadiumFit.ProtEntry['VanadiumProperties']:
        for hiL in np.arange(1,VanadiumFit.ProtEntry['VanadiumProperties']['NumberOfLorentz']):
            CurScal=VanadiumFit.Data[-1].Fits[-1].Result[currq].params['Lscaling%i'%hiL]
            CurCenter=VanadiumFit.Data[-1].Fits[-1].Result[currq].params['Lcenter%i'%hiL]
            HWHM=VanadiumFit.Data[-1].Fits[-1].Result[currq].params['LHWHM%i'%hiL]
            values+=CurScal*np.real( wofz( (x -CurCenter + 1j*HWHM) / sigma / np.sqrt(2.) ) ) / sigma / np.sqrt( 2.*np.pi )
    if 'NumberVoigt' in VanadiumFit.ProtEntry['VanadiumProperties']:
        for hiG in np.arange(1,VanadiumFit.ProtEntry['VanadiumProperties']['NumberVoigt']):
            CurScal=VanadiumFit.Data[-1].Fits[-1].Result[currq].params['Voigtscaling%i'%hiG]
            CurCenter=VanadiumFit.Data[-1].Fits[-1].Result[currq].params['Voigtcenter%i'%hiG]
            CurSigma=VanadiumFit.Data[-1].Fits[-1].Result[currq].params['VoigtSigma%i'%hiG]+sigma
            Curgamma=VanadiumFit.Data[-1].Fits[-1].Result[currq].params['Voigtgamma%i'%hiG]
            values+=CurScal*np.real( wofz( (x -CurCenter + 1j*Curgamma) / CurSigma / np.sqrt(2.) ) ) / CurSigma / np.sqrt( 2.*np.pi )
    return values

def D(x):
    #%%% convoluted Dirac
    return G(0,x)