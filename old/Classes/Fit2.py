#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 17:52:57 2020

@author: beck
"""
class Fit:
    """ 
    this is a class which contains the fits, this class is normally part of
    the class data or experiment. If it is part of "Data", the fits performed are either "q-dependent" or take also some "q"-dependence into account.
    If the class is part of Experiment, global fits, taking control-parameter independent parameters into account.
     
     """
    from lmfit import Parameters, minimize, report_fit
    
    def __init__(self, FitName='', ProteinModel = '',Model = '', ConvolutedModel='', Parameters = [], Comments = [],Subfit=[]):
        # initialising all the variables
        self.FitName = FitName
        self.ProteinModel = ProteinModel #string
        self.FitModel = FitModel #string
        self.ConvolutedModel = ConvolutedModel  #string
        self.Parameters = Parameters()
        self.Comments = Comments
        self.Subfit=Subfit
    #%%static methods
    @staticmethod
    #%% FIT of D2O with fixed linewidth
    def fit_QENS(data,FitName):
        if FitName=="Vanadium":
           Model='b'
           Fit.Parameters=Parameters()
           for hi=
        if FitName=="D2O":
           Fit.Model='betaD2O*L(gammaD2O,omega)'
        # if FitName=="twolorentz_free_q":
        #     Fit.ProteinModel='beta*(A*L(gamma,omega)+(1-A)*L(gamma+Gamma,omega))',
        #     Fit.Parameters=:['beta','A','gamma','Gamma'],
        #                 'Limits_q_min':[0,0,0,0],
        #                 'Limits_q_max':[Inf,1,Inf,Inf],
        #                 'Startval_q':[1,0.5,5,10],
        #     Fit.Parameters_Units_q' =['','','mueV','mueV']})
        # #%% fixed Dq^2, free internal
        # if FitName=="borwnian_diffision_internal_free":
        #     Fitsmodels.append({'ProteinModel':'beta*(A*L(D*q^2,omega)+(1-A)*L(D*q^2+Gamma,omega))',
        #                 'Parameters_Sample':['D'],
        #                 'Limits_Sample_min':[0],
        #                 'Limits_Sample_max':[Inf],
        #                 'Startval_Sample':[5],
        #                 'Units_Sample':['AA^2mueV'],
        #                 'Parameters_q':['beta','A','Gamma'],
        #                 'Limits_q_min':[0,0,0],
        #                 'Limits_q_max':[Inf,1,Inf],
        #                 'Startval_q':[1,0.5,10],
        #                 'Units_q' :['','','mueV']})   
        #%% FIT fixed Dq^2, jump diffusion for internal diffusion
        if FitName=="borwnian_diffision_internal_jump":
            Fit.ProteinModel='beta*(A*L(D*q^2,omega)+(1-A)*L(D*q^2+Dint*q^2/(tauint*Dint*q^2+1),omega))'
            # define parameters and limits
            Fit.Parameters.add('D' % (iy+1), value=5, min=0.01)# in AA^2mueV
            Fit.Parameters.add('Dint' % (iy+1), value=5, min=0.01)#in AA^2mueV
            Fit.Parameters.add('tauint' % (iy+1), value=5, min=0.01)#in 1/mueV
            for hiq,_ in enumerate(data.q):
                Fit.Parameters.add('beta_%i'%hiq, value=1, min=0.0, max=200)
                Fit.Parameters.add('A_%i'%hiq, value=0.5, min=0.0, max=1)
            #treat solvent contribution
            Fit.Model,Params=TreatSolventContribution(Fit.ProteinModel,Params,data)
            Fit.Parameters+=Params
            #convolute mofel
            Fit.ConvolutedModel=ConvoluteModel(Fit.Model)
            
            
        # #%% fixed jump diffusion, jump diffusion for internal diffusion
        # if FitName=="jump_diffision_internal_jump":
        #    Fitsmodels.append({'ProteinModel':'beta*(A*L(D*q^2/(tau*D*q^2+1),omega)+(1-A)*L(D*q^2/(tau*D*q^2+1)+Dint*q^2/(tauint*Dint*q^2+1),omega))',
        #                 'Parameters_Sample':['D','tau','Dint','tauint'],
        #                 'Limits_Sample_min':[0,0,0,0],
        #                 'Limits_Sample_max':[Inf,Inf,Inf,Inf],
        #                 'Startval_Sample':[5,0.1,10,0.1],
        #                 'Units_Sample':['AA^2mueV','1/mueV','AA^2mueV','1/mueV'],
        #                 'Parameters_q':['beta','A'],
        #                 'Limits_q_min':[0,0],
        #                 'Limits_q_max':[Inf,1],
        #                 'Startval_q':[1,0.5],
        #                 'Units_q' :['','']})
        # #%% immobile fraction fixed jump diffusion, jump diffusion for internal diffusion
        # if FitName=="immobile_jump_diffision_internal_jump":
        #    Fitsmodels.append({'ProteinModel':'beta*(Imm*Delta(q,omega)+(1-Imm))*(A*L(D*q^2/(tau*D*q^2+1),omega)+(1-A)*L(D*q^2/(tau*D*q^2+1)+Dint*q^2/(tauint*Dint*q^2+1),omega))',
        #                 'Parameters_Sample':['Imm','D','tau','Dint','tauint'],
        #                 'Limits_Sample_min':[0,0,0,0,0],
        #                 'Limits_Sample_max':[1,Inf,Inf,Inf,Inf],
        #                 'Startval_Sample':[0.5,5,0.1,10,0.1],
        #                 'Units_Sample':['','AA^2mueV','1/mueV','AA^2mueV','1/mueV'],
        #                 'Parameters_q':['beta','A'],
        #                 'Limits_q_min':[0,0],
        #                 'Limits_q_max':[Inf,1],
        #                 'Startval_q':[1,0.5],
        #                 'Units_q' :['','']})
    return Fit
    #%% Treat Solvent contribution
    def TreatSolventContribution(Fitmodel):
        if =='subtractsolvent':
            Model=Fitmodel
            params=Parameters()
        elif =='modelsolvent':
            Model=Fitmodel+'SFSolvent*L(gammasolvent,omega)'
            params=Parameters()
            for hiq,_ in enumerate(D2OFit):
        return Model,params

