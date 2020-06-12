#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 17:52:57 2020
Version 06 Mai 2020 19:15
@author: beck
"""

from lmfit import Parameters, minimize, report_fit
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
import time
import importlib
import builtins
from tqdm import tqdm
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{upgreek}']
class Fit:
    """
    this is a class which contains the fits, this class is normally part of
    the class data or experiment. If it is part of "Data", the fits performed
    are either "q-dependent" or take also some "q"-dependence into account.
    If the class is part of Experiment, global fits, taking control-parameter
    independent parameters into account.
     
     """
    # TODO : hand over index of sample for saving the figures for q and S fits "Sampleidx"
    # TODO : write fit on Experimental level
    # TODO : define global flag "solventtreatment"
    # TODO : finish solvent treatment
    # TODO : finish fit of D2O
     
    def __init__(self, FitName='', ProteinModel='', Model='',
                 ConvolutedModel='', Parameter=[], Comments=[], Subfit=[]):
        # initialising all the variables
        self.FitName = FitName
        self.ProteinModel = ProteinModel  # string
        self.Parameter = Parameters()
        self.Comments = Comments
        self.Subfit = Subfit
        self.Parameterproperties = {}
        self.Model = ''
        self.Result = []
    # %% Function fit_QENS
    def fit_QENS(self,data,FitName,Vanadium):
        builtins.VanadiumFit=Vanadium
        VanadiumData=Vanadium.ProtEntry['VanadiumProperties']
        import fitfunctions as ff
        #%%% get fitfunction
        self=ff.getfitfunction(self,FitName,VanadiumData)
    #%%% prepare fit
    #%%%% obtain Fit level, which is relevant for further fit
        if 'E' in self.Parameterproperties['Level']:
            self.Level='E'
        elif 'S' in self.Parameterproperties['Level']:
            self.Level='S'
        else:
            self.Level='q'
    #%%%% generate list of parameters
        self.Parameter=Parameters()
        if self.Level=='E':
            raise NameError('Fits on experimental level have to be implemented!')
        elif self.Level=='S':#=================SAMPLE FITS=================================================================
            # defines parameters for Sample wise Fits
            for hi,_ in enumerate(self.Parameterproperties['Name']):
                if self.Parameterproperties['Level'][hi]=='S': # "q" independent parameters
                    self.Parameter.add(self.Parameterproperties['Name'][hi],
                                       value=self.Parameterproperties['Strtvl'][hi],
                                       min=self.Parameterproperties['Lower'][hi],
                                       max=self.Parameterproperties['Upper'][hi])
                else:# q dependent parameters
                    for hiq,_ in enumerate(data.q):
                        self.Parameter.add(self.Parameterproperties['Name'][hi]+'_%i'%hiq,
                                       value=self.Parameterproperties['Strtvl'][hi],
                                       min=self.Parameterproperties['Lower'][hi],
                                       max=self.Parameterproperties['Upper'][hi])
        else:#=========================== Q DEPENDENT FITS=======================================================
            for hi,_ in enumerate(self.Parameterproperties['Name']):# add Parameters
                self.Parameter.add(self.Parameterproperties['Name'][hi],
                                    value=self.Parameterproperties['Strtvl'][hi],
                                    min=self.Parameterproperties['Lower'][hi],
                                    max=self.Parameterproperties['Upper'][hi])
        #================================PLOTTING FLAGS==================================================
        #Separators used for plotting later
        for Sep in self.Plotutilities['Separators']:
            self.Parameter.add(Sep[0], value=1, vary=False)
        #%%%%treat solvent contribution
        self.TreatSolventContribution(data)
        #%%%%create fit function 
        self.createfitfunction()
        import Fitfunction_temp as Fitfun
        del Fitfun
        del sys.modules['Fitfunction_temp']
        import Fitfunction_temp as Fitfun
        #%%%% define objective functions
        if self.Level=='S' or self.Level=='E':#define objectove function for Sample or Experiment fit
            def objective(params, x,q, data,errors):
                """Calculate total residual for fits of Lorentzian to several data sets."""
                _, ndata = data.shape
                resid = 0.0*data[:]
                # make residual per data set
                for i in range(ndata):
                    currq=i
# =============================================================================
#                     if len(data[:,i])==len(q):
#                         resid[:,i] = 100000.0*resid[:,i-1]
#                     else:
# =============================================================================
                    resid[:,i] = (data[:,i] - Fitfun.model_dataset(params, i, x,q[i]))/errors[:,i]**2
                # now flatten this to a 1D array, as minimize() needs
                return resid.flatten()
        elif self.Level=='q':# define objective function for q fit
            def objective(params, x, data,errors):
                resid = (data - Fitfun.model_dataset(params, x))/errors**2
                return resid
        #%%% perform fit
        if self.Level=='E':
            raise NameError("Fits for Experiment Level have to be written!")
        elif self.Level=='S':
            self.Result = minimize(objective, self.Parameter, args=(data.hw,data.q,data.sqw,data.dsqw))
            tqdm.write(self.Result.message)
        elif self.Level=="q":
            with tqdm(total=len(data.q), file=sys.stdout,desc='Fit on "q" level',position=1) as pbar:
                for hiq , _ in enumerate(data.q):
                    builtins.currq = hiq
                    parameter = self.Parameter
                    self.Result.append(minimize(objective, parameter, 
                                                args=(data.hw, data.sqw[:,hiq],
                                                      data.dsqw[:,hiq])))
                    tqdm.write(self.Result[-1].message)
                    pbar.update(1)
        del objective
        #%%% plot fit Results 
        if flags['plotresults']=='showfits' or flags['plotresults']=='all' or flags['plotresults']=='save' or flags['plotresults']=='showsave':
            if self.Level=='E':
                # TODO write plot routine for Experimental fit
                #=========================SAMPLE FIT PLOTS============================
                raise NameError("to be written")
            elif self.Level=='S':
                builtins.data=data
                for hiq,qval in enumerate(data.q):
                    plt.figure()
                    plt.errorbar(data.hw, data.sqw[:,hiq],yerr=data.dsqw[:,hiq],label='experimental data',zorder=0)
                    zo=0
                    #plot different contributions by setting the plotflags ==0
                    for _, Seplab in enumerate(self.Plotutilities['Separators']):
                        for _, Seplabzero in enumerate(self.Plotutilities['Separators']):
                            self.Result.params[Seplabzero[0]].value = 0
                        self.Result.params[Seplab[0]].value = 1
                        y_fit = Fitfun.model_dataset(self.Result.params, hiq,data.q,data.hw)
                        zo+=5
                        plt.plot(data.hw, y_fit, '-',label=Seplab[1],zorder=zo)
                        #plothandels=+(p1)
                        for _, Seplabzero in enumerate(self.Plotutilities['Separators']):# set all flags back to 1!
                            self.Result.params[Seplabzero[0]].value = 1
                    zo+=5
                    y_fit = Fitfun.model_dataset(self.Result.params, hiq,data.q,data.hw)
                    plt.plot(data.hw, y_fit, '-',label='Total fit',zorder=zo)
# =============================================================================
#                     for _, Seplab in enumerate(self.Plotutilities['Separators']):
#                         for _, Seplabzero in enumerate(self.Plotutilities['Separators']):
#                             self.Result.params[Seplabzero[0]]=0
#                         self.Result.params[Seplab[0]]=1
#                         y_fit = Fitfun.model_dataset(self.Result.params, hiq, hw,q)
#                         plt.plot(hw, y_fit, '-',label=Seplab[1])
#                         for _, Seplabzero in enumerate(self.Plotutilities['Separators']):# set all flags back to 1!
#                             self.Result.params[Seplabzero[0]]=1
# =============================================================================
                    plt.yscale("log")
                    plt.legend()
                    plt.xlabel(r"$\hbar\omega$ [$\upmu$eV]")
                    plt.ylabel(r"$S(q,\omega)$ [arb. units]")
                    plt.title(r"$q=$" + str(qval) +"\AA$^{-1}$")
                    if flags['plotresults']=='save' or flags['plotresults']=='showsave':
                        print(flags['savepath'] + "Fitresult_" + self.FitName + "_S_" + str(Sampleidx) + "_q_" + str(qval).replace('.','p') + "." + flags['pictureextension'])
                        plt.savefig(flags['savepath'] + "Fitresult_" + self.FitName + "_S_" + str(Sampleidx) + "_q_" + str(qval).replace('.','p') + "." + flags['pictureextension'], format = flags['pictureextension'], dpi=1000)
                    if  flags['plotresults']=='showfits' or  flags['plotresults']=='showsave':
                        plt.show()
#=========================SAMPLE FIT PLOTS============================
            elif self.Level=='q':
                for hiq,qval in enumerate(data.q):
                    plt.figure()
                    errors=data.dsqw[:,hiq]
                    values=data.sqw[:,hiq]
                    values[np.isinf(errors)]=0
                    errors[np.isinf(errors)]=0
                    builtins.currq=hiq
                    p1=plt.errorbar(data.hw, values,yerr=errors,label='experimental data',zorder=0)
                    zo=0
                    #plot different contributions by setting the plotflags ==0
                    for _, Seplab in enumerate(self.Plotutilities['Separators']):
                        for _, Seplabzero in enumerate(self.Plotutilities['Separators']):
                            self.Result[hiq].params[Seplabzero[0]].value = 0
                        self.Result[hiq].params[Seplab[0]].value = 1
                        y_fit = Fitfun.model_dataset(self.Result[hiq].params,data.hw)
                        zo+=5
                        plt.plot(data.hw, y_fit, '-',label=Seplab[1],zorder=zo)
                        #plothandels=+(p1)
                        for _, Seplabzero in enumerate(self.Plotutilities['Separators']):# set all flags back to 1!
                            self.Result[hiq].params[Seplabzero[0]].value = 1
                    zo+=5
                    y_fit = Fitfun.model_dataset(self.Result[hiq].params, data.hw)
                    plt.plot(data.hw, y_fit, '-',label='Total fit',zorder=zo)
                    plt.yscale("log")
                    plt.legend()
                    plt.ylim(bottom=1e-6,top=1)
                    plt.xlabel(r"$\hbar\omega$ [$\upmu$eV]")
                    plt.ylabel(r"$S(q,\omega)$ [arb. units]")
                    plt.title(r"$q=$" + str(qval) +"\AA$^{-1}$")
                    if flags['plotresults']=='save' or flags['plotresults']=='showsave':
                        plt.savefig(flags['savepath'] + "Fitresult_" +
                                    self.FitName + "_S_" + str(Sampleidx) +
                                    "_q_" + str(qval).replace('.','p') +
                                    "." + flags['pictureextension'],
                                    format = flags['pictureextension'], dpi=1000)
                    if flags['plotresults']=='showfits' or flags['plotresults']=='showsave':
                        plt.show()
            plt.close("all")
        #%%%save results in global variables for Vanadium and D2O
# =============================================================================
#         if self.FitName=='Vanadium':
#             del self.Parameter
#             self.Parameter=Parameters()
#             for hiq,Res in enumerate(self.Result):
#                 for _ , parameter in enumerate(Res.params.keys()):
#                     self.Parameter.add(parameter + '_q_' + str(hiq),
#                                           value = Res.params[parameter],
#                                           vary=False)
# =============================================================================
        elif self.FitName=='D2O':
            D2OData.update({'Fit':Fit})
        return Fit
    #%% Function TreatSolventContribution
    def TreatSolventContribution(self,data):
        if self.FitName=='SolventOneLorentz':
            self.Model='SFSolvent*L(gammasolvent,omega)'
            self.Parameter.add('SFSolvent',value=1,min=0)
            self.Parameter.add('gammasolvent',value=10,min=0)
            self.Parameterproperties['Name']=['SFSolvent','gammasolvent']
            self.createfitfunction(self)
            import Fitfunction_temp as Fitfun
            for hiq,qval in enumerate(data.q):
                def objective(params, x, data,errors):
                    ndata, _ = data.shape
                    for i in range(ndata):
                        resid = (data - Fitfun.model_dataset(params, x))/errors^2
                    return resid
                for hiq,_ in enumerate(data.q):
                    global currq
                    currq=hiq
                    self.Result[hiq]=minimize(objective, self.Parameter, args=(data.hw,data.sqw[hiq,:],data.dsqw[hiq,:]))
        elif self.FitName=='D2O':
            print(',,,')
        elif self.FitName=='Vanadium':
                self.Model=self.Model
        else:
            if data.ProtEntry['Solvent']['treatment']=='subtract':
                self.Model=self.ProteinModel
            elif data.ProtEntry['Solvent']['treatment']=='modelsolvent':
                self.Model=self.ProteinModel + 'solcont*SFSolvent*L(gammasolvent,omega)'
                self.Parameter.add('solcont', value=1, vary=False)
                self.Parameter.add('SFSolvent', value=1, vary=False)# TODO: change the scaling value!!!
                self.Plotutilities['Separators'].append(['solcont','solvent contribution'])
                for hiq,_ in enumerate(data.q):
                    self.Parameter.add('gammasolvent_%i'%hiq, value=1, vary=False)# TODO: change the width!!!
                self.Constants={'Name':['SFSolvent','gammasolvent'],
                               'Unit':['','mueV'],
                               'Level':['S','q']}
            else:
                self.Model=self.ProteinModel
    #%% Function createfitfunction
    def createfitfunction(self):
        space='    '
        #define fitfunction
        Parameterstring=''
        for param in self.Parameterproperties['Name']:
            Parameterstring+= param + ','
        for param in self.Plotutilities['Separators']:
            Parameterstring+= param[0] + ','
        for param in self.Constants['Name']:
            Parameterstring+= param + ','
        #%%% intelligent code Function fitfunction
        f = open("Fitfunction_temp.py", "w")
        f.write("import numpy as np\n")
        f.write("import convmod as CM\n")
        if self.Level=='q':
            f.write("def fitfunction(" + Parameterstring + "x):\n")
            stringtocallfitfun="fitfunction(" + Parameterstring + "x)"
        else:
            f.write("def fitfunction(" + Parameterstring + "q,x):\n")
            stringtocallfitfun="fitfunction(" + Parameterstring + "q,x)"
        Model=self.Model
        Model=Model.replace('L(','CM.L(')
        Model=Model.replace('G(','CM.G(')
        Model=Model.replace('D(','CM.D(')
        f.write(space + "Sqw=" + Model + "\n")
        f.write(space + "return Sqw\n")
        #%%% intelligent code Function model_dataset
        f.write("\n\n")
        if self.Level=="E":
            raise NameError("Still to be written")# TODO: write intelligent code model data set function Experiment level
        #======================================================================================
        elif self.Level=='S':
            f.write("def model_dataset(params,hiq,q,x):\n")
            f.write(space + "Q=q[hiq]\n")
            for hi,_ in enumerate(self.Parameterproperties['Name']):# assign fit parameters
                if self.Parameterproperties['Level'][hi]=="S":
                    f.write(space + self.Parameterproperties['Name'][hi] + "=params['" + self.Parameterproperties['Name'][hi] + "']\n")
                elif self.Parameterproperties['Level'][hi]=="q":
                    f.write(space + self.Parameterproperties['Name'][hi] + "=params['" + self.Parameterproperties['Name'][hi]+"_%i'%hiq]\n")
                else:
                    raise NameError("Something is wrong here!!!")
            for hi,_ in enumerate(self.Constants['Name']):# assign constants
                if self.Constants['Level']=='S':
                    f.write(space + self.Constants['Name'][hi] + "=params['" + self.Constants['Name'][hi] + "']\n")
                elif self.Constants['Level']=='q':
                    f.write(space + self.Constants['Name'][hi] + "=params['" + self.Constants['Name'][hi]+"'_%i'%hiq]\n")
                else:
                    raise NameError("Something is wrong here!!!")
        #======================================================================================
        elif self.Level=='q':
            f.write("def model_dataset(params,x):\n")
            for hi,_ in enumerate(self.Parameterproperties['Name']):# assign fit parameters
                f.write(space + self.Parameterproperties['Name'][hi] + "=params['" + self.Parameterproperties['Name'][hi] + "']\n")
            for hi,_ in enumerate(self.Constants['Name']):# assign constants
                f.write(space + self.Constants['Name'][hi] + "=params['" + self.Constants['Name'][hi] + "']\n")
        #======================================================================================
        #independent from level: plotutilities
        for hi,_ in enumerate(self.Plotutilities['Separators']):# assign plotparameters
                f.write(space + self.Plotutilities['Separators'][hi][0] + "=params['" + self.Plotutilities['Separators'][hi][0] + "']\n")
        f.write(space + "return " + stringtocallfitfun)
        f.close()
    
