# -*- coding: utf-8
"""
This will be the main skript for the analysis for the QENS analysis

Code written by Christian Beck
christian.beck@uni-tuebingen.de
ORCID: 0000-0001-7214-3447
"""
import PyQENSLib as PQL
import PyQENSAnalysis as PQA
import pickle
import os
import numpy as np
import builtins
import matplotlib.pyplot as plt
from lmfit.models import ExpressionModel
#%% flags
BackgroundFileName='EC_QENS_q.nxs'
CalibrationFileName='Vanadium_QENS_q.nxs'
autoupdate=0
mergeQENSFWS=1
#builtins are accessible within other functions
builtins.backgroundtreatment='subtract'#fit
builtins.backgrounduseage='data'#model
builtins.backgroundscaling=1#'Paalman-Pings'# scalar,  list of scalars
builtins.NumberOfGaussians=2
#builtins.plot_fitresults=0
builtins.hbar=0.6582119569
#builtins.save_fitresult='pdf'#eps,png and combinations of these
#%% specify paths
experiment_folder='/home/christian/Documents/ILLFS/Documents/Copy9-13-876_new/'
data_folder=experiment_folder + 'processed/reduced_data/'
builtins.figure_folder=experiment_folder + 'processed/fit_results/'
if not os.path.isdir(figure_folder):
    os.mkdir(figure_folder)
#%%
Data=PQL.load_data(data_folder)
builtins.Calibration=PQL.FindEntry(Data,CalibrationFileName)
builtins.Background=PQL.FindEntry(Data,BackgroundFileName)
Data=PQL.IntegrateSQW(Data,suffix='aftersubtraction')
Data=PQL.MergeData(Data,mergeQENSFWS)
#Data=PQL.BackgroundTreatment(Data,exclude=['Vanadiumrod_QENS_q.nxs'])
Data=PQA.QENS_Fit(Data,Type='Calibration')
#%% subtract EC from water
# 0bar
# =============================================================================
# water=Data[17]
# plt.figure()
# fig,ax=plt.subplots(nrows=4,ncols=5)
# axs=PQL.trim_axs(ax,20)
# for hiq,q in enumerate(water[0]['q']):
#     axs[hiq].errorbar(water[0]['hw'],water[0]['sqw'],water[0]['dsqw'])
#     axs[hiq].errorbar(Background[0]['hw'],Background[0]['sqw'],Background[0]['dsqw'])
# =============================================================================
#%%
for hid,dat in enumerate(Data):
    if 'BSA' in dat[0]['Name'] and 'QENS' in dat[0]['Name']:
        if '1000bar' in dat[0]['Name']:
            pressure='1000bar'
            dat[0].update({'pressure':1000})
        elif '1500bar' in dat[0]['Name']:
            pressure='1500bar'
            dat[0].update({'pressure':1500})
        elif '2000bar' in dat[0]['Name']:
            pressure='2000bar'
            dat[0].update({'pressure':2000})
        elif '2500bar' in dat[0]['Name']:
            pressure='2500bar'
            dat[0].update({'pressure':2500})
        elif '3000bar' in dat[0]['Name']:
            pressure='3000bar'
            dat[0].update({'pressure':3000})
        elif '500bar' in dat[0]['Name']:
            pressure='500bar'
            dat[0].update({'pressure':500})
        else:
            pressure='0bar'
            dat[0].update({'pressure':0})
        datsub=[]
        for hid1,dat1 in enumerate(Data):
            if 'BSA' not in dat1[0]['Name'] and pressure in dat1[0]['Name'] :
                datsub=dat1
        if datsub==[]:
            raise NameError(dat[0]['Name'] + ': No Corresponding D2O found')
        phi_dry,dphi_dry,phi_hydr,dphi_hydr,n,dn=PQL.VolumeFraction(Protein='BSA',cpreal=200,T=295)
        scalingparameter=1-phi_dry
        dat=PQL.SubtractSolvent(dat,datsub,scalingparameter)#
# =============================================================================
# #%%
# plt.figure()
# leg=()
# for _,Dat in enumerate(Data):
#     for _,dat in enumerate(Dat):
#         if dat['Type']=='QENS' and 'BSA' not in dat['Name'] and 'D2O' in dat['Name'] and 'D2O_2000bar_QENS_q' not in dat['Name'] and 'D2O_2500bar_QENS_q' not in dat['Name'] or 'EC_QENS' in dat['Name']:
#             leg=leg+(dat['Name'],)
#             plt.plot(dat['q'],np.array(dat['integrSqaftersubtraction'])/np.array(Calibration[0][0]['integrSqaftersubtraction']))
# plt.legend(leg)
# plt.xlabel(r'q[$\AA^{-1}$]')
# plt.ylabel(r'$\sum_{\hbar\omega} S(q,\omega)$')
# =============================================================================
#%%
tobedeleted=[]
for hid,dat in enumerate(Data):
    if 'BSA' not in dat[0]['Name'] or 'failior' in dat[0]['Name'] or 'FWS' in dat[0]['Name'] \
        or 'defekt' in dat[0]['Name']:
        tobedeleted.append(hid)
for hid in range(len(tobedeleted),0,-1):
    del Data[tobedeleted[hid-1]]
#%%
# =============================================================================
# plt.figure()
# for _,Dat in enumerate(Data):
#     for _,dat in enumerate(Dat):
#         if dat['Type']=='QENS':
#             plt.plot(dat['q'],np.array(dat['integrSqaftersubtraction'])/np.array(Calibration[0][0]['integrSqaftersubtraction']))
# =============================================================================
#%%
Data=PQA.QENS_Fit(Data,Type='twofreeLorentz')
#Data=PQA.QENS_Fit(Data,Type='a+browndiff_jump',qrange=[0.3,1.5])
#Data=PQA.QENS_Fit(Data,Type='browndiff_free',qrange=[0.3,1.5])
#Data=PQA.QENS_Fit(Data,Type='browndiff_free+R',qrange=[0.3,1.5])
Data=PQA.QENS_Fit(Data,Type='browndiff_jump',qrange=[0.3,1.5])
#%%
for _,Dat in enumerate(Data):
    for _,DAT in enumerate(Dat):
        gamma=[]
        dgamma=[]
        A0=[]
        dA0=[]
        Gamma=[]
        dGamma=[]
        FIT=DAT['Fit'][0]
        for hiq,_ in enumerate(DAT['q']):
            gamma.append(FIT[hiq].params['L1width'].value/ 0.6582119569)
            dgamma.append(FIT[hiq].params['L1width'].stderr/ 0.6582119569)
            A0.append(FIT[hiq].params['A0'].value)
            dA0.append(FIT[hiq].params['A0'].stderr)
            Gamma.append(FIT[hiq].params['Gamma'].value/ 0.6582119569)
            dGamma.append(FIT[hiq].params['Gamma'].stderr/ 0.6582119569)
        gamma=np.array(gamma[2:])
        dgamma=np.array(dgamma[2:])
        Gamma=np.array(Gamma[2:])
        dGamma=np.array(dGamma[2:])
        A0=np.array(A0)
        dA0=np.array(dA0)
        EISFFIT = ExpressionModel("p+(1-p)*exp(-(a*x)**2/5)*(s*(3/x/R*(sin(x*R)/x**2/R**2-cos(x*R)/x/R))**2+(1-s)/3*(1+2*sin(x*am)/x/am))")
        DAT.update({'EISF':EISFFIT.fit(A0,x=DAT['q'],weights=1/dA0**2,p=1,a=1,s=1,R=10,am=5)})
        gammaFIT = ExpressionModel("a+D * x")
        DAT.update({'gamma':gammaFIT.fit(gamma,x=DAT['q'][2:]**2,weights=1/dgamma**2,D=1,a=0)})
        GammaFIT = ExpressionModel("D/ 0.6582119569 * x/(1+D/ 0.6582119569 * x*tau)")
        DAT.update({'Gamma':GammaFIT.fit(Gamma,x=DAT['q'][2:]**2,weights=1/dGamma**2,D=1,tau=0.1)})
        DAT['EISF'].plot()
        plt.xlim(left=0)
        plt.ylim(bottom=0,top=1)
        plt.savefig(figure_folder + DAT['Name'] + '/EISF_fit.eps')
        plt.close()
        DAT['gamma'].plot()
        plt.xlim(left=0)
        plt.ylim(bottom=0)
        plt.savefig(figure_folder + DAT['Name'] + '/gamma_fit.eps')
        plt.close()
        DAT['Gamma'].plot()
        plt.xlim(left=0)
        plt.ylim(bottom=0)
        plt.savefig(figure_folder + DAT['Name'] + '/Gamma_fit.eps')
        plt.close()
#%%
Dglob=[]
dDglob=[]    
D=[]
dD=[]
salt=[]
pressure=[]
for  _,Dat in enumerate(Data):
    for _,DAT in enumerate(Dat):
        if 'YCl3' in DAT['Name']:
            salt.append(18)
        else:
            salt.append(0)
        D.append(DAT['gamma'].params['D'].value)
        dD.append(DAT['gamma'].params['D'].stderr)
        pressure.append(DAT['pressure'])
        Dglob.append(DAT['Fit'][1].params['D'].value)
        dDglob.append(DAT['Fit'][1].params['D'].stderr)
plt.figure()
salt=np.array(salt)
D=np.array(D)
dD=np.array(dD)
Dglob=np.array(Dglob)/hbar
dDglob=np.array(dDglob)/hbar
pressure=np.array(pressure)
for hip in [0,500,1500,2000,2500,3000]:
    plt.errorbar(salt[hip==pressure],D[hip==pressure],dD[hip==pressure],label=str(hip),ls='None',marker='o')
plt.legend()
plt.figure()
for hiS in [0,18]:
    plt.errorbar(pressure[hiS==salt],D[hiS==salt],dD[hiS==salt],label=str(hiS)+ ' mM',ls='None',marker='o')
    plt.errorbar(pressure[hiS==salt],Dglob[hiS==salt],dDglob[hiS==salt],label=str(hiS)+ ' mM',ls='None',marker='o')
plt.legend()
plt.xlabel('pressure [bar]')
plt.ylabel(r'D[$\mathrm{\AA}^2$/ns]')
plt.savefig(figure_folder  + 'pressuredependence.eps')
plt.figure()
plt.errorbar(pressure[salt==0],D[salt==18]/D[salt==0],np.power(np.power(dD[salt==18]/D[salt==0],2)+np.power(D[salt==18]/D[salt==0]/D[salt==0]*dD[salt==0],2),0.5),label=str(hiS)+ ' mM',ls='None',marker='o')
plt.errorbar(pressure[salt==0],Dglob[salt==18]/Dglob[salt==0],np.power(np.power(dDglob[salt==18]/Dglob[salt==0],2)+np.power(Dglob[salt==18]/Dglob[salt==0]/Dglob[salt==0]*dDglob[salt==0],2),0.5),label=str(hiS)+ ' mM',ls='None',marker='o')
plt.legend()
plt.xlabel('pressure [bar]')
plt.ylabel(r'D/D(c_s=0)')
plt.savefig(figure_folder  + 'normalizedpressuredependence.eps')

#%%
#PQA.FWS_MSD(Data,exclude=[BackgroundFileName,CalibrationFileName])
#PQA.FWS_Ratio(Data,exclude=[BackgroundFileName,CalibrationFileName])
#%% 
#PQA.QENS_Fit(Data,exclude=[BackgroundFileName,CalibrationFileName])
