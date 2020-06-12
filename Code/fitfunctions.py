import numpy as np
def getfitfunction(Fit,FitName,VanadiumData):
    if FitName=="Vanadium":
        Fit.FitName=FitName
        if VanadiumData['bkg']=="flat":
            Fit.Model='b'
            Fit.Parameterproperties={'Name':['b'],
                                     'Lower':[0 ],
                                     'Upper':[np.inf],
                                     'Strtvl':[0.00005],
                                     'Level':['q']}
        elif VanadiumData['bkg']=="slope":
            Fit.Model='a*x+b'
            Fit.Parameterproperties={'Name':['a','b'],
                                     'Lower':[-np.inf,0 ],
                                     'Upper':[np.inf,np.inf],
                                     'Strtvl':[0,0.00005],
                                     'Level':['q','q']}
        elif VanadiumData['bkg']=="none":
            Fit.Model='0'
            Fit.Parameterproperties={'Name':[],
                                     'Lower':[ ],
                                     'Upper':[],
                                     'Strtvl':[],
                                     'Level':[]}
        Fit.Plotutilities={'Separators':[]}
        for hiG in np.arange(VanadiumData['NumberOfGaussians']):
            Fit.Model+='+Gauss%i*Gscaling%i/Gsigma%i*(2*np.pi)**(-0.5)*np.exp(-0.5*(x-Gcenter%i)**2/Gsigma%i**2)' % (hiG,hiG,hiG,hiG,hiG)
            Fit.Plotutilities['Separators'].append(['Gauss%i'%hiG,'Gauss %i'%hiG])
            Fit.Parameterproperties['Name']+=['Gscaling%i'%hiG,'Gsigma%i'%hiG,'Gcenter%i'%hiG]
            Fit.Parameterproperties['Lower']+=[0,0,-2]
            Fit.Parameterproperties['Upper']+=[np.inf,np.inf,2]
            Fit.Parameterproperties['Strtvl']+=[1,1,0]
            Fit.Parameterproperties['Level']+=['q','q','q']
        for hiL in np.arange(VanadiumData['NumberOfLorentzians']):
            Fit.Model+='+Lorentz%i*Lscaling%i/np.pi*LHWHM%i/((x-Lcenter%i)**2+(LHWHM%i)**2)' % (hiL,hiL,hiL,hiL,hiL)
            Fit.Plotutilities['Separators'].append(['Lorentz%i'%hiL,':Lorentz %i'%hiL])
            Fit.Parameterproperties['Name']+=['Lscaling%i'%hiL,'LHWHM%i'%hiL,'Lcenter%i'%hiL]
            Fit.Parameterproperties['Lower']+=[0,0,-np.inf]
            Fit.Parameterproperties['Upper']+=[np.inf,np.inf,np.inf]
            Fit.Parameterproperties['Strtvl']+=[1,10,0]
            Fit.Parameterproperties['Level']+=['q','q','q']
        for hiV in np.arange(VanadiumData['NumberOfVoigts']):
            Fit.Model+='+Voigt%i*Voigtscaling%i*np.real( wofz( (x -Voigtcenter%i + 1j*Voigtgamma%i) / VoigtSigma%i/ np.sqrt(2.) ) ) / VoigtSigma%i / np.sqrt( 2.*np.pi )' % (hiV,hiV,hiV,hiV,hiV,hiV)
            Fit.Plotutilities['Separators'].append(['Voigt%i'%hiV,':Voigt %i'%hiV])
            Fit.Parameterproperties['Name']+=['Voigtscaling%i'%hiV,'Voigtgamma%i'%hiV,'VoigtSigma%i'%hiV,'Voigtcenter%i'%hiV]
            Fit.Parameterproperties['Lower']+=[0,0,0,-np.inf]
            Fit.Parameterproperties['Upper']+=[np.inf,np.inf,np.inf,np.inf]
            Fit.Parameterproperties['Strtvl']+=[1,10,10,0]
            Fit.Parameterproperties['Level']+=['q','q','q','q']
        Fit.Constants={'Name':[],
                       'Level':[],
                       'Value':[]}
    #%%%% FIT of D2O with fixed linewidth
    elif FitName=="D2O":
        Fit.FitName=FitName
        Fit.Model='betaD2O*L(gammaD2O,omega)'
        Fit.Parameterproperties={'Name':['betaD2O'],
                                 'Lower':[0 ],
                                 'Upper':[np.inf],
                                 'Strtvl':[1],
                                 'Unit':[''],
                                 'Level':['q']}
        Fit.Plotutilities={'Separators':[['internalcont','internal contribution'],
                                         ['globalcont','global contribution']]}
        Fit.Constants={'Name':[],
                       'Level':[],
                       'Value':[]}
    #%%%% Fit of two free lorentzians
    elif FitName=="twolorentz_free_q":
        Fit.FitName=FitName
        Fit.ProteinModel='beta*(globalcont*A*L(gamma,x)+internalcont*(1-A)*L(gamma+Gamma,x))'
        Fit.Parameterproperties={'Name':['beta',       'A',     'gamma', 'Gamma'],
                                 'Lower':[0 ,        0,          0,           0],
                                 'Upper':[np.inf,    1,         np.inf,  np.inf],
                                 'Strtvl':[100,         1,         5,       15],
                                 'Unit':['','','mueV','mueV'],
                                 'Level':['q','q','q','q','q']}
        Fit.Plotutilities={'Separators':[['internalcont','internal contribution'],
                                         ['globalcont','global contribution']]}
        Fit.Constants={'Name':[],
                       'Level':[],
                       'Value':[]}
    elif FitName=="twolorentz_free_dirac_q":
        Fit.FitName=FitName
        Fit.ProteinModel='beta*dirac*C*D(x)+beta*(1-C)*(globalcont*A*L(gamma,x)+internalcont*(1-A)*L(gamma+Gamma,x))'
        Fit.Parameterproperties={'Name':['beta',  'C',     'A',     'gamma', 'Gamma'],
                                 'Lower':[0 ,      0,  0,          0,           0],
                                 'Upper':[np.inf,   1, 1,         np.inf,  np.inf],
                                 'Strtvl':[100,      0.5,   1,         5,       15],
                                 'Unit':['','','mueV','mueV'],
                                 'Level':['q','q','q','q','q','q']}
        Fit.Plotutilities={'Separators':[['internalcont','internal contribution'],
                                         ['globalcont','global contribution'],
                                         ['dirac','immobile contribution']]}
        Fit.Constants={'Name':[],
                       'Level':[],
                       'Value':[]}
    #%%%% fixed Dq^2, free internal
    elif FitName=="brownian_diffusion_internal_free":
        Fit.FitName=FitName
        Fit.ProteinModel='beta*(globalcont*A*L(D*q**2,omega)+internalcont*(1-A)*L(D*q**2+Gamma,omega))'
        Fit.Parameterproperties={'Name':['D', 'Gamma','beta',       'A'],
                                 'Lower':[0 ,        0,          0,           0],
                                 'Upper':[np.inf,   np.inf,  np.inf,    1],
                                 'Strtvl':[1,       15,         1,       0.5],
                                 'Unit':['\AA^2mueV','mueV','',''],
                                 'Level':['S','q','q','q']}
        Fit.Plotutilities={'Separators':[['internalcont','internal contribution'],
                                         ['globalcont','global contribution']]}
        Fit.Constants={'Name':[],
                       'Level':[],
                       'Value':[]}
    #%%%% fixed jump diffusion, jump diffusion for internal diffusion
    elif FitName=="jump_diffusion_internal_jump":
        Fit.FitName=FitName
        Fit.ProteinModel='beta*(A*L(D*q**2/(tau*D*q**2+1),omega)+(1-A)*L(D*q**2/(tau*D*q**2+1)+Dint*q**2/(tauint*Dint*q**2+1),omega))'
        Fit.Parameterproperties={'Name':['D', 'tau','Dint','tauint','beta','A'],
                                 'Lower':[0 ,   0    , 0,    0,      0,     0],
                                 'Upper':[np.inf, np.inf,np.inf,np.inf,np.inf,    1],
                                 'Strtvl':[1,   0.1, 15,0.1,1,0.5],
                                 'Unit':['\AA^2mueV','m1/mueV','\AA^2mueV','m1/mueV','',''],
                                 'Level':['S','S','S','S','q','q']}
        Fit.Plotutilities={'Separators':[['internalcont','internal contribution'],
                                         ['globalcont','global contribution']]}
        Fit.Constants={'Name':[],
                       'Level':[],
                       'Value':[]}
    #%%%% immobile fraction fixed jump diffusion, jump diffusion for internal diffusion
    elif FitName=="immobile_jump_diffusion_internal_jump":
        Fit.FitName=FitName
        Fit.ProteinModel='beta*((immobile*Imm*D(x)+globalcont*(1-Imm)*(A*L(D*q**2/(tau*D*q**2+1),x)))+internalcont*(Imm*(L(Dint*q**2/(tauint*Dint*q**2+1),x))+(1-Imm)*(1-A)*L(D*q**2/(tau*D*q**2+1)+Dint*q**2/(tauint*Dint*q**2+1),x)))'
        Fit.Parameterproperties={'Name':['Imm','D', 'tau','Dint','tauint','beta','A'],
                                 'Lower':[0 ,   0    , 0,    0,      0,     0,     0],
                                 'Upper':[1, np.inf,np.inf,np.inf,np.inf,np.inf,    1],
                                 'Strtvl':[0.1,  5, 0.1,15,0.1,1,0.5],
                                 'Unit':['','\AA^2mueV','m1/mueV','\AA^2mueV','m1/mueV','',''],
                                 'Level':['S','S','S','S','S','q','q']}
        Fit.Plotutilities={'Separators':[['internalcont','internal contribution'],
                                         ['globalcont','global contribution'],
                                         ['immobile','immobile contribution']]}
        Fit.Constants={'Name':[],
                       'Level':[],
                       'Value':[]}
    #%%%% FIT fixed Dq^2, jump diffusion for internal diffusion
    elif FitName=="brownian_diffusion_internal_jump":
        Fit.FitName=FitName
        Fit.ProteinModel="beta*(globalcont*A*L(D*q**2,omega)+internalcont*(1-A)*L(D*q**2+Dint*q**2/(tauint*Dint*q**2+1),omega))"
        # define parameters and limits
        Fit.Parameterproperties={'Name':['D',       'Dint',     'tauint', 'beta','A'],
                                 'Lower':[0 ,        0,          0,          0,   0],
                                 'Upper':[np.inf,    np.inf,     np.inf,  np.inf, 1],
                                 'Strtvl':[5,         5,         0.1,       1,   0.5],
                                 'Unit':['\AA^2mueV','\AA^2mueV','1/mueV','',    ''],
                                 'Level':['S','S','S','q','q']}
        Fit.Plotutilities={'Separators':[['internalcont','internal contribution'],
                                         ['globalcont','global contribution']]}
        Fit.Constants={'Name':[],
                       'Level':[],
                       'Value':[]}
    #%%%% FIT fixed Dq^2, switching diffusion for internal diffusion
    elif FitName=="brownian_diffusion_internal_switching":
        Fit.FitName=FitName
        Gamma1='(D1*q**2)'
        Gamma2='(D2*q**2)'
        tau1inv='(1/tau1)'
        tau2inv='(1/tau2)'
        Lambda='(((' + Gamma1 + '-' + Gamma2 + '+' + tau1inv + '-' + tau2inv + ')**2+4*' + tau1inv + '*' + tau2inv + ')**0.5)'
        lambda1='(0.5*(' + Gamma1 + '+' + tau1inv + '+' + Gamma2 + '+' + tau2inv + '+' + Lambda + '))'
        lambda2='(0.5*(' + Gamma1 + '+' + tau1inv + '+' + Gamma2 + '+' + tau2inv + '-' + Lambda + '))'
        alpha='((1/' + Lambda + ')(tau1/(tau1+tau2)*(' + Gamma2 + '+' + tau1inv + '+' + tau2inv + '-'  + lambda1 + ')+tau2/(tau1+tau2)*(' + Gamma1 + '+' + tau1inv + '+' + tau2inv + '-'  + lambda1 + '))'
        internalcontribution = alpha + '*L(D*q**2+' + lambda1 +',omega)+(1-' + alpha + ')*L(D*q**2+' + lambda2 +')'
        Fit.ProteinModel="beta*(globalcont*A*L(D*q**2,omega)+internalcont*(1-A)" + internalcontribution
        # define parameters and limits
        Fit.Parameterproperties={'Name':['D',       'Dint',     'tauint', 'beta','A'],
                                 'Lower':[0 ,        0,          0,          0,   0],
                                 'Upper':[np.inf,    np.inf,     np.inf,  np.inf, 1],
                                 'Strtvl':[5,         5,         0.1,       1,   0.5],
                                 'Unit':['\AA^2mueV','\AA^2mueV','1/mueV','',    ''],
                                 'Level':['S','S','S','q','q']}
        Fit.Plotutilities={'Separators':[['internalcont','internal contribution'],
                                         ['globalcont','global contribution']]}
        Fit.Constants={'Name':[],
                       'Level':[],
                       'Value':[]}
    return Fit