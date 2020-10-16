"""
Fit Multiple Data Sets
======================

Fitting multiple (simulated) Gaussian data sets simultaneously.

All minimizers require the residual array to be one-dimensional. Therefore, in
the ``objective`` we need to ```flatten``` the array before returning it.

TODO: this should be using the Model interface / built-in models!

"""
import matplotlib.pyplot as plt
import numpy as np
from lmfit import Parameters, minimize, report_fit,Model

global q
q=np.linspace(0.1, 2, 20)

def gauss(x, amp, cen, sigma):
    """Gaussian lineshape."""
    return amp * np.exp(-(x-cen)**2 / (2.*sigma**2))

def objective(x,**parameterlist):
    yf=[]
    for hiq,qv in enumerate(q):
        amp=parameterlist['amp%i'%hiq]
        cen=parameterlist['cen%i'%hiq]
        sigma=parameterlist['sig%i'%hiq]
        bkg=parameterlist['bkg%i'%hiq]
        yf.append(gauss(x,amp,cen,sigma)+bkg)
    # now flatten this to a 1D array, as minimize() needs
    return [y for x in yf for y in x]
params = Parameters()
xval = np.linspace(-10, 10, 1000)
xval=np.array(xval)
for hiq,qv in enumerate(q):
    params.add('amp%i' % hiq,value=1, min=0)
    params.add('cen%i' % hiq,value=0)
    params.add('bkg%i' % hiq,value=0)
    if hiq!=0:
        params.add('sig%i' % hiq,expr='sig0')
    else:
        params.add('sig%i' % hiq,value=0.5,min=1e-7)
lambdastr='lambda x,' + ','.join(params.keys()) + ':objective(x' 
for PK in params.keys():
    lambdastr+=',' + PK + '=' + PK 
lambdastr+=')'
LambdaModel=Model(eval(lambdastr))
ObMod=Model(objective,prefix='test')
print(params)
print(ObMod)
###############################################################################
# Create five simulated Gaussian data sets
data = []
for i,qv in enumerate(q):
    amp = 0.60 + 9.50*np.random.rand()
    cen = -0.20 + 1.20*np.random.rand()
    sig = 0.25 + 0.03*np.random.rand()
    dat = gauss(xval, amp, cen, sig) + np.random.normal(size=xval.size, scale=0.1)
    data.append(dat)
data=np.array(data)
data = data.flatten()

out = ObMod.fit(data, params, x=xval)
outlambd=LambdaModel.fit(data, params,x=xval)
print(out.fit_report())
result=np.array(out.eval())
result.resize(q.shape[0],xval.shape[0])
data.resize(q.shape[0],xval.shape[0])
for i,qv in enumerate(q):
    plt.figure()
    plt.plot(xval,data[i,:])
    plt.plot(xval,result[i,:])
