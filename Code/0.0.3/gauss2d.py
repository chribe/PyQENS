#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 11:03:33 2020
https://stackoverflow.com/questions/52148141/2d-gaussian-fit-using-lmfit
@author: christian
"""


import numpy as np
from lmfit import Parameters, minimize, report_fit, Model

x, y = np.meshgrid(np.linspace(-1, 1, 10), np.linspace(-1, 1, 10))

def gaussian2D(x, y, cen_x, cen_y, sig_x, sig_y, offset):
    return np.exp(-(((cen_x-x)/sig_x)**2 + ((cen_y-y)/sig_y)**2)/2.0) + offset
def gaussian2D_1(x, y, height,cen_x, cen_y, sig_x, sig_y, offset):
    return np.exp(-(((cen_x-x)/sig_x)**2 + ((cen_y-y)/sig_y)**2)/2.0) + offset

def residuals(p, x, y, z):
    height = p["height"].value
    cen_x = p["centroid_x"].value
    cen_y = p["centroid_y"].value
    sigma_x = p["sigma_x"].value
    sigma_y = p["sigma_y"].value
    offset = p["background"].value
    return (z - height*gaussian2D(x,y, cen_x, cen_y, sigma_x, sigma_y, offset))

# test data
g = gaussian2D(x, y, 1.2, 2.1, 0.5, 0.7, 1.1)+0.1*np.random.random(x.shape)


initial = Parameters()
initial.add("height",value=1.)
initial.add("centroid_x",value=0.)
initial.add("centroid_y",value=0.)
initial.add("sigma_x",value=1.)
initial.add("sigma_y",value=3.)
initial.add("background",value=0.)

params=Parameters()
params.add("height",value=5)
params.add("cen_x",value=0.)
params.add("cen_y",value=0.)
params.add("sig_x",value=1.)
params.add("sig_y",value=3.)
params.add("offset",value=0.)

fit = minimize(residuals, initial, args=(x, y, g))
print(report_fit(fit))

GModel=Model(gaussian2D_1,independent_vars=['x','y'])
FT=GModel.fit(g,params,x=x,y=y)
print(report_fit(FT))