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
    def __init__(self, FitName='', FitModel = '', Parameters = [], Comments = [],Subfit=[],Background=''):
        # initialising all the variables
        self.FitName = FitName
        self.FitModel = FitModel #string
        self.Parameters = Parameters
        self.Comments = Comments
        self.Subfit=Subfit
    @staticmethod
    def fit_QENS(data,FitName):
        if FitName=='whatever':


