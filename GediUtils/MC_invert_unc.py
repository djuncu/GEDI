#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 13:57:29 2019

MCMC inversion for noisy data

@author: Daniel Juncu
"""
import numpy as np

def MC_invert_unc(nMC, C, alldata, GTG, GT, T, ndata, ngrid):
    vcmc = np.zeros([ngrid,nMC])
    #noisydata = np.empty([ndata])
    L = np.linalg.cholesky(C)
    
    noisydata = alldata + np.random.normal(size=[nMC,ndata]) @ L
    noisydata = noisydata.T
    
    for i in range(int(nMC)):
        T = np.concatenate((noisydata[:,i],np.zeros(ngrid)),axis=0)
        vcmc[:,i] = np.linalg.solve(GTG, GT @ T)
        
    return vcmc
