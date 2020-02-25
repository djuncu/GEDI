#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 12:24:17 2019

@author: eardj
"""

import numpy as np

def calc_R(G, ndata):

    ########################################################################
    ## singular value decomposition of G
    # compare to matlab https://stackoverflow.com/questions/50930899/svd-command-in-python-v-s-matlab
    U, lamb, VH = np.linalg.svd(G, full_matrices=True)
    #L = np.diag(lamb)
    
    lambinv = lamb**-1
    
    L = np.zeros((U.shape[0], lamb.shape[0]))
    np.fill_diagonal(L, lamb)
    
    Linv = np.zeros((U.shape[0], lamb.shape[0]))
    np.fill_diagonal(Linv, lambinv)
    
    #Linv = np.diag(lamb**-1)   
    
    # p ~= 228?
    p=650 # 691 found by mossop script
    
       
        ##########################################
        
    #    print(p)
        
    # conjugate transpose of VH
    V = VH.T.conj()
    Vp = V[:,0:p]    
    
    Lp = np.diag(lamb[0:p])
    
    lambinvp = lamb[0:p]**-1
    Linvp = np.diag(lambinvp)
    #Linvp= np.zeros((U.shape[0], lambinvp.shape[0]))
    #np.fill_diagonal(Linvp, lambinvp)
    
    # truncated inverse -OR- natural generalized inverse
#    Gtr = Vp @ Linvp @ UTp
    
    # was this from mossop?
#    H = G @ Gtr
#    h = np.diag(H)    
#    r = ( np.identity(H.shape[0]) - H ) @ T          
        
    U2 = U[ndata*3::,:]
    U2p = U2[:,0:p]
    U2pT = np.transpose(U2p)
    
    palme = np.identity(p) - Linvp @ U2pT @ U2p @ Lp
    
    
    
    R = Vp @ palme @ np.transpose(Vp)
#    R2 = Vp @ np.transpose(Vp)
#    R3 = Gtr @ G
#    R4 = Gtr @ g


    return R

