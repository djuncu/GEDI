#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 17:02:51 2019

Leave-one-out algorithm

@author: Daniel Juncu
"""

import numpy as np
def find_Wr(g,P,d):
    Wr_start = 1e-8
    
    Wr_opt = Wr_start
    for w in range(1,101):
        Wr = w*Wr_start
        
        G = np.concatenate((g,P*np.sqrt(Wr)))
        GT = np.transpose(G)
        
        H = G @ np.linalg.inv(GT @ G) @ GT
        h = np.diag(H)
        r = ( np.identity(H.shape[0]) - H ) @ d  
        
        CVSS = np.sum((r / (1. - h[w]) ) ** 2.)
        
        if w == 1:
            CVSS_min = CVSS
            
        if CVSS < CVSS_min:
            CVSS_min = CVSS
            Wr_opt = Wr
        print(w)   
     
    return Wr_opt


def find_p(G,d):
    ########################################################################
    ## singular value decomposition of G
    # compare to matlab https://stackoverflow.com/questions/50930899/svd-command-in-python-v-s-matlab
    U, lamb, VH = np.linalg.svd(G, full_matrices=False)
    #L = np.diag(lamb)
    
    lambinv = lamb**-1
    
    L = np.zeros((U.shape[0], lamb.shape[0]))
    np.fill_diagonal(L, lamb)
    
    Linv = np.zeros((U.shape[0], lamb.shape[0]))
    np.fill_diagonal(Linv, lambinv)
    
    
    # find optimal value for p
    pmin = 1
    for p in range(1,lamb.shape[0]):
        
        ##########################################                
        # conjugate transpose of VH
        V = VH.T.conj()
        Vp = V[:,0:p]
        #VHp = VH[:,0:p]
        Up =   U[:,0:p]        

        UTp = np.transpose(Up)
        
        lambinvp = lamb[0:p]**-1
        Linvp = np.diag(lambinvp)

        # truncated inverse -OR- natural generalized inverse
        Gtr = Vp @ Linvp @ UTp
        
        H = G @ Gtr
        h = np.diag(H)
        
        r = ( np.identity(H.shape[0]) - H ) @ d  
        
        CVSS = np.sum((r / (1. - h[p]) ) ** 2.)
                
        if p == 1:
            CVSS_min = CVSS
            
        if CVSS < CVSS_min:
            CVSS_min = CVSS
            pmin = p
        print(p)      
    
        ##################################################    
     
    return pmin
