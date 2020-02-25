#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 14:22:00 2019

Wrapper function for the quadtree algorithm

@author: Daniel Juncu
"""
import numpy as np
from scipy.interpolate import griddata
import GediUtils as gu

def downsample(xy, data, quadpixsize=1e2,tolerance=2,fittype=1,startlevel=3,maxdim=6):
    leftLim = np.min(xy[:,0])
    rightLim = np.max(xy[:,0])
    botLim = np.min(xy[:,1])
    topLim = np.max(xy[:,1])

    [xquad, yquad] = np.meshgrid(np.linspace(leftLim,rightLim,round((rightLim - leftLim) / quadpixsize)),\
        np.linspace(botLim,topLim,round((topLim - botLim) / quadpixsize)));
    
    # try bivar too?
    phgrid = griddata(xy, data, (xquad, yquad), method='linear')
        
    indmat, subdata, cx, cy, cntp, _ = gu.quadtree(phgrid,tolerance,fittype,startlevel,maxdim)
    
    # remove NaN quads
    nanquads = np.isnan(subdata)
    subdata = subdata[~nanquads]
    cx      = cx[:,~nanquads]
    cy      = cy[:,~nanquads]
    cntp    = cntp[~nanquads,:]       
        
    cntp *= quadpixsize
    cntp[:,0] = cntp[:,0] + leftLim
    cntp[:,1] = cntp[:,1] + botLim
    
    cx *= quadpixsize
    cx += leftLim
        
    cy *= quadpixsize
    cy += botLim
    
    # remove empty quads
    emptyquads = np.empty((subdata.size,), dtype=bool)  
    
    for i in range(subdata.size):
        boo = gu.inpolygon(xy[:,0:2], cx[:,i], cy[:,i])
        emptyquads[i] = ~np.any(boo)    
        
    subdata = subdata[~emptyquads]
    cx      = cx[:,~emptyquads]
    cy      = cy[:,~emptyquads]
    cntp    = cntp[~emptyquads,:]    
    
    nquads = subdata.size
    
    return indmat, subdata, cx, cy, cntp, nquads
