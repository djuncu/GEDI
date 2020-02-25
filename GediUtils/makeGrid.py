#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 16:47:10 2019

GDIS module for creating the inversion grid

contains:
makeGrid
makeSubGrid
testGridSurf

@author: Daniel Juncu
"""
import numpy as np
import GediUtils as gu

from matplotlib import path

def makeGrid(nx,ny,nz,lon0,lat0,z0,lonExtent,latExtent,zMax,projOrigin):   

    dz = (zMax-z0) / nz
    
    if projOrigin[0] == lon0 and projOrigin[1] == lat0:
        xsrc = np.linspace(-lonExtent/2.,lonExtent/2., nx)
        ysrc = np.linspace(-latExtent/2.,latExtent/2., ny)
        zsrc = -z0 - dz/2. + -1. *np.linspace(0., zMax-z0-dz, nz)
    else: # not currently included
        # this is required when Origin of grid and origin of map projection are different
        xy0 = gu.wgs2local(np.array(([lon0, lat0])),projOrigin)        
        xsrc = np.linspace(xy0[0],xy0[0]+lonExtent, nx)
        ysrc = np.linspace(xy0[1],xy0[1]+latExtent, ny)
        zsrc = np.linspace(z0+dz/2., z0+zMax, nz)
    
    zgrid, ygrid, xgrid = np.meshgrid(zsrc, ysrc, xsrc, sparse=False, indexing='ij')

    return xgrid, ygrid, zgrid, dz

def makeSubGrid(xgrid, ygrid, zgrid, nx, ny, nz, dz, ngrid, nlvl):
    
    dx = np.diff(xgrid,axis=2)[0][0][0]
    dy = np.diff(ygrid,axis=1)[0][0][0]
    
    xsubgrid = np.zeros((ngrid,(nlvl+1)**3))
    ysubgrid = np.zeros((ngrid,(nlvl+1)**3))
    zsubgrid = np.zeros((ngrid,(nlvl+1)**3))

    index = 0
    for j in range(0,nx):
        for i in range(0,ny):
            for k in range(0,nz): 
                subindex=0
                for l in range(nlvl+1):
                    for m in range(nlvl+1):
                        for n in range(nlvl+1):
                            xsubgrid[index,subindex] = xgrid[k,i,j] - dx/2. + 0.5 * dx/(nlvl+1) + l * dx/(nlvl+1)
                            ysubgrid[index,subindex] = ygrid[k,i,j] - dy/2. + 0.5 * dy/(nlvl+1) + m * dy/(nlvl+1)
                            zsubgrid[index,subindex] = zgrid[k,i,j] - dz/2. + 0.5 * dz/(nlvl+1) + n * dz/(nlvl+1)
                            subindex+=1
                index += 1
                
    return xsubgrid, ysubgrid, zsubgrid
    

def testGridSurf(data_xyz, xgrid, ygrid, zgrid, dz):
    # test if Grid and terrain surface intersect    

    nx = xgrid.shape[1]
    ny = xgrid.shape[2]   
    
    xExtent = np.max(xgrid) - np.min(xgrid)
    yExtent = np.max(ygrid) - np.min(ygrid)
    
    xSize = xExtent / nx
    ySize = yExtent / ny
    
    xgrid = xgrid[0,:,:]
    ygrid = ygrid[0,:,:]
    
    gridshape = xgrid.shape
    
    px1 = xgrid.reshape(-1) - xSize/2.
    px2 = xgrid.reshape(-1) + xSize/2.
    py1 = ygrid.reshape(-1) - ySize/2.
    py2 = ygrid.reshape(-1) + ySize/2.
    
    ndata = data_xyz.shape[0]
    
    inElm = np.zeros((nx*ny,ndata), dtype=bool)
    inElmLow = np.zeros((nx*ny), dtype=int)
    topo=[]    
    heightFlag = data_xyz[:,2] < (np.max(zgrid) + dz/2.)
    
    for i in range(nx*ny):
        p = path.Path([(px1[i],py1[i]), (px2[i], py1[i]),
                            (px2[i], py2[i]), (px1[i], py2[i])]) 
        inElm[i,:] = p.contains_points(data_xyz[:,0:2])        
        topo.append(data_xyz[inElm[i,:],2])             
        #no of points in element and lower than element
        inElmLow[i] = np.sum(np.logical_and(inElm[i,:], heightFlag))
        

#    heightFlag = data_xyz[:,2] < 0 # datapoint lower than element surface   
#    dataFlag =  np.sum(inElm,axis=0).astype(bool)  # flagged each datapoint if in an element    
#    inflags_heights = np.logical_and(dataFlag, heightFlag) # datapoint in element and lower than element surface    
#    elmNumbers = np.sum(inElm,axis=1) #.reshape(gridshape) # number of data per element
#    elmFlag = np.sum(inElm,axis=1).astype(bool).astype(int).reshape(gridshape)    

#    lowRatio = elmNumbers / inElmLow.astype(float)        
#    # IS THIS REASONABLE??
#    lowRatio_boo = lowRatio < 2. 

    elmAboveSurf = inElmLow.astype(bool)
    return elmAboveSurf.reshape(gridshape)
    
#    cmap_redgreen = cm.get_cmap('RdYlGn_r')
#    fig = plt.figure()
#    im = plt.pcolormesh(xgrid, ygrid, elmFlag,cmap = cmap_redgreen, edgecolors='black')    
    
#    cmap_redgreen = cm.get_cmap('RdYlGn_r')
#    fig = plt.figure()
#    im = plt.pcolormesh(xgrid, ygrid, lowRatio_boo.reshape(gridshape),cmap = cmap_redgreen, edgecolors='black')       
#    
#    cmap_redgreen = cm.get_cmap('RdYlGn_r')
#    plt.figure()
#    plt.pcolormesh(xgrid, ygrid, lowRatio_boo_all.reshape(gridshape),cmap = cmap_redgreen, edgecolors='black')       
    
