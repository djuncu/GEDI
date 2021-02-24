#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 12:02:25 2019

This script runs GEDI, a geodetic direct inversion software
Argument is the parameter file, see example

@author: Daniel Juncu
"""
import sys
#parmfile =

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
#from mpl_toolkits.mplot3d import Axes3D

from pathlib import Path
import time

#import parms_aluto_blob as parms
#import parmfile as parms


import GediUtils as gu
from scipy import optimize

#import parameters
parms=__import__(sys.argv[1]) 

#import warnings
#warnings.filterwarnings('error')
np.seterr(all='warn') # or 'raise'

savedir = Path('data/' + parms.prjname)
savedir.mkdir(parents=True, exist_ok=True)

# obtain parameters
origin         = np.asarray(parms.origin)
nx             = int(parms.nx)
ny             = int(parms.ny)
nz             = int(parms.nz)
grid_z0        = float(parms.grid_z0)
grid_lonExtent = float(parms.grid_lonExtent)
grid_latExtent = float(parms.grid_latExtent)
grid_zMax  = float(parms.grid_zMax)
nMC            = int(parms.nMC)
nu = parms.nu
Wr = parms.Wr
quadpixsize = parms.quadpixsize
do_dwnsmpl = bool(parms.do_dwnsmpl)
load_dwnsmpl = bool(parms.load_dwnsmpl)
regularization = bool(parms.regularization)
datatype = parms.datatype
subgridlvl = parms.subgridlvl
determine_pmin = bool(parms.determine_pmin)
#datashift = parms.datashift
#east_look = parms.east_look
#north_look = parms.north_look

quad_cosE = []
quad_cosN = []
quad_cosU = []

ngrid = nx*ny*nz

print('Loading data...')
nInsarFrames = len(parms.datalocation)
if datatype == 'los':
    insar_data, insar_xyz, insar_ll, cosN, cosE, cosU, ndata, origin, topo = gu.loadInsarData(parms.datalocation, origin, estimateOrigin=False)
elif datatype == 'eu':
    # this for E-U insar data
    e_data, u_data, insar_xyz, insar_ll, ndata, origin, topo = gu.loadInsarData_eu(parms.datalocation, origin, estimateOrigin=False)
elif datatype == 'enu':
    e_data, n_data, u_data, insar_xyz, insar_ll, ndata, origin, topo = gu.loadInsarData_enu(parms.datalocation, origin, estimateOrigin=False)            

# cropping
# make this better in the future and remove hardcoding
aoi_lim = 38.5

crop_id = np.argwhere(insar_ll[:,0] < aoi_lim) # 38.5

ndata_all = np.sum(ndata)

# just for counting
for i in range(nInsarFrames):
    crop_id_i = np.argwhere(insar_ll[ndata[i-1]*i:ndata[i]+i*ndata[i-1],0] < aoi_lim)
    ndata[i] -= len(crop_id_i)   

if datatype == 'los': 
    insar_data = np.delete(insar_data, crop_id)
    cosE = np.delete(cosE, crop_id)
    cosN = np.delete(cosN, crop_id)
    cosU = np.delete(cosU, crop_id)
elif datatype == 'eu':  
    e_data = np.delete(e_data, crop_id)
    u_data = np.delete(u_data, crop_id)
elif datatype == 'enu':
    e_data = np.delete(e_data, crop_id)
    n_data = np.delete(n_data, crop_id)
    u_data = np.delete(u_data, crop_id)
    
insar_xyz = np.delete(insar_xyz, crop_id, axis=0)
insar_ll = np.delete(insar_ll, crop_id, axis=0)
ndata_all -= len(crop_id)

if do_dwnsmpl == False:
    
    if datatype == 'eu':
        alldata = np.concatenate((e_data,u_data), axis=0)
        full_xyz = np.tile(insar_xyz, (2,1))  
    elif datatype == 'enu':
        alldata = np.concatenate((e_data,n_data, u_data), axis=0)
        full_xyz = np.tile(insar_xyz, (3,1))  
    elif datatype =='los':
        alldata = insar_data
    
else:
    if load_dwnsmpl:
        print('Loading downsampled data...')
    #    cx = np.loadtxt(str(savedir / 'quad_cx'))
    #    cy = np.loadtxt(str(savedir / 'quad_cy'))
        
        quad_data = np.load(str(savedir / 'quad_data.npy'), allow_pickle=False )
        quad_xy = np.load(str(savedir / 'quad_xy.npy'),  allow_pickle=False )
        cx = np.load(str(savedir / 'cx.npy'), allow_pickle=False )
        cy = np.load(str(savedir / 'cy.npy'), allow_pickle=False )
        if datatype == 'los':
            quad_cosE = np.load(str(savedir / 'quad_cosE.npy'), allow_pickle=False )
            quad_cosN = np.load(str(savedir / 'quad_cosN.npy'), allow_pickle=False )
            quad_cosU = np.load(str(savedir / 'quad_cosU.npy'), allow_pickle=False )
        
    #    quad_data = np.loadtxt(str(savedir / 'quad_data'))
    #    quad_xy = quad_data[:,0:2]
    #    quad_data = np.delete(quad_data, [0,1], axis=1)
    #    quad_data = quad_data.flatten()
        nquads = cx.shape[1]
        
        alldata = quad_data
        
    else:    
        nquads = []
        quad_xy = np.empty((0,2))   
        cx = np.empty((4,0))
        cy = np.empty((4,0))
        
        quad_cosN = []
        quad_cosE = []
        quad_cosU = []
        
        alldata = []
                
        for i in range(nInsarFrames):
        #alldata = np.concatenate((data[:,3],data[:,4],data[:,5]),axis=0)
            print('Downsampling data...')
            if datatype == 'los':
                indmat, quad_data, cx_i, cy_i, quad_xy_i, nquads_i = gu.downsample(insar_xyz[ndata[i-1]*i:ndata[i]+i*ndata[i-1],0:2],
                                                                           insar_data[ndata[i-1]*i:ndata[i]+i*ndata[i-1]], quadpixsize=parms.quadpixsize ,tolerance=parms.quadtolerance,fittype=parms.quadfittype,startlevel=parms.quadstartlevel,maxdim=parms.quadmaxdim)
                
                nquads.append(nquads_i) 
                cx = np.append(cx,cx_i, axis=1)
                cy = np.append(cy,cy_i, axis=1)
                quad_xy = np.append(quad_xy,quad_xy_i, axis=0)
                
                print('Reprojecting onto quadtree grid...')
                quad_cosE = np.append(quad_cosE, gu.projectOnQuads(insar_xyz[ndata[i-1]*i:ndata[i]+i*ndata[i-1],0:2],cx_i,cy_i,cosE[ndata[i-1]*i:ndata[i]+i*ndata[i-1]] ))
                quad_cosN = np.append(quad_cosN, gu.projectOnQuads(insar_xyz[ndata[i-1]*i:ndata[i]+i*ndata[i-1],0:2],cx_i,cy_i,cosN[ndata[i-1]*i:ndata[i]+i*ndata[i-1]] ))
                quad_cosU = np.append(quad_cosU, gu.projectOnQuads(insar_xyz[ndata[i-1]*i:ndata[i]+i*ndata[i-1],0:2],cx_i,cy_i,cosU[ndata[i-1]*i:ndata[i]+i*ndata[i-1]] ))
                
                alldata = np.append(alldata, quad_data)
           
            elif datatype == 'eu':
                indmat, quad_e, cx, cy, quad_xy, nquads = gu.downsample(insar_xyz[:,0:2], 
                                                                           e_data, quadpixsize=parms.quadpixsize,
                                                                           tolerance=parms.quadtolerance,fittype=parms.quadfittype,
                                                                           startlevel=parms.quadstartlevel,maxdim=parms.quadmaxdim)
        #        nquads *= 2
                print('Reprojecting onto quadtree grid...')
                quad_u = gu.projectOnQuads(insar_xyz[:,0:2],cx,cy,u_data)
                quad_data = np.concatenate((quad_e,quad_u), axis=0)
    #            quad_xy = np.tile(quad_xy, (2,1))
                alldata = quad_data
                
            elif datatype == 'enu':
                indmat, quad_u, cx, cy, quad_xy, nquads = gu.downsample(insar_xyz[:,0:2], 
                                                                           u_data, quadpixsize=parms.quadpixsize,
                                                                           tolerance=parms.quadtolerance,fittype=parms.quadfittype,
                                                                           startlevel=parms.quadstartlevel,maxdim=parms.quadmaxdim)
        #        nquads *= 2
                print('Reprojecting onto quadtree grid...')
                quad_e = gu.projectOnQuads(insar_xyz[:,0:2],cx,cy,e_data)
                quad_n = gu.projectOnQuads(insar_xyz[:,0:2],cx,cy,n_data)
                quad_data = np.concatenate((quad_e,quad_n,quad_u), axis=0)
                
                alldata = quad_data
    #            quad_xy = np.tile(quad_xy, (3,1))
                
    
        if datatype == 'los':
            np.save(str(savedir / 'quad_cosE'), quad_cosE, allow_pickle=False )
            np.save(str(savedir / 'quad_cosN'), quad_cosN, allow_pickle=False )
            np.save(str(savedir / 'quad_cosU'), quad_cosU, allow_pickle=False )
    
        np.save(str(savedir / 'quad_data'), quad_data, allow_pickle=False )
        np.save(str(savedir / 'quad_xy'), quad_xy, allow_pickle=False )
        np.save(str(savedir / 'cx'), cx, allow_pickle=False )
        np.save(str(savedir / 'cy'), cy, allow_pickle=False )         
        
        

if topo:
    # adjust heights: heights are referenced to grid_z0. grid_z0 becomes the vertical
    # distance between the highest (now: mean?) surface altitude and the grid
    
#    quad_heights = gu.projectOnQuads(insar_xyz[:,0:2],cx,cy,insar_xyz[:,2])
#    #quad_heights = -1 * ( quad_heights - ( np.max(quad_heights) - grid_z0 ) )
#    quad_heights =  np.max(quad_heights) - quad_heights 
#    quad_xyz = np.hstack((quad_xy,quad_heights[:,None]))  
   
    #insar_xyz[:,2] = -1 * ( insar_xyz[:,2] - ( np.max(insar_xyz[:,2]) - grid_z0 ) )
    refHeight = np.mean(insar_xyz[:,2])
    print('Reference surface elevation is ' + str(refHeight.astype(int)) + 'm.' )
    insar_xyz[:,2] = -1. * (( refHeight - insar_xyz[:,2] ) ) # - grid_z0  )
    #insar_xyz[:,2] -= np.max(insar_xyz[:,2])
    
    if do_dwnsmpl:
        quad_heights = np.empty((np.sum(nquads)))
        for i in range(nInsarFrames):
            quad_heights[nquads[i-1]*i:nquads[i]+i*nquads[i-1]] = gu.projectOnQuads(insar_xyz[ndata[i-1]*i:ndata[i]+i*ndata[i-1],0:2],
                        cx[:,nquads[i-1]*i:nquads[i]+i*nquads[i-1]],
                        cy[:,nquads[i-1]*i:nquads[i]+i*nquads[i-1]],
                        insar_xyz[ndata[i-1]*i:ndata[i]+i*ndata[i-1],2])
        
        if datatype == 'eu': quad_heights = np.tile(quad_heights, (2,1))
        if datatype == 'enu': quad_heights = np.tile(quad_heights, (3,1))
        
        quad_xyz = np.hstack((quad_xy,quad_heights[:,None]))   
else:
    refHeight = 0    
    if do_dwnsmpl == True:
        quad_xyz = np.hstack((quad_xy,np.zeros((nquads,1))))
#        if datatype == 'los': quad_xyz = np.hstack((quad_xy,np.zeros((nquads,1))))
#        if datatype == 'eu':  quad_xyz = np.hstack((quad_xy,np.zeros((2*nquads,1))))
#        if datatype == 'enu': quad_xyz = np.hstack((quad_xy,np.zeros((3*nquads,1))))

# need to average bunch of quantities for each quad here (e.g. surface height, if applicable)
#quad_inc = gu.projectOnQuads(insar_xyz[:,0:2],cx,cy,insar_inc)
# 

# combine all datasets here
alldata_xyz = insar_xyz



if do_dwnsmpl:
    xyz_inv = quad_xyz
    ndata_inv = np.sum(nquads)
    cosE = quad_cosE
    cosN = quad_cosN
    cosU = quad_cosU
    
    # shift data if necessary
    for i in range(nInsarFrames):        
        parms.datashift[i]
        alldata[nquads[i-1]*i:nquads[i]+i*nquads[i-1]] += parms.datashift[i]
else:
#    xyz_inv = full_xyz
    xyz_inv = insar_xyz
    ndata_inv = np.sum(ndata)
    
    # shift data if necessary
    for i in range(nInsarFrames):        
        alldata[ndata[i-1]*i:ndata[i]+i*ndata[i-1]] += parms.datashift[i]

#nalldata = nquads

#gu.plotQuads(cx/1e3,cy/1e3,quad_data[:nquads], np.min(quad_data), np.max(quad_data), 'Greens')

if datatype == 'los':
    for i in range(nInsarFrames):
        fig = plt.figure()
        h2 = plt.scatter(insar_ll[ndata[i-1]*i:ndata[i]+i*ndata[i-1],0],insar_ll[ndata[i-1]*i:ndata[i]+i*ndata[i-1],1],15,
                         insar_data[ndata[i-1]*i:ndata[i]+i*ndata[i-1]])
        fig.colorbar(h2)
elif (datatype == 'eu') or (datatype == 'enu'):
    fig = plt.figure()
    h2 = plt.scatter(insar_xyz[:,0],insar_xyz[:,1],15,e_data,vmin=np.min(e_data),vmax=np.max(e_data))
    fig.colorbar(h2)

    
print('Creating model grid...')
xgrid, ygrid, zgrid, dz = gu.makeGrid(nx,ny,nz,origin[0], origin[1], grid_z0, grid_lonExtent, grid_latExtent, grid_zMax, origin)

xsubgrid, ysubgrid, zsubgrid = gu.makeSubGrid(xgrid, ygrid, zgrid, nx, ny, nz, dz, ngrid, subgridlvl)



print('Creating G matrix...')
G, g, d, P  = gu.createGmatrix(nu,xyz_inv,xgrid,ygrid,zgrid, xsubgrid, ysubgrid, zsubgrid, subgridlvl, nx,ny,nz,dz,alldata,ndata_inv,datatype,Wr, cosE=cosE, cosN=cosN, cosU=cosU, regularization=regularization)
#G, g, d, P  = gu.createGmatrix_nosubgrid(nu,xyz_inv,xgrid,ygrid,zgrid,nx,ny,nz,alldata,ndata_inv,datatype,Wr, cosE=quad_cosE, cosN=quad_cosN, cosU=quad_cosU, regularization=regularization)

#print('Find optimal regularization weight factor...')
#Wr = gu.find_Wr(g,P,d)
#print('Optimal Wr value found: ' + str(Wr))
#G = np.concatenate((g,P*np.sqrt(Wr)))
print('Testing wether grid and topography intersect...')
# only test on first dataset if there's 2
elmAboveSurf = gu.testGridSurf(insar_xyz[0:ndata[0]], xgrid, ygrid, zgrid, dz)
if np.any(elmAboveSurf):
    print('WARNING: Grid and topography intersect. Appending constraints to G matrix..')
    cmap_redgreen = cm.get_cmap('RdYlGn_r')
    plt.figure()
    plt.pcolormesh(xgrid[0,:,:], ygrid[0,:,:], elmAboveSurf,cmap = cmap_redgreen, edgecolors='black') 
    
    G, d = gu.appendGmatrix(G,d,elmAboveSurf,ny, nz)      
else:
    print('No intersection between grid and topography.')

if determine_pmin:
    print('Determining p value....')
    pmin = gu.find_p(G,d)
    print('P cut-off value found: ' + str(pmin))
else:
    pmin = G.shape[1] - parms.cutoff
    print('Using cut-off value of ' + str(parms.cutoff))
    
print('Setting up singular value decomposition...')
### SVD
U, lamb, VH = np.linalg.svd(G, full_matrices=False)
#L = np.diag(lamb)
lambinv = lamb**-1
L = np.zeros((U.shape[0], lamb.shape[0]))
np.fill_diagonal(L, lamb)
Linv = np.zeros((U.shape[0], lamb.shape[0]))
np.fill_diagonal(Linv, lambinv)

#V = VH.T.conj()
V = np.transpose(VH)

Vp = V[:,0:pmin]
#VHp = VH[:,0:p]
Up =   U[:,0:pmin]
#Vp = np.transpose(VHp)
UTp = np.transpose(Up)
lambinvp = lamb[0:pmin]**-1
Linvp = np.diag(lambinvp)
#Linvp= np.zeros((U.shape[0], lambinvp.shape[0]))
#np.fill_diagonal(Linvp, lambinvp)


### resolution matrix   
if regularization:
    Lp = np.diag(lamb[0:pmin]) 
    if datatype == 'los':  U2 = U[ndata_inv::,:]
    if datatype == 'eu':   U2 = U[ndata_inv*2::,:]
    if datatype == 'enu':  U2 = U[ndata_inv*3::,:]
    U2p = U2[:,0:pmin]
    U2pT = np.transpose(U2p)
    
    palme = np.identity(pmin) - Linvp @ U2pT @ U2p @ Lp
    R = Vp @ palme @ np.transpose(Vp)
else:
    R = Vp @ np.transpose(Vp)


# truncated inverse -OR- natural generalized inverse
#Ginvtr = Vp @ Linvp @ UTp
#vc = Ginvtr @ d


print('Setting up non-negative inversion...')
### testarea for nonneg
# trans1
H = np.diag(1.*np.ones(ngrid))
h = np.zeros(ngrid)
Hp = -H @ Vp @ Linvp
hp = h + Hp @ UTp @ d
#trans2
Gp = np.transpose(np.concatenate((Hp, hp[:,None]), axis=1))
dp  = np.transpose(np.append(np.zeros(Hp.shape[1]), 1.))
mpp,_ = optimize.nnls(Gp, dp)
#mpp,_ = gu.nnls(Gp, dp,1e5)
ep = dp - Gp @ mpp
mp = -ep[0:-1] / ep[-1]
#take mp back to m


print('Solving...')
# non neg
vc = Vp @ Linvp @ (UTp @ d - mp)
# 'regular' SVD
#vc = Vp @ Linvp @ UTp @ d
# regular least squares
# this is using the least sqaures general inverse: (G.T * G)^-1 * G.T | G.T: transpose of G
#GT = np.transpose(G)
#vc = np.linalg.solve(GT @ G, GT @ d)

#Pdiag = np.diag(P)
Pdiag = np.sum(P,axis=1)


print('Reprojecting...')
vcc = np.zeros((nz,ny,nx))
Rc = np.zeros((nz,ny,nx))
#Pc = np.zeros((nz,nx,nx))
testcount = np.zeros((nz,nx,nx))

index=0
for i in range(0,nx):
    for j in range(0,ny):
        for k in range(0,nz): 
            Rc[k,j,i] = R[index,index]
            vcc[k,j,i] = vc[index]
#            Pc[k,j,i] =Pdiag[index]
            testcount[k,j,i] = index
            index+=1

# calculate forward model
dpred = G @ vc

if datatype == 'enu':
    residuals = alldata-dpred[0:3*ndata]
elif datatype == 'eu':
    residuals = alldata-dpred[0:2*ndata]
elif datatype == 'los':
    residuals = alldata-dpred[0:ndata_inv]
RSS = np.dot(residuals, residuals)
print('Residual sum of squares: % 6.4f' % RSS)

print('Plotting...')

# make grids for pcolormesh plotting...
dx = np.diff(xgrid,axis=2)[0][0][0]
dy = np.diff(ygrid,axis=1)[0][0][0]
xgplot = np.concatenate((xgrid[0,:,:]-dy/2, xgrid[0,:,-1][:,None]+dy/2),axis=1)
xgplot = np.concatenate((xgplot,xgplot[-1,:][None,:]),axis=0)    
ygplot = np.concatenate((ygrid[0,:,:]-dy/2, ygrid[0,-1,:][None,:]+dy/2),axis=0)
ygplot = np.concatenate((ygplot, ygplot[:,-1][:,None]), axis=1)

## PLOT QUADS
#
#if datatype == 'los':
#    ulos = gu.forwardmodel_los(nu,insar_xyz,ngrid,xsubgrid,ysubgrid,zsubgrid,subgridlvl,vc,nx,ny,nz, cosE, cosN, cosU)
#    
#    if do_dwnsmpl == True:
#        ulos_quad = gu.forwardmodel_los(nu,quad_xyz,ngrid,xsubgrid,ysubgrid,zsubgrid,subgridlvl,vc,nx,ny,nz, quad_cosE, quad_cosN, quad_cosU)    
#        gu.plotQuads(cx/1e3,cy/1e3,ulos_quad, np.min(quad_data), np.max(quad_data), 'Greens')
#        gu.plotQuads(cx/1e3,cy/1e3,quad_xyz[:,2], np.min(quad_xyz[:,2]), np.max(quad_xyz[:,2]), 'terrain')
#elif (datatype == 'eu') or (datatype == 'enu'):
#    ue, un, uu = gu.forwardmodel_enu(nu,insar_xyz,ngrid,xsubgrid,ysubgrid,zsubgrid,subgridlvl,vc,nx,ny,nz)
#    ue, un, uu = gu.forwardmodel_enu_nosubgrid(nu,insar_xyz,xgrid,ygrid,zgrid,vcc,nx,ny,nz)
#    

#    if do_dwnsmpl == True:
#        ue_quad, un_quad, uu_quad = gu.forwardmodel_enu(nu,quad_xyz[:nquads,:],ngrid,xsubgrid,ysubgrid,zsubgrid,subgridlvl,vc,nx,ny,nz)    
#        gu.plotQuads(cx/1e3,cy/1e3,ue_quad, np.min(quad_data), np.max(quad_data), 'Greens')
#        gu.plotQuads(cx/1e3,cy/1e3,uu_quad, np.min(quad_data), np.max(quad_data), 'Greens')
    
    
#    gu.plotQuads(cx/1e3,cy/1e3,quad_xyz[:,2], np.min(quad_xyz[:,2]), np.max(quad_xyz[:,2]), 'terrain')

cmap_reversed = cm.get_cmap('viridis')

plotsize = np.ceil(np.sqrt(nz)).astype(int)
if plotsize == 1:
    plotsize = 2

fig, ax = plt.subplots(nrows=plotsize, ncols=plotsize)
#plt.figure()
depthid=0
for row in ax:
    for col in row:
       if depthid < nz: 
           im = col.pcolormesh(xgplot/1e3, ygplot/1e3, vcc[depthid,:,:],vmin=np.min(vc),vmax=np.max(vc),cmap = cmap_reversed)
           col.text(.85*np.min(xgplot/1e3), .75*np.max(ygplot/1e3), str(int(refHeight+zgrid[depthid,0,0])) + ' m', fontsize=12, color='white')
#           col.imshow(vcc[depthid,:,:],interpolation='bilinear',
#                        cmap = cmap_reversed,
#                        origin='lower',vmin=np.min(vc),vmax=np.max(vc))
           depthid+=1

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im, cax=cbar_ax)


fig, ax = plt.subplots(nrows=plotsize, ncols=plotsize)
#plt.figure()
depthid=0
for row in ax:
    for col in row:
       if depthid < nz: 
#           im = col.matshow(Rc[depthid,:,:],extent=(np.min(xgrid[0,:,:]/1e3), np.max(xgrid[0,:,:]/1e3), np.min(ygrid[0,:,:]/1e3), np.max(ygrid[0,:,:]/1e3)))
           im = col.pcolormesh(xgplot/1e3, ygplot/1e3, Rc[depthid,:,:],vmin=np.min(Rc),vmax=np.max(Rc),cmap = 'Greys')
           col.text(-6, 6, str(int(refHeight+zgrid[depthid,0,0])) + ' m', fontsize=12, color='white')

#           col.imshow(vcc[depthid,:,:],interpolation='bilinear',
#                        cmap = cmap_reversed,
#                        origin='lower',vmin=np.min(vc),vmax=np.max(vc))
           depthid+=1

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im, cax=cbar_ax)

### plot P diagonal
#fig, ax = plt.subplots(nrows=plotsize, ncols=plotsize)
##plt.figure()
#depthid=0
#for row in ax:
#    for col in row:
#       if depthid < nz: 
##           im = col.matshow(Rc[depthid,:,:],extent=(np.min(xgrid[0,:,:]/1e3), np.max(xgrid[0,:,:]/1e3), np.min(ygrid[0,:,:]/1e3), np.max(ygrid[0,:,:]/1e3)))
#           im = col.pcolormesh(xgplot/1e3, ygplot/1e3, Pc[depthid,:,:],cmap = 'RdYlGn')
#           col.text(-6, 6, str(int(refHeight+zgrid[depthid,0,0])) + ' m', fontsize=12, color='white')
#
##           col.imshow(vcc[depthid,:,:],interpolation='bilinear',
##                        cmap = cmap_reversed,
##                        origin='lower',vmin=np.min(vc),vmax=np.max(vc))
#           depthid+=1
#
#fig.subplots_adjust(right=0.8)
#cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
#fig.colorbar(im, cax=cbar_ax)


#cmap_reversed = cm.get_cmap('viridis_r')
#fig, ax = plt.subplots(nrows=2, ncols=3)
##plt.figure()
#depthid=0
#for row in ax:
#    for col in row:
#       if depthid < nz: 
#           im = col.pcolormesh(xgrid[0,:,:], ygrid[0,:,:], Rc[depthid,:,:],vmin=0,vmax=1,cmap = cmap_reversed)
##           col.imshow(vcc[depthid,:,:],interpolation='bilinear',
##                        cmap = cmap_reversed,
##                        origin='lower',vmin=np.min(vc),vmax=np.max(vc))
#           depthid+=1





#fig = plt.figure()
#h2 = plt.scatter(insar_xyz[:,0],insar_xyz[:,1],15,insar_data)
#fig.colorbar(h2)

fig = plt.figure()
h2 = plt.scatter(insar_xyz[:,0],insar_xyz[:,1],15,insar_xyz[:,2])
fig.colorbar(h2)

if not do_dwnsmpl:   
    
#        else:
#            fig = plt.figure()
#            h2 = plt.scatter(insar_xyz[:,0],insar_xyz[:,1],15,dpred[0:ndata_inv])
#            fig.colorbar(h2)
    if (datatype == 'enu'):
        fig = plt.figure()
    #    h2 = plt.scatter(insar_xyz[:,0],insar_xyz[:,1],15,dpred[0:ndata],vmin=np.min(e_data),vmax=np.max(e_data))
    #    fig.colorbar(h2)
        
        plt.subplot(3,3,1)
        he = plt.scatter(insar_xyz[:,0]/3e4,insar_xyz[:,1]/3e4,15,e_data,vmin=np.min(e_data),vmax=np.max(e_data))  
        fig.colorbar(he)
        plt.subplot(3,3,2)
        hn = plt.scatter(insar_xyz[:,0]/3e4,insar_xyz[:,1]/3e4,15,n_data,vmin=np.min(n_data),vmax=np.max(n_data))   
        fig.colorbar(hn)     
        plt.subplot(3,3,3)
        hu = plt.scatter(insar_xyz[:,0]/3e4,insar_xyz[:,1]/3e4,15,u_data,vmin=np.min(u_data),vmax=np.max(u_data)) 
        fig.colorbar(hu)
        
        plt.subplot(3,3,4)
        he = plt.scatter(insar_xyz[:,0]/3e4,insar_xyz[:,1]/3e4,15,dpred[0:ndata],vmin=np.min(e_data),vmax=np.max(e_data))  
        fig.colorbar(he)
        plt.subplot(3,3,5)
        hn = plt.scatter(insar_xyz[:,0]/3e4,insar_xyz[:,1]/3e4,15,dpred[ndata:2*ndata],vmin=np.min(n_data),vmax=np.max(n_data))   
        fig.colorbar(hn)     
        plt.subplot(3,3,6)
        hu = plt.scatter(insar_xyz[:,0]/3e4,insar_xyz[:,1]/3e4,15,dpred[2*ndata:3*ndata],vmin=np.min(u_data),vmax=np.max(u_data)) 
        fig.colorbar(hu)
        
        plt.subplot(3,3,7)
        hr1 = plt.scatter(insar_xyz[:,0]/3e4,insar_xyz[:,1]/3e4,15,e_data-dpred[0:ndata])  
        fig.colorbar(hr1)
        plt.subplot(3,3,8)
        hr2 = plt.scatter(insar_xyz[:,0]/3e4,insar_xyz[:,1]/3e4,15,n_data-dpred[ndata:2*ndata])   
        fig.colorbar(hr2)     
        plt.subplot(3,3,9)
        hr3 = plt.scatter(insar_xyz[:,0]/3e4,insar_xyz[:,1]/3e4,15,u_data-dpred[2*ndata:3*ndata]) 
        fig.colorbar(hr3)
        
    elif (datatype == 'eu'):  
        fig = plt.figure()
    #    h2 = plt.scatter(insar_xyz[:,0],insar_xyz[:,1],15,dpred[0:ndata],vmin=np.min(e_data),vmax=np.max(e_data))
    #    fig.colorbar(h2)
        
        plt.subplot(3,2,1)
        he = plt.scatter(insar_xyz[:,0]/3e4,insar_xyz[:,1]/3e4,15,e_data,vmin=np.min(e_data),vmax=np.max(e_data))  
        fig.colorbar(he)
        plt.subplot(3,2,2)
        hu = plt.scatter(insar_xyz[:,0]/3e4,insar_xyz[:,1]/3e4,15,u_data,vmin=np.min(u_data),vmax=np.max(u_data)) 
        fig.colorbar(hu)
        
        plt.subplot(3,2,3)
        he = plt.scatter(insar_xyz[:,0]/3e4,insar_xyz[:,1]/3e4,15,dpred[0:ndata],vmin=np.min(e_data),vmax=np.max(e_data))  
        fig.colorbar(he)
        plt.subplot(3,2,4)
        hu = plt.scatter(insar_xyz[:,0]/3e4,insar_xyz[:,1]/3e4,15,dpred[ndata:2*ndata],vmin=np.min(u_data),vmax=np.max(u_data)) 
        fig.colorbar(hu)
        
        plt.subplot(3,2,5)
        hr1 = plt.scatter(insar_xyz[:,0]/3e4,insar_xyz[:,1]/3e4,15,e_data-dpred[0:ndata])  
        fig.colorbar(hr1)  
        plt.subplot(3,2,6)
        hr3 = plt.scatter(insar_xyz[:,0]/3e4,insar_xyz[:,1]/3e4,15,u_data-dpred[ndata:2*ndata])    

else:
    if datatype == 'los':
        for i in range(nInsarFrames):
            gu.plotQuads(cx[:,nquads[i-1]*i:nquads[i]+i*nquads[i-1]]/1e3,
                         cy[:,nquads[i-1]*i:nquads[i]+i*nquads[i-1]]/1e3,
                         alldata[nquads[i-1]*i:nquads[i]+i*nquads[i-1]],
                         np.min(alldata), np.max(alldata), 'Greens')
            gu.plotQuads(cx[:,nquads[i-1]*i:nquads[i]+i*nquads[i-1]]/1e3,
                         cy[:,nquads[i-1]*i:nquads[i]+i*nquads[i-1]]/1e3,
                         dpred[nquads[i-1]*i:nquads[i]+i*nquads[i-1]], np.min(alldata), np.max(alldata), 'Greens')

#fig = plt.figure()
#h2 = plt.scatter(insar_xyz[:,0],insar_xyz[:,1],15,e_data-dpred[0:ndata])
#fig.colorbar(h2)


#    fig = plt.figure()
#    h2 = plt.scatter(insar_xyz[:,0],insar_xyz[:,1],15,ue)
#    fig.colorbar(h2)


plt.show()

#fig = plt.figure()
#ax = Axes3D(fig)
#h2 = ax.scatter(xsubgrid.reshape((xsubgrid.size)),ysubgrid.reshape((ysubgrid.size)),zsubgrid.reshape((zsubgrid.size)))
