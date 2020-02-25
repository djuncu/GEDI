#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 16:44:29 2019

GDIS module for loading geodetic data

@author: Daniel Juncu
"""
import numpy as np
import GediUtils as gu

def loadInsarData(datalocation,projOrigin,estimateOrigin=False):
    
    disp = []
    cosN = []
    cosE = []
    cosU = []
    lonlat = np.empty((0, 2))
    xyz = np.empty((0, 3))
    #ndata_all = 0
    ndata = []
    
    for i in range(len(datalocation)):
    
        data = np.loadtxt(datalocation[i])   
        
        ndata.append(data.shape[0])        
        #ndata = data.shape[0] 
        ncol = data.shape[1]
        
        # optional: automatically place origin in scene center
        if estimateOrigin: 
            westLim = np.min(data[:,0])
            eastLim = np.max(data[:,0])
            southLim = np.min(data[:,1])
            northLim = np.max(data[:,1])    
            
            projOrigin = np.array(([np.mean([westLim, eastLim]), np.mean([southLim, northLim]) ]))       
            
        if ncol == 7:
            xyz_i = np.hstack((gu.wgs2local(data[:,0:2],projOrigin),data[:,6][:,None]))
            topo = True       
        elif ncol == 6:
            xyz_i = np.hstack((gu.wgs2local(data[:,0:2],projOrigin), np.zeros((ndata,1)) ))
            topo = False       
        else:
            raise ValueError('Check input file format')

        disp = np.append(disp,data[:,2])
            
        # need to switch order in input file....
        cosN = np.append(cosN,data[:,3])
        cosE = np.append(cosE,data[:,4])
        cosU = np.append(cosU,data[:,5])   
        
        lonlat = np.append(lonlat,data[:,0:2],axis=0)
        xyz    = np.append(xyz,xyz_i,axis=0)
        

    
    return disp, xyz, lonlat, cosN, cosE, cosU, ndata, projOrigin, topo
    
def loadInsarData_eu(datalocation,projOrigin,estimateOrigin=False):
    data = np.loadtxt(datalocation)   
    
    ndata = data.shape[0]  
    ncol = data.shape[1]
    
    # optional: automatically place origin in scene center
    if estimateOrigin: 
        westLim = np.min(data[:,0])
        eastLim = np.max(data[:,0])
        southLim = np.min(data[:,1])
        northLim = np.max(data[:,1])    
        
        projOrigin = np.array(([np.mean([westLim, eastLim]), np.mean([southLim, northLim]) ]))       
        
    if ncol == 5:
        xyz = np.hstack((gu.wgs2local(data[:,0:2],projOrigin),data[:,6][:,None]))
        topo = True
        e_disp = data[:,3]
        u_disp = data[:,4]       
        
    elif ncol == 4:
        xyz = np.hstack((gu.wgs2local(data[:,0:2],projOrigin), np.zeros((ndata,1)) ))
        topo = False
        e_disp = data[:,2]
        u_disp = data[:,3]
    else:
        raise ValueError('Check input file format')    

    return e_disp, u_disp, xyz, data[:,0:2], ndata, projOrigin, topo

def loadInsarData_enu(datalocation,projOrigin,estimateOrigin=False):
    data = np.loadtxt(datalocation)   
    
    ndata = data.shape[0]  
    ncol = data.shape[1]
    
    # optional: automatically place origin in scene center
    if estimateOrigin: 
        westLim = np.min(data[:,0])
        eastLim = np.max(data[:,0])
        southLim = np.min(data[:,1])
        northLim = np.max(data[:,1])    
        
        projOrigin = np.array(([np.mean([westLim, eastLim]), np.mean([southLim, northLim]) ]))       
        
    if ncol == 6:
        xyz = np.hstack((gu.wgs2local(data[:,0:2],projOrigin),data[:,2][:,None]))
        topo = True
        e_disp = data[:,3]
        n_disp = data[:,4]
        u_disp = data[:,5]       
        
    elif ncol == 5:
        xyz = np.hstack((gu.wgs2local(data[:,0:2],projOrigin), np.zeros((ndata,1)) ))
        topo = False
        e_disp = data[:,2]
        n_disp = data[:,3]
        u_disp = data[:,4]
    else:
        raise ValueError('Check input file format')
    
   
    return e_disp, n_disp, u_disp, xyz, data[:,0:2], ndata, projOrigin, topo

def loadInsarData_old(datalocation,inc_file,projOrigin,estimateOrigin=False):
    data = np.loadtxt(datalocation)
    inc = np.loadtxt(inc_file)
    
    ndata = data.shape[0]  
    ncol = data.shape[1]
    
    if estimateOrigin: # optional: automatically place origin in scene center
        westLim = np.min(data[:,0])
        eastLim = np.max(data[:,0])
        southLim = np.min(data[:,1])
        northLim = np.max(data[:,1])    
        
        projOrigin = np.array(([np.mean([westLim, eastLim]), np.mean([southLim, northLim]) ]))       
        
    
    if ncol == 4:
        xyz = gu.wgs2local(data[:,0:3],projOrigin) 
        disp = data[:,3]
    elif ncol == 3:
        xyz = np.hstack((gu.wgs2local(data[:,0:2],projOrigin), np.zeros((ndata,1)) ))
        disp = data[:,2] 
    
    
    return data, xyz, inc, ndata, projOrigin
    
