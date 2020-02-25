#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 16:30:00 2019

Calculate forward deformation model

@author: Daniel Juncu
"""

import numpy as np

def forwardmodel_los(nu,datalocs,ngrid, xsubgrid, ysubgrid, zsubgrid, subgridlvl,vc,nx,ny,nz, cosE, cosN, cosU, straintype='constrained'):

    ux = np.zeros((datalocs.shape[0]))
    uy = np.zeros((datalocs.shape[0]))
    uz = np.zeros((datalocs.shape[0]))   
    ulos = np.zeros((datalocs.shape[0]))   

    for i in range(ngrid):
                
                if straintype.lower() == 'constrained':
                    strainfac = (1 - nu) / np.pi  
                elif straintype.lower() == 'free':
                    strainfac = (nu + 1) / (3.*np.pi)
                else:
                    raise ValueError('invalid straintype')
                         
                    
                gx = np.zeros((datalocs.shape[0]))
                gy = np.zeros((datalocs.shape[0]))
                gz = np.zeros((datalocs.shape[0]))  
                        
                for l in range((subgridlvl+1)**3):
                    S = np.sqrt((datalocs[:,0] - xsubgrid[i,l])**2 + (datalocs[:,1] - ysubgrid[i,l])**2 + (datalocs[:,2] - zsubgrid[i,l])**2)**3
                                                            
                    gx +=  (( strainfac * (datalocs[:,0] - xsubgrid[i,l])) / S / (subgridlvl+1.)**3)                    
                    gy +=  (( strainfac * (datalocs[:,1] -  ysubgrid[i,l])) / S / (subgridlvl+1.)**3)                  
                    gz +=  (( strainfac * (datalocs[:,2] - zsubgrid[i,l])) / S/  (subgridlvl+1.)**3) 
               
                
                ux = gx*vc[i] 
                uy = gy*vc[i] 
                uz = gz*vc[i]  

                #ulos = (ux*east_look + uy*north_look) * np.sin(inc) + uz*np.cos(inc) +ulos 
                
                # new method
                ulos = ux*cosE + uy*cosN + uz*cosU  + ulos              
                
                         
    return ulos

def forwardmodel_enu(nu,datalocs,ngrid, xsubgrid, ysubgrid, zsubgrid, nlvl, vc,nx,ny,nz, straintype='constrained'):
    gx = np.zeros((datalocs.shape[0]))
    gy = np.zeros((datalocs.shape[0]))
    gz = np.zeros((datalocs.shape[0]))  

    ux = np.zeros((datalocs.shape[0]))
    uy = np.zeros((datalocs.shape[0]))
    uz = np.zeros((datalocs.shape[0]))   

    #for i in range(0,datalocs.shape[0]):
    for i in range(ngrid): 
                
                if straintype.lower() == 'constrained':
                    strainfac = (1 - nu) / np.pi  
                elif straintype.lower() == 'free':
                    strainfac = (nu + 1) / (3.*np.pi)
                else:
                    raise ValueError('invalid straintype')
                                         
                gx = np.zeros((datalocs.shape[0]))
                gy = np.zeros((datalocs.shape[0]))
                gz = np.zeros((datalocs.shape[0]))  

                S = np.sqrt((datalocs[:,0][:,None] - xsubgrid[i,:][None,:])**2 \
                                + (datalocs[:,1][:,None] - ysubgrid[i,:][None,:])**2 \
                                + (datalocs[:,2][:,None] - zsubgrid[i,:][None,:])**2)**3
                                
                gx =  np.sum(( strainfac* (datalocs[:,0][:,None] - xsubgrid[i,:][None,:]) / S / (nlvl+1)**3), axis=1)                       
                gy =  np.sum(( strainfac* (datalocs[:,1][:,None] - ysubgrid[i,:][None,:]) / S / (nlvl+1)**3), axis=1)                       
                gz =  np.sum(( strainfac* (datalocs[:,2][:,None] - zsubgrid[i,:][None,:]) / S / (nlvl+1)**3), axis=1)

               
                ux = gx*vc[i] + ux
                uy = gy*vc[i] + uy
                uz = gz*vc[i] + uz                             

    return ux, uy, uz

def forwardmodel_enu_nosubgrid(nu,datalocs,xgrid,ygrid,zgrid,vc,nx,ny,nz, straintype='constrained'):
    gx = np.zeros((datalocs.shape[0]))
    gy = np.zeros((datalocs.shape[0]))
    gz = np.zeros((datalocs.shape[0]))  

    ux = np.zeros((datalocs.shape[0]))
    uy = np.zeros((datalocs.shape[0]))
    uz = np.zeros((datalocs.shape[0]))   

    for j in range(0,nx):
        for i in range(0,ny):
            for k in range(0,nz):  
                
                if straintype.lower() == 'constrained':
                    strainfac = (1 - nu) / np.pi  
                elif straintype.lower() == 'free':
                    strainfac = (nu + 1) / (3.*np.pi)
                else:
                    raise ValueError('invalid straintype')
                                         

                gx =  (( strainfac * (datalocs[:,0] - xgrid[k,i,j])) / \
                  (np.sqrt((datalocs[:,0] - xgrid[k,i,j])**2 + (datalocs[:,1] - ygrid[k,i,j])**2 + (datalocs[:,2] - zgrid[k,i,j])**2)**3))  
                ux = gx*vc[k,i,j] + ux
                  
                gy =  (( strainfac* (datalocs[:,1] - ygrid[k,i,j])) / \
                  (np.sqrt((datalocs[:,0] - xgrid[k,i,j])**2 + (datalocs[:,1] - ygrid[k,i,j])**2 + (datalocs[:,2] - zgrid[k,i,j])**2)**3))  
                uy = gy*vc[k,i,j] + uy
                  
                gz =  (( strainfac * (datalocs[:,2] - zgrid[k,i,j])) / \
                  (np.sqrt((datalocs[:,0] - xgrid[k,i,j])**2 + (datalocs[:,1] - ygrid[k,i,j])**2 + (datalocs[:,2] - zgrid[k,i,j])**2)**3)) 
                uz = gz*vc[k,i,j] + uz                             

    return ux, uy, uz
