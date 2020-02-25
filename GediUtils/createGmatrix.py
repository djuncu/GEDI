#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 16:37:17 2019

GDIS module to create G matrix for inversion

@author: Daniel Juncu
"""

import numpy as np

def createPmatrix(ngrid,P,Ptype,gridbound,xleft,xright,yfront,yback,ztop,zbot,nx,ny,nz,dx,dy,dz,xgrid,ygrid,zgrid):                         
                    
    if Ptype == '2bounds-2D':
        for i in range(0,ngrid):
                P[i,i] = -2./dx**2. -2./dy**2. 
                if not gridbound[i]:  
                    P[i,i+nz*ny] = 1./dx**2.  
                    P[i,i-nz*ny] = 1./dx**2.  
                    P[i,i+nz] =    1./dy**2.
                    P[i,i-nz] =    1./dy**2.
                else:
                    if xleft[i] and yfront[i]:
                        P[i,i+nz*ny] = 1./dx**2.    
                        P[i,i+nz] =    1./dy**2.
                    elif xright[i] and yfront[i]:
                        P[i,i-nz*ny] = 1./dx**2.    
                        P[i,i+nz] =    1./dy**2.
                    elif xright[i] and yback[i]:
                        P[i,i-nz*ny] = 1./dx**2.    
                        P[i,i-nz] =    1./dy**2.
                    elif xleft[i] and yback[i]:
                        P[i,i+nz*ny] = 1./dx**2.    
                        P[i,i-nz] =    1./dy**2.                  
                        
                    elif xleft[i]:
                        P[i,i+nz*ny] = 1./dx**2.                    
                        P[i,i+nz] =    1./dy**2.
                        P[i,i-nz] =    1./dy**2.
                    elif xright[i]: 
                        P[i,i-nz*ny] = 1./dx**2.  
                        P[i,i+nz] =    1./dy**2.
                        P[i,i-nz] =    1./dy**2.
                    elif yfront[i]:
                        P[i,i+nz*ny] = 1./dx**2.  
                        P[i,i-nz*ny] = 1./dx**2.  
                        P[i,i+nz] =    1./dy**2.
                    elif yback[i]:
                        P[i,i+nz*ny] = 1./dx**2.  
                        P[i,i-nz*ny] = 1./dx**2.  
                        P[i,i-nz] =    1./dy**2.
        
        
    if Ptype == '2bounds':        
        for i in range(0,ngrid):
                P[i,i] = -2./dx**2. -2./dy**2. - 2./dz**2.
                if not gridbound[i]:  
                    P[i,i+nz*ny] = 1./dx**2.  
                    P[i,i-nz*ny] = 1./dx**2.  
                    P[i,i+nz] =    1./dy**2.
                    P[i,i-nz] =    1./dy**2.
                    P[i,i+1] =     1./dz**2.  
                    P[i,i-1] =     1./dz**2. 
                else:
                    if xleft[i] and yfront[i] and ztop[i]:
#                        P[i,i] = -1./dx**2. -1./dy**2. - 1./dz**2.
                        P[i,i+nz*ny] = 1./dx**2.    
                        P[i,i+nz] =    1./dy**2.
                        P[i,i+1] =     1./dz**2.  
                    elif xright[i] and yfront[i] and ztop[i]:
#                        P[i,i] = -1./dx**2. -1./dy**2. - 1./dz**2.
                        P[i,i-nz*ny] = 1./dx**2.    
                        P[i,i+nz] =    1./dy**2.
                        P[i,i+1] =     1./dz**2.      
                    elif xright[i] and yback[i] and ztop[i]:
#                        P[i,i] = -1./dx**2. -1./dy**2. - 1./dz**2.
                        P[i,i-nz*ny] = 1./dx**2.    
                        P[i,i-nz] =    1./dy**2.
                        P[i,i+1] =     1./dz**2.  
                    elif xleft[i] and yback[i] and ztop[i]:
#                        P[i,i] = -1./dx**2. -1./dy**2. - 1./dz**2.
                        P[i,i+nz*ny] = 1./dx**2.    
                        P[i,i-nz] =    1./dy**2.
                        P[i,i+1] =     1./dz**2.  
                    elif xleft[i] and yfront[i] and zbot[i]:
#                        P[i,i] = -1./dx**2. -1./dy**2. - 1./dz**2.
                        P[i,i+nz*ny] = 1./dx**2.    
                        P[i,i+nz] =    1./dy**2.
                        P[i,i-1] =     1./dz**2.  
                    elif xright[i] and yfront[i] and zbot[i]:
#                        P[i,i] = -1./dx**2. -1./dy**2. - 1./dz**2.
                        P[i,i-nz*ny] = 1./dx**2.    
                        P[i,i+nz] =    1./dy**2.
                        P[i,i-1] =     1./dz**2. 
                    elif xright[i] and yback[i] and zbot[i]:
#                        P[i,i] = -1./dx**2. -1./dy**2. - 1./dz**2.
                        P[i,i-nz*ny] = 1./dx**2.    
                        P[i,i-nz] =    1./dy**2.
                        P[i,i-1] =     1./dz**2. 
                    elif xleft[i] and yback[i] and zbot[i]:
#                        P[i,i] = -1./dx**2. -1./dy**2. - 1./dz**2.
                        P[i,i+nz*ny] = 1./dx**2.    
                        P[i,i-nz] =    1./dy**2.
                        P[i,i-1] =     1./dz**2. 
                    elif xleft[i] and yfront[i]:
                        P[i,i+nz*ny] = 1./dx**2.    
                        P[i,i+nz] =    1./dy**2.
                        P[i,i+1] =     1./dz**2.
                        P[i,i-1] =     1./dz**2. 
                    elif xright[i] and yfront[i]:
                        P[i,i-nz*ny] = 1./dx**2.    
                        P[i,i+nz] =    1./dy**2.
                        P[i,i+1] =     1./dz**2.
                        P[i,i-1] =     1./dz**2. 
                    elif xright[i] and yback[i]:
                        P[i,i-nz*ny] = 1./dx**2.    
                        P[i,i-nz] =    1./dy**2.
                        P[i,i+1] =     1./dz**2.
                        P[i,i-1] =     1./dz**2. 
                    elif xleft[i] and yback[i]:
                        P[i,i+nz*ny] = 1./dx**2.    
                        P[i,i-nz] =    1./dy**2.
                        P[i,i+1] =     1./dz**2.
                        P[i,i-1] =     1./dz**2. 
                    elif xleft[i] and ztop[i]:     
                        P[i,i+nz*ny] = 1./dx**2.    
                        P[i,i+nz] =    1./dy**2.
                        P[i,i-nz] =    1./dy**2.
                        P[i,i+1] =     1./dz**2. 
                    elif xright[i] and ztop[i]:     
                        P[i,i-nz*ny] = 1./dx**2.    
                        P[i,i+nz] =    1./dy**2.
                        P[i,i-nz] =    1./dy**2.
                        P[i,i+1] =     1./dz**2. 
                    elif yfront[i] and ztop[i]:     
                        P[i,i+nz*ny] = 1./dx**2.  
                        P[i,i-nz*ny] = 1./dx**2.
                        P[i,i+nz] =    1./dy**2.
                        P[i,i+1] =     1./dz**2. 
                    elif yback[i] and ztop[i]:     
                        P[i,i+nz*ny] = 1./dx**2.  
                        P[i,i-nz*ny] = 1./dx**2.
                        P[i,i-nz] =    1./dy**2.
                        P[i,i+1] =     1./dz**2.                     
                    elif xleft[i] and zbot[i]:     
                        P[i,i+nz*ny] = 1./dx**2.    
                        P[i,i+nz] =    1./dy**2.
                        P[i,i-nz] =    1./dy**2.
                        P[i,i-1] =     1./dz**2. 
                    elif xright[i] and zbot[i]:     
                        P[i,i-nz*ny] = 1./dx**2.    
                        P[i,i+nz] =    1./dy**2.
                        P[i,i-nz] =    1./dy**2.
                        P[i,i-1] =     1./dz**2. 
                    elif yfront[i] and zbot[i]:     
                        P[i,i+nz*ny] = 1./dx**2.  
                        P[i,i-nz*ny] = 1./dx**2.
                        P[i,i+nz] =    1./dy**2.
                        P[i,i-1] =     1./dz**2. 
                    elif yback[i] and zbot[i]:     
                        P[i,i+nz*ny] = 1./dx**2.  
                        P[i,i-nz*ny] = 1./dx**2.
                        P[i,i-nz] =    1./dy**2.
                        P[i,i-1] =     1./dz**2. 
                        
                    elif xleft[i]:
                        P[i,i+nz*ny] = 1./dx**2.                    
                        P[i,i+nz] =    1./dy**2.
                        P[i,i-nz] =    1./dy**2.
                        P[i,i+1] =     1./dz**2.  
                        P[i,i-1] =     1./dz**2. 
                    elif xright[i]: 
                        P[i,i-nz*ny] = 1./dx**2.  
                        P[i,i+nz] =    1./dy**2.
                        P[i,i-nz] =    1./dy**2.
                        P[i,i+1] =     1./dz**2.  
                        P[i,i-1] =     1./dz**2. 
                    elif yfront[i]:
                        P[i,i+nz*ny] = 1./dx**2.  
                        P[i,i-nz*ny] = 1./dx**2.  
                        P[i,i+nz] =    1./dy**2.
                        P[i,i+1] =     1./dz**2.  
                        P[i,i-1] =     1./dz**2. 
                    elif yback[i]:
                        P[i,i+nz*ny] = 1./dx**2.  
                        P[i,i-nz*ny] = 1./dx**2.  
                        P[i,i-nz] =    1./dy**2.
                        P[i,i+1] =     1./dz**2.  
                        P[i,i-1] =     1./dz**2. 
                    elif ztop[i]:
                        P[i,i+nz*ny] = 1./dx**2.  
                        P[i,i-nz*ny] = 1./dx**2.  
                        P[i,i+nz] =    1./dy**2.
                        P[i,i-nz] =    1./dy**2.
                        P[i,i+1] =     1./dz**2.  
                    elif zbot[i]:
                        P[i,i+nz*ny] = 1./dx**2.  
                        P[i,i-nz*ny] = 1./dx**2.  
                        P[i,i+nz] =    1./dy**2.
                        P[i,i-nz] =    1./dy**2.
                        P[i,i-1] =     1./dz**2. 

#        G = np.concatenate((g,P*np.sqrt(Wr)))
##            d = np.concatenate((d,np.zeros(ngrid)),axis=0)
#        d = np.concatenate((d,np.zeros(ngrid-outerCells)),axis=0)        

    if Ptype == '2bounds-b':
        j=0        
        for i in range(0,ngrid):                
                if not gridbound[i]:  
                    P[j,i] = -2./dx**2. -2./dy**2. - 2./dz**2.
                    P[j,i+nz*ny] = 1./dx**2.  
                    P[j,i-nz*ny] = 1./dx**2.  
                    P[j,i+nz] =    1./dy**2.
                    P[j,i-nz] =    1./dy**2.
                    P[j,i+1] =     1./dz**2.  
                    P[j,i-1] =     1./dz**2. 
                    j+=1
                else:
                    if xleft[i] and yfront[i] and ztop[i]:
                        continue
                    elif xright[i] and yfront[i] and ztop[i]:
                        continue   
                    elif xright[i] and yback[i] and ztop[i]:
                        continue
                    elif xleft[i] and yback[i] and ztop[i]:
                        continue 
                    elif xleft[i] and yfront[i] and zbot[i]:
                        continue  
                    elif xright[i] and yfront[i] and zbot[i]:
                        continue
                    elif xright[i] and yback[i] and zbot[i]:
                        continue
                    elif xleft[i] and yback[i] and zbot[i]:
                        continue 
                    elif (xleft[i] and yfront[i]) or (xright[i] and yfront[i]) \
                         or (xright[i] and yback[i]) or (xleft[i] and yfront[i]):
                        P[j,i] = - 2./dz**2.
                        P[j,i+1] =     1./dz**2.
                        P[j,i-1] =     1./dz**2. 
                        j+=1                  
                    elif (xleft[i] and ztop[i]) or (xright[i] and ztop[i]) \
                         or (xleft[i] and zbot[i]) or (xright[i] and zbot[i]):
                        P[j,i] = - 2./dy**2.  
                        P[j,i+nz] =    1./dy**2.
                        P[j,i-nz] =    1./dy**2.
                        j+=1          
                        
                    elif (yfront[i] and ztop[i]) or (yback[i] and ztop[i]) \
                         or (yfront[i] and zbot[i]) or (yback[i] and zbot[i]):     
                        P[j,i] = - 2./dx**2.  
                        P[j,i+nz*ny] = 1./dx**2.  
                        P[j,i-nz*ny] = 1./dx**2.
                        j+=1                                    
                    elif xleft[i] or xright[i]:
                        P[j,i] = - 2./dy**2. - 2./dz**2.                  
                        P[j,i+nz] =    1./dy**2.
                        P[j,i-nz] =    1./dy**2.
                        P[j,i+1] =     1./dz**2.  
                        P[j,i-1] =     1./dz**2. 
                        j+=1
                    elif yfront[i] or yback[i]:
                        P[j,i] = - 2./dx**2. - 2./dz**2. 
                        P[j,i+nz*ny] = 1./dx**2.  
                        P[j,i-nz*ny] = 1./dx**2.   
                        P[j,i+1] =     1./dz**2.  
                        P[j,i-1] =     1./dz**2. 
                        j+=1
                    elif ztop[i] or zbot[i]:
                        P[j,i] = - 2./dx**2. - 2./dy**2. 
                        P[j,i+nz*ny] = 1./dx**2.  
                        P[j,i-nz*ny] = 1./dx**2.  
                        P[j,i+nz] =    1./dy**2.
                        P[j,i-nz] =    1./dy**2. 
                        j+=1

    if Ptype == '1bounds-b':
        j=0        
        for i in range(0,ngrid):                
                if not gridbound[i]:  
                    P[j,i] = -1./dx -1./dy - 1./dz
                    P[j,i+nz*ny] = 1./dx   
                    P[j,i+nz] =    1./dy
                    P[j,i+1] =     1./dz 
                    j+=1
                else:
                    if xleft[i] and yfront[i] and ztop[i]:
                        P[j,i] = -1./dx -1./dy - 1./dz
                        P[j,i+nz*ny] = 1./dx   
                        P[j,i+nz] =    1./dy
                        P[j,i+1] =     1./dz 
                        j+=1
                    elif xright[i] and yfront[i] and ztop[i]:
                        P[j,i] =  -1./dy - 1./dz  
                        P[j,i+nz] =    1./dy
                        P[j,i+1] =     1./dz 
                        j+=1  
                    elif xright[i] and yback[i] and ztop[i]:
                        P[j,i] = - 1./dz
                        P[j,i+1] =     1./dz 
                        j+=1
                    elif xleft[i] and yback[i] and ztop[i]:
                        P[j,i] = -1./dx - 1./dz
                        P[j,i+nz*ny] = 1./dx 
                        P[j,i+1] =     1./dz 
                        j+=1
                    elif xleft[i] and yfront[i] and zbot[i]:
                        P[j,i] = -1./dx -1./dy 
                        P[j,i+nz*ny] = 1./dx   
                        P[j,i+nz] =    1./dy
                        j+=1
                    elif xright[i] and yfront[i] and zbot[i]:
                        P[j,i] = -1./dy  
                        P[j,i+nz] =    1./dy   
                        j+=1
                    elif xright[i] and yback[i] and zbot[i]:
                        continue
                    elif xleft[i] and yback[i] and zbot[i]:
                        P[j,i] = -1./dx 
                        P[j,i+nz*ny] = 1./dx          
                        j+=1
                    elif (xleft[i] and yfront[i]):
                        P[j,i] = -1./dx -1./dy - 1./dz
                        P[j,i+nz*ny] = 1./dx   
                        P[j,i+nz] =    1./dy
                        P[j,i+1] =     1./dz 
                        j+=1                        
                    elif (xright[i] and yfront[i]):
                        P[j,i] =  -1./dy - 1./dz      
                        P[j,i+nz] =    1./dy
                        P[j,i+1] =     1./dz 
                        j+=1                         
                    elif (xright[i] and yback[i]):
                        P[j,i] =  - 1./dz
                        P[j,i+1] =     1./dz 
                        j+=1                           
                    elif (xleft[i] and yfront[i]):
                        P[j,i] = -1./dx - 1./dz
                        P[j,i+nz*ny] = 1./dx   
                        P[j,i+1] =     1./dz 
                        j+=1                   
                    elif (xleft[i] and ztop[i]):
                        P[j,i] = -1./dx -1./dy - 1./dz
                        P[j,i+nz*ny] = 1./dx   
                        P[j,i+nz] =    1./dy
                        P[j,i+1] =     1./dz 
                        j+=1                            
                    elif (xright[i] and ztop[i]):
                        P[j,i] = -1./dy - 1./dz 
                        P[j,i+nz] =    1./dy
                        P[j,i+1] =     1./dz 
                        j+=1    
                    elif (xleft[i] and zbot[i]):
                        P[j,i] = -1./dx -1./dy 
                        P[j,i+nz*ny] = 1./dx   
                        P[j,i+nz] =    1./dy 
                        j+=1    
                    elif (xright[i] and zbot[i]):
                        P[j,i] =-1./dy
                        P[j,i+nz] =    1./dy
                        j+=1                            
                    elif (yfront[i] and ztop[i]):
                        P[j,i] = -1./dx -1./dy - 1./dz
                        P[j,i+nz*ny] = 1./dx   
                        P[j,i+nz] =    1./dy
                        P[j,i+1] =     1./dz 
                        j+=1                         
                    elif (yback[i] and ztop[i]):
                        P[j,i] = -1./dx - 1./dz
                        P[j,i+nz*ny] = 1./dx                       
                        P[j,i+1] =     1./dz 
                        j+=1                           
                    elif (yfront[i] and zbot[i]):
                        P[j,i] = -1./dx -1./dy
                        P[j,i+nz*ny] = 1./dx   
                        P[j,i+nz] =    1./dy            
                        j+=1    
                    elif (yback[i] and zbot[i]):     
                        P[j,i] = -1./dx
                        P[j,i+nz*ny] = 1./dx  
                        j+=1                              
                    elif xright[i]:
                        P[j,i] = -1./dy - 1./dz                
                        P[j,i+nz] =    1./dy
                        P[j,i+1] =     1./dz 
                        j+=1                        
                    elif yback[i]:
                        P[j,i] = -1./dx  - 1./dz
                        P[j,i+nz*ny] = 1./dx   
                        P[j,i+1] =     1./dz 
                        j+=1    
                    elif zbot[i]:
                        P[j,i] = -1./dx -1./dy - 1./dz
                        P[j,i+nz*ny] = 1./dx   
                        P[j,i+nz] =    1./dy
                        j+=1    
                    elif xleft[i] or yfront[i] or ztop[i]:
                        P[j,i] = -1./dx -1./dy - 1./dz
                        P[j,i+nz*ny] = 1./dx   
                        P[j,i+nz] =    1./dy
                        P[j,i+1] =     1./dz 
                        j+=1    
                        

    ## OLD VERSION UNTIL 22.10.19
    ## REGULARIZATION
    if Ptype == '2simple':
        index = 0
        # second derivativePtype
        for i in range(0,ngrid):
            #print(gridbound[i])
            if not gridbound[i]:  
                P[index,i] = -2./dx**2. -2./dy**2. - 2./dz**2.   
                #print(P)
                P[index,i+nz*ny] = 1./dx**2.  
                P[index,i-nz*ny] = 1./dx**2.  
                P[index,i+nz] =    1./dy**2.
                P[index,i-nz] =    1./dy**2.
                P[index,i+1] =     1./dz**2.  
                P[index,i-1] =     1./dz**2. 
                index+=1

            else:
                if xleft[i]:
#                    P[i,i+nz*ny] = 1./dx
                    P[i,i] -= 0 #1./dx
                elif xright[i]:
#                    P[i,i-nz*ny] = 1./dx
                    P[i,i] -= 0#1./dx
                else:                                     
                    P[i,i+nz*ny] = 1./dx**2.  
                    P[i,i-nz*ny] = 1./dx**2. 
                    P[i,i] += -2./dx**2.              
#                    print('hello')
                
                if yfront[i]:
#                    P[i,i+nz] = 1./dy
                    P[i,i] -= 0#1./dy
                elif yback[i]:
#                    P[i,i-nz] = 1./dy
                    P[i,i] -= 0#1./dy
                else:
                    P[i,i+nz] = 1./dy**2.
                    P[i,i-nz] = 1./dy**2.
                    P[i,i]   += -2./dy**2. 
                    
                if ztop[i]:
#                    P[i,i+1] = 1./dz
                    P[i,i] -= 0#1./dz
                elif zbot[i]:
#                    P[i,i-1] = 1./dz
                    P[i,i] -= 0#1./dz    
                else:
                    P[i,i+1] = 1./dz**2.  
                    P[i,i-1] = 1./dz**2. 
                    P[i,i]   += -2./dz**2. 
                         
    return P                    


def createGmatrix(nu,datalocs,xgrid,ygrid,zgrid, xsubgrid, ysubgrid, zsubgrid, nlvl, nx,ny,nz,dz, d,ndata,datatype,Wr =0,cosE=0, cosN=0, cosU=0, regularization=False, straintype='constrained'):
    ngrid = nx*ny*nz
    datatype = datatype.lower()
    
    SM = np.zeros((xgrid.shape))
    
    dx = xgrid[0,0,1]-xgrid[0,0,0]
    dy = ygrid[0,1,0]-ygrid[0,0,0]
#    dz = zgrid[1,0,0]-zgrid[0,0,0]
    
    if datatype == 'enu':        
        gx = np.zeros((ngrid,ndata), dtype=np.single)
        gy = np.zeros((ngrid,ndata), dtype=np.single)
        gz = np.zeros((ngrid,ndata), dtype=np.single)
    elif datatype == 'eu':
        gx = np.zeros((ngrid,ndata), dtype=np.single)
        gz = np.zeros((ngrid,ndata), dtype=np.single)
    elif datatype == 'los':
        glos = np.zeros((ngrid,ndata))

#    xbound = np.zeros((ncells))
#    ybound = np.zeros((ncells))
#    zbound = np.zeros((ncells))

    xleft = np.zeros((ngrid), dtype=bool)
    xright = np.zeros((ngrid), dtype=bool)
    yfront = np.zeros((ngrid), dtype=bool)
    yback = np.zeros((ngrid), dtype=bool)
    zbot = np.zeros((ngrid), dtype=bool)
    ztop = np.zeros((ngrid), dtype=bool)
    
    gridbound = np.zeros((ngrid), dtype=bool)
    
    xyz = np.zeros((ngrid,3))
    
    if datatype == 'los':  g = np.zeros((ndata,ngrid), dtype=np.single)
    if datatype == 'eu':   g = np.zeros((2*ndata,ngrid), dtype=np.single)
    if datatype == 'enu':  g = np.zeros((3*ndata,ngrid), dtype=np.single)

#    P = np.zeros((ngrid,ngrid), dtype=np.single)
    Ptype = '2bounds-2D'
    
    if Ptype == 'nobounds':
        ## 2nd derivative with out specific boundary terms 
        outerCells = 2*nx*ny+ 2*nx*(nz-2) + 2*(ny-2)*(nz-2)
    elif Ptype == 'nobounds-1':
        ## first derivative without specific boundary terms
        outerCells = nx*ny + nx*(nz-1) + (ny-1)*(nz-1)
    elif Ptype == '2bounds':
        ## square matrix
        outerCells=0
    elif Ptype == '2bounds-b':
        outerCells = 8      
    elif Ptype == '1bounds-b':
        outerCells = 1
    elif Ptype == 'priorlocs':
        outerCells=0
    elif Ptype == '2bounds-2D':
        # no vertical regularization
        outerCells = 0
        
    P = np.zeros((ngrid-outerCells,ngrid), dtype=np.single)
    
    
    # index system for all dimensions combined
    index=0
    
    for j in range(0,nx):
        for i in range(0,ny):
            for k in range(0,nz):                
                if straintype.lower() == 'constrained':
                    strainfac = (1 - nu) / np.pi  
                elif straintype.lower() == 'free':
                    strainfac = (nu + 1) / (3.*np.pi)
                else:
                    raise ValueError('invalid straintype')
                
                # boundaries with Dirichlet BC (dV=0)
#                if i == nx-1 or i == 0:
#                    xbound[index] = 1
#                if j == ny-1 or j == 0:
#                    ybound[index] = 1
#                if k == nz-1 or :
#                    zbound[index] = 1                         

#                if i == nx-1 or i == 0 or j == ny-1 or j == 0 or k == nz-1:
#                    gridbound[index] = 1   
#                #elif i == nx-2 or i==1 or j == ny-2 or j==1 or k == nz-2 or k == 1: # for 2nd derivative
#                elif k == 1:
#                    topbound[index] = 1
                    
                if j == nx-1 or j == 0 or i == ny-1 or i == 0 or k == nz-1 or k == 0:
                    gridbound[index] = True   
                    
                    if j == 0:
                        xleft[index] = True
                    if j == nx-1:
                        xright[index] = True
                    if i == 0:
                        yfront[index] = True
                    if i == ny-1:
                        yback[index] = True
                    if k == 0:
                        ztop[index] = True
                    if k == nz-1:
                        zbot[index] = True                

#                if i == nx-1  or j == ny-1 or k == nz-1:
#                    gridbound[index] = 1                                         

  

         
                gxtmp = np.zeros((ndata))
                gytmp = np.zeros((ndata))
                gztmp = np.zeros((ndata))


                S = np.sqrt((datalocs[:,0][:,None] - xsubgrid[index,:][None,:])**2 \
                                + (datalocs[:,1][:,None] - ysubgrid[index,:][None,:])**2 \
                                + (datalocs[:,2][:,None] - zsubgrid[index,:][None,:])**2)**3
                            
                SM[k,i,j] = S[0,0]
        
                gxtmp =  np.sum(( strainfac* (datalocs[:,0][:,None] - xsubgrid[index,:][None,:]) / S / (nlvl+1)**3), axis=1)                       
                gytmp =  np.sum(( strainfac* (datalocs[:,1][:,None] - ysubgrid[index,:][None,:]) / S / (nlvl+1)**3), axis=1)                       
                gztmp =  np.sum(( strainfac* (datalocs[:,2][:,None] - zsubgrid[index,:][None,:]) / S / (nlvl+1)**3), axis=1)
                

   
                if datatype == 'los':
                    #glos[index,:] = (gx*east_look + gy*north_look) * np.sin(inc) + gz*np.cos(inc) 
                    glos[index,:] = gxtmp*cosE + gytmp*cosN + gztmp*cosU
                elif datatype == 'eu':
                    gx[index,:] = gxtmp
                    gz[index,:] = gztmp
                elif datatype == 'enu':
                    gx[index,:] = gxtmp
                    gy[index,:] = gytmp
                    gz[index,:] = gztmp               
                
        
                xyz[index,0] =  xgrid[k,i,j]
                xyz[index,1] =  ygrid[k,i,j]
                xyz[index,2] =  zgrid[k,i,j]               
                
                index+=1
  
    if datatype == 'enu':    
        g = np.transpose(np.concatenate((gx,gy,gz),axis=1))
    elif datatype == 'eu':
        g = np.transpose(np.concatenate((gx,gz),axis=1))
    elif datatype == 'los':
        g = np.transpose(glos)
        

    # maybe for this case just: g = np.transpose[glos] ?
    
    
#    ## Regularization, smoothing
#    # first derivative
#    i=0
#    j=0
#    for i in range(0,ngrid-nx-ny-nz):
#        if not gridbound[i]:
#            P[j,i] = -1./dx - 1./dy - 1./dz   
#            P[j,i+nz] = 1./dx  
#            P[j,i+nz*nx] = 1./dy
#            P[j,i+1] = 1./dz
#            j+=1
    

    ## TEST VERSION: GHOST BOUNDARIES = 0  / 22.10.19
    if regularization:
        P = createPmatrix(ngrid,P,Ptype,gridbound,xleft,xright,yfront,yback,ztop,zbot,nx,ny,nz,dx,dy,dz,xgrid,ygrid,zgrid)  
        G = np.concatenate((g,P*np.sqrt(Wr)))
##        d = np.concatenate((d,np.zeros(ngrid)),axis=0)
        d = np.concatenate((d,np.zeros(ngrid-outerCells)),axis=0)        
    else:
        G = g
     
    return G, g, d, P


def appendGmatrix(G,d,elmAboveSurf, ny, nz):
    nConstraints = np.sum(elmAboveSurf)
    GC = np.zeros((nConstraints,G.shape[1]))
    
    index = 0
    for i in range(elmAboveSurf.shape[1]): # x-direction
        for j in range(elmAboveSurf.shape[0]): # y-direction
            if elmAboveSurf[j,i]:
                GC[index,nz*j+nz*ny*i] = 1
                index += 1
     
    G = np.concatenate((G,GC),axis=0)
    d = np.concatenate((d,np.zeros((nConstraints))))
    
    return G,d

## this did not really work too well...
## intended to force more cells to be zero
def appendGmatrix2(G,d, ny, nz, ngrid):
    nConstraints = ngrid
    GC = np.zeros((nConstraints,G.shape[1]))
    
    for i in range(ngrid):
        GC[i,i] = 1
     
    G = np.concatenate((G,GC),axis=0)
    d = np.concatenate((d,np.zeros((nConstraints))))
    
    return G,d

def createGmatrix_nosubgrid(nu,datalocs,xgrid,ygrid,zgrid,nx,ny,nz,d,ndata,datatype,Wr =0,cosE=0, cosN=0, cosU=0, regularization=False, straintype='constrained'):
    ngrid = nx*ny*nz
    datatype = datatype.lower()
    
    dx = xgrid[0,0,1]-xgrid[0,0,0]
    dy = ygrid[0,1,0]-ygrid[0,0,0]
    dz = zgrid[1,0,0]-zgrid[0,0,0]
    
    if datatype == 'enu':
        gxtmp = np.zeros((3*ndata))
        gytmp = np.zeros((3*ndata))
        gztmp = np.zeros((3*ndata))
        
        gx = np.zeros((ngrid,3*ndata))
        gy = np.zeros((ngrid,3*ndata))
        gz = np.zeros((ngrid,3*ndata))
    elif datatype == 'eu':
        gx = np.zeros((ngrid,2*ndata))
        gz = np.zeros((ngrid,2*ndata))
        gxtmp = np.zeros((3*ndata))
        gytmp = np.zeros((3*ndata))
        gztmp = np.zeros((3*ndata))
    elif datatype == 'los':
        gxtmp = np.zeros((ndata))
        gytmp = np.zeros((ndata))
        gztmp = np.zeros((ndata))
        glos = np.zeros((ngrid,ndata))

#    xbound = np.zeros((ncells))
#    ybound = np.zeros((ncells))
#    zbound = np.zeros((ncells))

    xleft = np.zeros((ngrid), dtype=bool)
    xright = np.zeros((ngrid), dtype=bool)
    yfront = np.zeros((ngrid), dtype=bool)
    yback = np.zeros((ngrid), dtype=bool)
    zbot = np.zeros((ngrid), dtype=bool)
    ztop = np.zeros((ngrid), dtype=bool)
    
    gridbound = np.zeros((ngrid), dtype=bool)
    
    xyz = np.zeros((ngrid,3))
    
    if datatype == 'los':  g = np.zeros((ndata,ngrid))
    if datatype == 'eu':   g = np.zeros((2*ndata,ngrid))
    if datatype == 'enu':  g = np.zeros((3*ndata,ngrid))

    P = np.zeros((ngrid,ngrid))
    
    # index system for all dimensions combined
    index=0
    
    for j in range(0,nx):
        for i in range(0,ny):
            for k in range(0,nz):                
                if straintype.lower() == 'constrained':
                    strainfac = (1 - nu) / np.pi  
                elif straintype.lower() == 'free':
                    strainfac = (nu + 1) / (3.*np.pi)
                else:
                    raise ValueError('invalid straintype')
                
                # boundaries with Dirichlet BC (dV=0)
#                if i == nx-1 or i == 0:
#                    xbound[index] = 1
#                if j == ny-1 or j == 0:
#                    ybound[index] = 1
#                if k == nz-1 or :
#                    zbound[index] = 1                         

#                if i == nx-1 or i == 0 or j == ny-1P or j == 0 or k == nz-1:
#                    gridbound[index] = 1   
#                #elif i == nx-2 or i==1 or j == ny-2 or j==1 or k == nz-2 or k == 1: # for 2nd derivative
#                elif k == 1:
#                    topbound[index] = 1
                    
                if j == nx-1 or j == 0 or i == ny-1 or i == 0 or k == nz-1 or k == 0:
                    gridbound[index] = True   
                    
                    if j == 0:
                        xleft[index] = True
                    if j == nx-1:
                        xright[index] = True
                    if i == 0:
                        yfront[index] = True
                    if i == ny-1:
                        yback[index] = True
                    if k == 0:
                        ztop[index] = True
                    if k == nz-1:
                        zbot[index] = True                

#                if i == nx-1  or j == ny-1 or k == nz-1:
#                    gridbound[index] = 1                                         

                gxtmp =  (( strainfac * (datalocs[:,0] - xgrid[k,i,j])) / \
                  (np.sqrt((datalocs[:,0] - xgrid[k,i,j])**2 + (datalocs[:,1] - ygrid[k,i,j])**2 + (datalocs[:,2] - zgrid[k,i,j])**2)**3)) 
                                    
                gytmp =  (( strainfac* (datalocs[:,1] - ygrid[k,i,j])) / \
                  (np.sqrt((datalocs[:,0] - xgrid[k,i,j])**2 + (datalocs[:,1] - ygrid[k,i,j])**2 + (datalocs[:,2] - zgrid[k,i,j])**2)**3)) 
                
                  
                gztmp =  (( strainfac * (datalocs[:,2] - zgrid[k,i,j])) / \
                  (np.sqrt((datalocs[:,0] - xgrid[k,i,j])**2 + (datalocs[:,1] - ygrid[k,i,j])**2 + (datalocs[:,2] - zgrid[k,i,j])**2)**3)) 
   
                
                if datatype == 'los':
                    #glos[index,:] = (gx*east_look + gy*north_look) * np.sin(inc) + gz*np.cos(inc) 
                    glos[index,:] = gxtmp*cosE + gytmp*cosN + gztmp*cosU
                elif datatype == 'eu':
                    gx[index,:] = gxtmp
                    gz[index,:] = gztmp
                elif datatype == 'enu':
                    gx[index,:] = gxtmp
                    gy[index,:] = gytmp
                    gz[index,:] = gztmp               
                
        
                xyz[index,0] =  xgrid[k,i,j]
                xyz[index,1] =  ygrid[k,i,j]
                xyz[index,2] =  zgrid[k,i,j]               
                
                index+=1
        
    for i in range(0,ndata):
        if datatype == 'enu':
            g[i,:]         = np.transpose(gx[:,i])
            g[i+ndata,:]   = np.transpose(gy[:,i])
            g[i+2*ndata,:] = np.transpose(gz[:,i])
        elif datatype == 'eu':
            g[i,:]         = np.transpose(gx[:,i])
            g[i+ndata,:]   = np.transpose(gz[:,i])
        elif datatype == 'los':
            g[i,:]         = np.transpose(glos[:,i])
            
    # maybe for this case just: g = np.transpose[gloPs] ?

    ## REGULARIZATION
        
    if regularization:
    
        # second derivative
        for i in range(0,ngrid):
            if not gridbound[i]:  
                P[i,i] = -2./dx**2. -2./dy**2. - 2./dz**2.   
                P[i,i+nz*ny] = 1./dx**2.  
                P[i,i-nz*ny] = 1./dx**2.  
                P[i,i+nz] = 1./dy**2.
                P[i,i-nz] = 1./dy**2.
                P[i,i+1] = 1./dz**2.  
                P[i,i-1] = 1./dz**2.                 

            else:
                if xleft[i]:
                    P[i,i+nz*ny] = -1/dx
                    P[i,i] += 1/dx
                if xright[i]:
                    P[i,i-nz*ny] = -1/dx
                    P[i,i] += 1/dx
                if yfront[i]:
                    P[i,i+nz] = -1/dy
                    P[i,i] += 1/dy
                if yback[i]:
                    P[i,i-nz] = -1/dy
                    P[i,i] += 1/dy
                if ztop[i]:
                    P[i,i+1] = -1/dz
                    P[i,i] += 1/dz
                if zbot[i]:
                    P[i,i-1] = -1/dz
                    P[i,i] += 1/dz           
    
                         
                    
        #P *= np.sqrt(Wr)    
        G = np.concatenate((g,P*np.sqrt(Wr)))
        d = np.concatenate((d,np.zeros(ngrid)),axis=0)
    else:
        G = g
        
    return G, g, d, P
