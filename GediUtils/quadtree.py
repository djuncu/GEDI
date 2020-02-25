#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 13:33:40 2019

Quadtree partioning of data with given tolerance
Based on Matlab code by Siguron Jonsson

@author: Daniel Juncu
"""

import numpy as np

def quadcoords(indexmatrix, data):
    lin, col = data.shape
    
    length = indexmatrix.shape[0]
    level = indexmatrix.shape[1] - 1 
    
    coord = np.empty((0,2))
    coordl = np.empty((0,2))
    
    cx = np.empty((4,0))
    cy = np.empty((4,0))
    
    for k in range(length):
        blcksz = lin
        lst = 0
        cst = 0
        
        for l in range(level):
            if not indexmatrix[k,l] == 0: 
                blcksz //= 2
                
                if indexmatrix[k,l] == 1:
                    lst = lst
                    cst = cst
                elif indexmatrix[k,l] == 2:
                    lst = lst
                    cst += blcksz
                elif indexmatrix[k,l] == 3:
                    lst += blcksz
                    cst += blcksz
                elif indexmatrix[k,l] == 4:
                    lst += blcksz
                    cst = cst
                    
        coord = np.vstack((coord, np.array(([cst-1+blcksz//2, lst-1+blcksz//2]))))
        coordl = np.vstack((coordl, np.array((np.nan, np.nan)), np.array((lst,cst)), \
                            np.array((lst, cst+blcksz)), np.array((lst+blcksz,cst+blcksz)), \
                            np.array((lst+blcksz,cst)), np.array((lst,cst)) ))
        
        cx = np.hstack((cx, np.vstack((cst,cst+blcksz, cst+blcksz, cst))))
        cy = np.hstack((cy, np.vstack((lst,lst,lst+blcksz,lst+blcksz))))
                    
    return coord, coordl, cx, cy

def getchunk(index, data):
    length = index.size - 1
    lin, col = data.shape
    
    blcksz = lin
    
    lst = 0
    cst = 0
    
   
    for k in range(length):
        blcksz //= 2        
        if index[k] == 1:
            lst = lst
            cst = cst
        elif index[k] == 2:
            lst = lst
            cst += blcksz
        elif index[k] == 3:
            lst += blcksz
            cst += blcksz
        elif index[k] == 4:
            lst += blcksz
            cst = cst               

    chunk = data[lst:lst+blcksz, cst:cst+blcksz]
    return chunk

def check_quad(oldindmat, data, tolerance, fittype):
    lin, col = oldindmat.shape
    mvalue = np.empty((lin,), dtype=float)
    
    for k in range(lin):
        chunk = getchunk(oldindmat[k,:], data)        
        chunk.reshape(-1)
        
        nonan = ~np.isnan(chunk) 
        chunk_nn = chunk[nonan]
        
        scn = chunk_nn.size

        if scn >= chunk.size / 2 and np.any(nonan):
            if fittype == 1:
                mvalue[k] = np.median(chunk_nn)           
            elif fittype == 2:
                mvalue[k] = np.median(chunk_nn)
                
            tmp = np.ones((scn)) * mvalue[k]
            diff = chunk_nn - tmp   
            rms =  np.sqrt(np.mean(diff**2.))        
                                   
            if rms <= tolerance:                
                oldindmat[k,-1] = 1
                
            else:
                oldindmat[k,-1] = 0
        
        elif scn < chunk.size / 2 and np.any(nonan):
            mvalue[k] = np.nan
            oldindmat[k,-1] = 0
        
        else:
            mvalue[k] = np.nan
            oldindmat[k,-1] = 1
        
            
    return oldindmat, mvalue


def quadtree_level(oldindmat):
    
    lin, col = oldindmat.shape    
    indexmatrix = np.empty((0,col+1), dtype=int)    
    
    for k in range(lin):    
        if oldindmat[k,-1] == 1:
            tmp1 = np.append(oldindmat[k,0:-1], 0)
            tmp2 = np.append(tmp1, 1)            
            indexmatrix = np.vstack((indexmatrix, tmp2))
        else:           
            indexmatrix = np.vstack((indexmatrix,np.hstack((np.ones((4,1))*oldindmat[k,0:-1], np.array(([1],[2],[3],[4],)), \
                                    np.zeros((4,1), dtype=int) ))))

    return indexmatrix.astype(int)


def quadtree(data,tolerance,fittype,startlevel=1,maxdim=13):
    # quadtree_part  - Quadtree partitioning of data with given tolerance
    #
    # function newindmat,cc,cx,cy = quadtree_part(data,tolerance);
    #
    # INPUT:
    #          data      - (MxN) data matrix (can have NaN's)
    #          tolerance - (1x1) tolerance value for each square
    #          fittype   - (1x1) type of fit (mean=0; median=1; bilinear=2)
    #          startlevel- (1x1) starting quadtree level (default =1)
    #          maxdim    - (1x1) maximum number of quadtree levels (default=13)
    #
    # OUTPUT:
    #          newindmat - (txs) quadtree index matrix of level s-2
    #                      and t squares, the last column gives the 
    #                      values for each square.
    #          sqval     - (1xt) row vector of quadtree squares values
    #          cx        - (4xt) x-coordinate of square corners
    #          cy        - (4xt) y-coordinate of square corners
    #          cntp      - (tx2) center point coordinates of squares
    #          nlin      - (1x1) # number of lines in working image (e.g. 512)
    #
    # This function need functions 'check', 'quadtree_level', 'plot_quadtree', and
    # 'get_chunck'.
    #
    # original Matlab code by Sigurjon Jonsson
    # translated to Python by Daniel Juncu
    
    
    lin, col = data.shape
    
    # add NaN's to change data size ti (2^k x 2^k)
    dim = 1
    dimax = np.max([lin, col])
    while dimax > 2.**dim:
        dim += 1
        
    nlin = 2**dim
    ncol = nlin
    
    dataexp = np.full((nlin,ncol), np.nan)   
    dataexp[0:lin,0:col] = data

    # This is the starting quadtree index matrix (values 
    # 10 are arbitrary).  The first number tells you the
    # quadrant (clockwise) and the second value is
    # a test index (zero = further partioning needed), and
    # the last three columns are the bilinear values a, b, 
    # and c in  a + bx + cy
    
    indmat = np.array(([1, 0], [2, 0], [3, 0], [4, 0]), dtype=int)
    
    # Here we add to the index matrix in case we don not want
    # to start at the top level.

    if startlevel > 1:
        for k in range(2,startlevel+1):
            indmat = quadtree_level(indmat)
            
    fini = False
    while not fini:
        newindmat, mvalue = check_quad(indmat, dataexp, tolerance, fittype)
        che = newindmat[:,indmat.shape[1] - 1]


        if np.prod(che) == 0:
            indmat = quadtree_level(newindmat)
        else:
            fini = True

    
    cntp, co2, cx, cy = quadcoords(newindmat,dataexp)
    
    return newindmat, mvalue, cx, cy, cntp, nlin
        
            
            
            
    
    
    
    
