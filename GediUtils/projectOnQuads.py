#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 14:47:19 2019

Projecting values of scattered pixels onto quads
using pixel mean value

@author: Daniel Juncu
"""

import numpy as np
import GediUtils as gu

def projectOnQuads(insar_xy,cx,cy,val):
    
    nquads = cx.shape[1]
    quadval = np.empty((nquads,))
    
    for i in range(nquads):
        boo = gu.inpolygon(insar_xy, cx[:,i], cy[:,i]) 
        quadval[i] = np.mean(val[boo])
        
    return quadval
        



    
