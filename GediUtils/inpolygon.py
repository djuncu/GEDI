#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 15:05:53 2019

Find points xyq located inside of polygon defined by xv, yv. 
Structure follows Matlab's inpolygon. 
modified (improved speed) from https://stackoverflow.com/a/49733403 

@author: Daniel Juncu
"""

import numpy as np
from matplotlib import path

def inpolygon(xyq, xv, yv):
    shape = xyq.shape[0]
    
    xv = xv.reshape(-1)
    yv = yv.reshape(-1)

    p = path.Path([(xv[i], yv[i]) for i in range(xv.shape[0])])
    return p.contains_points(xyq).reshape(shape)


