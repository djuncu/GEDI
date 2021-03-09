#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 15:59:44 2019

Plot quadtree downsampled deformation maps

@author: Daniel Juncu
"""

import geopandas as gpd
from shapely.geometry import Polygon 
import matplotlib.pyplot as plt

def plotQuads(cx,cy,subdata,minval,maxval,colormap, title=None, lengthUnits = 'km'):
    columns = [ 'geometry', 'defo']
    quads = gpd.GeoDataFrame(columns=columns)
    quads['defo']= subdata
    
    for index, row in quads.iterrows(): 
        quads.loc[index, 'geometry'] = Polygon(zip(cx[:,index],cy[:,index]))

    fig, ax = plt.subplots(1,1)
    quads.plot(column='defo', cmap=colormap, legend=True, alpha=1, vmin=minval, vmax=maxval, ax=ax)
    ax.set_title(title)
    ax.set_ylabel(f'Northing in {lengthUnits}')
    ax.set_xlabel(f'Easting in {lengthUnits}')
