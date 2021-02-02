#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 17:19:11 2019

Coordinate projection between WGS and local cartesian coordinates 
using Lambert azimuthal equal area projection, Snyder (1987) p.182pp

Reference:
J. P. Snyder, 1987. Map Projections-A Working Manual
USGS Professional Paper 1395

@author: Daniel Juncu
"""
import numpy as np

def wgs2local(lonlat, origin):
    # Lambert azimuthal equal area projection, Snyder (1987) p.182pp
    # suitable for projections of continental and regional scale
    #
    # Reference:
    # J. P. Snyder, 1987. Map Projections-A Working Manual
    # USGS Professional Paper 1395
    #
    # origin: sequence of lon, lat
    #
    # Daniel Juncu, 2019

    # WGS84 ellipsoid parameters
    a = 6378137.0
    e = 0.081819190842621
    
    origin = np.array(origin)

    # create tuples to access coordinates arrays. necessary in case only 1 coordinate pair is given
    # in that case lonlat must be of size 1x2
    dim = lonlat.ndim

    tpl0 = [slice(None)] * dim
    tpl1 = [slice(None)] * dim
    tpl0[dim-1] = 0
    tpl1[dim-1] = 1
    tpl0 = tuple(tpl0)
    tpl1 = tuple(tpl1)
    
    lon = lonlat[tpl0] * np.pi / 180.
    lat = lonlat[tpl1] * np.pi / 180.

    ori = origin * np.pi / 180.

    lambda0 = ori[0]
    phi1 = ori[1]
    
    #nonzero_lat = lat != 0.
    #d_lon = lon[nonzero_lat] - ori[0]
    
    q = (1.-e**2.) * (np.sin(lat) / (1.-e**2.*np.sin(lat)**2.) - 1./(2.*e)* np.log((1.-e*np.sin(lat))/ (1.+e*np.sin(lat))))
    qp = (1.-e**2.) * (np.sin(np.pi/2.) / (1.-e**2.*np.sin(np.pi/2.)**2.) - 1./(2.*e) * np.log((1.-e*np.sin(np.pi/2.))/ (1.+e*np.sin(np.pi/2.))))
    q1 =  (1.-e**2.) * (np.sin(phi1) / (1.-e**2.*np.sin(phi1)**2.) - 1./(2.*e) * np.log((1.-e*np.sin(phi1))/ (1.+e*np.sin(phi1))))

    beta = np.arcsin(q/qp)
    beta1 = np.arcsin(q1/qp)
    Rq = a * np.sqrt(qp/2.)
    m1 = np.cos(phi1) / np.sqrt(1.-e**2.*np.sin(phi1)**2.)
    
    D = a*m1 / (Rq * np.cos(beta1))
    B = Rq * np.sqrt(2. / (1. + np.sin(beta1) * np.sin(beta) + np.cos(beta1) * np.cos(beta) * np.cos(lon-lambda0)))
    
    x = np.zeros(shape=np.shape(lonlat[tpl0]))
    y = np.zeros(shape=np.shape(lonlat[tpl0]))
    xy = np.zeros(shape=np.shape(lonlat))
    
    x = B * D * np.cos(beta) *np.sin(lon-lambda0)
    y = (B/D) * (np.cos(beta1) * np.sin(beta) - np.sin(beta1) * np.cos(beta) * np.cos(lon-lambda0))
    
    xy[tpl0] = x
    xy[tpl1] = y
    
    return xy


def local2wgs(xy, origin):
    # Lambert azimuthal equal area projection, Snyder (1987) p.182pp
    # suitable for projections of continental and regional scale
    #
    # INPUT
    # xy - Nx2 coordinates as [easting, northing] in m
    # origin - 1x2 origin coordinates in decimal degrees
    # 
    # OUTPUT
    # lon - longitude in decimal degrees
    # lat - latitude in decimal degrees

    # Reference:
    # J. P. Snyder, 1987. Map Projections-A Working Manual
    # USGS Professional Paper 1395
    
    # Daniel Juncu, 2019

    # ellipsoid parameters
    a = 6378137.0
    e = 0.081819190842621

    # create tuples to access coordinates arrays. necessary in case only 1 coordinate pair is given
    dim = xy.ndim
    tpl0 = [slice(None)] * dim
    tpl1 = [slice(None)] * dim

    tpl0[dim-1] = 0
    tpl1[dim-1] = 1

    # convert to meters, radians
    x = xy[tpl0]
    y = xy[tpl1] * 1000.
    ori = origin * np.pi / 180.

    lambda0 = ori[0]
    phi1 = ori[1]
    
    qp = (1.-e**2.) * (np.sin(np.pi/2.) / (1.-e**2.*np.sin(np.pi/2.)**2.) - 1./(2.*e) * np.log((1.-e*np.sin(np.pi/2.))/ (1.+e*np.sin(np.pi/2.))))
    q1 =  (1.-e**2.) * (np.sin(phi1) / (1.-e**2.*np.sin(phi1)**2.) - 1./(2.*e) * np.log((1.-e*np.sin(phi1))/ (1.+e*np.sin(phi1))))

    beta1 = np.arcsin(q1/qp)
    Rq = a * np.sqrt(qp/2.)
    m1 = np.cos(phi1) / np.sqrt(1.-e**2.*np.sin(phi1)**2.)
    D = a*m1 / (Rq * np.cos(beta1))

    rho = np.sqrt( (x/D) ** 2. + (D*y) ** 2.)
    ce = 2. * np.arcsin(rho / (2. * Rq))

    q = qp * ( np.cos(ce) * np.sin(beta1) + (D*y*np.sin(ce)*np.cos(beta1)/rho) )

    lon = lambda0 + np.arctan( x*np.sin(ce) / ( D*rho*np.cos(beta1)*np.cos(ce) - D**2.*y*np.sin(beta1)*np.sin(ce) ))

    lat_trial = np.arcsin(q/2.)
    lat_diff = 1.
    lat_diff_old = 10.
    c=1

    while np.abs(1. - lat_diff/lat_diff_old) > 1e-12:
        lat_diff_old = lat_diff
        
        lat = lat_trial + (1. - e**2.*np.sin(lat_trial)**2.)**2. / (2.*np.cos(lat_trial)) * (q / (1.-e**2.) - np.sin(lat_trial) / \
              (1. - e**2.*np.sin(lat_trial)**2.) + 1./(2.*e) * np.log((1. - e*np.sin(lat_trial)) / (1. + e*np.sin(lat_trial) ))) 

        lat_diff = np.max(np.abs(lat-lat_trial))
        
        c+=1   
        lat_trial = lat

        if c>200:
            raise ValueError('Convergence Failure')

    lonlat[tpl0] = lon * 180. / np.pi
    lonlat[tpl1] = lat * 180. / np.pi

    return lonlat

