#import numpy as np

insar1 = 'data/al_envi.dat'
insar2 = 'data/al_alos.dat'

datalocation = [insar1,insar2]
dataNames = ['ENVISAT', 'ALOS']

prjname = 'example'

datatype = 'los'

# origin for parameter conversion
origin = [38.78799,  7.78826] 

# number of grid points
nx = 12
ny = 12
nz = 1

grid_z0 = 1000.
grid_lonExtent = 15000.
grid_latExtent = 15000.
grid_zMax   = 2000.

# subgrid level
subgridlvl = 1

# downsampling
quadpixsize    = 1e2
quadtolerance  = .0035
quadfittype    = 1
quadstartlevel = 3
quadmaxdim     = 12

# Poisson's ratio
nu = 0.25

# InSAR covariance parameters
# now in m and m2
c_nugget = 3.8e-6
c_range = 3159
c_sill = 1.3e-5

# Monte Carlo simulation for estimation of uncertainties
nMC = 1000

# skip downsampling
do_dwnsmpl   = 1 
load_dwnsmpl = 0 # 1 or 0

# regularization
regularization = 1
Wr = 0.002 

# shift data
datashift = [0,0.01]  

# determine pmin / predefine cut-off 
determine_pmin = 0
cutoff = 0 # needs to be integer
