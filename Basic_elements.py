#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 16:13:30 2019

@author: chenjh
"""
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords)
#####################################################################
dirin='/Volumes/DATA1/Models/WRF/WK01/Test01/'
dirout='/Volumes/MacHD/Work/Papers/CASE_MLYR/'
fname='wrfout_d01_2014-07-28_09~00~00'
fpath=dirin+fname
f=Dataset(fpath,'r')
#for dims in f.dimensions.values():
#    print dims       # this is a test to check the dimensions in the file
lon=f.variables['XLONG'][:]
lat=f.variables['XLAT'][:]
lev=f.variables['ZNW'][:]
time=f.variables['Times'][:]
#nx=len(lon)
#ny=len(lat)
#nt=len(time)
#nz=len(lev)
#print 'nx= ',nx
#print 'ny= ',ny
#print 'nt= ',nt
#print 'nz= ',nz
fout=open(dirout+'vars.txt','w') 
allvar=[] # 
for var in f.variables:
    print var       # 
    fout.write(var)
    fout.write('\n')
    allvar.append(var) # 
fout.close()   
#tmp=f.variables[allvar[3]][:]
slp=getvar(f,"refl_10cm")
print(slp)