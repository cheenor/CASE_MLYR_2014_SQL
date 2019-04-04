#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 12:23:26 2019

@author: chenjh
"""
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import datetime
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import from_levels_and_colors
from wrf import (getvar, to_np, ALL_TIMES, latlon_coords, vertcross,
                 cartopy_xlim, cartopy_ylim, interpline, CoordPair)
#####################################################################
def getdata(wrflist,tarn,un,region,fileout):
    wdir=0
    print tarn
    if "wspd_wdir" in tarn:
        tmp1 =  getvar(wrflist, tarn, units=un,timeidx=ALL_TIMES,method="cat")[0,:] # wspd
        tmp2 =  getvar(wrflist, tarn, units=un,timeidx=ALL_TIMES,method="cat")[1,:] # wdir
        wdir=1 
    else:
        tmp1 = getvar(wrflist, tarn,  units=un,timeidx=ALL_TIMES,method="cat")
    ht = getvar(wrflist, "z",  units="km",timeidx=ALL_TIMES,method="cat")
    ter = getvar(wrflist, "ter", timeidx=-1) # Model Terrain Height        
    lats, lons = latlon_coords(tmp1)
    nx=len(lats[0,:])
    ny=len(lats[:,0])
    ### in my understand
    fdat=open(fileout,'w')
    fdat.write(tarn)
    fdat.write(' '+un)
    fdat.write("\n")
    size=tmp1.shape
    ndim=len(size)
    if ndim ==4: # t,z, lat,lon
        ntt=size[0]
        nz=size[1]
    elif ndim==3: # t lat lon 
        ntt=size[0]
        nz=0
    for it in range(0,ntt):
        at=np.zeros(shape=(nz),dtype=float)
        if wdir==1:
            bt=np.zeros(shape=(nz),dtype=float)  
        cont=0.0
        for ix in range(0,nx):
            for iy in range(0,ny):
                if lats[iy,ix] in range(region[2],region[3]) and lons[iy,ix] in range(region[0],region[1]) :                
                    cross_start = CoordPair(lat=lats[iy,ix], lon=lons[iy,ix])
                    cross_end = CoordPair(lat=lats[iy,ix], lon=lons[iy,ix])
                    # Compute the vertical cross-section interpolation.  Also, include the
                    # lat/lon points along the cross-section in the metadata by setting latlon
                    # to True.
                    if nz>0:
                        tmp_cross = vertcross(tmp1, ht, wrfin=wrflist,
                            start_point=cross_start,
                            end_point=cross_end,
                            latlon=True, meta=True)
                        for iz in range(0,nz):
                            at[iz]=at[iz]+tmp_cross[iz]
                            if wdir==1:
                                tmp2_cross = vertcross(tmp2, ht, wrfin=wrflist,
                                    start_point=cross_start,
                                    end_point=cross_end,
                                    latlon=True, meta=True) 
                                bt[iz]=bt[iz]+tmp2_cross[iz]
                    elif nz==0:
                        at[0]=at[0]+tmp1[iy,ix]
                        if wdir==1:
                            bt[0]=bt[0]+tmp2[iy,ix]    
                    cont=cont+1.0                           
        if nz>0:
            for iz in range(0,nz):
                if cont>0:
                    at[iz]=at[iz]/cont
                    itm="%9.2e"%at[iz]
                    fdat.write(itm)
                    if wdir==1:
                        bt[iz]=bt[iz]/cont
                        itm="%9.2e"%at[iz]
                        fdat.write(" ")
                        fdat.write(itm)
                    fdat.write("\n")
        else:
            if cont>0:
                at[0]=at[0]/cont
                itm="%9.2e"%at[0]
                fdat.write(itm)
                if wdir==1:
                    bt[0]=bt[0]/cont
                    itm="%9.2e"%at[0]
                    fdat.write(" ")
                    fdat.write(itm)
                fdat.write("\n")

    fdat.close()
    return 
######################################################################
dirinSIM='../'
diroutF='../Pics/'
diroutD='../Data/'
sdate=[2014,7,26,12]
edate=[2014,7,28,0]
Cases=['M01','M02','M03','M04','M05','M06','M07','M08'] #,'M09','M10']
nc=len(Cases)
bnd=[110,115,30,35]
region=[29.32,34.0,119.6,122.1]
varname=["z","ua","va","uvmet_wspd_wdir","dbz","pw","uvmet10_wspd_wdir",
    "cape3d_only","cin3d_only","lcl","lfc","p","XLAND"]
varunits=["km","m s-1","m s-1","m s-1","dBz","kg m-2","m s-1","J kg-1","J kg-1","J kg-1","J kg-1","hPa",""]
nv=len(varname)
for i in range(0,nc):
    fold="M"+"%2.2d"%(i+1)
    casename='WRF_'+fold
    dirin=dirinSIM+fold+'/'
    tmp1=[]
    tmp2=[]
    id1=sdate[2]
    id2=edate[2]
    for dm in range(3,4):
        dmstr="%2.2d"%(dm+1)
        wrflist=[]
        for iday in range(id1,id2+1):
            ih1=0
            ih2=24
            im1=0
            im2=60
            if iday==id1:
                ih1=sdate[3]+6
            if iday==id2:
                ih2=edate[3]+1
            daystr="%2.2d"%iday
            for ih in range(ih1,ih2):
                hourstr="%2.2d"%ih
                minstr="%2.2d"%(ih*10)
                fname='wrfout_d'+dmstr+'_2014-07-'+daystr+'_'+hourstr+':00:00'
                fpath=dirin+fname
                #=========get var and plotting=========================
                f=Dataset(fpath)
                wrflist.append(f)
        for iv in range(0,nv):
            tarn=varname[iv]
            un=varunits[iv]
            fileout=diroutD+casename+"_ave_"+tarn+'.dat'
            getdata(wrflist,tarn,un,region,fileout)















