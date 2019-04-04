#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 12:23:26 2019

@author: chenjh
"""
import matplotlib
matplotlib.use('Agg')
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
#import cartopy.crs as crs
#from cartopy.feature import NaturalEarthFeature
from wrf import (to_np, getvar, smooth2d, get_basemap, latlon_coords)
#from matplotlib.cm import get_cmap
from mpl_toolkits.basemap import Basemap
#####################################################################
def  plot1(f,dirout,casename,datestr,varname,tilstr,imm,dm):                     
    temp1=getvar(f,varname,timeidx=imm) # Grid rainfall
    if dm+1>1:
        if imm==-1:
            temp=temp1
            if 'dbz' in varname:
                tmpe=temp1[0,:,:]
        else:
            if 'dbz' in varname:
                temp=temp1[0,:,:]
            else:
                temp=temp1[:,:]
    else:
        if imm==-1:
            temp=temp1
            if 'dbz' in varname:
                tmpe=temp1[0,:,:]
        else:
            if 'dbz' in varname:
                temp=temp1[0,:,:]
            else:
                temp=temp1[:,:]
    lats, lons = latlon_coords(temp)
    #cart_proj = get_cartopy(temp)
    bm=get_basemap(temp)
    fig= plt.figure(figsize=(12,6))
    # Add geographic outlines
    bm.drawcoastlines(linewidth=0.25)
    bm.drawstates(linewidth=0.25)
    bm.drawcountries(linewidth=0.25)
    # Download and add the states and coastlines
    x, y = bm(to_np(lons), to_np(lats))
    bm.contourf(x, y, to_np(temp), 10,
        cmap=get_cmap("Greens"))
    #    Add a color bar
    plt.colorbar(shrink=.98)
    # Set the map bounds
    #ax.set_xlim(cartopy_xlim(temp))
    #ax.set_ylim(cartopy_ylim(temp))
    # Add the gridlines
    #ax.gridlines(color="black", linestyle="dotted")
    plt.title(tilstr)
    plt.savefig(dirout+casename+varname+"_"+datestr+'.png',dpi=300)
    plt.close()
    return
######################################################################
dirinSIM='../'
diroutF='../Pics/'
diroutD='../Data/'
Cases=["M01","M02","M03","M04","M05","M06","M07"]
nc=len(Cases)
for i in range(0,nc):
    fold="M"+"%2.2d"%(i+1)
    casename='WRF_'+fold
    dirin=dirinSIM+fold+'/'
    for dm in range(0,3):
        dmstr="%2.2d"%(dm+1)
        for iday in range(26,28):
            ih1=0
            ih2=24
            if dm+1==3:
                im1=0
                im2=6
            elif dm+1==2:
                im1=0
                im2=2
            else:
                im1=0
                im2=1
            if iday==26:
                ih1=12
            if iday==28:
                ih2=1            
            daystr="%2.2d"%iday
            for ih in range(ih1,ih2):
                hourstr="%2.2d"%ih                
                imm=-1
                if dm+1==3:
                    imm=0
                for im in range(im1,im2):
                    if dm+1==3:
                        minstr="%2.2d"%(im*10)
                    elif dm+1==2:
                        minstr="%2.2d"%(im*30)
                    else:
                        minstr="%2.2d"%im
                    fname='wrfout_d'+dmstr+'_2014-07-'+daystr+'_'+hourstr+':00:00'
                    datestr='d'+dmstr+'_2014-07-'+daystr+'_'+hourstr+'_'+minstr+'_00'
                    fpath=dirin+fname
                    f=Dataset(fpath,'r')
                    #=========get var and plotting=========================
                    #Grid Rain
                    # plot 1
                    dirout=diroutF
                    varname='RAINNC'
                    tilstr='Grid Scale Rainfall'
                    plot1(f,dirout,casename,datestr,varname,tilstr,imm)
                    #
                    dirout=diroutF
                    varname='RAINC'
                    tilstr='Cumulus Rainfall'
                    plot1(f,dirout,casename,datestr,varname,tilstr,imm) 
                    # 
                    dirout=diroutF
                    varname='dbz'
                    tilstr='Radar Reflectivity'
                    plot1(f,dirout,casename,datestr,varname,tilstr,imm)
                    if dm+1==3:
                        imm=imm+1