#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 12:23:26 2019

@author: chenjh
"""
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                   cartopy_ylim, latlon_coords)
#from matplotlib.cm import get_cmap
#from mpl_toolkits.basemap import Basemap
#####################################################################
def  plot1(f,dirout,casename,datestr,varname,tilstr,imm):                     
    temp1=getvar(f,varname) # Grid rainfall
    if imm==-1:
        temp=temp1
    else:
        temp=temp1[imm,:,:]
    lats, lons = latlon_coords(temp)
    cart_proj = get_cartopy(temp)
    fig = plt.figure(figsize=(12,6))
    ax = plt.axes(projection=cart_proj)
    # Download and add the states and coastlines
    states = NaturalEarthFeature(category="cultural", scale="50m",
        facecolor="none",
        name="admin_1_states_provinces_shp")
    ax.add_feature(states, linewidth=.5, edgecolor="black")
    ax.coastlines('50m', linewidth=0.8)
    plt.contourf(to_np(lons), to_np(lats), to_np(temp), 10,
        transform=crs.PlateCarree(),
        cmap=get_cmap("Greens"))
    #    Add a color bar
    plt.colorbar(ax=ax, shrink=.98)
    # Set the map bounds
    ax.set_xlim(cartopy_xlim(temp))
    ax.set_ylim(cartopy_ylim(temp))
    # Add the gridlines
    ax.gridlines(color="black", linestyle="dotted")
    plt.title(tilstr)
    plt.savefig(dirout+casename+varname+"_"+datestr+'.png',dpi=300)
    plt.close()
    return
######################################################################
dirinSIM='../'
diroutF='../Pics/'
diroutD='../Data/'
Cases=[]
nc=len(Cases)
for i in range(0,nc):
    fold="M"+"%2.2d"%(i+1)
    casename='WRF_'+fold
    dirin=dirinSIM+fold+'/'
    for dm in range(0,3):
        dmstr="%2.2d"%(dm+1)
        for iday in range(8,9):
            ih1=0
            ih2=24
            if dm+1==3:
                im1=0
                im2=6
            else:
                im1=0
                im2=1
            if iday==8:
                ih1=12
            if iday==9:
                ih2=1
            
            daystr="%2.2d"%iday
            for ih in range(ih1,ih2):
                hourstr="%2.2d"%ih                
                imm=-1
                if dm+1==3:
                    imm=0
                for im in range(im1,im2):
                    minstr="%2.2d"%(im*10)
                    fname='wrfout_d'+dmstr+'_2014-06-'+daystr+'_'+hourstr+':00:00'
                    datestr='d'+dmstr+'_2014-06-'+daystr+'_'+hourstr+'_'+minstr+'_00'
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