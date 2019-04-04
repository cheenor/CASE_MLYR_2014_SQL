#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 08:27:20 2019

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
dirin='/Volumes/DATA1/Models/WRF/WK01/Test01/'
dirout='/Volumes/MacHD/Work/Papers/CASE_MLYR/PIC/RAW/'
fname='wrfout_d01_2014-07-28_09~00~00'

for dm in range(0,3):
    dmstr="%2.2d"%(dm+1)
    for iday in range(26,29):
        ih1=0
        ih2=24
        if iday==26:
            ih1=13
        if iday==28:
            ih2=13
        daystr="%2.2d"%iday
        for ih in range(ih1,ih2):
            hourstr="%2.2d"%ih
            fname='wrfout_d'+dmstr+'_2014-07-'+daystr+'_'+hourstr+'~00~00'
            datestr='d'+dmstr+'_2014-07-'+daystr+'_'+hourstr+'_00_00'        
            fpath=dirin+fname
            f=Dataset(fpath,'r')
            #=========get var and plotting=========================
            #Grid Rain
            temp=getvar(f,"RAINNC") # Grid rainfall
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
                         cmap=get_cmap("jet"))
            #    Add a color bar
            plt.colorbar(ax=ax, shrink=.98)
            # Set the map bounds
            ax.set_xlim(cartopy_xlim(temp))
            ax.set_ylim(cartopy_ylim(temp))
            # Add the gridlines
            ax.gridlines(color="black", linestyle="dotted")            
            plt.title("Grid scale rainfall")
            plt.savefig(dirout+"RAINNC"+datestr+'.png',dpi=300)
            plt.close()
            # cumulus rain 
            temp=getvar(f,"RAINC") # Grid rainfall
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
                         cmap=get_cmap("jet"))
            #    Add a color bar
            plt.colorbar(ax=ax, shrink=.98)
            # Set the map bounds
            ax.set_xlim(cartopy_xlim(temp))
            ax.set_ylim(cartopy_ylim(temp))
            # Add the gridlines
            ax.gridlines(color="black", linestyle="dotted")            
            plt.title("Cumulus rainfall")
            plt.savefig(dirout+"RAINC"+datestr+'.png',dpi=300)
            plt.close() 
            ############3
            temp=getvar(f,"dbz") # Grid rainfall
            temp1=temp[0,:,:]
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
            plt.contourf(to_np(lons), to_np(lats), to_np(temp[0,:,:]), 10,
                         transform=crs.PlateCarree(),
                         cmap=get_cmap("jet"))
            #    Add a color bar
            plt.colorbar(ax=ax, shrink=.98)
            # Set the map bounds
            ax.set_xlim(cartopy_xlim(temp))
            ax.set_ylim(cartopy_ylim(temp))
            # Add the gridlines
            ax.gridlines(color="black", linestyle="dotted")            
            plt.title("Reflectivity")
            plt.savefig(dirout+"Reflectivity"+datestr+'.png',dpi=300)
            plt.close()              
            f.close()
