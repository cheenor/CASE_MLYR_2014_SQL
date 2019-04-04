#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 07:46:04 2019

@author: chenjh
"""
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import numpy as np
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords,interplevel)
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
            ncfile=Dataset(fpath,'r')
            #=========get var and plotting=========================
            #Grid Rain
            p = getvar(ncfile, "pressure")
            z = getvar(ncfile, "z", units="dm")
            ua = getvar(ncfile, "ua", units="kt")
            va = getvar(ncfile, "va", units="kt")
            wspd = getvar(ncfile, "wspd_wdir", units="kts")[0,:]
            tmp=getvar(ncfile, "tc")
            # Interpolate geopotential height, u, and v winds to 500 hPa
            ht_500 = interplevel(z, p, 500)
            u_500 = interplevel(ua, p, 500)
            v_500 = interplevel(va, p, 500)
            wspd_500 = interplevel(wspd, p, 500)
            t_500=interplevel(tmp, p, 500)
            #
            cart_proj = get_cartopy(ht_500)
            lats, lons = latlon_coords(ht_500)
            fig = plt.figure(figsize=(12,6))
            ax = plt.axes(projection=cart_proj)
            # Download and add the states and coastlines
            states = NaturalEarthFeature(category="cultural", scale="50m",
                             facecolor="none",
                             name="admin_1_states_provinces_shp")
            ax.add_feature(states, linewidth=.5, edgecolor="black")
            ax.coastlines('50m', linewidth=0.8)           
            #Add the 500 hPa geopotential height contours
            levels = np.arange(520., 580., 6.)
            smooth_ht = smooth2d(ht_500, 3, cenweight=4)
            contours = plt.contour(to_np(lons), to_np(lats), to_np(smooth_ht),
                       colors="black",
                       transform=crs.PlateCarree())
            plt.clabel(contours, inline=1, fontsize=10, fmt="%i")
            # Add the wind speed contours
            levels = [25, 30, 35, 40, 50, 60, 70, 80, 90, 100, 110, 120]
            temp_contours = plt.contourf(to_np(lons), to_np(lats), to_np(t_500),
                             cmap=get_cmap("Reds"),
                             transform=crs.PlateCarree())
            #plt.colorbar(temp_contours, ax=ax, orientation="horizontal", pad=.05)
            plt.colorbar(temp_contours)
            # Add the 500 hPa wind barbs, only plotting every 125th data point.
            plt.barbs(to_np(lons[::15,::15]), to_np(lats[::15,::15]),
                      to_np(u_500[::15,::15]), to_np(v_500[::15,::15]),
                      transform=crs.PlateCarree(), color='blue',length=6)
            # Set the map bounds
            ax.set_xlim(cartopy_xlim(ht_500))
            ax.set_ylim(cartopy_ylim(ht_500))
            # Add the gridlines
            ax.gridlines(color="black", linestyle="dotted")            
            plt.title("500 hpa")
            plt.savefig(dirout+"500pha"+datestr+'.png',dpi=300)
            plt.close() 

