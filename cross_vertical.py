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
import datetime
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import from_levels_and_colors
from cartopy.feature import NaturalEarthFeature, COLORS
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,
                 cartopy_xlim, cartopy_ylim, interpline, CoordPair)
#####################################################################
def Cross_Section_with_Mountains(wrfncfile,lat1,lat2,lon1,lon2,tvar,colors,titlname,figurefilename):
    wrf_file = Dataset(wrfncname)
    # Define the cross section start and end points
    cross_start = CoordPair(lat=lat1, lon=lon1)
    cross_end = CoordPair(lat=lat2, lon=lon2)
    # Get the WRF variables
    ht = getvar(wrf_file, "z", timeidx=-1)
    wspd =  getvar(wrf_file, "uvmet_wspd_wdir", units="m s-1")[0,:]
    wdir =  getvar(wrf_file, "uvmet_wspd_wdir", units="m s-1")[1,:]
    ter = getvar(wrf_file, "ter", timeidx=-1) # Model Terrain Height
    #dbz = getvar(wrf_file, "dbz", timeidx=-1)
    #max_dbz = getvar(wrf_file, "mdbz", timeidx=-1)
    #Z = 10**(dbz/10.) # Use linear Z for interpolation
    # Compute the vertical cross-section interpolation.  Also, include the
    # lat/lon points along the cross-section in the metadata by setting latlon
    # to True.
    z_cross = vertcross(wspd, ht, wrfin=wrf_file,
                    start_point=cross_start,
                    end_point=cross_end,
                    latlon=True, meta=True)
    # Convert back to dBz after interpolation
    dbz_cross = z_cross #10.0 * np.log10(z_cross)
    # Add back the attributes that xarray dropped from the operations above
    dbz_cross.attrs.update(z_cross.attrs)
    dbz_cross.attrs["description"] = titlname
    dbz_cross.attrs["units"] = "m s-1"
    # To remove the slight gap between the dbz contours and terrain due to the
    # contouring of gridded data, a new vertical grid spacing, and model grid
    # staggering, fill in the lower grid cells with the first non-missing value
    # for each column.
    # Make a copy of the z cross data. Let's use regular numpy arrays for this.
    dbz_cross_filled = np.ma.copy(to_np(dbz_cross))
    # For each cross section column, find the first index with non-missing
    #    values and copy these to the missing elements below.
    for i in range(dbz_cross_filled.shape[-1]):
        column_vals = dbz_cross_filled[:,i]
        # Let's find the lowest index that isn't filled. The nonzero function
        # finds all unmasked values greater than 0. Since 0 is a valid value
        # for dBZ, let's change that threshold to be -200 dBZ instead.
        first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
        dbz_cross_filled[0:first_idx, i] = dbz_cross_filled[first_idx, i]
        # Get the terrain heights along the cross section line
    ter_line = interpline(ter, wrfin=wrf_file, start_point=cross_start,
                      end_point=cross_end)
    # Get the lat/lon points
    lats, lons = latlon_coords(dbz)
    # Get the cartopy projection object
    cart_proj = get_cartopy(dbz)
    # Create the figure
    fig = pyplot.figure(figsize=(8,6))
    ax_cross = pyplot.axes()
    dbz_levels = np.arange(2., 30., 2.)
    # Create the color table found on NWS pages.
    nlv=len(dbz_levels)
    dbz_rgb = colors[0:nlv]
    dbz_map, dbz_norm = from_levels_and_colors(dbz_levels, dbz_rgb, extend="max")
    # Make the cross section plot for dbz
    dbz_levels = np.arange(2., 30., 2.)
    xs = np.arange(0, dbz_cross.shape[-1], 1)
    ys = to_np(dbz_cross.coords["vertical"])
    dbz_contours = ax_cross.contourf(xs,
                                 ys,
                                 to_np(dbz_cross_filled),
    #                             levels=dbz_levels,
                                 cmap=dbz_map,
                                 norm=dbz_norm,
                                 extend="max")
    # Add the color bar
    cb_dbz = fig.colorbar(dbz_contours, ax=ax_cross)
    cb_dbz.ax.tick_params(labelsize=8)
    # Fill in the mountain area
    ht_fill = ax_cross.fill_between(xs, 0, to_np(ter_line),
                                facecolor="saddlebrown")
    # Set the x-ticks to use latitude and longitude labels
    coord_pairs = to_np(dbz_cross.coords["xy_loc"])
    x_ticks = np.arange(coord_pairs.shape[0])
    x_labels = [pair.latlon_str() for pair in to_np(coord_pairs)]
    # Set the desired number of x ticks below
    num_ticks = 5
    thin = int((len(x_ticks) / num_ticks) + .5)
    ax_cross.set_xticks(x_ticks[::thin])
    ax_cross.set_xticklabels(x_labels[::thin], rotation=45, fontsize=8)
    # Set the x-axis and  y-axis labels
    ax_cross.set_xlabel("Latitude, Longitude", fontsize=12)
    ax_cross.set_ylabel("Height (m)", fontsize=12)
    # Add a title
    ax_cross.set_title(titlname, {"fontsize" : 14})
    plt.savefig(figurefilename+'.png',dpi=300)
    return
######################################################################
dirinSIM='../'
diroutF='../Pics/'
diroutD='../Data/'
sdate=[2014,7,26,12]
edate=[2014,7,28,0]
Cases=['M01','M02','M03','M04','M05','M06','M07','M08','M09','M10']
nc=len(Cases)
bnd=[110,115,30,35]
allrain=[]
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
                minstr="%2.2d"%ih
                fname='wrfout_d'+dmstr+'_2014-07-'+daystr+'_'+hourstr+':00:00'
                fpath=dirin+fname
                #=========get var and plotting=========================
                colors=["#001889",
                        "#4D008A",
                        "#78008D",
                        "#98008C",
                        "#B32086",
                        "#C93F79",
                        "#DA5D66",
                        "#E67B48",
                        "#EC9A06",
                        "#ECBA00",
                        "#E6DC00",
                        "#DAFF47"]
                lat1 = 31;lat2=31; lon1=110;lon2=113
                tvar="uvmet_wspd_wdir"
                titlname= 'Wind'
                figurefilename=diroutF+casename+'_vertical_cross_wind'
                Cross_Section_with_Mountains(fpath,lat1,lat2,lon1,lon2,
                    tvar,colors,titlname,figurefilename)







