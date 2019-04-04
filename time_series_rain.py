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
import datetime
#from matplotlib.cm import get_cmap
#from mpl_toolkits.basemap import Basemap
#####################################################################
def  getdata(f,varname,bnd):                     
    temp=getvar(f,varname) # Grid rainfall
    lats, lons = latlon_coords(temp)
    nx=len(lat[0,:])
    ny=len(lat[:,0])
    e=bnd[1]
    w=bnd[0]
    s=bnd[2]
    n=bnd[3]
    b=[]
    for it in range(0,nt):
        a=0.0
        for ix in range(0,nx):
            for iy in range(0,ny):
                if lons[iy,ix]>w and lons[iy,ix]<e :
                    if lats[iy,ix]>s and lats[iy,ix]<n and temp[it,iy,ix]>0. :
                        a=a+temp[iy,ix]
                        c=c+1
        if c>0:
            a=a/c
        else:
            a=0.0
        b.append(a)
    return b
######################################################################
dirinSIM='../'
diroutF='../Pics/'
diroutD='../Data/'
sdate=[2014,7,26,12]
edate=[2014,7,28,0]
Cases=[]
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
                f=Dataset(fpath,'r')
                #=========get var and plotting=========================
                #Grid Rain
                # plot 1
                varname='RAINNC'
                a=getdata(f,varname,bnd)
                tmp1.append(a)
                varname='RAINC'
                a=getdata(f,varname,bnd)
                tmp2.append(a)
    allrain.append(tmp1)
    allrain.append(tmp2)
    nx=len(tmp1)
# plot
datestart=datetime.datetime(sdate[0],sdate[1],sdate[2],sdate[3]+6,0,0)
det=datetime.timedelta(minutes=10)            
dateiso=[]            
for dt in range(0,nx):
    dateiso.append(datestart+dt*det)
xdate=[]              
for tm in dateiso:
    xdate.append(datetime.datetime.strftime(tm,"%H:%M"))
#
colors=["#00138F",
    "#57008B",
    "#86008C",
    "#A90F87",
    "#C5397A",
    "#DA5F62",
    "#E78439",
    "#ECAB00",
    "#E7D400",
    "#DAFF47"]  # hclwizard.org 
wd=[1,1.5,1,1.5,1,1.5,1,1.5,1,1.5,1,1.5,1,1.5]
lss=['-','-']
if nc> len(colors):
    print 'The colors are not enough, Figure may be not clear enough'
fig,axs = plt.subplots(nrows=3,ncols=1,figsize=(14,15))
for i in range(0,nc):
    rc=allrain[i+1]
    rnc=allrain[i]
    nx=len(rc)
    xxx=range(0,nx)
    zdat=rnc
    ax[0].plot(xxx[0:nx],zdat,color=colors[i],lw=wd[i],ls=lss[1],label=casename[i])
    zdat=rc
    ax[1].plot(xxx[0:nx],zdat,color=colors[i],lw=wd[i],ls=lss[1],label=casename[i])
    zdat=rnc+rc
    ax[2].plot(xxx[0:nx],zdat,color=colors[i],lw=wd[i],ls=lss[1],label=casename[i])
ax[0].set_title(r'(a) Grid precipitation')
ax[1].set_title(r'(b) Cumulus precipitation')
ax[2].set_title(r'(c) Total precipitation')
plt.legend()
fig.subplots_adjust(left=0.1,bottom=0.06,right=1-0.08,top=1-0.06,wspace=0.3)                    
plt.savefig(diroutF+'AllCase_timeserious_rain.png',dpi=300)          
plt.show()
plt.close()                 


