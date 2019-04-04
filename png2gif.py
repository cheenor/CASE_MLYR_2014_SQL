#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 13:29:16 2019

@author: chenjh
"""
import imageio
#from PIL import Image 
#from images2gif import writeGif
import numpy as np
############################################################
dirpng='/Volumes/DATA02/ModelOUTPUT/MLYR_SQL/RawPics/'
dirgif='/Volumes/DATA02/ModelOUTPUT/MLYR_SQL/Gifs/'
Cases=['M01','M02','M03','M04','M05','M06','M07','M08']
nc=len(Cases)
for i in range(0,2):
    gifname=dirgif+Cases[i]+'_d03_dbz.gif'
    frames = []
    frames2 = []
    for ih in range(6,18):
        hourstr="%2.2d"%ih
        for im in range(0,6):
            minstr="%2.2d"%(im*10)
            fname='WRF_'+Cases[i]+'dbz_d03_2014-07-27_'+hourstr+'_'+minstr+'_00.png'
            fpath=dirpng+fname
            img = imageio.imread(fpath)
            #img=img.convert("RGB")
            frames.append(img)
            #img2=Image.open(fpath)
            #img2=img2.convert("RGB")
            #img2=np.array(img2)
            #frames2.append(img2)
    imageio.mimsave(gifname,frames,'GIF',duration=0.1)
    gifname=gifname=dirgif+Cases[i]+'_d03_dbz_v2.gif'                      
    #writeGif(gifname, frames2, duration=0.1, subRectangles=False)