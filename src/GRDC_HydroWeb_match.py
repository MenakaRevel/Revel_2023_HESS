#!/opt/local/bin/python
# -*- coding: utf-8 -*-
'''
Find GRDC gauges on rivers Virtual Stations exsits
'''
import numpy as np
import matplotlib.pyplot as plt
import datetime
from matplotlib.cbook import boxplot_stats
# from matplotlib.colors import LogNorm,Normalize,ListedColormap
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# from mpl_toolkits.basemap import Basemap
# import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import sys
import os
import calendar
from multiprocessing import Pool
# from multiprocessing import Process
from multiprocessing import sharedctypes
from numpy import ma
import re
# import my_colorbar as mbar
# import cartopy.crs as ccrs
# import cartopy
# from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
# import cartopy.feature as cfeature
import os
import seaborn as sns
import pandas as pd


# import params as pm
import read_grdc as grdc
import read_hydroweb as hweb
# import cal_stat as stat
#import plot_colors as pc
def along_river(ix,iy):
    xlist=[]
    ylist=[]
    while (ix>=0): #!= -9 or ix !=-10):
        xlist.append(ix)
        ylist.append(iy)
        iix=ix
        iiy=iy
        ix=nextxy[0,iiy,iix]
        iy=nextxy[1,iiy,iix]
        if (ix==-9999):
            break
        ix=ix-1
        iy=iy-1
    return np.array(xlist), np.array(ylist)
#====================================================================
# argv=sys.argv
# syear=int(argv[1])
# eyear=int(argv[2])
# CaMa_dir=argv[3]
# mapname=argv[4]
# figname=argv[5]
# ncpus=int(sys.argv[6])
ens_mem=49
seaborn_map=True
# seaborn_map=False
#==================================
# DA_dir="/cluster/data6/menaka/HydroDA"
syear=2009
eyear=2014
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
mapname="amz_06min"
# ens_mem=21
# lexp=7
# seaborn_map=False
#==================================
#===================
fname=CaMa_dir+"/map/"+mapname+"/params.txt"
with open(fname,"r") as f:
    lines=f.readlines()
#-------
nx     = int(filter(None, re.split(" ",lines[0]))[0])
ny     = int(filter(None, re.split(" ",lines[1]))[0])
gsize  = float(filter(None, re.split(" ",lines[3]))[0])
lon0   = float(filter(None, re.split(" ",lines[4]))[0])
lat0   = float(filter(None, re.split(" ",lines[7]))[0])
west   = float(filter(None, re.split(" ",lines[4]))[0])
east   = float(filter(None, re.split(" ",lines[5]))[0])
south  = float(filter(None, re.split(" ",lines[6]))[0])
north  = float(filter(None, re.split(" ",lines[7]))[0])
#----
nextxy = CaMa_dir+"/map/"+mapname+"/nextxy.bin"
rivwth = CaMa_dir+"/map/"+mapname+"/rivwth_gwdlr.bin"
rivhgt = CaMa_dir+"/map/"+mapname+"/rivhgt.bin"
rivlen = CaMa_dir+"/map/"+mapname+"/rivlen.bin"
elevtn = CaMa_dir+"/map/"+mapname+"/elevtn.bin"
lonlat = CaMa_dir+"/map/"+mapname+"/lonlat.bin"
uparea = CaMa_dir+"/map/"+mapname+"/uparea.bin"
nextxy = np.fromfile(nextxy,np.int32).reshape(2,ny,nx)
rivwth = np.fromfile(rivwth,np.float32).reshape(ny,nx)
rivhgt = np.fromfile(rivhgt,np.float32).reshape(ny,nx)
rivlen = np.fromfile(rivlen,np.float32).reshape(ny,nx)
elevtn = np.fromfile(elevtn,np.float32).reshape(ny,nx)
lonlat = np.fromfile(lonlat,np.float32).reshape(2,ny,nx)
uparea = np.fromfile(uparea,np.float32).reshape(ny,nx)
#===================
rivnum="/cluster/data6/menaka/HydroDA/dat/rivnum_"+mapname+".bin"
rivnum=np.fromfile(rivnum,np.int32).reshape(ny,nx)
rivermap=((nextxy[0]>0)*(rivnum==1))*1.0
#===================
obslist="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+mapname+"_QC0.txt"
with open(obslist,"r") as f:
    lines=f.readlines()
streams=[]
for line in lines[1::]:
    line    = re.split(" ",line)
    line    = list(filter(None, line))
    num     = line[0].strip()
    station = line[1].strip()
    ix      = int(line[4])
    iy      = int(line[5])
    stream  = station[0:-7]
    # distto  = float(station[-4::])
    #-------------------------
    if rivermap[iy-1,ix-1] !=1.0:
        continue
    #-------------------------
    if stream in streams:
        continue
    streams.append(stream)
print streams
pname=[]
dista=[]
xlist=[]
ylist=[]
for stream0 in streams:
    # print stream0
    dist0=0.0
    for line in lines[1::]:
        line    = re.split(" ",line)
        line    = list(filter(None, line))
        # print line
        num     = line[0].strip()
        station = line[1].strip()
        lon     = float(line[2])
        lat     = float(line[3])
        ix      = int(line[4])
        iy      = int(line[5])
        ele     = float(line[6])
        ele_dif = float(line[7])
        EGM08   = float(line[8])
        EGM96   = float(line[9])
        sat     = line[10].strip()
        stream  = station[0:-7]
        dstnow  = float(station[-4::])
        # print station, stream, dstnow
        if stream==stream0:
            # print stream, dstnow
            if dstnow>=dist0:
                dist0=dstnow
                loc0 =station
                ix0  =ix-1
                iy0  =iy-1
    # print loc0, dist0,ix0,iy0
    pname.append(loc0)
    dista.append(dist0)
    xlist.append(ix0)
    ylist.append(iy0)
#=========================
pnum=len(pname)
vscover=np.ones([ny,nx],np.float32)*-9999.0
for point in np.arange(pnum):
    print pname[point], dista[point], xlist[point], ylist[point]
    lx, ly = along_river(xlist[point], ylist[point])
    vscover[ly,lx]=1.0

#==== save file ===
vscover.tofile("./dat/satellite_coverage.bin")