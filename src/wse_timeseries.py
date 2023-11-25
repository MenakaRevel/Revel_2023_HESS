#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import datetime
from matplotlib.colors import LogNorm,Normalize,ListedColormap,BoundaryNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import sys
import os
import math
# import calendar
# from multiprocessing import Pool
# from multiprocessing import Process
# from multiprocessing import sharedctypes
from numpy import ma
import re
import my_colorbar as mbar
import cartopy.crs as ccrs
import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
import os
import string
import warnings
import seaborn as sns
warnings.filterwarnings('ignore')

import read_grdc as grdc
import my_colorbar as mbar
import read_hydroweb as hweb
from read_sfcelv import read_sfcelv, read_sfcelv_multi

#====================================================================
def read_wse_multi(ix, iy, syear, eyear, runoff_folder, nx, ny):
    #print ix1,iy1
    wse = np.zeros( (len(ix), nbdays), 'f')
    for year in range(syear, eyear+1):
        #print year
        s_days = int( (datetime.date(year , 1,1) - datetime.date(syear, 1, 1)). days)
        e_days = int( (datetime.date(year+1, 1, 1) - datetime.date(syear, 1, 1)). days)

        f = runoff_folder + '/sfcelv'+str(year)+'.bin'
        #print f, len(ix), len(iy)
        tmp = read_sfcelv_multi( ix, iy, e_days-s_days, f, nx, ny)

        #print year, e_days - s_days, s_days, e_days, outflw.shape
        wse[:,s_days:e_days] = tmp

    return wse #, wse_max, wse_min, wse_max_loc, wse_min_loc
#====================================================================
argv=sys.argv
syear=int(argv[1])
eyear=int(argv[2])
namea=argv[3]
CaMa_dir=argv[4]
mapname=argv[5]
figname=argv[6]
indirs=argv[7]
numvar=2
print (namea, figname)
# obstxt="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+mapname+"_QC0.txt"
obstxt="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_"+mapname+".txt"
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
# #higher resolution data
# catmxy = pm.CaMa_dir()+"/map/"+pm.mapname()+"/1min/1min.catmxy.bin"
# catmxy = np.fromfile(catmxy,np.int16).reshape(2,ny*60,nx*60)
#----
# rivnum="/cluster/data6/menaka/HydroDA/dat/rivnum_"+mapname+".bin"
# rivnum=np.fromfile(rivnum,np.int32).reshape(ny,nx)
# rivermap=((nextxy[0]>0)*(rivnum==1))*1.0
#----
# satcov="./dat/satellite_coverage.bin"
# satcov=np.fromfile(satcov,np.float32).reshape(ny,nx)
#======================================================
start_dt=datetime.date(syear,1,1)
end_dt=datetime.date(eyear,12,31)
size=60

start=0
last=(end_dt-start_dt).days + 1
N=int(last)
nbdays=int(last)
#======================================================
# HydroWeb 
IX=[]; IY=[];legm08=[];legm96=[]
for station in [namea]:
    ix,iy,EGM08,EGM96=hweb.get_hydroweb0(station,fname=obstxt)
    IX.append(ix)
    IY.append(iy)
    legm08.append(EGM08)
    legm96.append(EGM96)
#===================
# CaMa-Flood
wse_cmf = read_wse_multi(IX, IY, syear, eyear, indirs, nx, ny)
wse_cmf = np.array(wse_cmf)
wse_cmf = wse_cmf.T
#===================
#hydroweb
locs, org = hweb.HydroWeb_WSE(namea,syear,eyear)
#==================
def transform_data(cmf,method="dir"):
    if method=="dir":
        trans_wse=cmf
    if method=="ano":
        trans_wse=cmf-np.mean(cmf)
    if method=="nom":
        trans_wse=(cmf-np.mean(cmf))/(np.std(cmf)+1e-20)
    return trans_wse
#==================
def mk_plot(locs, org, cmf, method="dir", ax=None):
    ax=ax or plt.gca()
    labels=["CaMa-Flood","HydroWeb"]
    colors=["red","k"] # #fe0000
    #cmf
    trans_wse=transform_data(cmf,method)
    lines=[ax.plot(np.arange(0,len(cmf)),trans_wse,color=colors[0],label=labels[0],linewidth=1.0)[0]]
    #observation
    # print locss
    orga      = np.array(org)+np.array(legm08[0])-np.array(legm96[0])
    trans_org = transform_data(orga,method)
    lines.append(ax.plot(locs,trans_org,color=colors[1],label=labels[1],linestyle='None',linewidth=0,marker="o",fillstyle="none",markersize=8)[0])
    if eyear-syear > 1:
        dtt=1
        dt=int(math.ceil(((eyear-syear)+1)/dtt))
    else:
        dtt=1
        dt=(eyear-syear)+1
    xxlist=np.linspace(0,N,dt,endpoint=True)
    xxlab=np.arange(syear,eyear+1,dtt)
    ax.set_xlim(xmin=0,xmax=last+1)
    ax.set_xticks(xxlist)
    ax.set_xticklabels(xxlab,fontsize=6)
    ax.set_ylabel('WSE $(m)$', color='k',fontsize=12)
    ax.tick_params('y',labelsize=6, colors='k')
    # ax.set_xlabel('Years', color='k',fontsize=8)
    ax.tick_params('x',labelsize=6, colors='k')
    # legend
    # plt.legend(lines,labels,loc="upper right")
    #====================
    # plt.show()
    return 0
#====================
if __name__ == "__main__":
    #==================
    # figure in A4 size
    va_margin= 0.0#1.38#inch 
    ho_margin= 0.0#1.18#inch
    hgt=(11.69 - 2*va_margin)*(3.0/3.0)
    wdt=(8.27 - 2*ho_margin)
    fig=plt.figure(figsize=(wdt, hgt))
    G = gridspec.GridSpec(3,1)
    ############################
    ax1 = fig.add_subplot(G[0,0])
    mk_plot(locs, org, wse_cmf, method="dir")
    ax2 = fig.add_subplot(G[1,0])
    mk_plot(locs, org, wse_cmf, method="ano")
    ax3 = fig.add_subplot(G[2,0])
    mk_plot(locs, org, wse_cmf, method="nom")
    #============================
    ax3.set_xlabel('Years', color='k',fontsize=12)
    #============================
    # legend
    labels=["CaMa-Flood","Observation"]
    colors=["red","k"] 
    lines=[plt.plot([],[],color=colors[0],label=labels[0],linewidth=0.5)[0]]
    lines.append(plt.plot([],[],color=colors[1],label=labels[1],linestyle='None',linewidth=0,marker="o",fillstyle="none",markersize=5)[0])
    plt.legend(lines, labels, loc="upper right") #loc="lower center", bbox_to_anchor=(0.5,-.05))
    #============================
    plt.savefig("./pptimage/"+figname+".jpg",dpi=500,bbox_inches="tight", pad_inches=0.05)