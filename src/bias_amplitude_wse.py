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
nameb=argv[4]
CaMa_dir=argv[5]
mapname=argv[6]
figname=argv[7]
indirs=argv[8]
numvar=2
print (namea, nameb, figname)
# obstxt="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+mapname+"_QC0.txt"
obstxt="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_"+mapname+".txt"
#===================
fname=CaMa_dir+"/map/"+mapname+"/params.txt"
f=open(fname,"r")
lines=f.readlines()
f.close()
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
for station in [namea, nameb]:
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
#==================
labels=["CaMa-Flood","HydroWeb"]
#colors=["#ff8021","#34495e"]
#colors=["#f7022a","#0804f9"]
# colors=["#3174a1","#34495e"]
# colors=["#1600fe","#34495e"]
#colors=["#ff8021","#0804f9"]
colors=["red","k"] # #fe0000
#==================
# figure in A4 size
va_margin= 0.0#1.38#inch 
ho_margin= 0.0#1.18#inch
hgt=(11.69 - 2*va_margin)*(1.5/3.0)
wdt=(8.27 - 2*ho_margin)
fig=plt.figure(figsize=(wdt, hgt))
G = gridspec.GridSpec(7,4)
# ax0 = fig.add_subplot(G[0,:])
###########
ax0 = fig.add_subplot(G[0:3,2::])
#cmf
lines=[ax0.plot(np.arange(0,len(wse_cmf[:,0])),wse_cmf[:,0],color=colors[0],label=labels[0],linewidth=0.5)[0]]
#hydroweb
locs, org = hweb.HydroWeb_WSE(namea,syear,eyear)
# print locss
orga = np.array(org)+np.array(legm08[0])-np.array(legm96[0])
lines.append(ax0.plot(locs,orga,color=colors[1],label=labels[1],linestyle='None',linewidth=0,marker="o",fillstyle="none",markersize=5)[0])
if eyear-syear > 1:
    dtt=1
    dt=int(math.ceil(((eyear-syear)+1)/dtt))
else:
    dtt=1
    dt=(eyear-syear)+1
xxlist=np.linspace(0,N,dt,endpoint=True)
xxlab=np.arange(syear,eyear+1,dtt)
ax0.set_xlim(xmin=0,xmax=last+1)
ax0.set_xticks(xxlist)
ax0.set_xticklabels(xxlab,fontsize=6)
ax0.set_ylabel('WSE $(m)$', color='k',fontsize=8)
ax0.tick_params('y',labelsize=6, colors='k')
ax0.set_xlabel('Years', color='k',fontsize=8)
ax0.tick_params('x',labelsize=6, colors='k')
# legend
plt.legend(lines,labels,loc="upper right")
###########
ax1 = fig.add_subplot(G[4::,2::])
#cmf
lines=[ax1.plot(np.arange(0,len(wse_cmf[:,1])),wse_cmf[:,1],color=colors[0],label=labels[0],linewidth=0.5)[0]]
#hydroweb
locs, org = hweb.HydroWeb_WSE(nameb,syear,eyear)
# print locs
orgb = np.array(org)+np.array(legm08[1])-np.array(legm96[1])
lines.append(ax1.plot(locs,orgb,color=colors[1],label=labels[1],linestyle='None',linewidth=0,marker="o",fillstyle="none",markersize=5)[0])
if eyear-syear > 1:
    dtt=1
    dt=int(math.ceil(((eyear-syear)+1)/dtt))
else:
    dtt=1
    dt=(eyear-syear)+1
xxlist=np.linspace(0,N,dt,endpoint=True)
xxlab=np.arange(syear,eyear+1,dtt)
ax1.set_xlim(xmin=0,xmax=last+1)
ax1.set_ylim(ymin=84.0,ymax=98.0+1)
ax1.set_xticks(xxlist)
ax1.set_xticklabels(xxlab,fontsize=6)
ax1.set_ylabel('WSE $(m)$', color='k',fontsize=8)
ax1.tick_params('y',labelsize=6, colors='k')
ax1.set_xlabel('Years', color='k',fontsize=8)
ax1.tick_params('x',labelsize=6, colors='k')
# legend
plt.legend(lines,labels,loc="upper right")
# ax1.text(-0.05,1.05,"%s)"%(string.ascii_lowercase[1]),ha="left",va="center",transform=ax1.transAxes,fontsize=10)
##########
#boxplot
ax2 = fig.add_subplot(G[4::,0])
flierprops = dict(marker='o', markerfacecolor='none', markersize=12,linestyle='none', markeredgecolor='k')
boxprops = dict(color='grey')#facecolor='none'
whiskerprops = dict(color='grey',linestyle="--")
capprops = dict(color='grey')
medianprops = dict(color='grey')
box=ax2.boxplot([wse_cmf[:,1],orgb],labels=labels,boxprops=boxprops,showfliers=False, 
                whiskerprops=whiskerprops,capprops=capprops,medianprops=medianprops, 
                notch=False, sym=None, vert=True, whis=1.5,positions=None, widths=None, 
                patch_artist=True,bootstrap=None, usermedians=None, conf_intervals=None)#flierprops=flierprops,
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
ax2.set_ylabel('WSE $(m)$', color='k',fontsize=8)
ax2.tick_params('y',labelsize=6, colors='k')
ax2.tick_params('x',labelsize=6, colors='k')#,labelrotation=45)
ax2.set_xticklabels(labels,rotation=0)
##########
#pdf
ax3 = fig.add_subplot(G[4::,1])
sns.distplot(wse_cmf[:,1], ax=ax3, hist=True, color=colors[0], label="CaMa-Flood") #ax=ax3,
sns.distplot(orgb, ax=ax3, hist=True, color=colors[1], label="HydroWeb") #ax=ax3,
ax3.set_yticks([0.0,0.1,0.2])
ax3.set_yticklabels([0.0,0.1,0.2],fontsize=6)
ax3.set_ylabel('density', color='k',fontsize=8)
ax3.tick_params('y',labelsize=6, colors='k')
ax3.set_xlabel('WSE $(m)$', color='k',fontsize=8)
ax3.tick_params('x',labelsize=6, colors='k')
#ax3.set_title("Histogram of Bias",fontsize=8)
#ax3.set_xlim(xmin=-20.0,xmax=20.0)
#ax3.text(0.01,0.95,"b",transform=ax3.transAxes,fontsize=8)
# ax3.legend(ncol=2,bbox_to_anchor=(1.0, -0.12), loc=1)
#=======
patches=[]
for i,label in enumerate(labels):
    patches.append(mpatches.Patch(color=colors[i], label=label))
plt.legend(handles=patches,ncol=2,bbox_to_anchor=(-1.32, -0.2), loc="lower left",fontsize=6)
#
# plt.subplots_adjust(wspace=0, hspace=0)
plt.subplots_adjust(wspace=0.5, hspace=0.5)
# plt.tight_layout()
#plt.title(stitle)
plt.savefig("./figures/"+figname+".jpg",dpi=500,bbox_inches="tight", pad_inches=0.05)
plt.savefig("./figures/"+figname+".png",dpi=500,bbox_inches="tight", pad_inches=0.05)
plt.savefig("./figures/"+figname+".pdf",dpi=800,bbox_inches="tight", pad_inches=0.05)