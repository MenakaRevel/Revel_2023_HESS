#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import statistics
import numpy as np
import matplotlib.pyplot as plt
# import datetime
from matplotlib.colors import LogNorm,Normalize,ListedColormap,BoundaryNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import sys
import os
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
import seaborn as sns
import pandas as pd
import warnings;warnings.filterwarnings('ignore')

import read_grdc as grdc
import my_colorbar as mbar
import read_hydroweb as hweb
#========================================
#====  functions for making figures  ====
#========================================
def filter_nan(s,o):
    """
    this functions removed the data  from simulated and observed data
    where ever the observed data contains nan
    """
    data = np.array([s.flatten(),o.flatten()])
    data = np.transpose(data)
    data = data[~np.isnan(data).any(1)]

    return data[:,0],data[:,1]
#====================================================================
def vec_par(LEVEL,ax=None):
    ax=ax or plt.gca()
    txt=re.split("-",figname)[0]+"_%02d.txt"%(LEVEL)
    os.system("./bin/print_rivvec "+tmp0+" 1 "+str(LEVEL)+" > "+txt)
    width=(float(LEVEL)**sup)*w
    #print LEVEL, width#, lon1,lat1,lon2-lon1,lat2-lat1#x1[0],y1[0],x1[1]-x1[0],y1[1]-y1[0]
    # open tmp2.txt
    f = open(txt,"r")
    lines = f.readlines()
    f.close()
    #print LEVEL, width, lines, txt
    #---
    for line in lines:
        line = filter(None, re.split(" ",line))
        lon1 = float(line[0])
        lat1 = float(line[1])
        lon2 = float(line[3])
        lat2 = float(line[4])

        ix = int((lon1 - west)*(1/gsize))
        iy = int((-lat1 + north)*(1/gsize))

        if rivermap[iy-1,ix-1] == 0:
            continue

        if lon1-lon2 > 180.0:
            print (lon1, lon2)
            lon2=180.0
        elif lon2-lon1> 180.0:
            print (lon1,lon2)
            lon2=-180.0
        #--------
        colorVal="grey"#"k" 
        #print (lon1,lon2,lat1,lat2,width)
        plot_ax(lon1,lon2,lat1,lat2,width,colorVal,ax)
#====================================================================
def plot_ax(lon1,lon2,lat1,lat2,width,colorVal,ax=None):
    ax=ax or plt.gca()
    return ax.plot([lon1,lon2],[lat1,lat2],color=colorVal,linewidth=width,zorder=105,alpha=alpha)
#====================================================================
def mk_dir(sdir):
  try:
    os.makedirs(sdir)
  except:
    pass
#====================================================================
def NS(s,o):
    """
    Nash Sutcliffe efficiency coefficient
    input:
        s: simulated
        o: observed
    output:
        ns: Nash Sutcliffe efficient coefficient
    """
    s,o = filter_nan(s,o)
    o=ma.masked_where(o<=0.0,o).filled(0.0)
    s=ma.masked_where(o<=0.0,s).filled(0.0)
    o=np.compress(o>0.0,o)
    s=np.compress(o>0.0,s) 
    return 1 - sum((s-o)**2)/(sum((o-np.mean(o))**2)+1e-20)
#====================================================================
def correlation(s,o):
    """
    correlation coefficient
    input:
        s: simulated
        o: observed
    output:
        correlation: correlation coefficient
    """
    s,o = filter_nan(s,o)
    o=ma.masked_where(o<=0.0,o).filled(0.0)
    s=ma.masked_where(o<=0.0,s).filled(0.0)
    o=np.compress(o>0.0,o)
    s=np.compress(o>0.0,s)
    if s.size == 0:
        corr = 0.0 #np.NaN
    else:
        corr = np.corrcoef(o, s)[0,1]
        
    return corr
#========================================
def KGE(s,o):
    """
	Kling Gupta Efficiency (Kling et al., 2012, http://dx.doi.org/10.1016/j.jhydrol.2012.01.011)
	input:
        s: simulated
        o: observed
    output:
        KGE: Kling Gupta Efficiency
    """
    o=ma.masked_where(o<=0.0,o).filled(0.0)
    s=ma.masked_where(o<=0.0,s).filled(0.0)
    o=np.compress(o>0.0,o)
    s=np.compress(o>0.0,s)
    s,o = filter_nan(s,o)
    B = np.mean(s) / np.mean(o)
    y = (np.std(s) / np.mean(s)) / (np.std(o) / np.mean(o))
    r = np.corrcoef(o, s)[0,1]
    return 1 - np.sqrt((r - 1) ** 2 + (B - 1) ** 2 + (y - 1) ** 2)
#==========================================================
def RMSE(s,o):
    """
    Root Mean Squre Error
    input:
        s: simulated
        o: observed
    output:
        RMSE: Root Mean Squre Error
    """
    o=ma.masked_where(o==-9999.0,o).filled(0.0)
    s=ma.masked_where(o==-9999.0,s).filled(0.0)
    o=np.compress(o>0.0,o)
    s=np.compress(o>0.0,s)
    s,o = filter_nan(s,o)
    # return np.sqrt(np.mean((s-o)**2))
    return np.sqrt(np.ma.mean(np.ma.masked_where(o<=0.0,(s-o)**2)))
#====================================================================
def read_dis(experiment,station):
    asm=[]
    opn=[]
    fname="./txt/"+experiment+"/outflow/"+station+".txt"
    with open(fname,"r") as f:
        lines=f.readlines()
    for line in lines:
        line    = re.split(" ",line)
        line    = list(filter(None, line))
        try:
            asmval  = float(line[1])
        except:
            asmval  = 0.0
        # asmval  = float(line[1])
        opnval  = float(line[2])
        asm.append(asmval)
        opn.append(opnval)
    return np.array(asm), np.array(opn)
#====================================================================
def get_GRDClist(fname,satellite=False):
    # obslist="/cluster/data6/menaka/HydroDA/dat/grdc_"+mapname+".txt"
    with open(fname,"r") as f:
        lines=f.readlines()
    stationlist=[]
    for line in lines[1::]:
        line    = re.split(";",line)
        line    = list(filter(None, line))
        # print (line)
        num     = line[0].strip()
        # basin   = line[1].strip()
        # stream  = line[2].strip()
        ix1     = int(line[3])
        iy1     = int(line[4])
        # ix2     = int(line[5])
        # iy2     = int(line[6])
        # staid   = int(num)
        #-------------------------
        if uparea[iy1-1,ix1-1] < 1.0e9:
            continue
        if rivermap[iy1-1,ix1-1] !=1.0:
            continue
        if satellite:
            if satcov[iy1-1,ix1-1] !=1.0:
                continue
        org=grdc.grdc_dis(num,syear,eyear)
        org=np.array(org)
        if np.sum((org!=-9999.0)*1.0) < 365*1:
            continue
        stationlist.append(num)
    return np.array(stationlist)
#===================
argv=sys.argv
syear=int(argv[1])
eyear=int(argv[2])
CaMa_dir=argv[3]
mapname=argv[4]
figname=argv[5]
# numvar=2
# print (lexp, figname)
#=== read experiment names ===
exlist="./Fig08-experiment_list.nam"
with open(exlist,"r") as exf:
    linesexf=exf.readlines()
experiments=[]
labels=[]
for lineexf in linesexf:
    lineexf = re.split(":",lineexf)
    lineexf = list(filter(None, lineexf))
    labels.append(lineexf[0])
    experiments.append(lineexf[1].strip())
lexp=len(experiments)
statistics=["r","NSE","KGE"]
statistics0=["$r$","$NSE$","$KGE$"]
lstat=len(statistics)
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
rivnum="/cluster/data6/menaka/HydroDA/dat/rivnum_"+mapname+".bin"
rivnum=np.fromfile(rivnum,np.int32).reshape(ny,nx)
rivermap=((nextxy[0]>0)*(rivnum==1))*1.0
#----
satcov="./dat/satellite_coverage.bin"
satcov=np.fromfile(satcov,np.float32).reshape(ny,nx)
#======================================================
#files for making river network
tmp0=re.split("-",figname)[0]+".txt"

land="#C0C0C0"
water="#FFFFFF"

londiff=(east-west)*4
latdiff=(north-south)*4

npix=(90-north)*4
spix=(90-south)*4
wpix=(180+west)*4
epix=(180+east)*4

cmap=mbar.colormap("H01")
vmin=0.0
vmax=1.0
# norm=Normalize(vmin=vmin,vmax=vmax)
norm=BoundaryNorm(np.arange(vmin,vmax+0.1,0.1),cmap.N)

#------
# river width
sup=3 #2
w=0.002 #0.02
alpha=1
width=0.5
#------
west0=west
east0=east
south0=south-2
north0=north+2
#------
print (west,east,south,north)
resol=1
# figure in A4 size
va_margin= 0.0#1.38#inch 
ho_margin= 0.0#1.18#inch
hgt=(11.69 - 2*va_margin)*(3.0/5.0)
wdt=(8.27 - 2*ho_margin)
fig=plt.figure(figsize=(wdt,hgt))
G = gridspec.GridSpec(nrows=lexp,ncols=lstat)
G.update(wspace=0.0, hspace=0.0) # set the spacing between axes. 
#=========================
ax=[]
metric=[]
for i,exp in enumerate(experiments):
    metric1=[]
    for j,stat in enumerate(statistics):
        ax.append(fig.add_subplot(G[i,j]))
        m = Basemap(projection='cyl',llcrnrlat=south0,urcrnrlat=north0,llcrnrlon=west0,urcrnrlon=east0, lat_ts=0,resolution='c',ax=plt.gca())
        # amazon shapefile
        m.readshapefile('./etc/Amazon_basin/Amazon_boundry', 'Amazon_boundry', drawbounds = False)
        for info, shape in zip(m.Amazon_boundry, m.Amazon_boundry):
            # if info['nombre'] == 'Amazon_boundry':
            x, y = zip(*shape) 
            m.plot(x, y, marker=None,color='grey',linewidth=1.0)
        box="%f %f %f %f"%(west,east,north,south) 
        os.system("./bin/txt_vector "+box+" "+CaMa_dir+" "+mapname+" > "+tmp0)  
        #map(vec_par,np.arange(1,10+1,1))
        map(vec_par,np.arange(2,10+1,1))
        # map(vec_par,np.arange(5,10+1,1))
        #--
        obslist="/cluster/data6/menaka/HydroDA/dat/grdc_"+mapname+".txt"
        with open(obslist,"r") as f:
            lines=f.readlines()
        metric_frag=[]
        for line in lines[1::]:
            line    = re.split(";",line)
            line    = list(filter(None, line))
            # print (line)
            num     = line[0].strip()
            basin   = line[1].strip()
            stream  = line[2].strip()
            ix1     = int(line[3])
            iy1     = int(line[4])
            ix2     = int(line[5])
            iy2     = int(line[6])
            staid   = int(num)
            # print (num, basin ,stream)
            org=grdc.grdc_dis(num,syear,eyear)
            org=np.array(org)
            asm, opn = read_dis(exp,num)
            # CORopn=correlation(opn,org)
            # DCOR=CORasm-CORopn
            #-------------------------
            if uparea[iy1-1,ix1-1] < 1.0e9:
                continue
            #-------------------------
            if rivermap[iy1-1,ix1-1] !=1.0:
                continue
            # #-------------------------
            # if satcov[iy1-1,ix1-1] !=1.0:
            #     continue
            #--------------
            if np.sum((org!=-9999.0)*1.0) < 365*1:
                # print ("no obs: ",stream, np.sum((org!=-9999.0)*1.0), "< 365 days")
                continue
            if stat=="r":
                val=correlation(asm,org)
            elif stat=="NSE":
                val=NS(asm,org)
            elif stat=="KGE":
                val=KGE(asm,org)
            metric_frag.append(val)
            #lon=lon0+ix*gsize
            #lat=lat0-iy*gsize
            lon=lonlat[0,iy1-1,ix1-1]
            lat=lonlat[1,iy1-1,ix1-1]
            c=cmap(norm(val))
            # plt.gca().scatter(lon,lat,s=15,marker="o",edgecolors="k",linewidth=0.5, facecolors=c,zorder=106)
            if satcov[iy1-1,ix1-1] !=1.0:
                plt.gca().scatter(lon,lat,s=15,marker="d",edgecolors="k",linewidth=0.3, facecolors=c,zorder=106)
            else:
                plt.gca().scatter(lon,lat,s=15,marker="o",edgecolors="k",linewidth=0.3, facecolors=c,zorder=106)
        #===
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['bottom'].set_visible(False)
        plt.gca().spines['left'].set_visible(False)
        if i==0:
            plt.gca().text(0.5,1.05,"%s"%(statistics0[j]),ha="center",va="center",transform=plt.gca().transAxes,fontsize=10)
        # ax[lstat*i+j].text(-0.05,1.05,"%s)"%(string.ascii_lowercase[2*i]),ha="left",va="center",transform=ax[2*i].transAxes,fontsize=10)
        if j==0:
            plt.gca().text(-0.05,1.05,"%s)"%(string.ascii_lowercase[i]),ha="left",va="center",transform=plt.gca().transAxes,fontsize=10)
            plt.gca().text(-0.05,0.60,"%s"%(labels[i]),ha="center",va="center",transform=plt.gca().transAxes,fontsize=10,rotation=90)
        metric1.append(metric_frag)
    metric.append(metric1)
#======================
im=plt.scatter([],[],c=[],cmap=cmap,s=0.1,vmin=vmin,vmax=vmax,norm=norm,zorder=101)
im.set_visible(False)
cax=fig.add_axes([0.9,0.1,0.01,0.8])
cbar=plt.colorbar(im,cax=cax,ticks=np.arange(vmin,vmax+0.1,0.2),orientation='vertical',extend="min")#np.arange(-0.5,0.5+0.001,0.25)
# cbarQ.set_label("$%s$") #$(r_{assim} - r_{open})$") correlation
cbar.ax.tick_params(labelsize=6.0)
#======================
# legend 
feature1=mlines.Line2D([], [], color='grey', marker='d',
                          markersize=5, label='Outside Satellite Coverage', linewidth=0.0)
feature2=mlines.Line2D([], [], color='grey', marker='o',
                          markersize=5, label='Inside Satellite Coverage',linewidth=0.0)
legend=plt.legend(handles=[feature1,feature2],bbox_to_anchor=(0.5,0.10), loc="lower center",
           bbox_transform=fig.transFigure, ncol=2,  borderaxespad=0.0, frameon=False)#
#======================
plt.subplots_adjust(wspace=0, hspace=0)
#plt.title(stitle)
plt.savefig("./figures/"+figname+".jpg",dpi=800,bbox_inches="tight", pad_inches=0.0)
plt.savefig("./figures/"+figname+".png",dpi=800,bbox_inches="tight", pad_inches=0.0)
plt.savefig("./figures/"+figname+".pdf",dpi=800,bbox_inches="tight", pad_inches=0.0)
os.system("rm -r "+re.split("-",figname)[0]+"*.txt")
#=======================
# # calculate statistics
# metric=np.array(metric)
# # All
# obslist="/cluster/data6/menaka/HydroDA/dat/grdc_"+mapname+".txt"
# stationlist=get_GRDClist(obslist)
# df1=pd.DataFrame(metric[0].T,columns=statistics,index=stationlist)
# df2=pd.DataFrame(metric[1].T,columns=statistics,index=stationlist)
# pf =pd.Panel({labels[0]:df1,labels[1]:df2})
# print(labels[0])
# print(df1.describe())
# print(df1[df1>0.0].count())
# print(labels[1])
# print(df2.describe())
# print(df2[df2>0.0].count())
# # print(pf.describe())
# # satellite coverage
# obslist="/cluster/data6/menaka/HydroDA/dat/grdc_"+mapname+".txt"
# stationlist1=get_GRDClist(obslist,satellite=True)
# print(labels[0])
# print(df1.loc[stationlist1].describe())
# print(df1.loc[stationlist1][df1.loc[stationlist1]>0.0].count())
# print(labels[1])
# print(df2.loc[stationlist1].describe())
# print(df2.loc[stationlist1][df2.loc[stationlist1]>0.0].count())