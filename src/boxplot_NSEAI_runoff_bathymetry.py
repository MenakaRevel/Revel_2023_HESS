#!/opt/local/bin/python
# -*- coding: utf-8 -*-

from cProfile import label
from math import expm1
from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
# import datetime
from matplotlib.colors import LogNorm,Normalize,ListedColormap,BoundaryNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
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
import warnings;warnings.filterwarnings('ignore')
from adjustText import adjust_text
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
import seaborn as sns
import pandas as pd
import warnings;warnings.filterwarnings('ignore')

from statistics import *
import read_grdc as grdc
import my_colorbar as mbar
import read_hydroweb as hweb
from river_patch import patchms
#========================================
#====    functions for read files    ====
#========================================
def mk_dir(sdir):
  try:
    os.makedirs(sdir)
  except:
    pass
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
def read_wse(experiment,station,type0):
    asm=[]
    opn=[]
    fname="./txt/"+experiment+"/wse."+type0+"/"+station+".txt"
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
def read_ul_all(experiment,station):
    ua=[]
    la=[]
    uo=[]
    lo=[]
    # fname="./txt/"+experiment+"/Q_interval/"+station+".txt"
    fname="./txt/"+experiment+"/Q_percentile/"+station+".txt"
    with open(fname,"r") as f:
        lines=f.readlines()
    for line in lines:
        line = re.split(" ",line)
        line = list(filter(None, line))
        try:
            uvala=float(line[1])
            lvala=float(line[2])
            uvalo=float(line[3])
            lvalo=float(line[4])
        except:
            uvala=0.0
            lvala=0.0
            uvalo=0.0
            lvalo=0.0
        ua.append(uvala)
        la.append(lvala)
        uo.append(uvalo)
        lo.append(lvalo)
    return np.array(ua), np.array(la), np.array(uo), np.array(lo)
#====================================================================
def read_expname(exlist):
    "read experiment names"
    with open(exlist,"r") as exf:
        linesexf=exf.readlines()
    indirs=[]
    labels=[]
    for lineexf in linesexf:
        lineexf = re.split(":",lineexf)
        lineexf = list(filter(None, lineexf))
        labels.append(lineexf[0])
        indirs.append(lineexf[1].strip())
    return labels, indirs
#====================================================================
def ret_metric(exp,mapname="amz_06min",sim="assm"):
    lNSEAI=[]
    # GRDC
    obslist="/cluster/data6/menaka/HydroDA/dat/grdc_"+mapname+".txt"
    # print obslist
    with open(obslist,"r") as f:
        lines=f.readlines()
    for line in lines[1::]:
        line    = re.split(";",line)
        line    = list(filter(None, line))
        # print (line)
        num     = line[0].strip()
        basin   = line[1].split(",")[0].strip()
        stream  = line[2].strip()
        ix1     = int(line[3])
        iy1     = int(line[4])
        ix2     = int(line[5])
        iy2     = int(line[6])
        staid   = int(num)
        # print (num, basin ,stream)
        #-------------------------
        if uparea[iy1-1,ix1-1] < 1.0e9:
            continue
        #-------------------------
        if rivermap[iy1-1,ix1-1] !=1.0:
            continue
        #-------------------------
        if satcov[iy1-1,ix1-1] !=1.0:
            continue
        #--------------
        org=grdc.grdc_dis(num,syear,eyear)
        org=np.array(org)
        #--------------
        if np.sum((org!=-9999.0)*1.0) < 365*1:
            # print ("no obs: ",stream, np.sum((org!=-9999.0)*1.0), "< 365 days")
            continue
        #--------------
        # uncalibrated
        asm, opn = read_dis(exp,num)
        # NSEAI
        NSEasm=NS(asm,org)
        NSEopn=NS(opn,org)
        NAI=NSEasm #(NSEasm-NSEopn)/(1.0-NSEopn+1.0e-20)
        # print sim
        if sim=="open":
            # print sim,NSEopn 
            NAI=NSEopn
        lNSEAI.append(NAI)
    return lNSEAI
#====================================================================
def plot_boxplot(lstat0, lstat1, lstat2, lstat3, color1, color2, color3, numch, ylabel, xlabel,ax=None):
    # plot boxplot using seaborn
    ax=ax or plt.gca()
    #
    flierprops = dict(marker='o', markerfacecolor='none', markersize=8,linestyle='none', markeredgecolor='k')
    boxprops = dict(color='grey')#facecolor='none'
    whiskerprops = dict(color='grey',linestyle="--")
    capprops = dict(color='grey')
    medianprops = dict(color='grey',linestyle="-",linewidth=1.0)
    meanprops = dict(marker='D', markeredgecolor='black',markerfacecolor='green',markersize=8)
    #==========
    box=sns.boxplot(ax=ax,data=[lstat0,lstat1,lstat2, lstat3], fliersize=0.0, palette=["w",color1,color2,color3], whis=1.5\
        ,meanline=True, width=0.8, linewidth=0.3, dodge=True\
        ,meanprops=meanprops,capprops=capprops,medianprops=medianprops)
    # median line
    ax.axhline(y=0.0, linestyle="--", color="k",linewidth=0.5)
    # ylabel
    ax.set_ylabel(ylabel,fontsize=8)
    # xlabel
    # ax.set_xlabel(xlabel,fontsize=8)
    ax.tick_params(axis='both', which='major', labelsize=6)
    ax.set_ylim(ymin=-10.2,ymax=1.2)
    ax.set_xticklabels(["Open-loop","Direct", "Anomaly", "Normalized"],rotation = 0)
    ax.text(0.0,1.05,"%s"%(numch),ha="left",va="center",transform=ax.transAxes,fontsize=10)
    # plot_line(lstat1, lstat2, lstat3, ax=ax)
    return 0
#====================================================================
def plot_line(lstato1, lstato2, lstato3, ax=None):
    ax=ax or plt.gca()
    data=[np.median(np.array(lstato)) for lstato in [lstato1, lstato2, lstato3]]
    ax.plot(np.arange(0,len(data)),data,marker="*",linewidth=0.5,linestyle="--",color="grey")
    return
#====================================================================
argv=sys.argv
syear=int(argv[1])
eyear=int(argv[2])
CaMa_dir=argv[3]
mapname=argv[4]
exlist1=argv[5]
exlist2=argv[6]
exlist3=argv[7]
exlist4=argv[8]
figname=argv[9]
ncpus=int(sys.argv[10])
ens_mem=21
seaborn_map=True
numvar=2
#=============================
# read CaMa-Flood files
#=============================
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
nxtdst = CaMa_dir+"/map/"+mapname+"/nxtdst.bin"
nextxy = np.fromfile(nextxy,np.int32).reshape(2,ny,nx)
rivwth = np.fromfile(rivwth,np.float32).reshape(ny,nx)
rivhgt = np.fromfile(rivhgt,np.float32).reshape(ny,nx)
rivlen = np.fromfile(rivlen,np.float32).reshape(ny,nx)
elevtn = np.fromfile(elevtn,np.float32).reshape(ny,nx)
lonlat = np.fromfile(lonlat,np.float32).reshape(2,ny,nx)
uparea = np.fromfile(uparea,np.float32).reshape(ny,nx)
nxtdst = np.fromfile(nxtdst,np.float32).reshape(ny,nx)
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
#==============
# making figure
#==============
# figure in A4 size
# colors=["#3174a1","#3f913a","#bf353c"]
# colors=["#004488","#ddaa33","#ba5566"]
colors=["#D81B60","#FFC107","#004D40"]
# cmap   = matplotlib.colors.ListedColormap(matplotlib.cm.get_cmap("viridis_r").colors[:len(basins)])
va_margin= 0.0#1.38#inch 
ho_margin= 0.0#1.18#inch
hgt=(11.69 - 2*va_margin)*(8.27/11.69)*(2.0/3.0)
wdt=(8.27 - 2*ho_margin)*(8.27/8.27)
#fig=plt.figure(figsize=(8.27,11.69))
fig=plt.figure(figsize=(wdt,hgt))
#fig.suptitle("auto-correalated area")#"maximum autocorrelation length")
#G = gridspec.GridSpec(2,1)
G = gridspec.GridSpec(ncols=2,nrows=2)
ax1 = fig.add_subplot(G[0,0])
ax2 = fig.add_subplot(G[0,1])
ax3 = fig.add_subplot(G[1,0])
ax4 = fig.add_subplot(G[1,1])
#====
# no runoff bias, orignal bathymetry
labels,indirs=read_expname(exlist1)
lNSEAI0= ret_metric(indirs[0],mapname=mapname,sim="open")
lNSEAI1= ret_metric(indirs[0],mapname=mapname,sim="assm")
lNSEAI2= ret_metric(indirs[1],mapname=mapname,sim="assm")
lNSEAI3= ret_metric(indirs[2],mapname=mapname,sim="assm")
plot_boxplot(lNSEAI0, lNSEAI1, lNSEAI2, lNSEAI3, colors[0], colors[1], colors[2], "a)", "$NSE$", "DA methods",ax=ax1)
print "a)", np.median(np.array(lNSEAI0)), np.median(np.array(lNSEAI1)), np.median(np.array(lNSEAI2)), np.median(np.array(lNSEAI3))
ax1.text(-0.20,0.50,"without runoff bias",rotation=90,ha="right",va="center",transform=ax1.transAxes,fontsize=10)
ax1.text(0.5,1.15,"without bathymetry error",rotation=0,ha="center",va="bottom",transform=ax1.transAxes,fontsize=10)
# no runoff bias, courrpted bathymetry
labels,indirs=read_expname(exlist2)
lNSEAI0= ret_metric(indirs[0],mapname=mapname,sim="open")
lNSEAI1= ret_metric(indirs[0],mapname=mapname,sim="assm")
lNSEAI2= ret_metric(indirs[1],mapname=mapname,sim="assm")
lNSEAI3= ret_metric(indirs[2],mapname=mapname,sim="assm")
plot_boxplot(lNSEAI0,lNSEAI1, lNSEAI2, lNSEAI3, colors[0], colors[1], colors[2], "b)", "$NSE$", "DA methods",ax=ax2)
print "b)", np.median(np.array(lNSEAI0)), np.median(np.array(lNSEAI1)), np.median(np.array(lNSEAI2)), np.median(np.array(lNSEAI3))
ax2.text(0.5,1.15,"with bathymetry error",rotation=0,ha="center",va="bottom",transform=ax2.transAxes,fontsize=10)
# runoff bias, orignal bathymetry
labels,indirs=read_expname(exlist3)
lNSEAI0= ret_metric(indirs[0],mapname=mapname,sim="open")
lNSEAI1= ret_metric(indirs[0],mapname=mapname,sim="assm")
lNSEAI2= ret_metric(indirs[1],mapname=mapname,sim="assm")
lNSEAI3= ret_metric(indirs[2],mapname=mapname,sim="assm")
plot_boxplot(lNSEAI0, lNSEAI1, lNSEAI2, lNSEAI3, colors[0], colors[1], colors[2], "c)", "$NSE$", "DA methods",ax=ax3)
print "c)", np.median(np.array(lNSEAI0)), np.median(np.array(lNSEAI1)), np.median(np.array(lNSEAI2)), np.median(np.array(lNSEAI3))
ax3.text(-0.20,0.50,"with runoff bias",rotation=90,ha="right",va="center",transform=ax3.transAxes,fontsize=10)
# runoff bias, courrpted bathymetry
labels,indirs=read_expname(exlist4)
lNSEAI0= ret_metric(indirs[0],mapname=mapname,sim="open")
lNSEAI1= ret_metric(indirs[0],mapname=mapname,sim="assm")
lNSEAI2= ret_metric(indirs[1],mapname=mapname,sim="assm")
lNSEAI3= ret_metric(indirs[2],mapname=mapname,sim="assm")
plot_boxplot(lNSEAI0, lNSEAI1, lNSEAI2, lNSEAI3, colors[0], colors[1], colors[2], "d)", "$NSE$", "DA methods",ax=ax4)
print "d)", np.median(np.array(lNSEAI0)), np.median(np.array(lNSEAI1)), np.median(np.array(lNSEAI2)), np.median(np.array(lNSEAI3))
#==============================
# patches=[]
# for point in np.arange(len(labels)):
#     patches.append(mpatches.Patch(color=colors[point],label=labels[point]))
# #--
# legend1=plt.legend(handles=patches,bbox_to_anchor=(0.5,0.05), loc="upper center",
#            bbox_transform=fig.transFigure, ncol=3,  borderaxespad=0.0, frameon=False)#
#===
plt.savefig("./figures/"+figname+".pdf",dpi=800,bbox_inches="tight", pad_inches=0.05)
plt.savefig("./figures/"+figname+".png",dpi=800,bbox_inches="tight", pad_inches=0.05)
plt.savefig("./figures/"+figname+".jpg",dpi=800,bbox_inches="tight", pad_inches=0.05)