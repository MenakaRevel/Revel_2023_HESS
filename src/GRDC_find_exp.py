#!/opt/local/bin/python
# -*- coding: utf-8 -*-
'''
find GRDC loactions with best statistics
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
# import cal_stat as stat
#import plot_colors as pc

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
#========================================
def ISS(l,u,o,alpha=0.05):
    """
    Interval Skill Score
    input:
        l: lower bound
        u: upper bound
    output:
        ISS: Interval Skill Score 
    """
    u,l = filter_nan(u,l)
    o=ma.masked_where(o<=0.0,o).filled(0.0)
    u=ma.masked_where(o<=0.0,u).filled(0.0)
    l=ma.masked_where(o<=0.0,l).filled(0.0)
    o=np.compress(o>0.0,o)
    u=np.compress(o>0.0,u)
    l=np.compress(o>0.0,l)
    pnum=len(o)
    issalpha=[]
    for i in np.arange(pnum):
        if l[i] <= o[i] <= u[i]:
            issalpha.append(u[i]-l[i])
        elif o[i] < l[i]:
            issalpha.append((u[i]-l[i])+(2.0/alpha)*(l[i]-o[i]))
        elif u[i] < o[i]:
            issalpha.append((u[i]-l[i])+(2.0/alpha)*(o[i]-u[i]))
        else:
            issalpha.append(0.0)
    issalpha=np.array(issalpha)
    return np.sum(issalpha)
#====================================================================
def read_ul(experiment,station):
    u=[]
    l=[]
    fname="./txt/"+experiment+"/Q_interval/"+station+".txt"
    with open(fname,"r") as f:
        lines=f.readlines()
    for line in lines:
        line = re.split(" ",line)
        line = list(filter(None, line))
        try:
            uval=float(line[3])
            lval=float(line[2])
        except:
            uval=0.0
            lval=0.0
        u.append(uval)
        l.append(lval)
    return np.array(u), np.array(l)
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
def read_dis(experiment,station):
    asm=[]
    opn=[]
    fname="./txt/"+experiment+"/outflow/"+station+".txt"
    with open(fname,"r") as f:
        lines=f.readlines()
    for line in lines:
        line = re.split(" ",line)
        line = list(filter(None, line))
        asm.append(float(line[1]))
        opn.append(float(line[2]))
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
        asm.append(float(line[1]))
        opn.append(float(line[2]))
    return np.array(asm), np.array(opn)
#====================================================================
def get_GRDClist(fname):
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
        if rivermap[iy1-1,ix1-1] !=1.0:
            continue
        stationlist.append(num)
    return np.array(stationlist)
#====================================================================
#=====================================
#==================================
# DA_dir="/cluster/data6/menaka/HydroDA"
syear=2009
eyear=2014
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
mapname="amz_06min"

experiments=["DIR_WSE_E2O_HWEB_001","DIR_WSE_E2O_HWEB_002"\
            ,"ANO_WSE_E2O_HWEB_004","ANO_WSE_E2O_HWEB_005"\
            ,"NOM_WSE_E2O_HWEB_004","NOM_WSE_E2O_HWEB_006"]
lexp=len(experiments)
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
#=======
satcov="./dat/satellite_coverage.bin"
satcov=np.fromfile(satcov,np.float32).reshape(ny,nx)
#===================
rivnum="/cluster/data6/menaka/HydroDA/dat/rivnum_"+mapname+".bin"
rivnum=np.fromfile(rivnum,np.int32).reshape(ny,nx)
rivermap=((nextxy[0]>0)*(rivnum==1))*1.0
#======================================================
grdc_satcov=[]
obslist="/cluster/data6/menaka/HydroDA/dat/grdc_"+mapname+".txt"
with open(obslist,"r") as f:
    lines=f.readlines()
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
    org=grdc.grdc_dis(num,syear,eyear)
    org=np.array(org)
    # asm, opn = read_dis(exp,num)
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
    if np.sum((org!=-9999.0)*1.0) < 365*1:
        # print ("no obs: ",stream, np.sum((org!=-9999.0)*1.0), "< 365 days")
        continue
    grdc_satcov.append(num)
#==========================
grdcs=[]
thr=3
for i,exp in enumerate(experiments[::-1]):
    # if i==2 or i==4:
    #    grdcs=[] 
    df0=pd.read_csv(filepath_or_buffer="./tables/realtive_metrics_"+exp+".csv")
    print ("++++++++++++++++++++++++++++++++")
    print (exp)
    # print (df["GRDC ID"][df["delta r"].idxmax()])
    # print (df["GRDC ID"]grdc_satcov)
    # print (df["GRDC ID"](df["GRDC ID"][df["delta r"].idxmax()]))
    # print grdc_satcov
    # print (df0["GRDC ID"])
    # df=df0[df0.index.isin(grdc_satcov)]
    # print (df0.loc[df0["GRDC ID"].isin(grdc_satcov)])
    df=df0.loc[df0["GRDC ID"].isin(grdc_satcov)]
    # print ("***********",df[df["delta r"]>0.01].sort_values(by=["delta r"])[["GRDC ID"]].values)
    # print ("***********",df[df["delta r"]>0.01].sort_values(by=["delta r"])[["GRDC ID"]].values[2,0])
    # print ("***********",df[df["delta r"]>0.01].sort_values(by=["delta r"])[["GRDC ID"]].values[3,0])
    # print ("+++++++++++",df[df["delta r"]>0.01].sort_values(by=["delta r"])[["GRDC ID"]].sample().values[0,0])
    #====================
    # grdc1=df["GRDC ID"][df["delta r"].idxmax()] #.values[0][0]
    grdc1=df[df["delta r"]>-0.01].sort_values(by=["delta r"])[["GRDC ID"]][0:10].sample().values[0,0]
    ii=0
    while grdc1 in grdcs:
        grdc1=df[df["delta r"]>-0.01].sort_values(by=["delta r"])[["GRDC ID"]].values[ii,0]
        ii=ii+1
        if ii > thr:
            grdc1=df["GRDC ID"][df["delta r"].idxmax()]
            break
        if grdc1 in grdcs:
            grdc1=df[df["delta r"]>-0.01].sort_values(by=["delta r"])[["GRDC ID"]].values[-1,0]
    grdcs.append(grdc1)
    #====================
    # grdc2=df["GRDC ID"][df["NSEAI"].idxmax()] #.values[0][0]
    grdc2=df[df["NSEAI"]>-1.0].sort_values(by=["NSEAI"])[["GRDC ID"]][0:10].sample().values[0,0]
    ii=0
    while grdc2 in grdcs:
        grdc2=df[df["NSEAI"]>-1.0].sort_values(by=["NSEAI"])[["GRDC ID"]].values[ii,0]
        ii=ii+1
        if ii > thr:
            grdc1=df["GRDC ID"][df["NSEAI"].idxmax()]
            if grdc1==grdc2:
                grdc2==df["GRDC ID"][df["KGEAI"].idxmax()]
            break
        if grdc2 in grdcs:
            grdc2=df[df["NSEAI"]>-1.0].sort_values(by=["NSEAI"])[["GRDC ID"]].values[-1,0]
    grdcs.append(grdc2)
    #====================
    # grdc3=df["GRDC ID"][df["rISS"].idxmin()] #.values[0][0]
    # if grdc1==grdc3:
    #     grdc3=df["GRDC ID"][df["Sharpness"].idxmin()] #.values[0][0]
    # if grdc2==grdc3:
    #     grdc3=df["GRDC ID"][df["Reliability"].idxmax()] #.values[0][0]
    grdc3=df[df["rISS"]<0.0].sort_values(by=["rISS"])[["GRDC ID"]][0:10].sample().values[0][0]
    ii=0
    while grdc3 in grdcs:
        grdc3=df[df["rISS"]<0.0].sort_values(by=["rISS"])[["GRDC ID"]].values[ii,0]
        ii=ii+1
        if ii > thr:
            grdc1=df["GRDC ID"][df["rISS"].idxmax()]
            if grdc1==grdc3:
                grdc3=df["GRDC ID"][df["Sharpness"].idxmin()] #.values[0][0]
            if grdc2==grdc3:
                grdc3=df["GRDC ID"][df["Reliability"].idxmax()] #.values[0][0]
            break
        if grdc3 in grdcs:
            grdc3=df[df["rISS"]<0.0].sort_values(by=["rISS"])[["GRDC ID"]].values[-1,0]
    grdcs.append(grdc3)
    # print ("****",df.sort_values(by=["rISS"])[["GRDC ID"]][0:10].values[:])
    # print (df.head())
    # # grdc1=df["GRDC ID"][df["delta r"].idxmax()]
    # # grdcs.append(grdc1)
    # # grdc2=df["GRDC ID"][df["NSEAI"].idxmax()]
    # # if grdc1==grdc2:
    # #     grdc2==df["GRDC ID"][df["KGEAI"].idxmax()]
    # # grdc3=df["GRDC ID"][df["rISS"].idxmin()]
    # # if grdc1==grdc3:
    # #     grdc3=df["GRDC ID"][df["Sharpness"].idxmin()]
    # # if grdc2==grdc3:
    # #     grdc3=df["GRDC ID"][df["Reliability"].idxmax()]
    # print ("max delta r                       :", grdc1)
    # print ("max NSEAI/KGEAI                   :", grdc2)
    # print ("better rISS/Sharpness/Reliability :", grdc3)
    print grdc1
    print grdc2
    print grdc3