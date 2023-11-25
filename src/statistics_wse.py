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
import datetime
import calendar
import string
import seaborn as sns
import pandas as pd
import warnings;warnings.filterwarnings('ignore')

import read_grdc as grdc
import my_colorbar as mbar
import read_hydroweb as hweb
from statistics import *
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
'''def NS(s,o):
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
#==========================================================
def NRMSE(s,o):
    """
    Normalised Root Mean Squre Error
    input:
        s: simulated
        o: observed
    output:
        NRMSE: Normalised Root Mean Squre Error
    """
    o=ma.masked_where(o==-9999.0,o).filled(0.0)
    s=ma.masked_where(o==-9999.0,s).filled(0.0)
    o=np.compress(o>0.0,o)
    s=np.compress(o>0.0,s)
    s,o = filter_nan(s,o)
    # return np.sqrt(np.mean((s-o)**2))
    return np.sqrt(np.ma.mean(np.ma.masked_where(o<=0.0,(s-o)**2)))/np.mean(o)
#==========================================================
def pBIAS(s,o):
    """
    Precentage Bias
    input:
        s: simulated
        o: observed
    output:
        pBias: Precentage Bias
    """
    s,o = filter_nan(s,o)
    o=ma.masked_where(o==-9999.0,o).filled(0.0)
    s=ma.masked_where(o==-9999.0,s).filled(0.0)
    o=np.compress(o>0.0,o)
    s=np.compress(o>0.0,s)
    return abs((np.mean(s)-np.mean(o))/np.mean(o))
#==========================================================
def BIAS(s,o):
    """
    Bias
    input:
        s: simulated
        o: observed
    output:
        Bias: Bias
    """
    s,o = filter_nan(s,o)
    o=ma.masked_where(o==-9999.0,o).filled(0.0)
    s=ma.masked_where(o==-9999.0,s).filled(0.0)
    o=np.compress(o>0.0,o)
    s=np.compress(o>0.0,s)
    # o=np.compress(o==-9999.0,o)
    # s=np.compress(o==-9999.0,s)
    return abs((np.mean(s)-np.mean(o)))
#====================================================================
def rAmplitude(s,o,syear=2009,eyear=2014):
    """
    Realtive Amplitude Difference
    input:
        s: simulated
        o: observed
    output:
        rAmplitude: Realtive Amplitude Difference
    """ 
    s,o = filter_nan(s,o)
    # o=ma.masked_where(o==-9999.0,o).filled(0.0)
    # s=ma.masked_where(o==-9999.0,s).filled(0.0)
    # o=np.compress(o>0.0,o)
    # s=np.compress(o>0.0,s)
    alti_diff=[]
    for year in np.arange(syear,eyear+1):
        date1 = datetime.date(year, 1, 1) # month and day are 1-base
        date2 = datetime.date(year, 12, 31)
        days  = (date2-date1).days + 1
        if year == syear:
            st_dt = 0
            ed_dt = days
        elif year == eyear:
            st_dt = ed_dt
            ed_dt = -1
        else:
            st_dt = ed_dt
            ed_dt = ed_dt + days
        maxloc = np.argmax(s[st_dt:ed_dt])
        minloc = np.argmin(s[st_dt:ed_dt])
        smax   = np.amax(s[st_dt:ed_dt])
        smin   = np.amin(s[st_dt:ed_dt])
        maxloc1=max(maxloc-15,0)
        maxloc2=min(maxloc+15,len(s)-1)
        minloc1=max(minloc-15,0)
        minloc2=min(minloc+15,len(s)-1)
        maxarry=ma.masked_equal(o[st_dt:ed_dt],-9999.0).filled(-9999.0)
        minarry=ma.masked_equal(o[st_dt:ed_dt],-9999.0).filled(9999.0)
        omax   = np.amax(maxarry[maxloc1:maxloc2])
        omin   = np.amin(minarry[minloc1:minloc2])
        if omax == -9999.0 or omin == 9999.0:
            continue
        alti_diff.append((smax-smin)-(omax-omin))
    return np.mean(alti_diff)
'''
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
#====================================================================
def get_HydroWeblist(fname):
    with open(fname,"r") as f:
        lines=f.readlines()
    stationlist=[]
    for line in lines[1::]:
        line    = re.split(" ",line)
        line    = list(filter(None, line))
        # print line
        num     = line[0].strip()
        station = line[1].strip()
        # lon     = float(line[2])
        # lat     = float(line[3])
        ix      = int(line[4])
        iy      = int(line[5])
        # ele     = float(line[6])
        # ele_dif = float(line[7])
        # EGM08   = float(line[8])
        # EGM96   = float(line[9])
        # sat     = line[10].strip()
        #-------------------------
        if rivermap[iy-1,ix-1] !=1.0:
            continue
        #--------------
        stationlist.append(station)
    return stationlist
#====================================================================
def get_HydroWebAll():
    stationlist=[]
    obslist="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+mapname+"_QC0_simulation.txt"
    with open(obslist,"r") as f:
        lines=f.readlines()
    for line in lines[1::]:
        line    = re.split(" ",line)
        line    = list(filter(None, line))
        # print line
        num     = line[0].strip()
        station = line[1].strip()
        # lon     = float(line[2])
        # lat     = float(line[3])
        ix      = int(line[4])
        iy      = int(line[5])
        # ele     = float(line[6])
        # ele_dif = float(line[7])
        # EGM08   = float(line[8])
        # EGM96   = float(line[9])
        # sat     = line[10].strip()
        #-------------------------
        if rivermap[iy-1,ix-1] !=1.0:
            continue
        #--------------
        stationlist.append(station)
    #==============================
    obslist="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+mapname+"_QC0_validation.txt"
    with open(obslist,"r") as f:
        lines=f.readlines()
    for line in lines[1::]:
        line    = re.split(" ",line)
        line    = list(filter(None, line))
        # print line
        num     = line[0].strip()
        station = line[1].strip()
        # lon     = float(line[2])
        # lat     = float(line[3])
        ix      = int(line[4])
        iy      = int(line[5])
        # ele     = float(line[6])
        # ele_dif = float(line[7])
        # EGM08   = float(line[8])
        # EGM96   = float(line[9])
        # sat     = line[10].strip()
        #-------------------------
        if rivermap[iy-1,ix-1] !=1.0:
            continue
        #--------------
        stationlist.append(station)
    return stationlist
#===================
# argv=sys.argv
# syear=int(argv[1])
# eyear=int(argv[2])
# CaMa_dir=argv[3]
# mapname=argv[4]
# figname=argv[5]
syear=2009
eyear=2014
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"
mapname="amz_06min"
# numvar=2
# print (lexp, figname)
# #=== read experiment names ===
# exlist="./Fig09-experiment_list.nam"
# with open(exlist,"r") as exf:
#     linesexf=exf.readlines()
# experiments=[]
# labels=[]
# for lineexf in linesexf:
#     lineexf = re.split(":",lineexf)
#     lineexf = list(filter(None, lineexf))
#     labels.append(lineexf[0])
#     experiments.append(lineexf[1].strip())
# experiments=["DIR_WSE_E2O_HWEB_001","DIR_WSE_E2O_HWEB_002"\
#             ,"ANO_WSE_E2O_HWEB_004","ANO_WSE_E2O_HWEB_005"\
#             ,"NOM_WSE_E2O_HWEB_004","NOM_WSE_E2O_HWEB_006"]
experiments=["DIR_WSE_E2O_HWEB_001","ANO_WSE_E2O_HWEB_004","NOM_WSE_E2O_HWEB_004"]
lexp=len(experiments)
# statistics=["NRMSE","pBIAS","rAmplitude"]
# statistics0=["$NRMSE$","$pBIAS$","$\Delta$$A$"]
statistics=["RMSE","BIAS","rAmplitude"]
statistics0=["$RMSE$","$BIAS$","$\Delta$$A$"]
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
for i,exp in enumerate(experiments):
    metric=[]
    for j,stat in enumerate(statistics):
        obslist="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+mapname+"_QC0_simulation.txt"
        with open(obslist,"r") as f:
            lines=f.readlines()
        metric_frag=[]
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
            #-------------------------
            if rivermap[iy-1,ix-1] !=1.0:
                continue
            #--------------
            org=hweb.HydroWeb_continous_WSE(station,syear,1,1,eyear,12,31,EGM08,EGM96)
            org=np.array(org)#+np.array(EGM08)-np.array(EGM96)
            asm, opn = read_wse(exp,station,"simulation")
            # if stat=="NRMSE":
            #     val=NRMSE(asm,org)
            # elif stat=="pBIAS":
            #     val=pBIAS(asm,org)
            # elif stat=="rAmplitude":
            #     val=rAmplitude(asm,org)
            if stat=="RMSE":
                val=RMSE(asm,org)
            elif stat=="BIAS":
                val=BIAS(asm,org)
            elif stat=="rAmplitude":
                # val=rAmplitude(asm,org)
                val=dAmplitude(asm,org)
            #======================
            metric_frag.append(val) 
        obslist="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+mapname+"_QC0_validation.txt"
        with open(obslist,"r") as f:
            lines=f.readlines()
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
            #-------------------------
            if rivermap[iy-1,ix-1] !=1.0:
                continue
            #--------------
            org=hweb.HydroWeb_continous_WSE(station,syear,1,1,eyear,12,31,EGM08,EGM96)
            org=np.array(org)#+np.array(EGM08)-np.array(EGM96)
            asm, opn = read_wse(exp,station,"validation")
            # if stat=="NRMSE":
            #     val=NRMSE(asm,org)
            # elif stat=="pBIAS":
            #     val=pBIAS(asm,org)
            # elif stat=="rAmplitude":
            #     val=rAmplitude(asm,org)
            if stat=="RMSE":
                val=RMSE(asm,org)
            elif stat=="BIAS":
                val=BIAS(asm,org)
            elif stat=="rAmplitude":
                # val=rAmplitude(asm,org)
                val=dAmplitude(asm,org)
            #======================
            metric_frag.append(val)
        metric.append(metric_frag)
    #=====================
    # calculate statistics
    metric=np.array(metric)
    # All
    #"RMSE","BIAS","rAmplitude"
    obslist="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+mapname+"_QC0.txt"
    # stationlist=get_HydroWeblist(obslist)
    stationlist=get_HydroWebAll()
    df=pd.DataFrame(metric.T,columns=statistics,index=stationlist)
    print(len(stationlist))
    print(exp)
    print("===============================")
    print("    All     ")
    # print(df.describe())
    # print(df[df>0.0].count())
    print(len(stationlist))
    print(df.median(axis=0))
    # print(df.mean(axis=0))
    # for i in [5.0,7.0,10.0,15.0]:
    #     ii="%3.1f"%(i)
    #     print("RMSE<"+ii,df["RMSE"][df["RMSE"]<i].count())
    # for i in [1.0,2.0,3.0,5.0,10.0]:
    #     ii="%3.1f"%(i)
    #     print("|BIAS|<"+ii,df["BIAS"][(df["BIAS"]>-1.0*i) & (df["BIAS"]<1.0*i)].count())
    #     print("|rAmplitude|<"+ii,df["rAmplitude"][(df["rAmplitude"]>-1.0*i) & (df["rAmplitude"]<1.0*i)].count())

    # for assimilation
    obslist="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+mapname+"_QC0_simulation.txt"
    stationlist1=get_HydroWeblist(obslist)
    print("===============================")
    print("    Assimilation     ")
    print(len(stationlist1))
    # print(df.loc[stationlist1].describe())
    print(df.loc[stationlist1].median(axis=0))
    # print(df.loc[stationlist1].mean(axis=0))
    # print(df.loc[stationlist1]["BIAS"].median(axis=0))
    # print(df.loc[stationlist1]["rAmplitude"].median(axis=0))
    # for i in [5.0,7.0,10.0,15.0]:
    #     ii="%3.1f"%(i)
    #     print("RMSE<"+ii,df.loc[stationlist1]["RMSE"][df.loc[stationlist1]["RMSE"]<i].count())
    # for i in [1.0,2.0,3.0,5.0,10.0]:
    #     ii="%3.1f"%(i)
    #     print("|BIAS|<"+ii,df.loc[stationlist1]["BIAS"][df.loc[stationlist1]["BIAS"][(df["BIAS"]>-1.0*i) & (df["BIAS"]<1.0*i)]].count())
    #     print("|rAmplitude|<"+ii,df.loc[stationlist1]["rAmplitude"][df.loc[stationlist1]["rAmplitude"][(df["rAmplitude"]>-1.0*i) & (df["rAmplitude"]<1.0*i)]].count())
    # for validation
    obslist="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+mapname+"_QC0_validation.txt"
    stationlist2=get_HydroWeblist(obslist)
    print("===============================")
    print("    Validation     ")
    print(len(stationlist2))
    # print(df.loc[stationlist1].describe())
    print(df.loc[stationlist2].median(axis=0))
    # print(df.loc[stationlist2].mean(axis=0))
    # print(df.loc[stationlist1]["BIAS"].median(axis=0))
    # print(df.loc[stationlist1]["rAmplitude"].median(axis=0))
    # for i in [5.0,7.0,10.0,15.0]:
    #     ii="%3.1f"%(i)
    #     print("RMSE<"+ii,df.loc[stationlist2]["RMSE"][df.loc[stationlist2]["RMSE"]<i].count())
    # for i in [1.0,2.0,3.0,5.0,10.0]:
    #     ii="%3.1f"%(i)
    #     print("|BIAS|<"+ii,df.loc[stationlist2]["BIAS"][df.loc[stationlist2]["BIAS"][(df["BIAS"]>-1.0*i) & (df["BIAS"]<1.0*i)]].count())
    #     print("|rAmplitude|<"+ii,df.loc[stationlist2]["rAmplitude"][df.loc[stationlist2]["rAmplitude"][(df["rAmplitude"]>-1.0*i) & (df["rAmplitude"]<1.0*i)]].count())
    # print("RMSE<5.0",df.loc[stationlist2]["RMSE"][df.loc[stationlist2]["RMSE"]<5.0].count())
    # print("|BIAS|<1.0",df.loc[stationlist2]["BIAS"][df.loc[stationlist2]["BIAS"][(df["BIAS"]>-1.0) & (df["BIAS"]<1.0)]].count())
    # print("|rAmplitude|<1.0",df.loc[stationlist2]["rAmplitude"][df.loc[stationlist2]["rAmplitude"][(df["rAmplitude"]>-1.0) & (df["rAmplitude"]<1.0)]].count())
    df.to_csv(path_or_buf="./tables/2.absolute_metrics_"+exp+".csv")
# for station in stationlist:
#     print station

# for station in stationlist1:
#     print station
# print stationlist, stationlist2