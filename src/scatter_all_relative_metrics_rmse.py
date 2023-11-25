#!/opt/local/bin/python
# -*- coding: utf-8 -*-

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

from statistics import *
import read_grdc as grdc
import my_colorbar as mbar
import read_hydroweb as hweb
from river_patch import patchms
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
    with open(txt,"r") as f:
        lines = f.readlines()

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
def func(x, a, b, c):
    return np.float64(a)*np.exp(-b*x) + np.float64(c)
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
#====================================================================
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
    return (np.mean(s)-np.mean(o))/np.mean(o)
#====================================================================
# def dAmplitude(s,o,dt=30,syear=2009,eyear=2014):
#     """
#     Amplitude Difference
#     input:
#         s: simulated
#         o: observed
#     output:
#         dAmplitude: Amplitude Difference
#     """ 
#     s,o = filter_nan(s,o)
#     # o=ma.masked_where(o==-9999.0,o).filled(0.0)
#     # s=ma.masked_where(o==-9999.0,s).filled(0.0)
#     # o=np.compress(o>0.0,o)
#     # s=np.compress(o>0.0,s)
#     alti_diff=[]
#     sim_amp=[]
#     obs_amp=[]
#     for year in np.arange(syear,eyear+1):
#         date1 = datetime.date(year, 1, 1) # month and day are 1-base
#         date2 = datetime.date(year, 12, 31)
#         days  = (date2-date1).days + 1
#         if year == syear:
#             st_dt = 0
#             ed_dt = days
#         elif year == eyear:
#             st_dt = ed_dt
#             ed_dt = -1
#         else:
#             st_dt = ed_dt
#             ed_dt = ed_dt + days
#         # print st_dt, ed_dt
#         # max min simuulated
#         smax   = np.amax(s[st_dt:ed_dt])
#         smin   = np.amin(s[st_dt:ed_dt])
#         sim_amp.append(smax-smin)
#         # max min observed
#         maxarry=ma.masked_equal(o[st_dt:ed_dt],-9999.0).filled(-9999.0)
#         minarry=ma.masked_equal(o[st_dt:ed_dt],-9999.0).filled(9999.0)
#         #===
#         omax   = np.amax(maxarry)
#         omin   = np.amin(minarry)
#         if omax == -9999.0 or omin == 9999.0:
#             continue
#         maxloc = np.argmax(maxarry)
#         minloc = np.argmin(minarry)
#         obs_amp.append(omax-omin)

#         # maxloc = np.argmax(s[st_dt:ed_dt])
#         # minloc = np.argmin(s[st_dt:ed_dt])
#         # smax   = np.amax(s[st_dt:ed_dt])
#         # smin   = np.amin(s[st_dt:ed_dt])
#         # maxloc1=max(maxloc-dt,0)
#         # maxloc2=min(maxloc+dt,len(s)-1)
#         # minloc1=max(minloc-dt,0)
#         # minloc2=min(minloc+dt,len(s)-1)
#         # maxarry=ma.masked_equal(o[st_dt:ed_dt],-9999.0).filled(-9999.0)
#         # minarry=ma.masked_equal(o[st_dt:ed_dt],-9999.0).filled(9999.0)
#         # omax   = np.amax(maxarry[maxloc1:maxloc2])
#         # omin   = np.amin(minarry[minloc1:minloc2])
#         # print omin, omax
#         # if omax == -9999.0 or omin == 9999.0:
#         #     continue
#         # alti_diff.append((smax-smin)-(omax-omin))
#     return np.mean(np.array(sim_amp))-np.mean(np.array(obs_amp)) #np.mean(alti_diff)
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
def plot_scatter(stat, lstato, color, numch, ylabel, xlabel,ax1=None,left=False,bottom=False):
    ax1=ax1 or plt.gca()
    # ax2=ax2 or plt.gca()
    ax1.plot(lstato, stat, color=color, linewidth=0.0,marker="o",markersize=5,fillstyle="none")
    # ax2.plot(lstata, stat, color=color, linewidth=0.0,marker="^",markersize=5,fillstyle="none")
    #=====
    linearR(lstato,stat,color,ax=ax1)
    # linearR(lstata,stat,color,ax=ax2)
    # print np.corrcoef(lstato, stat)[0,1]
    # print np.corrcoef(lstata, stat)[0,1]
    # ax1.set_xlim(xmin=-1.0,xmax=1.0)
    # ax2.set_xlim(xmin=-1.0,xmax=1.0)
    if left:
        ax1.set_ylabel(ylabel,fontsize=8)
    if bottom:
        ax1.set_xlabel(xlabel+" of WSE (Open-loop)",fontsize=8)
    # ax2.set_xlabel("$"+xlabel+"$ of WSE (Assimilation)",fontsize=8)
    ax1.tick_params(axis='both', which='major', labelsize=6)
    # ax2.tick_params(axis='both', which='major', labelsize=6)
    ax1.text(-0.10,1.05,"%s"%(numch),ha="left",va="center",transform=ax1.transAxes,fontsize=10)
    return 0
#====================================================================
def linearR(xx0,yy0,color,ax=None):
    xx0=np.array(xx0)
    yy0=np.array(yy0)
    ax=ax or plt.gca() 
    xx  = np.sort(xx0)[::-1]
    yy  = yy0[np.argsort(xx0)][::-1]
    # linear regression
    model = LinearRegression().fit(xx.reshape((-1, 1)), yy.reshape((-1, 1)))
    y_pred = model.predict(xx.reshape((-1, 1)))
    ax.plot(xx, y_pred, color=color, linewidth=0.5, linestyle="--")
    return 0
#====================================================================
def patchMS(ix,iy,ngrids,nextxy,uparea,nxtdst,nx,ny):
    outdir='/cluster/data6/menaka/Empirical_LocalPatch'
    nextX=nextxy[0]
    nextY=nextxy[1]
    # print "**********",ix,iy,nx,ny #outdir,
    xlist,ylist,dlist=patchms(ix+1,iy+1,ngrids,nextX.T,nextY.T,uparea.T,nxtdst.T) #,outdir)
    xlist=xlist[xlist>0]
    ylist=ylist[ylist>0]
    dlist=dlist[xlist>0]
    # print xlist-1, ylist-1, dlist
    return xlist-1, ylist-1, dlist
#====================================================================
def near_VS(ix0,iy0,ngrids,nextxy,uparea,nxtdst,nx,ny):
    # print ix0,iy0
    xlist,ylist,dlist=patchMS(ix0,iy0,ngrids,nextxy,uparea,nxtdst,nx,ny)
    # print xlist
    # print ylist
    # print dlist
    flag=0
    dist0=1.0e20
    station0="none"
    eg08=0.0
    eg96=0.0
    for iXX,iYY,iDD in zip(xlist,ylist,dlist):
        # print "***",iXX, iYY, iDD
        obslist="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+mapname+"_QC0_simulation.txt"
        with open(obslist,"r") as f:
            lines=f.readlines()
        # metric_frag=[]
        for line in lines[1::]:
            line    = re.split(" ",line)
            line    = list(filter(None, line))
            # print line
            num     = line[0].strip()
            station = line[1].strip()
            lon     = float(line[2])
            lat     = float(line[3])
            ix      = int(line[4]) - 1
            iy      = int(line[5]) - 1
            ele     = float(line[6])
            ele_dif = float(line[7])
            EGM08   = float(line[8])
            EGM96   = float(line[9])
            sat     = line[10].strip()
            #-------------------------
            if ix == iXX and iy == iYY:
                flag=1
                station0=station
                eg08=EGM08
                eg96=EGM96
                dist0=iDD
                # print "L390",iXX, iYY, station, iDD
                break
        #=============================
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
            ix      = int(line[4]) - 1
            iy      = int(line[5]) - 1
            ele     = float(line[6])
            ele_dif = float(line[7])
            EGM08   = float(line[8])
            EGM96   = float(line[9])
            sat     = line[10].strip()
            #-------------------------
            if ix == iXX and iy == iYY:
                if flag==0 and dist0>iDD:
                    flag=1
                    station0=station
                    eg08=EGM08
                    eg96=EGM96
                    dist0=iDD
                    # print "L419", iXX, iYY, station, iDD
                    break
                elif flag==0:
                    flag=1
                    station0=station
                    eg08=EGM08
                    eg96=EGM96
                    dist0=iDD
                    # print "L426",iXX, iYY, station, iDD
                    break
        if flag==1:
            break
    return station0, flag, eg08, eg96, dist0
#====================================================================
def get_local_patch(ix,iy):
    xlist=[]
    ylist=[]
    fname="/cluster/data6/menaka/Empirical_LocalPatch/local_patch/amz_06min_S14FD_60/patch%04d%04d.txt"%(ix+1,iy+1)
    with open(fname,"r") as f:
        lines=f.readlines()
    for line in lines:
        line    = re.split(" ",line)
        line    = list(filter(None, line))
        xx      = int(line[0])
        yy      = int(line[1])
        xlist.append(xx)
        ylist.append(yy)
    return xlist, ylist
#====================================================================
def local_patch_VS(ix,iy):
    xlist, ylist=get_local_patch(ix,iy)
    lstation=[]
    legm08=[]
    legm96=[]
    for iXX,iYY in zip(xlist,ylist):
        # print "***",iXX, iYY, iDD
        obslist="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+mapname+"_QC0_simulation.txt"
        with open(obslist,"r") as f:
            lines=f.readlines()
        # metric_frag=[]
        for line in lines[1::]:
            line    = re.split(" ",line)
            line    = list(filter(None, line))
            # print line
            num     = line[0].strip()
            station = line[1].strip()
            lon     = float(line[2])
            lat     = float(line[3])
            ix      = int(line[4]) - 1
            iy      = int(line[5]) - 1
            ele     = float(line[6])
            ele_dif = float(line[7])
            EGM08   = float(line[8])
            EGM96   = float(line[9])
            sat     = line[10].strip()
            #-------------------------
            if ix == iXX and iy == iYY:
                lstation.append(station)
                legm08.append(EGM08)
                legm96.append(EGM96)
        #=============================
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
            ix      = int(line[4]) - 1
            iy      = int(line[5]) - 1
            ele     = float(line[6])
            ele_dif = float(line[7])
            EGM08   = float(line[8])
            EGM96   = float(line[9])
            sat     = line[10].strip()
            #-------------------------
            if ix == iXX and iy == iYY:
                lstation.append(station)
                legm08.append(EGM08)
                legm96.append(EGM96)
    return lstation, legm08, legm96
#====================================================================
def calc_rmse(exp,lstation, legm08, legm96):
    lrmseo=[]
    for i in np.arange(len(lstation)):
        station=lstation[i]
        egm08=legm08[i]
        egm96=legm96[i]
        #=============
        try:
            asm, opn = read_wse(exp,station,"simulation")
        except:
            asm, opn = read_wse(exp,station,"validation")
        #--------------
        org=hweb.HydroWeb_continous_WSE(station,syear,1,1,eyear,12,31,egm08,egm96)
        org=np.array(org)#+np.array(EGM08)-np.array(EGM96)
        rmsea=RMSE(asm,org)
        rmseo=RMSE(opn,org)
        lrmseo.append(rmseo)
    rmseo=np.mean(np.array(lrmseo))
    return rmseo
#====================================================================
#===================
# mk_dir("./figures")
# mk_dir("./figures/corr")
#===================
# syear=2009
# eyear=2014
# # CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
# CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"
# #map="glb_15min"
# # map="glb_06min"
# mapname="amz_06min"
# # expname="NOM_WSE_E2O_HWEB_004"
# lexp=["DIR_WSE_E2O_HWEB_001","ANO_WSE_E2O_HWEB_004","NOM_WSE_E2O_HWEB_004"]
# figname="fig11-all_relative_metrics"
#===================
# argv=sys.argv
# syear=int(argv[1])
# eyear=int(argv[2])
# numexp=int(argv[3])
# lexp=argv[4:4+numexp]
# CaMa_dir=argv[-3]
# mapname=argv[-2]
# figname=argv[-1]

argv=sys.argv
syear=int(argv[1])
eyear=int(argv[2])
CaMa_dir=argv[3]
mapname=argv[4]
exlist=argv[5]
dist_thr=argv[6]
figname=argv[7]
# gname=argv[7]
ncpus=int(sys.argv[8])
ens_mem=21
seaborn_map=True
numvar=2
#====================================================================
#=== read experiment names ===
# exlist="./experiment_list.nam"
# exlist="./Fig10-experiment_list.nam"
with open(exlist,"r") as exf:
    linesexf=exf.readlines()
indirs=[]
labels=[]
for lineexf in linesexf:
    lineexf = re.split(":",lineexf)
    lineexf = list(filter(None, lineexf))
    labels.append(lineexf[0])
    indirs.append(lineexf[1].strip())
#===
lexp=indirs
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
#============
stations=["3629001","3618000","3626000","3625000"]
# read GRDC IDs
# gname="./Figs4-GRDClist.txt"
# with open(gname,"r") as f:
#     lineGRDCs=f.readlines()
# stations=[]
# annotate_grdc={}
# text_label={}
# for exp in lexp:
#     annotate_grdc[exp]=[]
# for i,lineGRDC in enumerate(lineGRDCs):
#     lineGRDC= re.split(":",lineGRDC)
#     lineGRDC = list(filter(None, lineGRDC))
#     stations.append(lineGRDC[0].split()[0])
#     if i < 3:
#         annotate_grdc[lexp[0]].append(lineGRDC[0].split()[0])
#     else:
#         annotate_grdc[lexp[1]].append(lineGRDC[0].split()[0])
#     #---------
#     text_label[lineGRDC[0].split()[0]]=string.ascii_lowercase[i]
# print stations
#======================================================
#files for making river network
tmp0=re.split("-",figname)[0]+".txt"
# exp=expname
#
# rivernames  = ["AMAZON","NEGRO, RIO","PURUS, RIO","MADEIRA, RIO","JURUA, RIO"
#               ,"TAPAJOS, RIO","XINGU, RIO","CURUA, RIO","JAPURA, RIO","BRANCO, RIO"
#               ,"JAVARI, RIO","IRIRI, RIO","JURUENA, RIO","ACRE, RIO","BENI, RIO"
#               ,"MAMORE, RIO","GUAPORE, RIO","ARINOS, RIO","TROMBETAS, RIO"]
basins = ["AMAZON","NEGRO","PURUS","MADEIRA","JURUA","JAPURA","JAVARI"
         ,"ARIPUANA","IRIRI","BRANCO","UAUPES"] # "XINGU", ,"CURUA"
#
def ret_metric(exp):
    lNSEAI=[]
    ldelr=[]
    lriss=[]
    lbiaso=[]
    lbiasa=[]
    ldelAo=[]
    ldelAa=[]
    ldpeako=[]
    ldpeaka=[]
    lrmse=[]
    i=1
    # GRDC
    obslist="/cluster/data6/menaka/HydroDA/dat/grdc_"+mapname+".txt"
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
        # if basin not in basins:
        #     if np.sum((org!=-9999.0)*1.0) >= 365*1:
        #         nolist.append(basin)
        #     continue
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
        asm, opn = read_dis(exp,num)
        # NSEAI
        NSEasm=NS(asm,org)
        NSEopn=NS(opn,org)
        NAI=(NSEasm-NSEopn)/(1.0-NSEopn+1.0e-20)
        # delat r
        CORasm=correlation(asm,org)
        CORopn=correlation(opn,org)
        DCOR=CORasm-CORopn
        # rISS
        ua, la, uo, lo = read_ul_all(exp,num)
        ISSasm=ISS(ua,la,org,0.05)
        ISSopn=ISS(uo,lo,org,0.05)
        rISS  =(ISSasm-ISSopn)/(ISSopn)
        #-----------------
        print stream
        VSstation,flag,egm08,egm96,dist=near_VS(ix1-1,iy1-1,20,nextxy,uparea,nxtdst,nx,ny)
        if flag==0:
            # print stream, flag
            continue
        if dist > dist_thr:
            continue
        print stream, "-", VSstation, dist
        # print i,stream,ix1-1, iy1-1, ix2-1,iy2-1, VSstation,flag
        i=i+1
        try:
            asm1, opn1 = read_wse(exp,VSstation,"simulation")
        except:
            asm1, opn1 = read_wse(exp,VSstation,"validation")
        #--------------
        org1=hweb.HydroWeb_continous_WSE(VSstation,syear,1,1,eyear,12,31,egm08,egm96)
        org1=np.array(org1)#+np.array(EGM08)-np.array(EGM96)
        rmsea=RMSE(asm1,org1)
        rmseo=RMSE(opn1,org1)
        # biasa=BIAS(asm1,org1)
        # biaso=BIAS(opn1,org1)
        # delAa=dAmplitude(asm1,org1)
        # delAo=dAmplitude(opn1,org1)
        # dpeaka=dpeak_time(asm1,org1)
        # dpeako=dpeak_time(opn1,org1)
        # #--
        # lstation, legm08, legm96 = local_patch_VS(ix1-1,iy1-1)
        # rmseo = calc_rmse(exp,lstation, legm08, legm96)
        #====
        ldelr.append(DCOR)
        lNSEAI.append(NAI)
        lriss.append(rISS)
        lrmse.append(rmseo)
        if NAI<0.0:
            print basin, stream, NAI, rmseo
        # lbiaso.append(biaso)
        # lbiasa.append(biasa)
        # ldelAo.append(delAo)
        # ldelAa.append(delAa)
        # ldpeako.append(dpeako)
        # ldpeaka.append(dpeaka)
    return ldelr, lNSEAI, lriss, lrmse # lbiaso, ldelAo, ldpeako
#==================
#------------------
# figure in A4 size
colors = plt.cm.rainbow_r(np.linspace(0,1,len(basins)))
vmin=1
vmax=len(basins)
# norm=norm=Normalize(vmin=vmin,vmax=vmax)
cmap=plt.cm.get_cmap("rainbow_r")
norm=BoundaryNorm(np.arange(vmin,vmax+0.1,1),cmap.N)
# cmap   = matplotlib.colors.ListedColormap(matplotlib.cm.get_cmap("viridis_r").colors[:len(basins)])
va_margin= 0.0#1.38#inch 
ho_margin= 0.0#1.18#inch
hgt=(11.69 - 2*va_margin)*(8.27/11.69)*(1.0/3.0)
wdt=(8.27 - 2*ho_margin)*(8.27/8.27)
#fig=plt.figure(figsize=(8.27,11.69))
fig=plt.figure(figsize=(wdt,hgt))
#fig.suptitle("auto-correalated area")#"maximum autocorrelation length")
#G = gridspec.GridSpec(2,1)
G = gridspec.GridSpec(ncols=3,nrows=1)
#--scatter
markers=["o","^"]
colors=["#3174a1","#3f913a","#bf353c"]
# colors=["#afccdc","#b5d294","#f4adae"]
# lNSEAI1,lRMSE1,lNSEopn1=ret_metric(lexp[0])
# lNSEAI1,lRMSE1,lNSEopn1=ret_metric(lexp[1])
ax1 = fig.add_subplot(G[0,0])
ax2 = fig.add_subplot(G[0,1])
ax3 = fig.add_subplot(G[0,2])

# ax1 = fig.add_subplot(G[0,0])
# ax2 = fig.add_subplot(G[1,0])
# ax3 = fig.add_subplot(G[2,0])
# ax4 = fig.add_subplot(G[0,1])
# ax5 = fig.add_subplot(G[1,1])
# ax6 = fig.add_subplot(G[2,1])
# ax7 = fig.add_subplot(G[0,2])
# ax8 = fig.add_subplot(G[1,2])
# ax9 = fig.add_subplot(G[2,2])
#============================
texts=[]
for j,exp in enumerate(lexp):
    print "=================="
    print exp
    ldelr, lNSEAI, lriss, lrmse=ret_metric(exp)
    # RMSE
    plot_scatter(ldelr, lrmse,colors[j],"a) $\Delta$$r$","$\Delta$$r$","$RMSE$",ax1=ax1,left=True,bottom=True)
    plot_scatter(lNSEAI, lrmse, colors[j],"b) $NSEAI$","$NSEAI$","$RMSE$",ax1=ax2,left=True,bottom=True)
    plot_scatter(lriss, lrmse, colors[j],"c) $rISS$","$rISS$","$RMSE$",ax1=ax3,left=True,bottom=True)
    
    # # for BIAS
    # plot_scatter(ldelr, lbiaso,colors[j],"a","$\Delta$$r$","$BIAS$",ax1=ax1,left=True,bottom=False)
    # plot_scatter(lNSEAI, lbiaso, colors[j],"d","$NSEAI$","$BIAS$",ax1=ax2,left=True,bottom=False)
    # plot_scatter(lriss, lbiaso, colors[j],"g","$rISS$","$BIAS$",ax1=ax3,left=True,bottom=True)
    # # for delat A
    # plot_scatter(ldelr, ldelAo,colors[j],"b","$\Delta$$r$","$\Delta$$A$",ax1=ax4,left=False,bottom=False)
    # plot_scatter(lNSEAI, ldelAo, colors[j],"e","$NSEAI$","$\Delta$$A$",ax1=ax5,left=False,bottom=False)
    # plot_scatter(lriss, ldelAo, colors[j],"h","$rISS$","$\Delta$$A$",ax1=ax6,left=False,bottom=True)
    # # for BIAS
    # plot_scatter(ldelr, ldpeako,colors[j],"c","$\Delta$$r$","$\Delta$$t$",ax1=ax7,left=False,bottom=False)
    # plot_scatter(lNSEAI, ldpeako, colors[j],"f","$NSEAI$","$\Delta$$t$",ax1=ax8,left=False,bottom=False)
    # plot_scatter(lriss, ldpeako, colors[j],"i","$rISS$","$\Delta$$t$",ax1=ax9,left=False,bottom=True)
#==============================
patches=[]
for point in np.arange(len(labels)):
    patches.append(mpatches.Patch(color=colors[point],label=labels[point]))
#--
legend1=plt.legend(handles=patches,bbox_to_anchor=(0.5,0.01), loc="upper center",
           bbox_transform=fig.transFigure, ncol=3,  borderaxespad=0.0, frameon=False)#
#===
plt.savefig("./figures/"+figname+".pdf",dpi=800,bbox_inches="tight", pad_inches=0.05)
plt.savefig("./figures/"+figname+".png",dpi=800,bbox_inches="tight", pad_inches=0.05)
plt.savefig("./figures/"+figname+".jpg",dpi=800,bbox_inches="tight", pad_inches=0.05)