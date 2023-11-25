#!/opt/local/bin/python
# -*- coding: utf-8 -*-
'''
making WSE graphs for given VSs compare among experiment
'''
from cProfile import label
from operator import mod
import numpy as np
import matplotlib.pyplot as plt
import datetime
from matplotlib.cbook import boxplot_stats
import statistics
# from matplotlib.colors import LogNorm,Normalize,ListedColormap
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# from mpl_toolkits.basemap import Basemap
# import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import sys
import os
import calendar
from multiprocessing import Pool
# from multiprocessing import Process
from multiprocessing import sharedctypes
from numpy import ma
import re
import math
# import my_colorbar as mbar
# import cartopy.crs as ccrs
# import cartopy
# from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
# import cartopy.feature as cfeature
import os
import seaborn as sns
import pandas as pd


# import params as pm
from statistics import *
import read_grdc as grdc
import read_hydroweb as hweb
from read_CMF import read_discharge, read_discharge_multi
from read_sfcelv import read_sfcelv, read_sfcelv_multi

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
    s,o = filter_nan(s,o)
    o=ma.masked_where(o==-9999.0,o).filled(0.0)
    s=ma.masked_where(o==-9999.0,s).filled(0.0)
    o=np.compress(o>0.0,o)
    s=np.compress(o>0.0,s)
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
        dt = 30
        maxloc1=max(maxloc-dt,0)
        maxloc2=min(maxloc+dt,len(s)-1)
        minloc1=max(minloc-dt,0)
        minloc2=min(minloc+dt,len(s)-1)
        maxarry=ma.masked_equal(o[st_dt:ed_dt],-9999.0).filled(-9999.0)
        minarry=ma.masked_equal(o[st_dt:ed_dt],-9999.0).filled(9999.0)
        omax   = np.amax(maxarry[maxloc1:maxloc2])
        omin   = np.amin(minarry[minloc1:minloc2])

        # maxarry=ma.masked_equal(o[st_dt:ed_dt],-9999.0).filled(-9999.0)
        # minarry=ma.masked_equal(o[st_dt:ed_dt],-9999.0).filled(9999.0)
        # omax   = np.amax(maxarry)
        # omin   = np.amin(minarry)
        # maxloc = np.argmax(maxarry)
        # minloc = np.argmin(minarry)
        # maxloc1=max(maxloc-15,0)
        # maxloc2=min(maxloc+15,len(s)-1)
        # minloc1=max(minloc-15,0)
        # minloc2=min(minloc+15,len(s)-1)
        # smax   = np.amax(s[st_dt:ed_dt][maxloc1:maxloc2])
        # smin   = np.amin(s[st_dt:ed_dt][minloc1:minloc2])
        # print year, omax, omin, smax, smin
        if omax == -9999.0 or omin == 9999.0:
            continue
        alti_diff.append((smax-smin)-(omax-omin))
        # print year, st_dt, ed_dt, smax, smin, omax, omin
    return np.mean(alti_diff)
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
def read_dis_multi(ix1, iy1, ix2, iy2, syear, eyear, runoff_folder, nx, ny):
    #print ix1,iy1
    dis = np.zeros( (len(ix1), nbdays), 'f')
    # dis_max = np.zeros( (len(ix1), nbyears), 'f')
    for year in range(syear, eyear+1):
        # print year
        s_days = int( (datetime.date(year , 1,1) - datetime.date(syear, 1, 1)). days)
        e_days = int( (datetime.date(year+1, 1, 1) - datetime.date(syear, 1, 1)). days)
        
        f = runoff_folder + '/outflw'+str(year)+'.bin'

        tmp = read_discharge_multi( ix1, iy1, ix2, iy2, e_days-s_days, f, nx, ny)

        #print year, e_days - s_days, s_days, e_days, outflw.shape
        dis[:,s_days:e_days] = tmp
        # dis_max[:,year-syear] = np.nanmax(tmp, axis=1)

    return dis #, dis_max
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
CaMa_dir=argv[3]
mapname=argv[4]
figname=argv[5]
exlist=argv[6]
hname=argv[7]
ncpus=int(sys.argv[8])
ens_mem=21
seaborn_map=True
#====================================================================
# station=argv[1]
# expname=argv[2]
# syear=2009
# eyear=2014
# CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
# mapname="amz_06min"
# expname="NOM_WSE_E2O_HWEB_006"
# expname="DIR_WSE_E2O_HWEB_001"
# expname="DIR_WSE_E2O_HWEB_002"
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
#=============================
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
satcov="./dat/satellite_coverage.bin"
satcov=np.fromfile(satcov,np.float32).reshape(ny,nx)
#======================================================
start_dt=datetime.date(syear,1,1)
end_dt=datetime.date(eyear,12,31)
size=60

start=0
last=(end_dt-start_dt).days + 1
N=int(last)
nbdays=int(last)
#====================================================================
# stations=["3629001","3618000","3626000","3625000"]
# read HydroWeb IDs
# hname="./Figs5-HydroWeblist.txt"
with open(hname,"r") as f:
    lineHWEBs=f.readlines()
stations=[]
for lineHWEB in lineHWEBs:
    lineHWEB= re.split(":",lineHWEB)
    lineHWEB = list(filter(None, lineHWEB))
    stations.append(lineHWEB[0].strip())
# print stations 
#================================
obslist="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+mapname+"_QC0_simulation.txt"
# obslist="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+mapname+"_QC0.txt"
with open(obslist,"r") as f:
    lines=f.readlines()
IX=[];IY=[];nums=[];names=[];EGM08=[];EGM96=[]
for station in stations:
    for line in lines[1::]:
        line    = re.split(" ",line)
        line    = list(filter(None, line))
        # print line
        num     = line[0].strip()
        station0= line[1].strip()
        lon     = float(line[2])
        lat     = float(line[3])
        ix      = int(line[4])
        iy      = int(line[5])
        ele     = float(line[6])
        ele_dif = float(line[7])
        egm08   = float(line[8])
        egm96   = float(line[9])
        sat     = line[10].strip()
        if station0 == station:
            stream=re.split("_",station)[2]+"-"+re.split("_",station)[-1]
            # print stream, ix, iy
            IX.append(ix-1)
            IY.append(iy-1)
            nums.append(num)
            names.append(stream)
            EGM08.append(egm08)
            EGM96.append(egm96)
metric=[]
for i,exp in enumerate(labels):
    exppre=exp.split(" ")[0]
    if exppre=="CaMa":
    # if i < 2:
        wse_cmf = read_wse_multi(IX, IY, syear, eyear, indirs[i], nx, ny) 
        # print np.shape(dis_cmf.T)
        metric_frag=[]
        for j,num in enumerate(nums):
            metric_frag.append(wse_cmf.T[:,j])
    else:
        metric_frag=[]
        # print len(nums)
        for station in stations:
            experiment=indirs[i]
            try:
                asm, opn = read_wse(experiment,station,"simulation")
            except:
                asm, opn = read_wse(experiment,station,"validation")
            metric_frag.append(asm)
    metric.append(metric_frag)
metric=np.array(metric)
print (np.shape(metric))
#=============
# make figure
#=============
plt.close()
colors=['#1f77b4','#9467bd','#2ca02c','#d62728','xkcd:goldenrod','xkcd:crimson']
lexp=len(labels)
print lexp/2
# figure in A4 size
va_margin= 0.0#1.38#inch 
ho_margin= 0.0#1.18#inch
hgt=(11.69 - 2*va_margin)*(1.0/4.0)*(lexp/2)
wdt=(8.27 - 2*ho_margin)*(2.0/2.0)
fig=plt.figure(figsize=(wdt,hgt))
G = gridspec.GridSpec(nrows=lexp/2,ncols=2)
# G.update(wspace=0.05, hspace=0.05) # set the spacing between axes. 
for j,num in enumerate(nums):
    x=int(j%2)
    y=int(j/2)
    ax = fig.add_subplot(G[y,x])
    time,org=hweb.HydroWeb_WSE(stations[j],syear,eyear)
    data=np.array(org)+np.array(EGM08[j])-np.array(EGM96[j])
    ax.plot(time,data,label="obs",marker="o",color="#34495e",markersize=5,linewidth=0.0,fillstyle="none",zorder=101)
    org0=hweb.HydroWeb_continous_WSE(stations[j],syear,1,1,eyear,12,31,EGM08[j],EGM96[j])
    org0=np.array(org0)#+np.array(EGM08)-np.array(EGM96)
    mxx=[]
    for i,exp in enumerate(labels):
        ax.plot(np.arange(start,last),metric[i,j,:],label=exp,color=colors[i],linewidth=0.3,alpha=1,zorder=104)
        RMSE1=RMSE(metric[i,j,:],org0)
        BIAS1=BIAS(metric[i,j,:],org0)
        # DeltA=rAmplitude(metric[i,j,:],org0)
        DeltA=dAmplitude(metric[i,j,:],org0)
        mxx.append([RMSE1,BIAS1,DeltA])
    # Make the y-axis label, ticks and tick labels match the line color.
    if x==0:
        ax.set_ylabel('WSE $(m)$', color='k')
    ax.set_xlim(xmin=0,xmax=last+1)
    ax.tick_params('y', colors='k')
    # # scentific notaion
    # ax.ticklabel_format(style="sci",axis="y",scilimits=(0,0))
    # ax.yaxis.major.formatter._useMathText=True 
    #
    #xxlist=np.linspace(0,N,(eyear-syear)+1)
    #xlab=np.arange(syear,eyear+1,1)
    #xxlab=[calendar.month_name[i][:3] for i in range(1,13)]
    if eyear-syear > 1:
        dtt=1
        dt=int(math.ceil(((eyear-syear)+1)/dtt))
    else:
        dtt=1
        dt=(eyear-syear)+1
    xxlist=np.linspace(0,N,dt,endpoint=True)
    #xxlab=[calendar.month_name[i][:3] for i in range(1,13)]
    xxlab=np.arange(syear,eyear+1,dtt)
    ax.set_xticks(xxlist)
    
    if y==(lexp/2) - 1:
        ax.set_xlabel('year', color='k')
        ax.set_xticklabels(xxlab,fontsize=10)
    else:
        ax.set_xticklabels([])
    # plt.legend(lines,labels,ncol=1,loc='upper right') #, bbox_to_anchor=(1.0, 1.0),transform=ax.transAxes)
    # station_loc_list=pname[point].split("/")
    # station_name="-".join(station_loc_list) 
    # print ('save',river[point] , station_name)
    # plt.savefig(assim_out+"/figures/disgraph/"+river[point]+"-"+station_name+".png",dpi=500)
    ax.set_title(names[j])
    #============
    # data frame
    #============
    print ("======================")
    print (names[j])
    mxx=np.array(mxx)
    df=pd.DataFrame(mxx,columns=["RMSE","BIAS","dAmplitude"],index=labels)
    print (df)
# legend
# labels0=["Observations"]
# labels0.extend(labels)
# colors0=["#34495e"]
# colors0.extend(colors)
# print labels0
features=[mlines.Line2D([],[],label="observations",marker="o",color="#34495e",fillstyle="none",linewidth=0.0,zorder=101)]
for i in np.arange(len(labels)):
    features.append(mlines.Line2D([], [], color=colors[i], marker=None,label=labels[i],linewidth=2.5)) #edgecolors="k",
legend=plt.legend(handles=features,bbox_to_anchor=(0.5,0.05), loc="upper center",
           bbox_transform=fig.transFigure, ncol=len(labels)+1,  borderaxespad=0.0, frameon=False)#
#======================
# fig.tight_layout()
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)#wspace=0, hspace=0)
# plt.subplots_adjust(wspace=0, hspace=0)
#plt.title(stitle)
plt.savefig("./figures/"+figname+".jpg",dpi=800,bbox_inches="tight", pad_inches=0.0)
plt.savefig("./figures/"+figname+".png",dpi=800,bbox_inches="tight", pad_inches=0.0)
plt.savefig("./figures/"+figname+".pdf",dpi=800,bbox_inches="tight", pad_inches=0.0)