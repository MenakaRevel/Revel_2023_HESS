#!/opt/local/bin/python
# -*- coding: utf-8 -*-
'''
making hydrographs for give GRDC ID compare among experiment
'''
from cProfile import label
from operator import mod
import numpy as np
import matplotlib
matplotlib.use('Agg')
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
import string
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
from read_CMF import read_discharge, read_discharge_multi
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
def sharpness(l,u,o):
    """
    Sharpness
    input:
        l: lower bound
        u: upper bound
    output:
        Sharpness 
    """
    u,l = filter_nan(u,l)
    o=ma.masked_where(o<=0.0,o).filled(0.0)
    u=ma.masked_where(o<=0.0,u).filled(0.0)
    l=ma.masked_where(o<=0.0,l).filled(0.0)
    o=np.compress(o>0.0,o)
    u=np.compress(o>0.0,u)
    l=np.compress(o>0.0,l)
    return np.sum(u-l)
#====================================================================
def reliability(l,u,o):
    """
    Reliability
    input:
        l: lower bound
        u: upper bound
    output:
        Reliability 
    """
    u,l = filter_nan(u,l)
    o=ma.masked_where(o<=0.0,o).filled(0.0)
    u=ma.masked_where(o<=0.0,u).filled(0.0)
    l=ma.masked_where(o<=0.0,l).filled(0.0)
    o=np.compress(o>0.0,o)
    u=np.compress(o>0.0,u)
    l=np.compress(o>0.0,l)
    pnum=len(o)
    rcount=0
    for i in np.arange(pnum):
        if l[i] <= o[i] <= u[i]:
            rcount=rcount+1
    return float(rcount)/float(pnum)
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
argv=sys.argv
syear=int(argv[1])
eyear=int(argv[2])
CaMa_dir=argv[3]
mapname=argv[4]
figname=argv[5]
exlist=argv[6]
gname=argv[7]
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
obslist="/cluster/data6/menaka/HydroDA/dat/grdc_"+mapname+".txt"
with open(obslist,"r") as f:
    lines=f.readlines()
pname=[]
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
pname.append(num)
#====================================================================
def make_fig(station,expname):
    # read data
    ua, la, uo, lo = read_ul_all(expname,station)
    asm, opn = read_dis(expname,station)
    org=grdc.grdc_dis(station,syear,eyear)
    org=np.array(org)
    #=============
    # make figure
    #=============
    plt.close()
    #labels=["GRDC","corrupted","assimilated"]
    labels=["GRDC","simulated","assimilated"]
    fig, ax1 = plt.subplots()
    lines=[ax1.plot(np.arange(start,last),ma.masked_less(org,0.0),label="GRDC",color="#34495e",linewidth=3.0,zorder=101)[0]] #,marker = "o",markevery=swt[point])
    # draw mean of ensembles
    lines.append(ax1.plot(np.arange(start,last),opn,label="corrupted",color="#4dc7ec",linewidth=1.0,alpha=1,zorder=104)[0])
    lines.append(ax1.plot(np.arange(start,last),asm,label="assimilated",color="#ff8021",linewidth=1.0,alpha=1,zorder=106)[0])
    # ax1.fill_between(np.arange(start,last),lo,uo,color="#4dc7ec",alpha=0.2,zorder=102)
    # ax1.fill_between(np.arange(start,last),la,ua,color="#ff8021",alpha=0.2,zorder=103)
    # print ua
    # print asm
    # print la
    #    plt.ylim(ymin=)
    # Make the y-axis label, ticks and tick labels match the line color.
    ax1.set_ylabel('discharge (m$^3$/s)', color='k')
    ax1.set_xlim(xmin=0,xmax=last+1)
    ax1.tick_params('y', colors='k')
    # scentific notaion
    ax1.ticklabel_format(style="sci",axis="y",scilimits=(0,0))
    ax1.yaxis.major.formatter._useMathText=True 
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
    ax1.set_xticks(xxlist)
    ax1.set_xticklabels(xxlab,fontsize=10)
    plt.legend(lines,labels,ncol=1,loc='upper right') #, bbox_to_anchor=(1.0, 1.0),transform=ax1.transAxes)
    # station_loc_list=pname[point].split("/")
    # station_name="-".join(station_loc_list) 
    # print ('save',river[point] , station_name)
    # plt.savefig(assim_out+"/figures/disgraph/"+river[point]+"-"+station_name+".png",dpi=500)
    plt.show()
    return 0
#============
stations=["3629001","3618000","3626000","3625000"]
# read GRDC IDs
# gname="./Figs4-GRDClist.txt"
with open(gname,"r") as f:
    lineGRDCs=f.readlines()
stations=[]
for lineGRDC in lineGRDCs:
    lineGRDC= re.split(":",lineGRDC)
    lineGRDC = list(filter(None, lineGRDC))
    stations.append(lineGRDC[0].split()[0])
print stations
#======================================================
obslist="/cluster/data6/menaka/HydroDA/dat/grdc_"+mapname+".txt"
with open(obslist,"r") as f:
    lines=f.readlines()
IX1=[];IY1=[];IX2=[];IY2=[];nums=[];names=[];basins=[]
for station in stations:
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
        # #-------------------------
        # if num not in stations:
        #     continue
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
        org=grdc.grdc_dis(num,syear,eyear)
        if np.sum((org!=-9999.0)*1.0) < 365*1:
            # print ("no obs: ",stream, np.sum((org!=-9999.0)*1.0), "< 365 days")
            continue
        #-------------------------
        if num == station:
            stream=re.split("-",stream)[0]
            basin=re.split(",",basin)[0]
            print num, station, stream
            IX1.append(ix1-1)
            IY1.append(iy1-1)
            IX2.append(ix2-1)
            IY2.append(iy2-1)
            nums.append(num)
            names.append(stream)
            basins.append(basin)
#===========================
# indirs=["/cluster/data6/menaka/ensemble_org/CaMa_out/AMZE2O003"\
#        ,"/cluster/data6/menaka/ensemble_org/CaMa_out/AMZCALecmwf001"\
#        ,"./txt/NOM_WSE_E2O_HWEB_003/outflw"]
# metric=[]
# lasm=[]
# lopn=[]
# for i,exp in enumerate(labels):
#     print exp
#     lasm_frag=[]
#     lopn_frag=[]
#     for j in np.arange(3):
#         num=nums[3*i+j]
#         print (num)
#         experiment=indirs[i]
#         asm, opn = read_dis(experiment,num)
#         lasm_frag.append(asm)
#         lopn_frag.append(opn)
#     lasm.append(lasm_frag)
#     lopn.append(lopn_frag)
# lasm=np.array(lasm)
# lopn=np.array(lopn)
# print (np.shape(lasm))
# print (np.shape(lopn))
#=============
# make figure
#=============
plt.close()
# lexp=len(labels)
lexp=len(stations)
if lexp==6:
    colors=['#1f77b4','#9467bd','#2ca02c','#d62728','xkcd:goldenrod','xkcd:crimson']
elif lexp==4:
    colors=['#1f77b4','#9467bd','#2ca02c','#d62728']
elif lexp==2:
    colors=["#4dc7ec","#ff8021"]
else:
    colors=["#4dc7ec","#ff8021"]
# figure in A4 size
va_margin= 0.0#1.38#inch 
ho_margin= 0.0#1.18#inch
hgt=(11.69 - 2*va_margin)*(1.0/3.0)
wdt=(8.27 - 2*ho_margin)*(2.0/2.0)
fig=plt.figure(figsize=(wdt,hgt))
# print lexp/2
cols=2
G = gridspec.GridSpec(nrows=2,ncols=3)
# G.update(wspace=0.0, hspace=0.0) # set the spacing between axes. 
for i,exp in enumerate(labels):
    for j in np.arange(3):
        num=nums[3*i+j]
        print i,j,exp,num
        experiment=indirs[i]
        station=num
        ua, la, uo, lo = read_ul_all(experiment,station)
        asm, opn = read_dis(experiment,station)
        org=grdc.grdc_dis(station,syear,eyear)
        org=np.array(org)
        # x=int(j%cols)
        # y=int(j/cols)
        # print x, y
        ax = fig.add_subplot(G[i,j])
        ax.plot(np.arange(start,last),ma.masked_less(org,0.0),label="GRDC",color="#34495e",linewidth=1.0,zorder=101) #,marker = "o",markevery=swt[point])
        mxx=[]
        linewidth=0.3 #1.0/float(len(labels))
        ax.plot(np.arange(start,last),asm,label=exp,color="#ff8021",linewidth=linewidth,alpha=1,zorder=104)
        ax.plot(np.arange(start,last),opn,label=exp,color="#4dc7ec",linewidth=linewidth,alpha=1,zorder=104)
        ax.fill_between(np.arange(start,last),la,ua,color="#ff8021",alpha=0.2,zorder=103)
        ax.fill_between(np.arange(start,last),lo,uo,color="#4dc7ec",alpha=0.2,zorder=102)
        #====================
        # r=correlation(metric[i,j,:],org)
        # NSE=NS(metric[i,j,:],org)
        # KGE1=KGE(metric[i,j,:],org)
        # mxx.append([r,NSE,KGE1])
        # Make the y-axis label, ticks and tick labels match the line color.
        if j==0:
            ax.set_ylabel('discharge (m$^3$/s)', color='k')
        ax.set_xlim(xmin=0,xmax=last+1)
        ax.tick_params('y', colors='k')
        # scentific notaion
        ax.ticklabel_format(style="sci",axis="y",scilimits=(0,0))
        ax.yaxis.major.formatter._useMathText=True
        ax.yaxis.offsetText.set_fontsize(6) 
        ax.yaxis.label.set_size(6)
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
        ax.set_xticklabels(xxlab,fontsize=6)
        # if i==1:
        #     ax.set_xlabel('year', color='k')
        #     ax.set_xticklabels(xxlab,fontsize=6)
        # else:
        #     ax.set_xticklabels([])
        ax.tick_params(axis='both', which='major', labelsize=6)
        ax.set_title(names[3*i+j]+", "+basins[3*i+j],fontsize=6)
        #====================
        # correlation
        rasm=correlation(asm,org)
        ropn=correlation(opn,org)
        deltar=rasm-ropn
        dr="$\Delta$$r$= %3.2f"%(deltar)
        ax.text(0.0,-0.2,dr,ha="left",va="center",transform=ax.transAxes,fontsize=6)
        #====================
        # NSEAI
        NSEasm=NS(asm,org)
        NSEopn=NS(opn,org)
        NSEAI=(NSEasm-NSEopn)/(1-NSEopn+1.0e-20)
        NAI="$NSEAI$= %3.2f"%(NSEAI)
        ax.text(0.35,-0.2,NAI,ha="left",va="center",transform=ax.transAxes,fontsize=6)
        #====================
        # rISS
        ISSasm=ISS(ua,la,org,0.05)
        ISSopn=ISS(uo,lo,org,0.05)
        rISS=(ISSasm-ISSopn)/(ISSopn)
        rS="$rISS$= %3.2f"%(rISS)
        ax.text(0.8,-0.2,rS,ha="left",va="center",transform=ax.transAxes,fontsize=6)
        #========
        # if j==0:
        plt.gca().text(-0.15,1.08,"%s)"%(string.ascii_lowercase[3*i+j]),ha="left",va="center",transform=plt.gca().transAxes,fontsize=10)
    #========
    # NSE="NSE: %3.2f"%(NS(s,o))
    #============
    # data frame
    #============
    # print ("======================")
    # print (names[j])
    # mxx=np.array(mxx)
    # df=pd.DataFrame(mxx,columns=["r","NSE","KGE"],index=labels)
    # print (df)
# legend
labels0=["Observations", "Open-loop", "Assimilated"]
# labels0.extend(labels)
colors0=["#34495e","#4dc7ec","#ff8021"]
# colors0.extend(colors)
# print labels0
features=[]
for i in np.arange(len(labels0)):
    label=labels0[i]
    features.append(mlines.Line2D([], [], color=colors0[i], marker=None,
                          label=labels0[i],linewidth=2.5)) #edgecolors="k",
legend=plt.legend(handles=features,bbox_to_anchor=(0.5,-0.10), loc="lower center",
           bbox_transform=fig.transFigure, ncol=len(labels0),  borderaxespad=0.0, frameon=False)#
#======================
# fig.tight_layout()
# plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.1)#wspace=0, hspace=0)
plt.subplots_adjust(wspace=0.2, hspace=0.5)
#plt.title(stitle)
print figname

plt.savefig("./figures/"+figname+".jpg",dpi=800,bbox_inches="tight", pad_inches=0.05)
plt.savefig("./figures/"+figname+".png",dpi=800,bbox_inches="tight", pad_inches=0.05)
plt.savefig("./figures/"+figname+".pdf",dpi=800,bbox_inches="tight", pad_inches=0.05)