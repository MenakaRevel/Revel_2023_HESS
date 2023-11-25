#!/opt/local/bin/python
# -*- coding: utf-8 -*-
'''
delta r cumilative distribution
NSEAI boxplot 
realtive sharpness boxplot
delta relatibility boxplot
rISS boxplot
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
import string
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
import warnings;warnings.filterwarnings('ignore')

# import params as pm
import read_grdc as grdc
# import cal_stat as stat
#import plot_colors as pc
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
def print_stat(data,labels):
    # get stats
    df=pd.DataFrame(data=data, columns=labels)
    stats = boxplot_stats(df.values)
    print ("Experiment","Q1 ","Q3","Mean","Median")
    for i in np.arange(0,len(labels)):
        print ("%10s%15.5f%15.5f%15.5f%15.5f")%(labels[i],stats[i]['q1'],stats[i]['q3'],stats[i]['mean'],stats[i]['med'])
    return 0
#====================================================================
def plot_boxplot(data,num,labels,xlabel,colors=["#004488","#ddaa33","#ba5566"],ax=None,ylim=[0,1]):
    ax=ax or plt.gca() 
    flierprops = dict(marker='o', markerfacecolor='none', markersize=8,linestyle='none', markeredgecolor='k')
    boxprops = dict(color='grey')#facecolor='none'
    whiskerprops = dict(color='grey',linestyle="--")
    capprops = dict(color='grey')
    medianprops = dict(color='grey',linestyle="-",linewidth=1.0)
    meanprops = dict(marker='D', markeredgecolor='black',markerfacecolor='green',markersize=8)
    #
    df=pd.DataFrame(data=data, columns=labels)#, index=stationlist)
    box=sns.boxplot(ax=ax,data=df, fliersize=0.0, palette=colors, whis=1.5\
        ,meanline=True, width=0.8, linewidth=0.3, dodge=True\
        ,meanprops=meanprops,capprops=capprops,medianprops=medianprops) #"Paired"
    ax.axhline(0.0,color="k",linestyle="--",linewidth=0.5)
    ax.set_ylabel(xlabel, color='k',fontsize=8)
    ax.set_xlabel('Experiment', color='k',fontsize=8)
    ax.set_ylim(ymin=ylim[0],ymax=ylim[1])
    ax.tick_params(labelsize=6)
    ax.set_xticklabels(labels,rotation = 90)
    ax.text(-0.05,1.05,"%s)"%(string.ascii_lowercase[num-1]),ha="left",va="center",transform=ax.transAxes,fontsize=8)
    return 0
#====================================================================
def plot_displot(data,num,labels,xlabel,colors=["#004488","#ddaa33","#ba5566"],ax=None):
    ax=ax or plt.gca()
    print (np.shape(data))
    for i,exp in enumerate(labels):
        sns.distplot(data[:,i],ax=ax, hist = False, kde = True,
                kde_kws = {'linewidth': 1,'linestyle':'-'},
                label = exp, color=colors[i],norm_hist=False)
    ax.set_xlabel(xlabel, color='k',fontsize=8) #$correlation$ $coefficient$ $(r)$
    ax.set_ylabel('Number of gauges', color='k',fontsize=8)
    ax.set_xlim(xmin=-0.5,xmax=0.5)
    ax.tick_params(labelsize=6)
    ax.text(-0.05,1.05,"%s)"%(string.ascii_lowercase[num-1]),ha="left",va="center",transform=ax.transAxes,fontsize=8)
    # plt.legend(ncol=1,loc="lower right",bbox_to_anchor=(-0.15, 0.0),borderaxespad=0.0,frameon=False)
    plt.legend(ncol=1,loc="best",borderaxespad=0.0,frameon=False)
    return 0
#====================================================================
def plot_cdfplot(data,num,labels,xlabel,colors=["#004488","#ddaa33","#ba5566"],ax=None):
    ax=ax or plt.gca()
    kwargs = {'cumulative': True}
    print (np.shape(data))
    for i,exp in enumerate(labels):
        sns.distplot(data[:,i],ax=ax, hist = False, kde = True,
                hist_kws = kwargs,
                kde_kws = {'cumulative': True, 'linewidth': 1,'linestyle':'-'},
                label = exp, color=colors[i],norm_hist=False)
    ax.set_xlabel(xlabel, color='k',fontsize=8) #$correlation$ $coefficient$ $(r)$
    ax.set_ylabel('Cumilative number of gauges', color='k',fontsize=8)
    ax.set_xlim(xmin=-1.0,xmax=0.2)
    ax.tick_params(labelsize=6)
    ax.text(-0.05,1.05,"%s)"%(string.ascii_lowercase[num-1]),ha="left",va="center",transform=ax.transAxes,fontsize=8)
    # plt.legend(ncol=1,loc="lower right",bbox_to_anchor=(-0.1, 0.0),borderaxespad=0.0,frameon=False)
    return 0
#====================================================================
argv=sys.argv
syear=int(argv[1])
eyear=int(argv[2])
CaMa_dir=argv[3]
mapname=argv[4]
exlist=argv[5]
figname=argv[6]
ncpus=int(sys.argv[7])
ens_mem=49
seaborn_map=True
# seaborn_map=False
#==================================
# DA_dir="/cluster/data6/menaka/HydroDA"
# CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
# mapname="amz_06min"
# ens_mem=21
# lexp=7
# seaborn_map=False
#==================================
#==================================
# experiment="DIR_WSE_E2O_HWEB_001"
# experiment="ANO_WSE_E2O_HWEB_001"
# experiment="NOM_WSE_E2O_HWEB_001"
# experiments=["DIR_WSE_E2O_HWEB_001","DIR_WSE_E2O_HWEB_002","ANO_WSE_E2O_HWEB_001"\
            # ,"ANO_WSE_E2O_HWEB_002","NOM_WSE_E2O_HWEB_001","NOM_WSE_E2O_HWEB_002"\
            # ,"NOM_WSE_E2O_HWEB_003","NOM_WSE_E2O_HWEB_004","NOM_WSE_E2O_HWEB_005"]
# experiments=["NOM_WSE_E2O_HWEB_001","NOM_WSE_E2O_HWEB_002"\
#             ,"NOM_WSE_E2O_HWEB_003","NOM_WSE_E2O_HWEB_004"]
# experiments=["DIR_WSE_E2O_HWEB_001","DIR_WSE_E2O_HWEB_002"\
#             ,"ANO_WSE_E2O_HWEB_001","ANO_WSE_E2O_HWEB_002","ANO_WSE_E2O_HWEB_003"\
#             ,"NOM_WSE_E2O_HWEB_003","NOM_WSE_E2O_HWEB_002","NOM_WSE_E2O_HWEB_003"]
# experiments=["DIR_WSE_E2O_HWEB_001","DIR_WSE_E2O_HWEB_002"\
#             ,"ANO_WSE_E2O_HWEB_001","ANO_WSE_E2O_HWEB_003"\
#             ,"NOM_WSE_E2O_HWEB_001","NOM_WSE_E2O_HWEB_003"]
# experiments=["DIR_WSE_E2O_HWEB_001","DIR_WSE_E2O_HWEB_002"\
#             ,"ANO_WSE_E2O_HWEB_001","ANO_WSE_E2O_HWEB_003"\
#             ,"NOM_WSE_E2O_HWEB_001","NOM_WSE_E2O_HWEB_003"\
#             ,"NOM_WSE_E2O_HWEB_004","NOM_WSE_E2O_HWEB_005"]
#=== read experiment names ===
# exlist="./experiment_list.nam"
# exlist="./Fig10-experiment_list.nam"
with open(exlist,"r") as exf:
    linesexf=exf.readlines()
experiments=[]
labels=[]
for lineexf in linesexf:
    lineexf = re.split(":",lineexf)
    lineexf = list(filter(None, lineexf))
    labels.append(lineexf[0])
    experiments.append(lineexf[1].strip())
#=============================
lexp=len(experiments)
# colors=['xkcd:aqua green','xkcd:pastel blue','xkcd:soft pink','xkcd:pastel blue','xkcd:aqua green','xkcd:soft pink']

#assim_out=pm.DA_dir()+"/out/"+pm.experiment()+"/assim_out"
#assim_out=pm.DA_dir()+"/out/"+experiment+"/assim_out"
# assim_out=DA_dir+"/out/"+experiment
# print (assim_out)
# mk_dir("./figures")
# mk_dir("/figures/NSEAI")
#===================
fname=CaMa_dir+"/map/"+mapname+"/params.txt"
with open(fname,"r") as f:
    lines=f.readlines()
#-------
nx     = int(list(filter(None, re.split(" ",lines[0])))[0])
ny     = int(list(filter(None, re.split(" ",lines[1])))[0])
gsize  = float(list(filter(None, re.split(" ",lines[3])))[0])
lon0   = float(list(filter(None, re.split(" ",lines[4])))[0])
lat0   = float(list(filter(None, re.split(" ",lines[7])))[0])
west   = float(list(filter(None, re.split(" ",lines[4])))[0])
east   = float(list(filter(None, re.split(" ",lines[5])))[0])
south  = float(list(filter(None, re.split(" ",lines[6])))[0])
north  = float(list(filter(None, re.split(" ",lines[7])))[0])
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
metric=[]
lKGEa=[]
lKGEo=[]
for i,exp in enumerate(experiments):
    print (exp)
    metric_frag=[]
    lKGEa_frag=[]
    lKGEo_frag=[]
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
        # read discharge
        asm, opn = read_dis(exp,num)
        # correlation
        CORasm=correlation(asm,org)
        CORopn=correlation(opn,org)
        DCOR=CORasm-CORopn
        # KGE
        KGEasm=KGE(asm,org)
        KGEopn=KGE(opn,org)
        # NSEAI
        NSEasm=NS(asm,org)
        NSEopn=NS(opn,org)
        NSEAI=(NSEasm-NSEopn)/(1.0-NSEopn+1.0e-20)
        # rISS
        ua, la, uo, lo = read_ul_all(exp,num)
        ISSasm=ISS(ua,la,org,0.05)
        ISSopn=ISS(uo,lo,org,0.05)
        rISS=(ISSasm-ISSopn)/(ISSopn)#+1.0e-20)
        # sharpness
        shrasm=sharpness(la,ua,org)
        shropn=sharpness(lo,uo,org)
        rshr=(shrasm-shropn)/shropn
        # reliability
        relasm=reliability(la,ua,org)
        relopn=reliability(lo,uo,org)
        rrel=relasm-relopn #/relopn
        #=====================
        metric_frag.append([DCOR,NSEAI,rshr,rrel,rISS])
    metric.append(metric_frag)
metric=np.array(metric)
#====================================
#---------- making fig --------------
#====================================
# data=np.transpose(metric)
# data=np.nan_to_num(data)
# print (np.shape(metric))
#====================================
# figure in A4 size
# colors=['xkcd:aqua green','xkcd:pastel blue','xkcd:soft pink','xkcd:pastel blue','xkcd:aqua green','xkcd:soft pink']
# colors=["#afccdc","#3174a1","#b5d294","#3f913a","#f4adae","#bf353c"]
# colors=["#3174a1","#3f913a","#bf353c"]
# colors=["#004488","#ddaa33","#ba5566"]
colors=["#D81B60","#FFC107","#004D40"]
va_margin= 0.0#1.38#inch 
ho_margin= 0.0#1.18#inch
hgt=(11.69 - 2*va_margin)*(1.0/3.0)
wdt=(8.27 - 2*ho_margin)*(2.0/2.0)
#fig=plt.figure(figsize=(8.27,11.69))
fig=plt.figure(figsize=(wdt,hgt))
#fig.suptitle("auto-correalated area")#"maximum autocorrelation length")
G = gridspec.GridSpec(ncols=3, nrows=2)
#--boxplot
stationlist=get_GRDClist("/cluster/data6/menaka/HydroDA/dat/grdc_"+mapname+".txt")
print ("making plot")
# delta r probabilty distribution plot
ax1 = fig.add_subplot(G[0,0])
data=metric[:,:,0]
data=np.transpose(data)
data=np.nan_to_num(data)
plot_displot(data,1,labels,'$\Delta$$r$',colors=colors,ax=ax1)
# NSEAI boxplot
ax2 = fig.add_subplot(G[0,1])
data=metric[:,:,1]
data=np.transpose(data)
data=np.nan_to_num(data)
plot_boxplot(data,2,labels,"$NSEAI$",colors=colors,ax=ax2,ylim=[-7.2,1.2])
print (np.shape(data))
for i in np.arange(0,3):
    print ((sum(data[:,i]>0.0)*1.0)/float(len(data[:,i]))*100.0)
# rISS boxplot
ax3 = fig.add_subplot(G[0,2])
data=metric[:,:,4]
data=np.transpose(data)
data=np.nan_to_num(data)
plot_boxplot(data,3,labels,"$rISS$",colors=colors,ax=ax3,ylim=[-1.0,0.5])

# # sharpness boxplot
# ax3 = fig.add_subplot(G[1,0:2])
# data=metric[:,:,2]
# data=np.transpose(data)
# data=np.nan_to_num(data)
# # plot_boxplot(data,3,labels,"$rSharpness$",colors=colors,ax=ax3)
# plot_cdfplot(data,3,labels,"$rSharpness$",colors=colors,ax=ax3)
# # reliability boxplot
# ax4 = fig.add_subplot(G[1,2:4])
# data=metric[:,:,3]
# data=np.transpose(data)
# data=np.nan_to_num(data)
# plot_boxplot(data,4,labels,"$\Delta$$Reliability$",colors=colors,ax=ax4)
# # rISS boxplot
# ax5 = fig.add_subplot(G[1,4::])
# data=metric[:,:,4]
# data=np.transpose(data)
# data=np.nan_to_num(data)
# plot_boxplot(data,5,labels,"$rISS$",colors=colors,ax=ax5)
#--
plt.subplots_adjust(wspace=0.05,hspace=0.1)
plt.tight_layout()
#--
print ("./figures/"+figname+".png")
plt.savefig("./figures/"+figname+".pdf",dpi=800,bbox_inches="tight", pad_inches=0.01)
plt.savefig("./figures/"+figname+".png",dpi=800,bbox_inches="tight", pad_inches=0.01)
plt.savefig("./figures/"+figname+".jpg",dpi=800,bbox_inches="tight", pad_inches=0.01)
#plt.show()