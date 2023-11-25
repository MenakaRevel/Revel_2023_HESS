#!/opt/local/bin/python
# -*- coding: utf-8 -*-
'''
making IIS boxplot for the experiments
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
argv=sys.argv
syear=int(argv[1])
eyear=int(argv[2])
CaMa_dir=argv[3]
mapname=argv[4]
figname=argv[5]
ncpus=int(sys.argv[6])
ens_mem=21
seaborn_map=True
# seaborn_map=False
print (syear, eyear, CaMa_dir, mapname, figname)
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
#====================================================================
argv=sys.argv
syear=int(argv[1])
eyear=int(argv[2])
CaMa_dir=argv[3]
mapname=argv[4]
figname=argv[5]
ncpus=int(sys.argv[6])
exlist=argv[7]
ens_mem=49
seaborn_map=True
# seaborn_map=False
#=== read experiment names ===
# exlist="./Figs3-experiment_list.nam"
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
#assim_out=pm.DA_dir()+"/out/"+pm.experiment()+"/assim_out"
#assim_out=pm.DA_dir()+"/out/"+experiment+"/assim_out"
# assim_out=DA_dir+"/out/"+experiment
# print (assim_out)
# mk_dir("./figures")
# mk_dir("/figures/NSEAI")
# colors=['xkcd:aqua green','xkcd:pastel blue','xkcd:soft pink','xkcd:pastel blue','xkcd:aqua green','xkcd:soft pink']
# colors=["#004488","#ddaa33","#ba5566"]
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
#===================
rivnum="/cluster/data6/menaka/HydroDA/dat/rivnum_"+mapname+".bin"
rivnum=np.fromfile(rivnum,np.int32).reshape(ny,nx)
rivermap=((nextxy[0]>0)*(rivnum==1))*1.0
#===================
satcov="./dat/satellite_coverage.bin"
satcov=np.fromfile(satcov,np.float32).reshape(ny,nx)
#======================================================
metric=[]
for i,exp in enumerate(experiments):
    print (exp)
    metric_frag=[]
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
        # if np.sum((org!=-9999.0)*1.0) < 365*2:
        #     print ("no obs: ",num, basin, stream, np.sum((org!=-9999.0)*1.0))
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
        if np.sum((org!=-9999.0)*1.0) < 365*1:
            # print ("no obs: ",stream, np.sum((org!=-9999.0)*1.0), "< 365 days")
            continue
        #====================================
        # rISS
        ua, la, uo, lo = read_ul_all(exp,num)
        ISSasm=ISS(ua,la,org,0.05)
        ISSopn=ISS(uo,lo,org,0.05)
        rISS=(ISSasm-ISSopn)/(ISSopn)#+1.0e-20)
        # print exp, num,ISSasm, ISSopn, rISS*100.0 #basin,stream,
        # print exp, num, "%5.2f"%(rISS*100.0)
        # print (exp,num,basin,stream,KGEv)
        # metric_frag.append(rISS)
        # sharpness
        shrasm=sharpness(la,ua,org)
        shropn=sharpness(lo,uo,org)
        rshr=(shrasm-shropn)/shropn
        # rshr_frag.append(rshr)
        # reliability
        relasm=reliability(la,ua,org)
        relopn=reliability(lo,uo,org)
        rrel=relasm-relopn #/relopn
        # rrel_frag.append(rrel)
        #=====================
        metric_frag.append([rISS,rshr,rrel])
    metric.append(metric_frag)
metric=np.array(metric)
#====================================
#---------- making fig --------------
#====================================
# colors=['xkcd:pastel blue','xkcd:aqua green','xkcd:pastel blue','xkcd:aqua green','xkcd:soft pink','xkcd:pastel blue','xkcd:aqua green','xkcd:soft pink']
# colors=['xkcd:aqua green','xkcd:pastel blue','xkcd:soft pink','xkcd:pastel blue','xkcd:aqua green','xkcd:soft pink']
# labels=["Exp 1a","Exp 1b","Exp 2a","Exp 2b","Exp 2c","Exp 3a","Exp 3b","Exp 3c"] 
# labels=["Exp 1a","Exp 1b","Exp 2a","Exp 2c","Exp 3a","Exp 3c"] 
# labels=["Exp 1a","Exp 1b","Exp 2a","Exp 2c","Exp 3a","Exp 3c","Exp 3d","Exp 3e"]
print (np.shape(metric))
# colors=["#afccdc","#3174a1","#b5d294","#3f913a","#f4adae","#bf353c"]
# colors=['#1f77b4','#9467bd','#2ca02c','#d62728']
# colors=["#004488","#ddaa33","#ba5566"]
colors=["#D81B60","#FFC107","#004D40"]
metric_name=["$rISS$", "$rSharpness$","$\Delta$$Reliability$"]
# data=ma.masked_equal(metric,-9999.0)
# data=metric[np.all(metric!=-9999.0,axis=1),:]
# data=ma.masked_less_equal(metric,-100.0)
# print (data)
# metric=np.array(metric)
# print (np.shape(metric))
# plt.clf()  
# plt.close() 
# figure in A4 size
va_margin= 0.0#1.38#inch 
ho_margin= 0.0#1.18#inch
hgt=(11.69 - 2*va_margin)*(1.0/5.0)
wdt=(8.27 - 2*ho_margin)*(2.0/3.0)
#fig=plt.figure(figsize=(8.27,11.69))
fig=plt.figure(figsize=(wdt,hgt))
#fig.suptitle("auto-correalated area")#"maximum autocorrelation length")
#G = gridspec.GridSpec(2,1)
G = gridspec.GridSpec(ncols=2, nrows=1)
#--boxplot
for ii in np.arange(1,2+1):
    #ax = fig.add_subplot(G[1,0])
    data=metric[:,:,ii]
    data=np.transpose(data)
    data=np.nan_to_num(data)
    print (np.shape(data))
    ax = fig.add_subplot(G[0,ii-1])
    #ax.boxplot(data_area, labels=labels, showmeans=True)
    flierprops = dict(marker='o', markerfacecolor='none', markersize=5,linestyle='none', markeredgecolor='k')
    boxprops = dict(color='grey')#facecolor='none'
    whiskerprops = dict(color='grey',linestyle="-",linewidth=1.0)
    capprops = dict(color='grey')
    medianprops = dict(color='grey',linestyle="-",linewidth=1.0)
    meanprops = dict(marker='D', markeredgecolor='black',markerfacecolor='green',markersize=5)
    #with plt.style.context('classic'):#'default'):
    print ("making box plot")
    # print (np.mean(data,axis=0))
    # box=ax.boxplot(data,labels=labels,meanline=True,patch_artist=True,showfliers=False)

    if seaborn_map==False:
        box=ax.boxplot(data,labels=labels,meanprops=meanprops,meanline=True\
            ,boxprops=boxprops,showfliers=False, flierprops=flierprops, whiskerprops=whiskerprops \
            ,capprops=capprops,medianprops=medianprops, notch=False, sym='+' \
            ,vert=None, whis=1.5,positions=None, widths=0.7, patch_artist=True \
            ,bootstrap=None, usermedians=None, conf_intervals=None)
        for patch, color in zip(box['boxes'], colors):
            patch.set_facecolor(color)
        ax.set_ylabel(metric_name[ii], color='k',fontsize=8)
        ax.set_xlabel('Experiment', color='k',fontsize=8)
        ax.tick_params('y',labelsize=6, colors='k')
        ax.tick_params('x',labelsize=6, colors='k')#,labelrotation=45)
        ax.set_xticklabels(labels,rotation=90)
        ax.plot(np.arange(1,lexp+1),np.mean(data,axis=0),linewidth=0,marker="D",markersize=5,color='green',zorder=110)
    #ax.tick_params(axis='x', rotation=45)
    #ax.xaxis.set_tick_params(direction=45)
    #ax.set_title("BoxPlot NSEAI Rivers",fontsize=9,loc="left")
    else:
        stationlist=get_GRDClist("/cluster/data6/menaka/HydroDA/dat/grdc_"+mapname+".txt")
        df=pd.DataFrame(data=data, columns=labels) #, index=stationlist)
        # df=pd.DataFrame(data=data, columns=["Anomaly", "Normalized"])
        # print (df)
        # print (pd.melt(df))
        ax=sns.boxplot(data=df, palette=colors, whis=1.5\
                ,meanline=True, meanprops=meanprops\
                ,capprops=capprops,medianprops=medianprops\
                ,width=0.8, linewidth=0.3, dodge=True\
                ,showfliers=False, saturation=0.75,flierprops=flierprops) #boxprops=boxprops\
        # ax=sns.violinplot(data=df,split=False,showfliers=False, fliersize=0.5, saturation=0.75,\
        #     meanline=True, meanprops=meanprops, width=0.8, whis=1.5, linewidth=0.3,\
        #     capprops=capprops, medianprops=medianprops,scale="count", \
        #     inner="box", bw=0.4, bw_method="silverman", palette="husl") #x=labels,y=stationlist, scale="count",
        # sns.swarmplot(data=df,color="0.25")
        # ax.set_ylabel(metric_name[ii], color='k',fontsize=8)
        ax.axhline(0.0,color="k",linestyle="--",linewidth=0.5)
        ax.set_xlabel('Experiment', color='k',fontsize=8)
        ax.tick_params(labelsize=5)
        ax.set_xticklabels(labels,rotation = 90)
        ax.plot(np.arange(0,lexp),np.mean(data,axis=0),linewidth=0,marker="D",markersize=0.75,color='grey',alpha=0.5,zorder=110)
        # ax.plot(np.arange(0,6),np.median(data,axis=0),linewidth=0,marker="D",markersize=5,color='xkcd:teal blue',zorder=110)
        ax.text(-0.05,1.05,"%s) %s"%(string.ascii_lowercase[ii-1], metric_name[ii]),ha="left",va="center",transform=ax.transAxes,fontsize=10)

        # # get stats
        # stats = boxplot_stats(df.values)
        # print "Experiment","Q1 ","Q3","Mean","Median"
        # for i in np.arange(0,len(labels)):
        #     print ("%10s%15.7f%15.7f%15.7f%15.7f")%(labels[i],stats[i]['q1'],stats[i]['q3'],stats[i]['mean'],stats[i]['med'])
        #     # print (labels[i],"Q1 :",stats[i]['q1'],"Q3 :",stats[i] ['q3'],"Mean :", stats[i]['mean'], "Median",np.median(df[labels[i]]))
        #     # print (labels[i],sum((df[labels[i]]<=-1.e-20)*1.0), len(df[labels[i]]),(1 - (sum(df[labels[i]]<=-1.e-20)*1.0)/float(len(df[labels[i]])))*100.0)
        #     # print labels[i],stats[i]['q1'],stats[i]['q3'],stats[i]['mean'],np.median(df[labels[i]])
    #--
# for xc in [5.5,10.5]:
#     ax.axvline(x=xc,color="k",linestyle="--")
# for i,text in enumerate(['low-latitudes','mid-latitudes','high-latitudes'],start=0):
#     ax.text(i*5.0+3.0,1.05,text,ha='center', va='bottom')
# ax.set_ylim(ymin=-1.2,ymax=1.5)
# figname="NSEAI_boxplot_20210901"
#--
plt.tight_layout(w_pad=1.5)
print ("./figures/"+figname+".png")
plt.savefig("./figures/"+figname+".pdf",dpi=800,bbox_inches="tight", pad_inches=0.01)
plt.savefig("./figures/"+figname+".png",dpi=800,bbox_inches="tight", pad_inches=0.01)
plt.savefig("./figures/"+figname+".jpg",dpi=800,bbox_inches="tight", pad_inches=0.01)
#plt.show()