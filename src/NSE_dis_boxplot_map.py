#!/usr/bin/env python
# coding=utf-8
# Extract CMF discharge at corresponding gauges 
#
import os, sys
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm,Normalize,ListedColormap,BoundaryNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.cbook import boxplot_stats
import numpy as np
from netCDF4 import Dataset
import datetime
import string
import re
import seaborn as sns
import pandas as pd
from numpy import ma
import warnings;warnings.filterwarnings('ignore')

from read_CMF import read_discharge, read_discharge_multi
import read_grdc as grdc
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
            # print (lon1, lon2)
            lon2=180.0
        elif lon2-lon1> 180.0:
            # print (lon1,lon2)
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
def global_coordinates(station,mapname="glb_06min",CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"):
    obslist=CaMa_dir+"/map/"+mapname+"/grdc_loc.txt"
    with open(obslist,"r") as f:
        lines=f.readlines()
    ixx1=-9999
    iyy1=-9999
    ixx2=-9999
    iyy2=-9999
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
        if staid == int(station):
            # print staid, stations
            ixx1=ix1
            iyy1=iy1
            ixx2=ix2
            iyy2=iy2
            break
    return ixx1, iyy1, ixx2, iyy2
#====================================================================
def global_list(lstation):
    IX1=[];IY1=[];IX2=[];IY2=[]
    for station in lstation:
        ix1,iy1,ix2,iy2=global_coordinates(station)
        # if ix1==-9999:
        #     print station
        IX1.append(ix1-1)
        IY1.append(iy1-1)
        IX2.append(ix2-1)
        IY2.append(iy2-1)
    return IX1,IY1,IX2,IY2
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
# === main code ====
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
glob=True
# glob=False
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
# exlist="./Fig09-experiment_list.nam"
with open(exlist,"r") as exf:
    linesexf=exf.readlines()
indirs=[]
labels=[]
for lineexf in linesexf:
    lineexf = re.split(":",lineexf)
    lineexf = list(filter(None, lineexf))
    labels.append(lineexf[0])
    indirs.append(lineexf[1].strip())
    # print lineexf
#=============================
lexp=len(labels)
print labels
print indirs
colors=['xkcd:aqua green','xkcd:pastel blue','xkcd:soft pink','xkcd:pastel blue','xkcd:aqua green','xkcd:soft pink']

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
nbdays=int(last)
#======================================================
obslist="/cluster/data6/menaka/HydroDA/dat/grdc_"+mapname+".txt"
with open(obslist,"r") as f:
    lines=f.readlines()
IX1=[];IY1=[];IX2=[];IY2=[];nums=[];lons=[];lats=[]
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
    # print num, stream
    IX1.append(ix1-1)
    IY1.append(iy1-1)
    IX2.append(ix2-1)
    IY2.append(iy2-1)
    nums.append(num)
    lons.append(lonlat[0,iy1-1,ix1-1])
    lats.append(lonlat[1,iy1-1,ix1-1])
#===========================
# indirs=["/cluster/data6/menaka/ensemble_org/CaMa_out/AMZE2O003"\
#        ,"/cluster/data6/menaka/ensemble_org/CaMa_out/AMZCALecmwf001"\
#        ,"./txt/NOM_WSE_E2O_HWEB_003/outflw"]
metric=[]
for i,exp in enumerate(labels):
    exppre=exp.split(" ")[0]
    # print exp, exppre
    # if i == 0:
    if exppre=="CaMa":
        if glob:
            IXX1,IYY1,IXX2,IYY2=global_list(nums)
            # print len(IX1), len(IXX1)
            nxx=3600
            nyy=1800
            dis_cmf = read_dis_multi(IXX1, IYY1, IXX2, IYY2, syear, eyear, indirs[i], nxx, nyy)
        else:
            dis_cmf = read_dis_multi(IX1, IY1, IX2, IY2, syear, eyear, indirs[i], nx, ny) 
        # print np.shape(dis_cmf.T)
        metric_frag=[]
        for j,num in enumerate(nums):
            org=grdc.grdc_dis(num,syear,eyear)
            # print len(dis_cmf.T[:,j])
            NSEasm=NS(dis_cmf.T[:,j],org)
            # if j==0:
            #     print "========================="
            #     for xx in np.arange(len(org)):
            #         print org[xx], dis_cmf.T[xx,j]
            # print num, NSEasm
            metric_frag.append(NSEasm)
    else:
        metric_frag=[]
        experiment=indirs[i]
        # print len(nums)
        # print experiment
        for num in nums:
            # experiment="NOM_WSE_E2O_HWEB_003"
            org=grdc.grdc_dis(num,syear,eyear)
            asm, opn = read_dis(experiment,num)
            NSEasm=NS(asm,org)
            # print experiment,num,NSEasm
            metric_frag.append(NSEasm)
    metric.append(metric_frag)
#====================================
#---------- making fig --------------
#====================================
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

# cmap=mbar.colormap("H01")
# vmin=0.0
# vmax=1.0
# # norm=Normalize(vmin=vmin,vmax=vmax)
# norm=BoundaryNorm(np.arange(vmin,vmax+0.1,0.1),cmap.N)

#------
# river width
sup=3 #2
w=0.002 #0.02
alpha=1
width=0.5
#------
west0=west-5
east0=east
south0=south-2
north0=north+2
#------
# colors=['xkcd:pastel blue','xkcd:aqua green','xkcd:pastel blue','xkcd:aqua green','xkcd:soft pink','xkcd:pastel blue','xkcd:aqua green','xkcd:soft pink']
# colors=['xkcd:aqua green','xkcd:pastel blue','xkcd:soft pink','xkcd:pastel blue','xkcd:aqua green','xkcd:soft pink']
# colors=['xkcd:sunflower','xkcd:deep lavender','xkcd:cerulean','xkcd:bluish green']
# colors=['xkcd:light urple','xkcd:turquoise blue','xkcd:green apple','xkcd:sandy yellow']
# colors=['#9467bd','#1f77b4','#2ca02c','#d62728']
markers=["o","D","s","^"]
# colors=['#9B62A7',"#004488","#ddaa33","#ba5566"]
colors=['#1E88E5',"#D81B60","#FFC107","#004D40"]
my_pallte={}
for i in np.arange(len(labels)):
    my_pallte[labels[i]]=colors[i]
# labels=["Exp 1a","Exp 1b","Exp 2a","Exp 2b","Exp 2c","Exp 3a","Exp 3b","Exp 3c"] 
# labels=["Exp 1a","Exp 1b","Exp 2a","Exp 2c","Exp 3a","Exp 3c"] 
# labels=["Exp 1a","Exp 1b","Exp 2a","Exp 2c","Exp 3a","Exp 3c","Exp 3d","Exp 3e"]
metric=np.array(metric)
data=np.transpose(metric)
# data=metric
data=np.nan_to_num(data)
print (np.shape(metric))
# figure in A4 size
va_margin= 0.0#1.38#inch 
ho_margin= 0.0#1.18#inch
hgt=(11.69 - 2*va_margin)*(1.0/3.0)
wdt=(8.27 - 2*ho_margin)*(2.0/2.0)
#fig=plt.figure(figsize=(8.27,11.69))
fig=plt.figure(figsize=(wdt,hgt))
#fig.suptitle("auto-correalated area")#"maximum autocorrelation length")
G = gridspec.GridSpec(ncols=2, nrows=1)
# G = gridspec.GridSpec(1,1)
#--boxplot
#ax = fig.add_subplot(G[1,0])
ax1 = fig.add_subplot(G[0,0])
#ax1.boxplot(data_area, labels=labels, showmeans=True)
flierprops = dict(marker='o', markerfacecolor='none', markersize=8,linestyle='none', markeredgecolor='k')
boxprops = dict(color='grey')#facecolor='none'
whiskerprops = dict(color='grey',linestyle="--")
capprops = dict(color='grey')
medianprops = dict(color='grey',linestyle="-",linewidth=1.0)
meanprops = dict(marker='D', markeredgecolor='black',markerfacecolor='green',markersize=8)
#with plt.style.context('classic'):#'default'):
print ("making box plot")
# print (np.mean(data,axis=0))
# box=ax.boxplot(data,labels=labels,meanline=True,patch_artist=True,showfliers=False)
# labels=["A","B","C","D"]
cols={}
for col,label in zip(colors,labels):
    cols[label]=col
maks={}
for mak,label in zip(markers,labels):
    maks[label]=mak
# colors={'CaMa NoCal':'#1f77b4','CaMa Cal':'#9467bd','Assim NoCal':'#2ca02c','Assim Cal':'#d62728'}
# colors={'A':'xkcd:light urple','B':'xkcd:turquoise blue','C':'xkcd:green apple','D':'xkcd:sandy yellow'}
df=pd.DataFrame(data=metric.T, columns=labels, index=nums)
print (df.head())
# df["lons"]=lons
# df["lats"]=lats
# print(df.idxmax(axis=1, skipna=True))
bestQ=df.idxmax(axis=1, skipna=True)
#====================================
if seaborn_map==False:
    box=ax1.boxplot(data,labels=labels,meanprops=meanprops,meanline=True\
        ,boxprops=boxprops,showfliers=False, flierprops=flierprops, whiskerprops=whiskerprops \
        ,capprops=capprops,medianprops=medianprops, notch=False, sym='+' \
        ,vert=None, whis=1.5,positions=None, widths=0.7, patch_artist=True \
        ,bootstrap=None, usermedians=None, conf_intervals=None)
    for patch, color in zip(box['boxes'], colors):
        patch.set_facecolor(color)
    ax1.set_ylabel('$NSE$', color='k',fontsize=8)
    # ax1.set_xlabel('Experiment', color='k',fontsize=8)
    ax1.tick_params('y',labelsize=6, colors='k')
    ax1.tick_params('x',labelsize=6, colors='k')#,labelrotation=45)
    ax1.set_xticklabels(labels,rotation=90)
    ax1.plot(np.arange(1,lexp+1),np.mean(data,axis=0),linewidth=0,marker="D",markersize=5,color='green',zorder=110)
#ax.tick_params(axis='x', rotation=45)
#ax.xaxis.set_tick_params(direction=45)
#ax.set_title("BoxPlot NSEAI Rivers",fontsize=9,loc="left")
else:
    # stationlist=get_GRDClist("/cluster/data6/menaka/HydroDA/dat/grdc_"+mapname+".txt")
    # df=pd.DataFrame(data=data, columns=labels)#, index=stationlist)
    # df=pd.DataFrame(data=data, columns=["Anomaly", "Normalized"])
    # print (df)
    # print (pd.melt(df))
    box=sns.boxplot(ax=ax1,data=df, fliersize=0.0, palette=my_pallte, whis=1.5\
        ,meanline=True, width=0.8, linewidth=0.3, dodge=True\
        ,meanprops=meanprops,capprops=capprops,medianprops=medianprops)
    # ax1=sns.violinplot(data=df,split=False, fliersize=0.5, saturation=0.75, \
    #     meanline=True, width=0.8, whis=1.5, linewidth=0.3, scale="count", \
    #     inner="point", bw=0.35, bw_method="silverman")#, palette="Set2") #x=labels,y=stationlist,
    # sns.swarmplot(data=df,color="0.25")
    ax1.set_ylabel('$NSE$', color='k',fontsize=12)
    # ax1.set_xlabel('Experiment', color='k',fontsize=12)
    ax1.tick_params(labelsize=8)
    ax1.set_xticklabels(labels,rotation = 90)
    # ax1.plot(np.arange(0,3),np.mean(data,ax1is=0),linewidth=0,marker="D",markersize=5,color='xkcd:teal blue',zorder=110)
    # ax.plot(np.arange(0,6),np.median(data,axis=0),linewidth=0,marker="D",markersize=5,color='xkcd:teal blue',zorder=110)

    # # get stats
    # stats = boxplot_stats(df.values)
    # print "Experiment","Q1 ","Q3","Mean","Median"
    # for i in np.arange(0,len(labels)):
    #     print ("%15s%15.7f%15.7f%15.7f%15.7f")%(labels[i],stats[i]['q1'],stats[i]['q3'],stats[i]['mean'],stats[i]['med'])
    #     # print (labels[i],"Q1 :",stats[i]['q1'],"Q3 :",stats[i] ['q3'],"Mean :", stats[i]['mean'], "Median",np.median(df[labels[i]]))
    #     # print (labels[i],sum((df[labels[i]]<=-1.e-20)*1.0), len(df[labels[i]]),(1 - (sum(df[labels[i]]<=-1.e-20)*1.0)/float(len(df[labels[i]])))*100.0)
    #     # print labels[i],stats[i]['q1'],stats[i]['q3'],stats[i]['mean'],np.median(df[labels[i]])
        
    # # satellite coverage
    # obslist="/cluster/data6/menaka/HydroDA/dat/grdc_"+mapname+".txt"
    # stationlist1=get_GRDClist(obslist,satellite=True)
    # df1=df.loc[stationlist1]
    # stats = boxplot_stats(df1.values)
    # print "Experiment","Q1 ","Q3","Mean","Median"
    # for i in np.arange(0,len(labels)):
    #     print ("%15s%15.7f%15.7f%15.7f%15.7f")%(labels[i],stats[i]['q1'],stats[i]['q3'],stats[i]['mean'],stats[i]['med'])       
#--
# for xc in [5.5,10.5]:
#     ax.axvline(x=xc,color="k",linestyle="--")
# for i,text in enumerate(['low-latitudes','mid-latitudes','high-latitudes'],start=0):
#     ax.text(i*5.0+3.0,1.05,text,ha='center', va='bottom')
ax1.set_ylim(ymin=-3.2,ymax=1.2)
#======================
ax1.text(-0.05,1.05,"%s)"%(string.ascii_lowercase[0]),ha="left",va="center",transform=ax1.transAxes,fontsize=10)
#================
# map of discharge product
ax2 = fig.add_subplot(G[0,1])
m = Basemap(projection='cyl',llcrnrlat=south0,urcrnrlat=north0,llcrnrlon=west0,urcrnrlon=east0, lat_ts=0,resolution='c',ax=ax2)
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
pnum=len(lons)
for i in np.arange(pnum):
    lon=lons[i]
    lat=lats[i]
    beQ=bestQ[i]
    col=cols[beQ]
    mak=maks[beQ]
    ax2.scatter(lon,lat,s=15,marker=mak,edgecolors="k",linewidth=0.5, facecolors=col,zorder=106)
#======================
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.spines['left'].set_visible(False)
#======================
ax2.text(0.05,1.20,"%s)"%(string.ascii_lowercase[1]),ha="left",va="center",transform=ax2.transAxes,fontsize=10)
#======================
# # legend 
# patches=[]
# pnum=len(colors)
# for i in np.arange(pnum):
#     label=labels[i]
#     patches.append(mpatches.Patch(color=cols[label],label=label))
# legend=plt.legend(handles=patches,bbox_to_anchor=(0.90,0.01), loc="lower right",
#            bbox_transform=fig.transFigure, ncol=2,  borderaxespad=0.0, frameon=False)#
#======================
# legend 
features=[]
pnum=len(colors)
for i in np.arange(pnum):
    label=labels[i]
    features.append(mlines.Line2D([], [], color=cols[label], marker=maks[label],
                          markeredgecolor="k",markersize=5, markeredgewidth=0.5,
                          label=label,linewidth=0.0)) #edgecolors="k",
legend=plt.legend(handles=features,bbox_to_anchor=(0.85,0.01), loc="lower right",
           bbox_transform=fig.transFigure, ncol=2,  borderaxespad=0.0, frameon=False)#
#======================
plt.subplots_adjust(wspace=0.01, hspace=0.01)
#======================
print ("./figures/"+figname+".png")
#======================
plt.savefig("./figures/"+figname+".pdf",dpi=800,bbox_inches="tight", pad_inches=0.0)
plt.savefig("./figures/"+figname+".png",dpi=800,bbox_inches="tight", pad_inches=0.0)
plt.savefig("./figures/"+figname+".jpg",dpi=800,bbox_inches="tight", pad_inches=0.0)
os.system("rm -r "+re.split("-",figname)[0]+"*.txt")