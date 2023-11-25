#!/usr/bin/env python
# coding=utf-8
# Extract CMF discharge at corresponding gauges 
#
from cProfile import label
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
def read_ul(experiment,station,assim=True):
    u=[]
    l=[]
    # fname="./txt/"+experiment+"/Q_interval/"+station+".txt"
    fname="./txt/"+experiment+"/Q_percentile/"+station+".txt"
    with open(fname,"r") as f:
        lines=f.readlines()
    for line in lines:
        line = re.split(" ",line)
        line = list(filter(None, line))
        if assim:
            try:
                uval=float(line[1])
                lval=float(line[2])
            except:
                uval=0.0
                lval=0.0
        else:
            try:
                uval=float(line[3])
                lval=float(line[4])
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
# argv=sys.argv
# syear=int(argv[1])
# eyear=int(argv[2])
# CaMa_dir=argv[3]
# mapname=argv[4]
# exlist=argv[5]
# figname=argv[6]
# ncpus=int(sys.argv[7])
# ens_mem=49
# seaborn_map=True
# glob=True
# # glob=False
# # seaborn_map=False
# #==================================
syear=2009
eyear=2014
DA_dir="/cluster/data6/menaka/HydroDA"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
mapname="amz_06min"
ens_mem=21
glob=True
seaborn_map=False
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
# with open(exlist,"r") as exf:
#     linesexf=exf.readlines()
# indirs=[]
# labels=[]
# for lineexf in linesexf:
#     lineexf = re.split(":",lineexf)
#     lineexf = list(filter(None, lineexf))
#     labels.append(lineexf[0])
#     indirs.append(lineexf[1].strip())
#     # print lineexf
labels=["CaMa VIC BC","open","DIR","ANO","NOM"]
indirs=["/home/yamadai/data/CaMa_v400_simulations/VIC_BC_3h_06min"\
       ,"DIR_WSE_E2O_HWEB_001","DIR_WSE_E2O_HWEB_001"\
       ,"ANO_WSE_E2O_HWEB_004","NOM_WSE_E2O_HWEB_004"]
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
IX1=[];IY1=[];IX2=[];IY2=[];nums=[];lons=[];lats=[];numsat=[]
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
    if satcov[iy1-1,ix1-1] ==1.0:
        numsat.append(num)
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
            NSEval=NS(dis_cmf.T[:,j],org)
            Rval=correlation(dis_cmf.T[:,j],org)
            KGEval=KGE(dis_cmf.T[:,j],org)
            metric_frag.append([Rval,NSEval,KGEval,-9999,-9999,-9999])
    elif exppre=="open":
        metric_frag=[]
        experiment=indirs[i]
        # print len(nums)
        # print experiment
        for num in nums:
            # experiment="NOM_WSE_E2O_HWEB_003"
            org=grdc.grdc_dis(num,syear,eyear)
            asm, opn = read_dis(experiment,num)
            u,l=read_ul(experiment,num,assim=False)
            NSEval=NS(opn,org)
            Rval=correlation(opn,org)
            KGEval=KGE(opn,org)
            Shapval=sharpness(l,u,org)
            Relval=reliability(l,u,org)
            ISSval=ISS(u,l,org,0.05)
            # print experiment,num,NSEasm
            metric_frag.append([Rval,NSEval,KGEval,Shapval,Relval,ISSval])
    else:
        metric_frag=[]
        experiment=indirs[i]
        # print len(nums)
        # print experiment
        for num in nums:
            # experiment="NOM_WSE_E2O_HWEB_003"
            org=grdc.grdc_dis(num,syear,eyear)
            asm, opn = read_dis(experiment,num)
            u,l=read_ul(experiment,num)
            NSEval=NS(asm,org)
            Rval=correlation(asm,org)
            KGEval=KGE(asm,org)
            Shapval=sharpness(l,u,org)
            Relval=reliability(l,u,org)
            ISSval=ISS(u,l,org,0.05)
            # print experiment,num,NSEasm
            metric_frag.append([Rval,NSEval,KGEval,Shapval,Relval,ISSval])
    metric=np.array(metric_frag)
    df=pd.DataFrame(data=metric, columns=["r","NSE","KGE","Sharpness","Reliability","ISS"], index=nums)
    df.index.name="GRDC ID"
    # print(df.head())
    print (exp)
    print ("===All===")
    print(df.median(axis=0))
    print ("===Sat Cov===")
    print(df.loc[numsat].median(axis=0))
    # df.to_csv(path_or_buf="./tables/absolute_metrics_"+exp+".csv")
