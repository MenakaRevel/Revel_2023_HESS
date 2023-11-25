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
IX1=[];IY1=[];IX2=[];IY2=[];nums=[];names=[]
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
            print num, station, stream
            IX1.append(ix1-1)
            IY1.append(iy1-1)
            IX2.append(ix2-1)
            IY2.append(iy2-1)
            nums.append(num)
            names.append(stream)
#===========================
# indirs=["/cluster/data6/menaka/ensemble_org/CaMa_out/AMZE2O003"\
#        ,"/cluster/data6/menaka/ensemble_org/CaMa_out/AMZCALecmwf001"\
#        ,"./txt/NOM_WSE_E2O_HWEB_003/outflw"]
metric=[]
for i,exp in enumerate(labels):
    exppre=exp.split(" ")[0]
    # if i < 2:
    if exppre=="CaMa":
        dis_cmf = read_dis_multi(IX1, IY1, IX2, IY2, syear, eyear, indirs[i], nx, ny) 
        # print np.shape(dis_cmf.T)
        metric_frag=[]
        for j,num in enumerate(nums):
            # org=grdc.grdc_dis(num,syear,eyear)
            # # print len(dis_cmf.T[:,j])
            # NSEasm=NS(dis_cmf.T[:,j],org)
            # if j==0:
            #     print "========================="
            #     for xx in np.arange(len(org)):
            #         print org[xx], dis_cmf.T[xx,j]
            # print num, NSEasm
            metric_frag.append(dis_cmf.T[:,j])
    elif exppre=="Open":
        metric_frag=[]
        # print len(nums)
        for num in nums:
            experiment=indirs[i]
            # experiment="NOM_WSE_E2O_HWEB_003"
            # org=grdc.grdc_dis(num,syear,eyear)
            asm, opn = read_dis(experiment,num)
            # NSEasm=NS(asm,org)
            # print experiment,num,NSEasm
            metric_frag.append(opn)
    else:
        metric_frag=[]
        # print len(nums)
        for num in nums:
            experiment=indirs[i]
            # experiment="NOM_WSE_E2O_HWEB_003"
            # org=grdc.grdc_dis(num,syear,eyear)
            asm, opn = read_dis(experiment,num)
            # NSEasm=NS(asm,org)
            # print experiment,num,NSEasm
            metric_frag.append(asm)
    metric.append(metric_frag)
metric=np.array(metric)
print (np.shape(metric))
#=============
# make figure
#=============
plt.close()
lexp=len(labels)
lexp0=len(stations)
# lexp=len(stations)
if lexp==6:
    colors=['#1f77b4','#9467bd','#2ca02c','#d62728','xkcd:goldenrod','xkcd:crimson']
elif lexp==4:
    colors=["grey","#004488","#ddaa33","#ba5566"]
    # colors=['#1f77b4','#9467bd','#2ca02c','#d62728']
elif lexp==3:
    colors=["#004488","#ddaa33","#ba5566"]
    # colors=["#4dc7ec","#ff8021",'#2ca02c']
elif lexp==2:
    colors=["#4dc7ec","#ff8021"]
else:
    colors=["#4dc7ec","#ff8021"]
# figure in A4 size
va_margin= 0.0#1.38#inch 
ho_margin= 0.0#1.18#inch
hgt=(11.69 - 2*va_margin)*(1.0/4.0)*(lexp0/2)
wdt=(8.27 - 2*ho_margin)*(2.0/2.0)
fig=plt.figure(figsize=(wdt,hgt))
# print lexp/2
cols=2
rows=len(stations)/cols
G = gridspec.GridSpec(nrows=rows,ncols=cols)
# G.update(wspace=0.0, hspace=0.0) # set the spacing between axes. 
for j,num in enumerate(nums):
    x=int(j%cols)
    y=int(j/cols)
    print ("x, y: ",x, y)
    ax = fig.add_subplot(G[y,x])
    station=num
    org=grdc.grdc_dis(station,syear,eyear)
    org=np.array(org)
    ax.plot(np.arange(start,last),ma.masked_less(org,0.0),label="GRDC",color="#34495e",linewidth=3.0,zorder=101) #,marker = "o",markevery=swt[point])
    mxx=[]
    for i,exp in enumerate(labels):
        linewidth=0.5 #1.0/float(len(labels))
        print ("i, j:",i, j, exp)
        print (np.shape(metric[i,j,:]), colors[i])
        exppre=exp.split(" ")[0]
        if exppre=="Open":
            linestyle="--"
        else:
            linestyle="-"
        ax.plot(np.arange(start,last),metric[i,j,:],label=exp,color=colors[i],linewidth=linewidth,linestyle=linestyle,alpha=1,zorder=104)
        r=correlation(metric[i,j,:],org)
        NSE=NS(metric[i,j,:],org)
        KGE1=KGE(metric[i,j,:],org)
        mxx.append([r,NSE,KGE1])
    # Make the y-axis label, ticks and tick labels match the line color.
    if x==0:
        ax.set_ylabel('discharge (m$^3$/s)', color='k')
    ax.set_xlim(xmin=0,xmax=last+1)
    ax.tick_params('y', colors='k')
    # scentific notaion
    ax.ticklabel_format(style="sci",axis="y",scilimits=(0,0))
    ax.yaxis.major.formatter._useMathText=True 
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
    if y==lexp0/2 - 1:
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
    #========
    # NSE="NSE: %3.2f"%(NS(s,o))
    #============
    # data frame
    #============
    print ("======================")
    print (names[j])
    mxx=np.array(mxx)
    df=pd.DataFrame(mxx,columns=["r","NSE","KGE"],index=labels)
    print (df)
# legend
labels0=["Observations"]
labels0.extend(labels)
colors0=["#34495e"]
colors0.extend(colors)
linestyles=["-","--","-","-","-"]
# print labels0
features=[]
for i in np.arange(len(labels0)):
    label=labels0[i]
    features.append(mlines.Line2D([], [], color=colors0[i], marker=None,
                          label=labels0[i],linewidth=2.5,linestyle=linestyles[i])) #edgecolors="k",
legend=plt.legend(handles=features,bbox_to_anchor=(0.5,0.0), loc="upper center",
           bbox_transform=fig.transFigure, ncol=len(labels0),  borderaxespad=0.0, frameon=False)#
#======================
# fig.tight_layout()
# plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.1)#wspace=0, hspace=0)
# plt.subplots_adjust(wspace=0, hspace=0)
#plt.title(stitle)
print figname
plt.savefig("./figures/"+figname+".jpg",dpi=800,bbox_inches="tight", pad_inches=0.05)
plt.savefig("./figures/"+figname+".png",dpi=800,bbox_inches="tight", pad_inches=0.05)
plt.savefig("./figures/"+figname+".pdf",dpi=800,bbox_inches="tight", pad_inches=0.05)