#!/opt/local/bin/python
# -*- coding: utf-8 -*-
'''
making NSEAI boxplot for the experiments
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
#==================================
DA_dir="/cluster/data6/menaka/HydroDA"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
mapname="amz_06min"
ens_mem=21
# lexp=7
colors=['xkcd:pastel blue','xkcd:aqua green','xkcd:soft pink']
seaborn_map=True
# seaborn_map=False
#==================================
ncpus=int(sys.argv[1])
#==================================
# experiment="DIR_WSE_E2O_HWEB_001"
# experiment="ANO_WSE_E2O_HWEB_001"
# experiment="NOM_WSE_E2O_HWEB_001"
# experiments=["DIR_WSE_E2O_HWEB_001","DIR_WSE_E2O_HWEB_002","ANO_WSE_E2O_HWEB_001"\
            # ,"ANO_WSE_E2O_HWEB_002","NOM_WSE_E2O_HWEB_001","NOM_WSE_E2O_HWEB_002"\
            # ,"NOM_WSE_E2O_HWEB_003","NOM_WSE_E2O_HWEB_004","NOM_WSE_E2O_HWEB_005"]
experiments=["NOM_WSE_E2O_HWEB_001","NOM_WSE_E2O_HWEB_002"\
            ,"NOM_WSE_E2O_HWEB_003","NOM_WSE_E2O_HWEB_004"]
lexp=len(experiments)
#assim_out=pm.DA_dir()+"/out/"+pm.experiment()+"/assim_out"
#assim_out=pm.DA_dir()+"/out/"+experiment+"/assim_out"
# assim_out=DA_dir+"/out/"+experiment
# print (assim_out)
#----
def mk_dir(sdir):
  try:
    os.makedirs(sdir)
  except:
    pass
#----
def NS(s,o):
    """
    Nash Sutcliffe efficiency coefficient
    input:
        s: simulated
        o: observed
    output:
        ns: Nash Sutcliffe efficient coefficient
    """
    #s,o = filter_nan(s,o)
    o=ma.masked_where(o<=0.0,o).filled(0.0)
    s=ma.masked_where(o<=0.0,s).filled(0.0)
    o=np.compress(o>0.0,o)
    s=np.compress(o>0.0,s) 
    return 1 - sum((s-o)**2)/(sum((o-np.mean(o))**2)+1e-20)
#----
mk_dir("./figures")
# mk_dir("/figures/NSEAI")
#----
fname=CaMa_dir+"/map/"+mapname+"/params.txt"
f=open(fname,"r")
lines=f.readlines()
f.close()
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
# #higher resolution data
# catmxy = CaMa_dir+"/map/"+mapname+"/1min/1min.catmxy.bin"
# catmxy = np.fromfile(catmxy,np.int16).reshape(2,ny*60,nx*60)
#----
rivnum="../dat/rivnum_"+mapname+".bin"
# rivnum=np.fromfile(rivnum,np.int32).reshape(ny,nx)
# rivermap=((nextxy[0]>0)*(rivnum==1))*1.0
#----
syear,smonth,sdate=2009,1,1 
eyear,emonth,edate=2011,1,1 
#month=1
#date=1
start_dt=datetime.date(syear,smonth,sdate)
end_dt=datetime.date(eyear,emonth,edate)
size=60

start=0
last=(end_dt-start_dt).days
N=int(last)

green2="greenyellow"
green ="green"
#-------
staid=[]
pname=[]
xlist=[]
ylist=[]
river=[]
rivernames = grdc.grdc_river_name_v396(CaMa_dir,mapname)
for rivername in rivernames:
    grdc_id,station_loc,x_list,y_list = grdc.get_grdc_loc_v396(rivername,CaMa_dir,mapname)
    # print (rivername, grdc_id,station_loc)
    river.append([rivername]*len(station_loc))
    staid.append(grdc_id)
    pname.append(station_loc)
    xlist.append(x_list)
    ylist.append(y_list)
#-----------------------
river=([flatten for inner in river for flatten in inner])
staid=([flatten for inner in staid for flatten in inner])
pname=([flatten for inner in pname for flatten in inner])
xlist=([flatten for inner in xlist for flatten in inner])
ylist=([flatten for inner in ylist for flatten in inner])
#-----------------------
# print (len(pname), len(xlist))
#########################################################
pnum=len(pname)
# print (N,ens_mem,pnum,lexp)
#=========================================================
opn=np.ctypeslib.as_ctypes(np.zeros([N,ens_mem,pnum],np.float32))
shared_array_opn  = sharedctypes.RawArray(opn._type_, opn)
asm=np.ctypeslib.as_ctypes(np.zeros([N,ens_mem,pnum],np.float32))
shared_array_asm  = sharedctypes.RawArray(asm._type_, asm)
#=========================================================
# for name in pname:
#     print (name)
# opn=np.ctypeslib.as_ctypes(np.zeros([N,ens_mem,pnum,lexp],np.float32))
# shared_array_opn  = sharedctypes.RawArray(opn._type_, opn)
# asm=np.ctypeslib.as_ctypes(np.zeros([N,ens_mem,pnum,lexp],np.float32))
# shared_array_asm  = sharedctypes.RawArray(asm._type_, asm)

# # for parallel calcualtion
# inputlist=[]
# for exp in np.arange(len(experiments)):
#     expch = "%d"%(exp)
#     for day in np.arange(start,last):
#         target_dt=start_dt+datetime.timedelta(days=day)
#         yyyy='%04d' % (target_dt.year)
#         mm='%02d' % (target_dt.month)
#         dd='%02d' % (target_dt.day)
#         for num in np.arange(1,ens_mem+1):
#             numch='%03d'%num
#             inputlist.append([yyyy,mm,dd,numch,expch])
#             print (yyyy,mm,dd,numch,expch)
#----------------------
# function to read data
def read_data(inputlist):
    yyyy   = inputlist[0]
    mm     = inputlist[1]
    dd     = inputlist[2]
    numch  = inputlist[3]
    exp    = int(inputlist[4])
    experiment = experiments[exp]
    #--
    tmp_opn  = np.ctypeslib.as_array(shared_array_opn)
    tmp_asm  = np.ctypeslib.as_array(shared_array_asm)

    # year, mon, day
    year=int(yyyy)
    mon=int(mm)
    day=int(dd)
    num=int(numch)-1
    #--
    target_dt=datetime.date(year,mon,day)
    dt=(target_dt-start_dt).days
    # corrpted
    fname=DA_dir+"/out/"+experiment+"/assim_out/outflw/open/outflw"+yyyy+mm+dd+"_"+numch+".bin"
    #fname=assim_out+"/assim_out/rivout/open/rivout"+yyyy+mm+dd+"_"+numch+".bin"
    opnfile=np.fromfile(fname,np.float32).reshape([ny,nx])
    # assimilated
    fname=DA_dir+"/out/"+experiment+"/assim_out/outflw/assim/outflw"+yyyy+mm+dd+"_"+numch+".bin"
    #fname=assim_out+"/assim_out/rivout/assim/rivout"+yyyy+mm+dd+"_"+numch+".bin"
    asmfile=np.fromfile(fname,np.float32).reshape([ny,nx])
    #-------------
    for point in np.arange(pnum):
        # print (pname[point],CaMa_dir,mapname)
        # print (grdc.get_grdc_station_v396(pname[point],CaMa_dir,mapname))
        ix1,iy1,ix2,iy2=grdc.get_grdc_station_v396(pname[point],CaMa_dir,mapname)
        if ix2 == -9999 or iy2 == -9999:
            tmp_opn[dt,num,point]=opnfile[iy1,ix1]
            tmp_asm[dt,num,point]=asmfile[iy1,ix1]
        else:
            tmp_opn[dt,num,point]=opnfile[iy1,ix1]+opnfile[iy2,ix2]
            tmp_asm[dt,num,point]=asmfile[iy1,ix1]+asmfile[iy2,ix2]
# #--------
# p   = Pool(20)
# res = p.map(read_data, inputlist)
# opn = np.ctypeslib.as_array(shared_array_opn)
# asm = np.ctypeslib.as_array(shared_array_asm)
# p.terminate()		
#--------------------------------------------
# para_flag = 0
def make_NSEAI():
    para_flag = 1
    metric=[]
    for exp in np.arange(len(experiments)):
        expch = "%d"%(exp)
        #================================================================
        # opn=np.ctypeslib.as_ctypes(np.zeros([N,ens_mem,pnum],np.float32))
        # shared_array_opn  = sharedctypes.RawArray(opn._type_, opn)
        # asm=np.ctypeslib.as_ctypes(np.zeros([N,ens_mem,pnum],np.float32))
        # shared_array_asm  = sharedctypes.RawArray(asm._type_, asm)
        #================================================================
        inputlist=[]
        for day in np.arange(start,last):
            target_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (target_dt.year)
            mm='%02d' % (target_dt.month)
            dd='%02d' % (target_dt.day)
            for num in np.arange(1,ens_mem+1):
                numch='%03d'%num
                inputlist.append([yyyy,mm,dd,numch,expch])
                # print (yyyy,mm,dd,numch,expch)
        #--------
        if para_flag == 1:
            p   = Pool(ncpus)
            res = p.map(read_data, inputlist)
            opn = np.ctypeslib.as_array(shared_array_opn)
            asm = np.ctypeslib.as_array(shared_array_asm)
            p.terminate()
        else:
            res = map(read_data, inputlist)
            opn = np.ctypeslib.as_array(shared_array_opn)
            asm = np.ctypeslib.as_array(shared_array_asm)
        #--------
        metric_frag=[]
        for point in np.arange(pnum):
            org=grdc.grdc_dis(staid[point],syear,eyear-1)
            org=np.array(org)
            obs=np.sum((org!=-9999.0)*1.0)
            #
            NSEasm=NS(np.mean(asm[:,:,point],axis=1),org)
            NSEopn=NS(np.mean(opn[:,:,point],axis=1),org)
            NSEAI=(NSEasm-NSEopn)/(1.0-NSEopn+1.0e-20)
            print (experiments[exp], pname[point], NSEAI)
            if obs==0.0:
                print ("no obs", obs)
                metric_frag.append(-9999.0)
            else:
                metric_frag.append(NSEAI)
        metric.append(metric_frag)
    return np.array(metric)
#====================================
#---------- making fig --------------
#====================================
# labels=["Direct", "Anomaly", "Normalized"\
#         ,"Normalized_AI60","Normalized_AI40","Normalized_AI20"]
# labels=["Direct", "Dir_RH","Anomaly", "Ano_RH", "Normalized"\
#         ,"Norm_AI60","Norm_AI40","Norm_AI20","Norm_RH"] #Restricted HydroWeb
labels=["Norm_90MS","Norm_92MS","Norm_90","Norm_Zero"] #Restricted HydroWeb
fname="./data/metric_20210901.bin"
if os.path.exists(fname):
    metric=np.transpose(np.fromfile(fname,np.float64).reshape(lexp,-1))
else:
    metric=make_NSEAI()
    metric.tofile(fname)
    metric=np.transpose(metric)
# print (np.shape(metric))
data=metric
# data=ma.masked_equal(metric,-9999.0)
data=metric[np.all(metric!=-9999.0,axis=1),:]
# data=ma.masked_less_equal(metric,-100.0)
# print (data)
# metric=np.array(metric)
# print (np.shape(metric))
# plt.clf()  
# plt.close() 
# figure in A4 size
va_margin= 0.0#1.38#inch 
ho_margin= 0.0#1.18#inch
hgt=(11.69 - 2*va_margin)*(1.0/3.0)
wdt=(8.27 - 2*ho_margin)*(1.0/2.0)
#fig=plt.figure(figsize=(8.27,11.69))
fig=plt.figure(figsize=(wdt,hgt))
#fig.suptitle("auto-correalated area")#"maximum autocorrelation length")
#G = gridspec.GridSpec(2,1)
G = gridspec.GridSpec(1,1)
#--boxplot
#ax = fig.add_subplot(G[1,0])
ax = fig.add_subplot(G[0,0])
#ax.boxplot(data_area, labels=labels, showmeans=True)
flierprops = dict(marker='o', markerfacecolor='none', markersize=8,linestyle='none', markeredgecolor='k')
boxprops = dict(color='grey')#facecolor='none'
whiskerprops = dict(color='grey',linestyle="--")
capprops = dict(color='grey')
medianprops = dict(color='r')
meanprops = dict(marker='D', markeredgecolor='black',markerfacecolor='green',markersize=8)
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
    ax.set_ylabel('NSEAI', color='k',fontsize=8)
    ax.set_xlabel('Experiment', color='k',fontsize=8)
    ax.tick_params('y',labelsize=6, colors='k')
    ax.tick_params('x',labelsize=6, colors='k')#,labelrotation=45)
    ax.set_xticklabels(labels,rotation=90)
    ax.plot(np.arange(1,lexp+1),np.mean(data,axis=0),linewidth=0,marker="D",markersize=5,color='green',zorder=110)
#ax.tick_params(axis='x', rotation=45)
#ax.xaxis.set_tick_params(direction=45)
#ax.set_title("BoxPlot NSEAI Rivers",fontsize=9,loc="left")
else:
    df=pd.DataFrame(data=data, columns=labels) #, rows=np.arange(0,993))
    # df=pd.DataFrame(data=data, columns=["Anomaly", "Normalized"])
    # print (df)
    # print (pd.melt(df))
    # ax=sns.boxplot(data=df, fliersize=0.0, palette="Set2", whis=1.5, dodge=...)
    ax=sns.violinplot(data=df,split=False, fliersize=0.5, saturation=0.75, \
        meanline=True, width=0.8, whis=1.5, linewidth=0.3, scale="count", \
        inner="point", bw=0.35, palette="Set2")
    # sns.swarmplot(data=df,color="0.25")
    ax.set_ylabel('NSEAI', color='k',fontsize=8)
    ax.set_xlabel('Experiment', color='k',fontsize=8)
    ax.tick_params(labelsize=5)
    ax.set_xticklabels(labels,rotation = 90)
    # ax.plot(np.arange(0,3),np.mean(data,axis=0),linewidth=0,marker="D",markersize=5,color='xkcd:teal blue',zorder=110)

    # get stats
    stats = boxplot_stats(df.values)
    for i in np.arange(0,len(labels)):
        # print ("Q1 :",stats[i]['q1'],"Q3 :",stats[i] ['q3'],"Mean :", stats[i]['mean'])
        print (labels[i],sum((df[labels[i]]<=-1.e-20)*1.0), len(df[labels[i]]),(1 - (sum(df[labels[i]]<=-1.e-20)*1.0)/float(len(df[labels[i]])))*100.0)
#--
# for xc in [5.5,10.5]:
#     ax.axvline(x=xc,color="k",linestyle="--")
# for i,text in enumerate(['low-latitudes','mid-latitudes','high-latitudes'],start=0):
#     ax.text(i*5.0+3.0,1.05,text,ha='center', va='bottom')
ax.set_ylim(ymin=-6,ymax=4)


figname="NSEAI_boxplot_20210901"
#--
print ("./figures/"+figname+".png")
plt.savefig("./figures/"+figname+".png",dpi=800,bbox_inches="tight", pad_inches=0.05)
# plt.savefig("./figures/"+figname+".jpg",dpi=800,bbox_inches="tight", pad_inches=0.0)
#plt.show()