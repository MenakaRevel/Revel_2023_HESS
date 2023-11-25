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
import string
import seaborn as sns
import pandas as pd
import datetime
import warnings;warnings.filterwarnings('ignore')

import read_grdc as grdc
import my_colorbar as mbar
import read_hydroweb as hweb
from read_CMF import read_discharge, read_discharge_multi
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
    txt="corrtmp_%02d.txt"%(LEVEL)
    os.system("./bin/print_rivvec corrtmp1.txt 1 "+str(LEVEL)+" > "+txt)
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
argv=sys.argv
syear=int(argv[1])
eyear=int(argv[2])
CaMa_dir=argv[3]
mapname=argv[4]
exlist=argv[5]
figname=argv[6]
ncpus=int(sys.argv[7])
ens_mem=21
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
#=============================
lexp=len(labels)
colors=['xkcd:aqua green','xkcd:pastel blue','xkcd:soft pink','xkcd:pastel blue','xkcd:aqua green','xkcd:soft pink']
labels0=labels
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
    if i == 0:
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
        # print len(nums)
        for num in nums:
            experiment=indirs[i]
            # experiment="NOM_WSE_E2O_HWEB_003"
            org=grdc.grdc_dis(num,syear,eyear)
            asm, opn = read_dis(experiment,num)
            NSEasm=NS(asm,org)
            # print experiment,num,NSEasm
            metric_frag.append(NSEasm)
    metric.append(metric_frag)
#=============================
metric=np.array(metric)
labels=["A","B","C","D"]
# colors=['xkcd:aqua green','xkcd:pastel blue','xkcd:soft pink','xkcd:brick']
# colors={'A':'xkcd:aqua green','B':'xkcd:pastel blue','C':'xkcd:soft pink','D':'xkcd:brick'}
# colors={'A':'xkcd:sunflower yellow','B':'xkcd:muted pink','C':'xkcd:dodger blue','D':'xkcd:green/blue'}
# colors={'A':'xkcd:sunflower','B':'xkcd:deep lavender','C':'xkcd:sky blue','D':'xkcd:tealish'}
# colors={'A':'xkcd:sunflower','B':'xkcd:deep lavender','C':'xkcd:cerulean','D':'xkcd:bluish green'}
# colors={'A':'xkcd:sandy yellow','B':'xkcd:green apple','C':'xkcd:turquoise blue','D':'xkcd:light urple'}
colors={'A':'#1f77b4','B':'#9467bd','C':'#2ca02c','D':'#d62728'}
# colors={'A':'xkcd:light urple','B':'xkcd:turquoise blue','C':'xkcd:green apple','D':'xkcd:sandy yellow'}
df=pd.DataFrame(data=metric.T, columns=labels, index=nums)
print (df.head())
# df["lons"]=lons
# df["lats"]=lats
# print(df.idxmax(axis=1, skipna=True))
bestQ=df.idxmax(axis=1, skipna=True)
for label in labels:
    print(label, np.sum((bestQ==label)*1.0), (np.sum((bestQ==label)*1.0)/len(lons))*100.0)
# # #====================================
# # #---------- making fig --------------
# # #====================================
# # land="#C0C0C0"
# # water="#FFFFFF"

# # londiff=(east-west)*4
# # latdiff=(north-south)*4

# # npix=(90-north)*4
# # spix=(90-south)*4
# # wpix=(180+west)*4
# # epix=(180+east)*4

# # cmap=mbar.colormap("H01")
# # vmin=0.0
# # vmax=1.0
# # # norm=Normalize(vmin=vmin,vmax=vmax)
# # norm=BoundaryNorm(np.arange(vmin,vmax+0.1,0.1),cmap.N)

# # #------
# # # river width
# # sup=3 #2
# # w=0.002 #0.02
# # alpha=1
# # width=0.5
# # #------
# # west0=west
# # east0=east
# # south0=south-2
# # north0=north+2
# # #------
# # # figure in A4 size
# # va_margin= 0.0#1.38#inch 
# # ho_margin= 0.0#1.18#inch
# # hgt=(11.69 - 2*va_margin)*(1.0/3.0)
# # wdt=(8.27 - 2*ho_margin)*(1.0/2.0)
# # fig=plt.figure(figsize=(wdt,hgt))
# # G = gridspec.GridSpec(1,1)
# # ax = fig.add_subplot(G[0,0])
# # m = Basemap(projection='cyl',llcrnrlat=south0,urcrnrlat=north0,llcrnrlon=west0,urcrnrlon=east0, lat_ts=0,resolution='c',ax=ax)
# # # amazon shapefile
# # m.readshapefile('./etc/Amazon_basin/Amazon_boundry', 'Amazon_boundry', drawbounds = False)
# # for info, shape in zip(m.Amazon_boundry, m.Amazon_boundry):
# #     # if info['nombre'] == 'Amazon_boundry':
# #     x, y = zip(*shape) 
# #     m.plot(x, y, marker=None,color='grey',linewidth=1.0)
# # box="%f %f %f %f"%(west,east,north,south) 
# # os.system("./bin/txt_vector "+box+" "+CaMa_dir+" "+mapname+" > corrtmp1.txt") 
# # #map(vec_par,np.arange(1,10+1,1))
# # map(vec_par,np.arange(2,10+1,1))
# # # map(vec_par,np.arange(5,10+1,1))
# # pnum=len(lons)
# # for i in np.arange(pnum):
# #     lon=lons[i]
# #     lat=lats[i]
# #     beQ=bestQ[i]
# #     col=colors[beQ]
# #     ax.scatter(lon,lat,s=15,marker="o",edgecolors="k",linewidth=0.5, facecolors=col,zorder=106)
# # #======================
# # ax.spines['top'].set_visible(False)
# # ax.spines['right'].set_visible(False)
# # ax.spines['bottom'].set_visible(False)
# # ax.spines['left'].set_visible(False)
# # #======================
# # # legend 
# # features=[]
# # # colors=['xkcd:aqua green','xkcd:pastel blue','xkcd:soft pink','xkcd:brick']
# # # colors=['xkcd:sunflower yellow','xkcd:muted pink','xkcd:dodger blue','xkcd:green/blue']
# # pnum=len(colors)
# # for i in np.arange(pnum):
# #     label=labels[i]
# #     features.append(mlines.Line2D([], [], color=colors[label], marker='o',
# #                           markeredgecolor="k",markersize=5, markeredgewidth=0.5,
# #                           label=labels0[i],linewidth=0.0)) #edgecolors="k",
# # # plt.axis('off')
# # # legend=plt.legend(handles=[feature1,feature2],bbox_to_anchor=(0.8, 0.001, 0.15, 0.1), loc=3,
# # #            ncol=1,  borderaxespad=0.0)#,fancybox=False) mode="expand",
# # legend=plt.legend(handles=features,bbox_to_anchor=(0.1,0.01,0.8,0.1), loc="lower left",
# #            bbox_transform=fig.transFigure, ncol=2,  borderaxespad=0.0, frameon=False)#
# # #======================
# # plt.subplots_adjust(wspace=0, hspace=0)
# # #plt.title(stitle)
# # plt.savefig("./figures/"+figname+".png",dpi=800,bbox_inches="tight", pad_inches=0.0)
# # plt.savefig("./figures/"+figname+".pdf",dpi=800,bbox_inches="tight", pad_inches=0.0)
# # os.system("rm -r corrtmp*.txt")