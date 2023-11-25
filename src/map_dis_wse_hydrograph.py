#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
# import datetime
from matplotlib.colors import LogNorm,Normalize,ListedColormap,BoundaryNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator,FormatStrFormatter
import sys
import os
import calendar
import datetime
import math
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

import read_grdc as grdc
import my_colorbar as mbar
import read_hydroweb as hweb
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
#===================
mk_dir("./figures")
# mk_dir("./figures/corr")
#===================
# # syear=2009
# # eyear=2014
# # # CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
# # CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"
# # #map="glb_15min"
# # # map="glb_06min"
# # mapname="amz_06min"
# # lexp=["DIR_WSE_E2O_HWEB_001","DIR_WSE_E2O_HWEB_002"]
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
figname=argv[5]
expname=argv[6]
gname=argv[7]
ncpus=int(sys.argv[8])
ens_mem=21
seaborn_map=True
numvar=2
#====================================================================
#=== read experiment names ===
# exlist="./experiment_list.nam"
# exlist="./Fig10-experiment_list.nam"
# # with open(exlist,"r") as exf:
# #     linesexf=exf.readlines()
# # indirs=[]
# # labels=[]
# # for lineexf in linesexf:
# #     lineexf = re.split(":",lineexf)
# #     lineexf = list(filter(None, lineexf))
# #     labels.append(lineexf[0])
# #     indirs.append(lineexf[1].strip())
# # #===
# # lexp=indirs
#=============================
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
with open(gname,"r") as f:
    lineGRDCs=f.readlines()
stations=[]
annotate_grdc=[]
text_label={}
# for exp in lexp:
#     annotate_grdc[exp]=[]
for i,lineGRDC in enumerate(lineGRDCs):
    lineGRDC = re.split(":",lineGRDC)
    lineGRDC = list(filter(None, lineGRDC))
    stations.append(lineGRDC[0].split()[0])
    annotate_grdc.append(lineGRDC[0].split()[0])
    #---------
    text_label[lineGRDC[0].split()[0]]=string.ascii_lowercase[i+2]
print stations
#======================================================
#files for making river network
tmp0=re.split("-",figname)[0]+".txt"

#  figure titles
titles=["Discharge", "WSE"]
# titles=["Uncalibrated model, Discharge", "Uncalibrated model, WSE"
        # , "Calibrated model, Discharge", "Calibrated model, WSE"]

land="#C0C0C0"
water="#FFFFFF"

londiff=(east-west)*4
latdiff=(north-south)*4

npix=(90-north)*4
spix=(90-south)*4
wpix=(180+west)*4
epix=(180+east)*4

#cmap=make_colormap(colors_list)
# cmap=mbar.colormap("H01")
# cmap=cm.get_cmap("bwr")
# cmap=cm.get_cmap("Spectral_r")
# cmap=cm.get_cmap("PiYG")
cmapQ=cm.get_cmap("PRGn")

# cmap.set_under("w",alpha=0)
# cmapL=cmap #cm.get_cmap("rainbow_r")
vmin=-0.50
vmax=0.50
norm=Normalize(vmin=vmin,vmax=vmax)
#------
# river width
sup=3 #2
w=0.002 #0.02
alpha=1
width=0.5
#------
west0=west
east0=east
south0=south-2
north0=north+2
#------
print (west,east,south,north)
resol=1
# figure in A4 size
va_margin= 0.0#1.38#inch 
ho_margin= 0.0#1.18#inch
numexp=1 #len(lexp)
#
exp=expname
hgt=(11.69 - 2*va_margin)*(6.0/10.0)
wdt=(8.27 - 2*ho_margin)
fig=plt.figure(figsize=(wdt,hgt),constrained_layout=True)
# G = gridspec.GridSpec(numexp,numvar)
G = gridspec.GridSpec(nrows=4,ncols=6)
# G.update(wspace=0.0, hspace=0.0) # set the spacing between axes. 
#=========================
# discharge correlation
# cmapQ=cm.get_cmap("PRGn")
# cmapQ=mbar.colormap("H07")
#colorlistL=['#7c0112','#b4791b','#f9c25c']
#colorlistR=['#81deae','#0288a6','#1f29a2']

# colorlistL=['#332288','#117733','#44aa99','#2fef10']
# colorlistR=['#fd411e','#cc6677','#aa4499','#882255']

# colorlistL=['#875f42','#882255','#aa4499','#cc6677','#fd411e']
# colorlistR=['#2fef10','#44aa99','#117733','#014600','#748500']

colorlistL=['#882255','#aa4499','#cc6677','#fd411e']
colorlistR=['#2fef10','#44aa99','#117733','#014600']
cmapQ=mbar.multiple_div_colormap(colorlistL,colorlistR)
vmin=-0.40
vmax=0.40
# norm=Normalize(vmin=vmin,vmax=vmax)
norm=BoundaryNorm(np.arange(vmin,vmax+0.1,0.1),cmapQ.N,clip=True)
ax1=fig.add_subplot(G[0:3,0:3])
m = Basemap(projection='cyl',llcrnrlat=south0,urcrnrlat=north0,llcrnrlon=west0,urcrnrlon=east0, lat_ts=0,resolution='c',ax=ax1)
# m.drawcoastlines( linewidth=0.1, color='k' )
# m.fillcontinents(color=land,lake_color=water,zorder=99)
# m.drawmapboundary(fill_color=water,zorder=100)
imQ=plt.scatter([],[],c=[],cmap=cmapQ,s=0.1,vmin=vmin,vmax=vmax,norm=norm,zorder=101)
imQ.set_visible(False)
# m.drawparallels(np.arange(south,north+0.1,5), labels = [1,0,0,0], fontsize=10,linewidth=0,zorder=102)
# m.drawmeridians(np.arange(west,east+0.1,5), labels = [0,0,0,1], fontsize=10,linewidth=0,zorder=102)
# Amazon basin boundry
m.readshapefile('./etc/Amazon_basin/Amazon_boundry', 'Amazon_boundry', drawbounds = False)
for info, shape in zip(m.Amazon_boundry, m.Amazon_boundry):
# if info['nombre'] == 'Amazon_boundry':
    x, y = zip(*shape) 
    m.plot(x, y, marker=None,color='grey',linewidth=1.0)
# plot river network
box="%f %f %f %f"%(west,east,north,south) 
os.system("./bin/txt_vector "+box+" "+CaMa_dir+" "+mapname+" > "+tmp0) 
#map(vec_par,np.arange(1,10+1,1))
map(vec_par,np.arange(2,10+1,1))
# map(vec_par,np.arange(7,10+1,1))
texts=[]
#--
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
    # print (num, basin ,stream)
    org=grdc.grdc_dis(num,syear,eyear)
    org=np.array(org)
    asm, opn = read_dis(exp,num)
    CORasm=correlation(asm,org)
    CORopn=correlation(opn,org)
    DCOR=CORasm-CORopn
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
    if np.sum((org!=-9999.0)*1.0) < 365*1:
        # print ("no obs: ",stream, np.sum((org!=-9999.0)*1.0), "< 365 days")
        continue
    #lon=lon0+ix*gsize
    #lat=lat0-iy*gsize
    lon=lonlat[0,iy1-1,ix1-1]
    lat=lonlat[1,iy1-1,ix1-1]
    c=cmapQ(norm(DCOR))
    ax1.scatter(lon,lat,s=15,marker="o",edgecolors="k",linewidth=0.5, facecolors=c,zorder=106)
    if num in annotate_grdc:
        texts.append(ax1.text(lon,lat,text_label[num], weight='bold',fontsize=10 ,zorder=110))
#===
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax1.spines['left'].set_visible(False)
# if i==0:
#     ax[2*i].text(0.5,1.05,"$Q$",ha="center",va="center",transform=ax[2*i].transAxes,fontsize=10)
ax1.text(-0.05,0.90,"%s) %s"%(string.ascii_lowercase[0],titles[0]),ha="left",va="center",transform=ax1.transAxes,fontsize=10)
#---
adjust_text(texts,ax=plt.gca(), expand_text=(1.0, 1.0), arrowprops=dict(arrowstyle='-', color='k', lw=0.5), zorder=110)
#======================
# wse rRMSE
#cmapW=mbar.colormap("H03")
# colorlistL=['#411900','#a90308','#ff1e02','#ff8000','#f9bc08']
# colorlistR=['#a2cffe','#4e518b','#1d5dec','#bf00fe','#720058']

colorlistL=['#a90308','#ff1e02','#ff8000','#f9bc08']
colorlistR=['#a2cffe','#0804f9','#6140ef','#bf00fe']
cmapW=mbar.multiple_div_colormap(colorlistL,colorlistR)
vmin=-0.8
vmax=0.8
# norm=Normalize(vmin=vmin,vmax=vmax)
norm=BoundaryNorm(np.arange(vmin,vmax+0.1,0.2),cmapW.N,clip=True)
ax2=fig.add_subplot(G[0:3,3::])
m = Basemap(projection='cyl',llcrnrlat=south0,urcrnrlat=north0,llcrnrlon=west0,urcrnrlon=east0, lat_ts=0,resolution='c',ax=ax2)
# m.drawcoastlines( linewidth=0.1, color='k' )
# m.fillcontinents(color=land,lake_color=water,zorder=99)
# m.drawmapboundary(fill_color=water,zorder=100)
imW=plt.scatter([],[],c=[],cmap=cmapW,s=0.1,vmin=vmin,vmax=vmax,norm=norm,zorder=101)
imW.set_visible(False)
# m.drawparallels(np.arange(south,north+0.1,5), labels = [1,0,0,0], fontsize=10,linewidth=0,zorder=102)
# m.drawmeridians(np.arange(west,east+0.1,5), labels = [0,0,0,1], fontsize=10,linewidth=0,zorder=102)
# Amzoon basin boundry
m.readshapefile('./etc/Amazon_basin/Amazon_boundry', 'Amazon_boundry', drawbounds = False)
for info, shape in zip(m.Amazon_boundry, m.Amazon_boundry):
# if info['nombre'] == 'Amazon_boundry':
    x, y = zip(*shape) 
    m.plot(x, y, marker=None,color='grey',linewidth=1.0)
# draw river network
box="%f %f %f %f"%(west,east,north,south) 
os.system("./bin/txt_vector "+box+" "+CaMa_dir+" "+mapname+" > corrtmp1.txt") 
#map(vec_par,np.arange(1,10+1,1))
map(vec_par,np.arange(2,10+1,1))
# map(vec_par,np.arange(7,10+1,1))
#--
obslist="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+mapname+"_QC0_simulation.txt"
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
    ix      = int(line[4])
    iy      = int(line[5])
    ele     = float(line[6])
    ele_dif = float(line[7])
    EGM08   = float(line[8])
    EGM96   = float(line[9])
    sat     = line[10].strip()
    org=hweb.HydroWeb_continous_WSE(station,syear,1,1,eyear,12,31,EGM08,EGM96)
    org=np.array(org)#+np.array(EGM08)-np.array(EGM96)
    asm, opn = read_wse(exp,station,"simulation")
    RMSEasm=RMSE(asm,org)
    RMSEopn=RMSE(opn,org)
    rRMSE=(RMSEasm-RMSEopn)/(RMSEopn+1.0e-20) 
    #-------------------------
    if rivermap[iy-1,ix-1] !=1.0:
        continue
    #--------------
    #lon=lon0+ix*gsize
    #lat=lat0-iy*gsize
    lon=lonlat[0,iy-1,ix-1]
    lat=lonlat[1,iy-1,ix-1]
    c=cmapW(norm(rRMSE))
    ax2.scatter(lon,lat,s=15,marker="s",edgecolors="k",linewidth=0.5, facecolors=c,zorder=106)
#===============================================================================
obslist="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+mapname+"_QC0_validation.txt"
with open(obslist,"r") as f:
    lines=f.readlines()
for line in lines[1::]:
    line    = re.split(" ",line)
    line    = list(filter(None, line))
    #print line
    num     = line[0].strip()
    station = line[1].strip()
    lon     = float(line[2])
    lat     = float(line[3])
    ix      = int(line[4])
    iy      = int(line[5])
    ele     = float(line[6])
    ele_dif = float(line[7])
    EGM08   = float(line[8])
    EGM96   = float(line[9])
    sat     = line[10].strip()
    org=hweb.HydroWeb_continous_WSE(station,syear,1,1,eyear,12,31,EGM08,EGM96)
    org=np.array(org)#+np.array(EGM08)-np.array(EGM96)
    asm, opn = read_wse(exp,station,"validation")
    RMSEasm=RMSE(asm,org)
    RMSEopn=RMSE(opn,org)
    rRMSE=(RMSEasm-RMSEopn)/(RMSEopn+1.0e-20) 
    #-------------------------
    if rivermap[iy-1,ix-1] !=1.0:
        continue
    #--------------
    #lon=lon0+ix*gsize
    #lat=lat0-iy*gsize
    lon=lonlat[0,iy-1,ix-1]
    lat=lonlat[1,iy-1,ix-1]
    c=cmapW(norm(rRMSE))
    ax2.scatter(lon,lat,s=20,marker="o",edgecolors="k",linewidth=0.5, facecolors=c,zorder=106)
#===
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.spines['left'].set_visible(False)
# if i==0:
#     ax[2*i+1].text(0.5,1.05,"$WSE$",ha="center",va="center",transform=ax[2*i+1].transAxes,fontsize=10)
ax2.text(-0.05,0.90,"%s) %s"%(string.ascii_lowercase[1],titles[1]),ha="left",va="center",transform=ax2.transAxes,fontsize=10)
# # #--
# # # add color bar below chart
# # # dividerQ = make_axes_locatable(ax1)
# # # caxQ = dividerQ.new_vertical(size='2%') #, pad=0.6, pack_start = True)
# # l,b,w,h = ax1.get_position().bounds
# # print l,b,w,h
# # # caxQ=fig.add_axes([0.1,0.42,0.35,0.01])
# # caxQ=fig.add_axes([l,b+0.15*h,w,0.01])
# # cbarQ=plt.colorbar(imQ,cax=caxQ,ticks=[-0.4,-0.2,0.0,0.2,0.4],orientation='horizontal',extend="both")#np.arange(-0.5,0.5+0.001,0.25)
# # cbarQ.set_label("$\Delta$$r$") #$(r_{assim} - r_{open})$") correlation
# # cbarQ.ax.tick_params(labelsize=6.0)
# # # cbar.set_label("$r $ $correlation$") #$(r_{assim} - r_{open})$")
# # # caxW=fig.add_axes([0.55,0.05,0.35,0.01])
# # # cbarW=plt.colorbar(imW,cax=caxW,ticks=[-0.8,-0.4,0.0,0.4,0.8],orientation='horizontal',extend="both") #np.arange(-1.0,1.0+0.001,0.2)
# # # cbarW.set_label("$rRMSE$") #$(r_{assim} - r_{open})$")
# # # cbarW.ax.tick_params(labelsize=6.0)

# # # dividerQ = make_axes_locatable(ax2)
# # # caxW = dividerQ.new_vertical(size='5%', pad=0.6, pack_start = True)
# # l,b,w,h=ax2.get_position().bounds
# # print l,b,w,h
# # # caxW=fig.add_axes([0.55,0.42,0.35,0.01])
# # caxW=fig.add_axes([l,b+0.15*h,w,0.01])
# # cbarW=plt.colorbar(imW,cax=caxW,ticks=[-0.8,-0.4,0.0,0.4,0.8],orientation='horizontal',extend="both") #np.arange(-1.0,1.0+0.001,0.2)
# # cbarW.set_label("$rRMSE$") #$(r_{assim} - r_{open})$")
# # cbarW.ax.tick_params(labelsize=6.0)

# # # legend 
# # feature1=mlines.Line2D([], [], color='grey', marker='s',
# #                           markersize=5, label='Assimilation',linewidth=0.0)
# # feature2=mlines.Line2D([], [], color='grey', marker='o',
# #                           markersize=5, label='Validation',linewidth=0.0)
# # # plt.axis('off')
# # # legend=plt.legend(handles=[feature1,feature2],bbox_to_anchor=(0.8, 0.001, 0.15, 0.1), loc=3,
# # #            ncol=1,  borderaxespad=0.0)#,fancybox=False) mode="expand",
# # l,b,w,h=ax2.get_position().bounds
# # legend=plt.legend(handles=[feature1,feature2],bbox_to_anchor=(l+0.99*w,b+0.2*h), loc="lower right",
# #            bbox_transform=fig.transFigure, ncol=1,  borderaxespad=0.0, frameon=False)#
# # # legend.get_frame().set_alpha(None)
# # # legend.get_frame().set_facecolor((0, 0, 1, 0.1))
# # # legend.set_facecolor('white')
# # # legend.set_edgecolor('white')
#=============================
# #============
# stations=["3629001","3618000","3626000","3625000"]
# # read GRDC IDs
# # gname="./Figs4-GRDClist.txt"
# with open(gname,"r") as f:
#     lineGRDCs=f.readlines()
# stations=[]
# for lineGRDC in lineGRDCs:
#     lineGRDC= re.split(":",lineGRDC)
#     lineGRDC = list(filter(None, lineGRDC))
#     stations.append(lineGRDC[0].split()[0])
# print stations
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
#======================================================
start_dt=datetime.date(syear,1,1)
end_dt=datetime.date(eyear,12,31)
size=60

start=0
last=(end_dt-start_dt).days + 1
N=int(last)
nbdays=int(last)
#-------------------------------
lexp=len(stations)
# if lexp==6:
#     colors=['#1f77b4','#9467bd','#2ca02c','#d62728','xkcd:goldenrod','xkcd:crimson']
# elif lexp==4:
#     colors=['#1f77b4','#9467bd','#2ca02c','#d62728']
# elif lexp==2:
#     colors=["#4dc7ec","#ff8021"]
# else:
#     colors=["#ff8021","#4dc7ec",]
colors=["#ff8021","#4dc7ec"]
#---
for j in np.arange(3):
    num=nums[j]
    print i,j,exp,num
    experiment=expname
    station=num
    ua, la, uo, lo = read_ul_all(experiment,station)
    asm, opn = read_dis(experiment,station)
    org=grdc.grdc_dis(station,syear,eyear)
    org=np.array(org)
    # x=int(j%cols)
    # y=int(j/cols)
    # print x, y
    if j==0:
        ax = fig.add_subplot(G[3,0:2])
    elif j==1:
        ax = fig.add_subplot(G[3,2:4])
    else:
        ax = fig.add_subplot(G[3,4::])
    ax.plot(np.arange(start,last),ma.masked_less(org,0.0),label="GRDC",color="#34495e",linewidth=1.0,zorder=101) #,marker = "o",markevery=swt[point])
    mxx=[]
    linewidth=0.3 #1.0/float(len(labels))
    ax.plot(np.arange(start,last),asm,label=exp,color="#ff8021",linewidth=linewidth,alpha=1,zorder=105)
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
        ax.set_ylabel('discharge ($m^3s^{-1}$)', color='k',fontsize=10)
        # ax.yaxis.label.set_size(10)
    ax.set_xlim(xmin=0,xmax=last+1)
    ax.tick_params('y', colors='k')
    # scentific notaion
    ax.ticklabel_format(style="sci",axis="y",scilimits=(0,0))
    ax.yaxis.major.formatter._useMathText=True
    # print (ax.get_ylim())
    # ymax=ax.get_ylim()
    # int(math.log(ymax))
    # ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    # ax.yaxis.set_major_formatter(FormatStrFormatter('%1.0e'))
    ax.yaxis.offsetText.set_fontsize(6) 
    ax.yaxis.get_offset_text().set_x(-0.10)
    ax.yaxis.get_offset_text().set_y(1.0)
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
    # ax.set_title("%s) "%(string.ascii_lowercase[j+2])+names[j]+", "+basins[j],ha="left",va="center",fontsize=10)
    # ax.text(0.0,1.15,"%s) %s%s"%(string.ascii_lowercase[j+2],basins[j][0],basins[j][1::].lower()),ha="left",va="center",transform=ax.transAxes,fontsize=10)
    # ax.text(0.0,1.15,"%s) GRDC ID: %s, lon: %4.1f, lat: %3.1f"%(string.ascii_lowercase[j+2],num[j],lonlat[0,IY1[j],IX1[j]],lonlat[1,IY1[j],IX1[j]])
    #     ,ha="left",va="center",transform=ax.transAxes,fontsize=10)
    ax.text(-0.10,1.20,"%s) %s%s, lon: %4.1f$\degree$, lat: %3.1f$\degree$"%(string.ascii_lowercase[j+2],basins[j][0],basins[j][1::].lower(),lonlat[0,IY1[j],IX1[j]],lonlat[1,IY1[j],IX1[j]])
        ,ha="left",va="center",transform=ax.transAxes,fontsize=10)
    #====================
    # correlation
    rasm=correlation(asm,org)
    ropn=correlation(opn,org)
    deltar=rasm-ropn
    dr="$\Delta$$r$= %3.2f"%(deltar)
    ax.text(0.0,-0.25,dr,ha="left",va="center",transform=ax.transAxes,fontsize=8)
    #====================
    # NSEAI
    NSEasm=NS(asm,org)
    NSEopn=NS(opn,org)
    NSEAI=(NSEasm-NSEopn)/(1-NSEopn+1.0e-20)
    NAI="$NSEAI$= %3.2f"%(NSEAI)
    ax.text(0.30,-0.25,NAI,ha="left",va="center",transform=ax.transAxes,fontsize=8)
    #====================
    # rISS
    ISSasm=ISS(ua,la,org,0.05)
    ISSopn=ISS(uo,lo,org,0.05)
    rISS=(ISSasm-ISSopn)/(ISSopn)
    rS="$rISS$= %3.2f"%(rISS)
    ax.text(0.75,-0.25,rS,ha="left",va="center",transform=ax.transAxes,fontsize=8)
    #========
    # if j==0:
    # plt.gca().text(-0.15,1.08,"%s)"%(string.ascii_lowercase[j+2]),ha="left",va="center",transform=plt.gca().transAxes,fontsize=10)
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
    #===
    ll, bb, ww, hh = ax.get_position().bounds 
    print (ll, bb, ww, hh)
    ax.set_position([ll, bb+hh*0.22, ww*0.95, hh*0.90])
    print (ax.get_position().bounds)
# add tight layout before colobar and legends
# fig.tight_layout()
#================================================
# add color bar below chart
# for river discharge
l,b,w,h = ax1.get_position().bounds
print ("colobar 1: ",l,b+0.15*h,w,0.01)
caxQ=fig.add_axes([l,b+0.15*h,w,0.01])
cbarQ=plt.colorbar(imQ,cax=caxQ,ticks=[-0.4,-0.2,0.0,0.2,0.4],orientation='horizontal',extend="both")#np.arange(-0.5,0.5+0.001,0.25)
cbarQ.set_label("$\Delta$$r$") #$(r_{assim} - r_{open})$") correlation
cbarQ.ax.tick_params(labelsize=6.0)

# for WSE
l,b,w,h=ax2.get_position().bounds
print ("colobar 1: ",l,b+0.15*h,w,0.01)
caxW=fig.add_axes([l,b+0.15*h,w,0.01])
cbarW=plt.colorbar(imW,cax=caxW,ticks=[-0.8,-0.4,0.0,0.4,0.8],orientation='horizontal',extend="both") #np.arange(-1.0,1.0+0.001,0.2)
cbarW.set_label("$rRMSE$") #$(r_{assim} - r_{open})$")
cbarW.ax.tick_params(labelsize=6.0)

# legend for scatter
feature1=mlines.Line2D([], [], color='grey', marker='s',
                          markersize=5, label='Assimilation',linewidth=0.0)
feature2=mlines.Line2D([], [], color='grey', marker='o',
                          markersize=5, label='Validation',linewidth=0.0)
l,b,w,h=ax2.get_position().bounds
legend0=plt.legend(handles=[feature1,feature2],bbox_to_anchor=(l+0.99*w,b+0.2*h), loc="lower right",
           bbox_transform=fig.transFigure, ncol=1,  borderaxespad=0.0, frameon=False)#
# Add the legend manually to the Axes.
plt.gca().add_artist(legend0)
#*********
# legend for time series
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
legend1=plt.legend(handles=features,bbox_to_anchor=(0.6,0.04), loc="lower right",
           bbox_transform=fig.transFigure, ncol=len(labels0),  borderaxespad=0.0, frameon=False, prop={"size":10})
# Add the legend manually to the Axes.
plt.gca().add_artist(legend1)
#
colors=["#ff8021","#4dc7ec"]
patches=[]
labelS=["95% Conf. Bound (Assimilation)","95% Conf. Bound (Open-loop)"]
for point in np.arange(2):
    patches.append(mpatches.Patch(color=colors[point],label=labelS[point],alpha=0.2))
#--
legend2=plt.legend(handles=patches,bbox_to_anchor=(0.6,0.025), loc="lower left",
           bbox_transform=fig.transFigure, ncol=1,  borderaxespad=0.0, frameon=False, prop={"size":10})#
#======================
# plt.subplots_adjust(wspace=0.05, hspace=0.05)
# fig.tight_layout()
#plt.title(stitle)
plt.savefig("./figures/"+figname+".jpg",dpi=500,bbox_inches="tight", pad_inches=0.0)
plt.savefig("./figures/"+figname+".png",dpi=500,bbox_inches="tight", pad_inches=0.0)
plt.savefig("./figures/"+figname+".pdf",dpi=800,bbox_inches="tight", pad_inches=0.0)
os.system("rm -r "+re.split("-",figname)[0]+"*.txt")