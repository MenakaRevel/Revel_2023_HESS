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
exlist=argv[6]
gname=argv[7]
ncpus=int(sys.argv[8])
ens_mem=21
seaborn_map=True
numvar=2
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
#===
lexp=indirs
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
annotate_grdc={}
text_label={}
for exp in lexp:
    annotate_grdc[exp]=[]
for i,lineGRDC in enumerate(lineGRDCs):
    lineGRDC= re.split(":",lineGRDC)
    lineGRDC = list(filter(None, lineGRDC))
    stations.append(lineGRDC[0].split()[0])
    if i < 3:
        annotate_grdc[lexp[0]].append(lineGRDC[0].split()[0])
    else:
        annotate_grdc[lexp[1]].append(lineGRDC[0].split()[0])
    #---------
    text_label[lineGRDC[0].split()[0]]=string.ascii_lowercase[i]
print stations
#======================================================
#files for making river network
tmp0=re.split("-",figname)[0]+".txt"

#  figure titles
titles=["Uncalibrated model, Discharge", "Uncalibrated model, WSE"
        , "Calibrated model, Discharge", "Calibrated model, WSE"]

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
numexp=len(lexp)
hgt=(11.69 - 2*va_margin)*(1.0/3.0)*numexp
wdt=(8.27 - 2*ho_margin)
fig=plt.figure(figsize=(wdt,hgt))
G = gridspec.GridSpec(numexp,numvar)
G.update(wspace=0.0, hspace=0.0) # set the spacing between axes. 
#=========================
ax=[]
for i,exp in enumerate(lexp):
    # discharge correlation
    # cmapQ=cm.get_cmap("PRGn")
    cmapQ=mbar.colormap("H07")
    vmin=-0.40
    vmax=0.40
    # norm=Normalize(vmin=vmin,vmax=vmax)
    norm=BoundaryNorm(np.arange(vmin,vmax+0.1,0.1),cmapQ.N)
    ax.append(fig.add_subplot(G[i,0]))
    m = Basemap(projection='cyl',llcrnrlat=south0,urcrnrlat=north0,llcrnrlon=west0,urcrnrlon=east0, lat_ts=0,resolution='c',ax=ax[2*i])
    # m.drawcoastlines( linewidth=0.1, color='k' )
    # m.fillcontinents(color=land,lake_color=water,zorder=99)
    # m.drawmapboundary(fill_color=water,zorder=100)
    imQ=plt.scatter([],[],c=[],cmap=cmapQ,s=0.1,vmin=vmin,vmax=vmax,norm=norm,zorder=101)
    imQ.set_visible(False)
    # m.drawparallels(np.arange(south,north+0.1,5), labels = [1,0,0,0], fontsize=10,linewidth=0,zorder=102)
    # m.drawmeridians(np.arange(west,east+0.1,5), labels = [0,0,0,1], fontsize=10,linewidth=0,zorder=102)
    #--
    m.readshapefile('./etc/Amazon_basin/Amazon_boundry', 'Amazon_boundry', drawbounds = False)
    for info, shape in zip(m.Amazon_boundry, m.Amazon_boundry):
        # if info['nombre'] == 'Amazon_boundry':
        x, y = zip(*shape) 
        m.plot(x, y, marker=None,color='grey',linewidth=1.0)
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
        ax[2*i].scatter(lon,lat,s=15,marker="o",edgecolors="k",linewidth=0.5, facecolors=c,zorder=106)
        if num in annotate_grdc[exp]:
            texts.append(ax[2*i].text(lon,lat,text_label[num], weight='bold',fontsize=10 ,zorder=110))
        #===
        ax[2*i].spines['top'].set_visible(False)
        ax[2*i].spines['right'].set_visible(False)
        ax[2*i].spines['bottom'].set_visible(False)
        ax[2*i].spines['left'].set_visible(False)
        # if i==0:
        #     ax[2*i].text(0.5,1.05,"$Q$",ha="center",va="center",transform=ax[2*i].transAxes,fontsize=10)
        ax[2*i].text(-0.05,1.05,"%s) %s"%(string.ascii_lowercase[2*i],titles[2*i]),ha="left",va="center",transform=ax[2*i].transAxes,fontsize=10)
    #---
    adjust_text(texts,ax=plt.gca(), expand_text=(1.0, 1.0), arrowprops=dict(arrowstyle='-', color='k', lw=0.5), zorder=110)
    #======================
    # wse rRMSE
    cmapW=mbar.colormap("H03")
    vmin=-0.8
    vmax=0.8
    # norm=Normalize(vmin=vmin,vmax=vmax)
    norm=BoundaryNorm(np.arange(vmin,vmax+0.1,0.2),cmapW.N)
    ax.append(fig.add_subplot(G[i,1]))
    m = Basemap(projection='cyl',llcrnrlat=south0,urcrnrlat=north0,llcrnrlon=west0,urcrnrlon=east0, lat_ts=0,resolution='c',ax=ax[2*i+1])
    # m.drawcoastlines( linewidth=0.1, color='k' )
    # m.fillcontinents(color=land,lake_color=water,zorder=99)
    # m.drawmapboundary(fill_color=water,zorder=100)
    imW=plt.scatter([],[],c=[],cmap=cmapW,s=0.1,vmin=vmin,vmax=vmax,norm=norm,zorder=101)
    imW.set_visible(False)
    # m.drawparallels(np.arange(south,north+0.1,5), labels = [1,0,0,0], fontsize=10,linewidth=0,zorder=102)
    # m.drawmeridians(np.arange(west,east+0.1,5), labels = [0,0,0,1], fontsize=10,linewidth=0,zorder=102)
    #--
    m.readshapefile('./etc/Amazon_basin/Amazon_boundry', 'Amazon_boundry', drawbounds = False)
    for info, shape in zip(m.Amazon_boundry, m.Amazon_boundry):
        # if info['nombre'] == 'Amazon_boundry':
        x, y = zip(*shape) 
        m.plot(x, y, marker=None,color='grey',linewidth=1.0)
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
        ax[2*i+1].scatter(lon,lat,s=15,marker="s",edgecolors="k",linewidth=0.5, facecolors=c,zorder=106)
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
        ax[2*i+1].scatter(lon,lat,s=20,marker="o",edgecolors="k",linewidth=0.5, facecolors=c,zorder=106)
        #===
        ax[2*i+1].spines['top'].set_visible(False)
        ax[2*i+1].spines['right'].set_visible(False)
        ax[2*i+1].spines['bottom'].set_visible(False)
        ax[2*i+1].spines['left'].set_visible(False)
        # if i==0:
        #     ax[2*i+1].text(0.5,1.05,"$WSE$",ha="center",va="center",transform=ax[2*i+1].transAxes,fontsize=10)
        ax[2*i+1].text(-0.05,1.05,"%s) %s"%(string.ascii_lowercase[2*i+1],titles[2*i+1]),ha="left",va="center",transform=ax[2*i+1].transAxes,fontsize=10)
#--
# add color bar below chart
# dividerQ = make_axes_locatable(ax[-2])
# caxQ = dividerQ.new_vertical(size='5%', pad=0.6, pack_start = True)
caxQ=fig.add_axes([0.1,0.1,0.35,0.01])
cbarQ=plt.colorbar(imQ,cax=caxQ,ticks=[-0.4,-0.2,0.0,0.2,0.4],orientation='horizontal',extend="both")#np.arange(-0.5,0.5+0.001,0.25)
cbarQ.set_label("$\Delta$$r$") #$(r_{assim} - r_{open})$") correlation
cbarQ.ax.tick_params(labelsize=6.0)
# cbar.set_label("$r $ $correlation$") #$(r_{assim} - r_{open})$")
# caxW=fig.add_axes([0.55,0.05,0.35,0.01])
# cbarW=plt.colorbar(imW,cax=caxW,ticks=[-0.8,-0.4,0.0,0.4,0.8],orientation='horizontal',extend="both") #np.arange(-1.0,1.0+0.001,0.2)
# cbarW.set_label("$rRMSE$") #$(r_{assim} - r_{open})$")
# cbarW.ax.tick_params(labelsize=6.0)

caxW=fig.add_axes([0.55,0.1,0.35,0.01])
cbarW=plt.colorbar(imW,cax=caxW,ticks=[-0.8,-0.4,0.0,0.4,0.8],orientation='horizontal',extend="both") #np.arange(-1.0,1.0+0.001,0.2)
cbarW.set_label("$rRMSE$") #$(r_{assim} - r_{open})$")
cbarW.ax.tick_params(labelsize=6.0)

# legend 
feature1=mlines.Line2D([], [], color='grey', marker='s',
                          markersize=5, label='Assimilation',linewidth=0.0)
feature2=mlines.Line2D([], [], color='grey', marker='o',
                          markersize=5, label='Validation',linewidth=0.0)
# plt.axis('off')
# legend=plt.legend(handles=[feature1,feature2],bbox_to_anchor=(0.8, 0.001, 0.15, 0.1), loc=3,
#            ncol=1,  borderaxespad=0.0)#,fancybox=False) mode="expand",
legend=plt.legend(handles=[feature1,feature2],bbox_to_anchor=(0.9,0.13), loc="lower right",
           bbox_transform=fig.transFigure, ncol=1,  borderaxespad=0.0, frameon=False)#
# legend.get_frame().set_alpha(None)
# legend.get_frame().set_facecolor((0, 0, 1, 0.1))
# legend.set_facecolor('white')
# legend.set_edgecolor('white')
plt.subplots_adjust(wspace=0, hspace=0)
#plt.title(stitle)
plt.savefig("./figures/"+figname+".jpg",dpi=500,bbox_inches="tight", pad_inches=0.0)
plt.savefig("./figures/"+figname+".png",dpi=500,bbox_inches="tight", pad_inches=0.0)
plt.savefig("./figures/"+figname+".pdf",dpi=800,bbox_inches="tight", pad_inches=0.0)
os.system("rm -r "+re.split("-",figname)[0]+"*.txt")