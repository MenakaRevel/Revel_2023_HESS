#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import datetime
from matplotlib.colors import LogNorm,BoundaryNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import sys
import os
import calendar
from multiprocessing import Pool
from multiprocessing import Process
from multiprocessing import sharedctypes
from numpy import ma
import re
import math
import string
from scipy import stats
import warnings;warnings.filterwarnings('ignore')

import read_grdc as grdc
import my_colorbar as mbar
import read_hydroweb as hweb
#========================================
DA_dir="/cluster/data6/menaka/HydroDA"
# CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"
# expname="NOM_WSE_E2O_HWEB_003"
# expname="ANO_WSE_E2O_HWEB_003"
expname="DIR_WSE_E2O_HWEB_002"
mapname="amz_06min"
ens_mem=49
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
#==========================================================
def NRMSE(s,o):
    """
    Normalised Root Mean Squre Error
    input:
        s: simulated
        o: observed
    output:
        NRMSE: Normalised Root Mean Squre Error
    """
    s,o = filter_nan(s,o)
    o=ma.masked_where(o==-9999.0,o).filled(0.0)
    s=ma.masked_where(o==-9999.0,s).filled(0.0)
    o=np.compress(o>0.0,o)
    s=np.compress(o>0.0,s)
    # return np.sqrt(np.mean((s-o)**2))
    # return np.sqrt(np.ma.mean(np.ma.masked_where(o<=0.0,(s-o)**2)))/np.mean(o)
    # return np.sqrt(np.ma.mean(np.ma.masked_where(o<=0.0,(s-o)**2)))/(np.amax(np.ma.masked_less(o,0.0))-np.amin(np.ma.masked_less(o,0.0)))
    return np.mean(s)/np.mean(o)
#========================================
def measures(y1,y2):
    corr_val = np.corrcoef(y1,y2)[0,1]
    r_squared = corr_val**2
    return corr_val, r_squared
#========================================
def mk_dir(sdir):
  try:
    os.makedirs(sdir)
  except:
    pass
#========================================
def read_data_txt(grdc_id,expname,syear,eyear,smon=1,emon=12,sday=1,eday=31):
    iid=grdc_id
    start_dt=datetime.date(syear,smon,sday)
    last_dt=datetime.date(eyear,emon,eday)
    #======================================
    start=0
    last=(last_dt-start_dt).days + 1
    #======================================
    fname="./txt/"+expname+"/outflow/"+iid+".txt"
    disassim={}
    disopen={}
    with open(fname,"r") as f:
        lines=f.readlines()
    for line in lines:
        line    = list(filter(None, re.split(" ",line)))
        year    = int(line[0][0:4])
        mon     = int(line[0][4:6])
        day     = int(line[0][6::])
        assim   = float(line[1])
        opens   = float(line[2].strip())
        if start_dt <= datetime.date(year,mon,day) and last_dt  >= datetime.date(year,mon,day):
            disassim[year,mon,day]=float(assim)
            disopen[year,mon,day]=float(opens)
        elif last_dt  < datetime.date(year,mon,day):
            break
    #==================================
    start=0
    last=(last_dt-start_dt).days + 1
    Qassim=[]
    Qopen=[]
    for day in np.arange(start,last):
        target_dt=start_dt+datetime.timedelta(days=day)
        if (target_dt.year,target_dt.month,target_dt.day) in disassim.keys():
            Qassim.append(disassim[target_dt.year,target_dt.month,target_dt.day])
            Qopen.append(disopen[target_dt.year,target_dt.month,target_dt.day])
        else:
            Qassim.append(-9999.0)
            Qopen.append(-9999.0)
    return Qassim, Qopen
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
        asm.append(float(line[1]))
        opn.append(float(line[2]))
    return np.array(asm), np.array(opn)
#====================================================================
def get_pnum(fname):
    # obslist="/cluster/data6/menaka/HydroDA/dat/grdc_"+mapname+".txt"
    with open(fname,"r") as f:
        lines=f.readlines()
    pnum=0
    for line in lines[1::]:
        line    = re.split(";",line)
        line    = list(filter(None, line))
        # print (line)
        num     = line[0].strip()
        ix1     = int(line[3])
        iy1     = int(line[4])
        #-------------------------
        if rivermap[iy1-1,ix1-1] !=1.0:
            continue
        #--------------
        org=grdc.grdc_dis(num,syear,eyear)
        if np.sum((org!=-9999.0)*1.0) < 365*1:
            continue
        pnum=pnum+1
    return pnum
#========================================
#====================================================================
argv=sys.argv
syear=int(argv[1])
eyear=int(argv[2])
CaMa_dir=argv[3]
mapname=argv[4]
figname=argv[5]
ncpus=int(sys.argv[6])
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
exlist="./Figs6-experiment_list.nam"
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
colors=['xkcd:aqua green','xkcd:pastel blue','xkcd:soft pink','xkcd:periwinkle','xkcd:peach','xkcd:brick']

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
# metric=[]
start_dt=datetime.date(syear,1,1)
end_dt=datetime.date(eyear,12,31)
size=60
start=0
last=(end_dt-start_dt).days + 1
N=int(last)
pnum=get_pnum("/cluster/data6/menaka/HydroDA/dat/grdc_"+mapname+".txt")
org_max=np.ones([lexp,N,pnum])*-9999.0
org_min=np.ones([lexp,N,pnum])*-9999.0
opn_max=np.ones([lexp,N,pnum])*-9999.0
opn_min=np.ones([lexp,N,pnum])*-9999.0
asm_max=np.ones([lexp,N,pnum])*-9999.0
asm_min=np.ones([lexp,N,pnum])*-9999.0
for i,exp in enumerate(experiments):
    print (exp)
    metric_frag=[]
    obslist="/cluster/data6/menaka/HydroDA/dat/grdc_"+mapname+".txt"
    with open(obslist,"r") as f:
        lines=f.readlines()
    point=-1
    IX1=[];IY1=[]
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
        if rivermap[iy1-1,ix1-1] !=1.0:
            continue
        # #-------------------------
        # if satcov[iy1-1,ix1-1] !=1.0:
        #     continue
        #--------------
        org=grdc.grdc_dis(num,syear,eyear)
        if np.sum((org!=-9999.0)*1.0) < 365:
            # print ("no obs: ",stream, np.sum((org!=-9999.0)*1.0), "< 365 days")
            continue
        point=point+1
        #--
        IX1.append(ix1)
        IY1.append(iy1)
        #==================================
        for year in np.arange(syear,eyear+1):
            asm,opn=read_data_txt(num,exp,year,year)
            org=grdc.grdc_dis(num,year,year)
            #-------------------------
            if rivermap[iy1-1,ix1-1] !=1.0:
                continue
            #--------------
            if np.sum((org!=-9999.0)*1.0) < 1:
                # print ("no obs: ",stream, np.sum((org!=-9999.0)*1.0), "< 365 days")
                continue
            asm_max[i,year-syear,point]=np.max(asm)
            asm_min[i,year-syear,point]=np.min(asm)
            max_index=np.argmax(asm)
            min_index=np.argmin(asm)
            #=======
            # org_max[i,year-syear,point]=np.max(org)
            # org_min[i,year-syear,point]=np.min(org)
            # max_index=np.argmax(org)
            # min_index=np.argmin(org)
            max1=max(max_index-15,0)
            max2=min(max_index+15,len(org)-1)
            min1=max(min_index-15,0)
            min2=min(min_index+15,len(org)-1)
            opn_max[i,year-syear,point]=np.max(opn[max1:max2])
            opn_min[i,year-syear,point]=np.min(opn[min1:min2])
            # asm_max[i,year-syear,point]=np.max(asm[max1:max2])
            # asm_min[i,year-syear,point]=np.min(asm[min1:min2])
            # maximum
            if np.sum((org[max1:max2]!=-9999.0)*1.0) < 1:
                org_max[i,year-syear,point]=-9999.0
            else:
                org_max[i,year-syear,point]=np.max(org[max1:max2])
            # minimum
            if np.sum((org[min1:min2]!=-9999.0)*1.0) < 1:
                org_min[i,year-syear,point]=-9999.0
            else:
                org_min[i,year-syear,point]=np.min(ma.masked_equal(org[min1:min2],-9999.0))
            # print np.max(asm), np.min(asm)
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

# cmap=cm.get_cmap("RdPu")
# cmap=cm.get_cmap("RdPu")
# cmap=mbar.colormap("H01")
cmap=mbar.colormap("H03")
vmin=0.0
vmax=2.0
# norm=Normalize(vmin=vmin,vmax=vmax)
norm=BoundaryNorm(np.arange(vmin,vmax+0.1,0.1),cmap.N)

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
va_margin= 0.0#1.38#inch 
ho_margin= 0.0#1.18#inch
hgt=(11.69 - 2*va_margin)*(1.0/3.0)
wdt=(8.27 - 2*ho_margin)*(2.0/3.0)
#fig=plt.figure(figsize=(8.27,11.69))
fig=plt.figure(figsize=(wdt,hgt))
#fig.suptitle("auto-correalated area")#"maximum autocorrelation length")
#G = gridspec.GridSpec(2,1)
G = gridspec.GridSpec(2,2)
#=========
opt0=["Annual Maximum","Annual Minimum"]
pnum=np.shape(org_max)[2]
for i,exp in enumerate(experiments):
    for j,opt in enumerate(["max","min"]):
        print i,j,exp,opt
        ax = fig.add_subplot(G[i,j])
        m = Basemap(projection='cyl',llcrnrlat=south0,urcrnrlat=north0,llcrnrlon=west0,urcrnrlon=east0, lat_ts=0,resolution='c',ax=ax)
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
        for point in np.arange(pnum):
            ix=IX1[point]
            iy=IY1[point]
            lon=lonlat[0,iy-1,ix-1]
            lat=lonlat[1,iy-1,ix-1]
            if opt=="max":
                # np.sum((org!=-9999.0)*1.0) < 365
                if np.sum((org_max[i,:,point]!=-9999.0)*1.0) < 1:
                    continue
                val=NRMSE(asm_max[i,:,point],org_max[i,:,point])
            else:
                if np.sum((org_min[i,:,point]!=-9999.0)*1.0) < 1:
                    continue
                val=NRMSE(asm_min[i,:,point],org_min[i,:,point])
            if np.isnan(val):
                continue
            # print val
            col=cmap(norm(val))
            ax.scatter(lon,lat,s=15,marker="o",edgecolors="k",linewidth=0.5, facecolors=col,zorder=106)
        if j==0:
            plt.gca().text(-0.05,1.05,"%s)"%(string.ascii_lowercase[i]),ha="left",va="center",transform=plt.gca().transAxes,fontsize=10)
        if i==0:
            plt.gca().text(0.5,1.05,"%s"%(opt0[j]),ha="center",va="center",transform=plt.gca().transAxes,fontsize=10)
        #======================
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
#=======================
im=plt.scatter([],[],c=[],cmap=cmap,s=0.1,vmin=vmin,vmax=vmax,norm=norm,zorder=101)
im.set_visible(False)
cax=fig.add_axes([0.20,0.05,0.60,0.01])
cbar=plt.colorbar(im,cax=cax,ticks=np.arange(vmin,vmax+0.1,0.2),orientation='horizontal',extend="max")#np.arange(-0.5,0.5+0.001,0.25)
cbar.set_label("$NRMSE$") #$(r_{assim} - r_{open})$") correlation
cbar.ax.tick_params(labelsize=6.0)
#======================
plt.subplots_adjust(wspace=0, hspace=0)
print ("./figures/"+figname+".pdf")
# print ("./figures/"+figname+".png")
plt.savefig("./figures/"+figname+".pdf",dpi=800,bbox_inches="tight", pad_inches=0.0)
plt.savefig("./figures/"+figname+".png",dpi=800,bbox_inches="tight", pad_inches=0.0)