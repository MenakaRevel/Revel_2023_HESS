#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import datetime
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
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
from scipy import stats

import read_grdc as grdc
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
experiments=[]
labels=[]
for lineexf in linesexf:
    lineexf = re.split(":",lineexf)
    lineexf = list(filter(None, lineexf))
    labels.append(lineexf[0])
    experiments.append(lineexf[1].strip())
#=============================
lexp=len(experiments)
# colors=['xkcd:aqua green','xkcd:pastel blue','xkcd:soft pink','xkcd:periwinkle','xkcd:peach','xkcd:brick']
print lexp, experiments
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
        #-------------------------
        if satcov[iy1-1,ix1-1] !=1.0:
            continue
        #--------------
        org=grdc.grdc_dis(num,syear,eyear)
        if np.sum((org!=-9999.0)*1.0) < 365:
            # print ("no obs: ",stream, np.sum((org!=-9999.0)*1.0), "< 365 days")
            continue
        point=point+1
        #==================================
        for year in np.arange(syear,eyear+1):
            asm,opn=read_data_txt(num,exp,year,year)
            org=grdc.grdc_dis(num,year,year)
            #-------------------------
            if rivermap[iy1-1,ix1-1] !=1.0:
                continue
            #--------------
            if np.sum((org!=-9999.0)*1.0) < 365:
                # print ("no obs: ",stream, np.sum((org!=-9999.0)*1.0), "< 365 days")
                continue
            org_max[i,year-syear,point]=np.max(org)
            org_min[i,year-syear,point]=np.min(org)
            max_index=np.argmax(org)
            min_index=np.argmin(org)
            max1=max(max_index-15,0)
            max2=min(max_index+15,len(org)-1)
            min1=max(min_index-15,0)
            min2=min(min_index+15,len(org)-1)
            opn_max[i,year-syear,point]=np.max(opn[max1:max2])
            opn_min[i,year-syear,point]=np.min(opn[min1:min2])
            asm_max[i,year-syear,point]=np.max(asm[max1:max2])
            asm_min[i,year-syear,point]=np.min(asm[min1:min2])
            # print np.max(asm), np.min(asm)
#====================================
#---------- making fig --------------
#====================================
# colorsA=["#3174a1","#3f913a","#bf353c"]
# colorsO=["#afccdc","#b5d294","#f4adae"]
colorsO=["xkcd:lime green","xkcd:salmon"]
colorsA=["xkcd:emerald","xkcd:burnt umber"]
va_margin= 0.0#1.38#inch 
ho_margin= 0.0#1.18#inch
hgt=(11.69 - 2*va_margin)*(1.0/3.0)
wdt=(8.27 - 2*ho_margin)*(1.0)
#fig=plt.figure(figsize=(8.27,11.69))
fig=plt.figure(figsize=(wdt,hgt))
#fig.suptitle("auto-correalated area")#"maximum autocorrelation length")
#G = gridspec.GridSpec(2,1)
G = gridspec.GridSpec(1,2)
#--scatter
ax1 = fig.add_subplot(G[0,0])
flag=(org_max!=-9999.0)*(asm_max!=opn_max)*1.0
# ax1.plot(org_max,asm_max,marker="o",linestyle="none",linewidth=0,color="#ff8021")
print "Annual Maximum"
print "%10s%10s%10s%10s%10s"%("Experiment", "r^2(assim)", "r^2(open)","slope(assim)","slope(open)")
for i in np.arange(lexp):
    # Assim
    # print ma.masked_less_equal(org_max[i],0.0).flatten(), colors[i]
    org_max_data=np.log10(ma.masked_less_equal(org_max[i],0.0).filled(1.0).flatten())
    asm_max_data=np.log10(ma.masked_where(org_max[i]<=0.0,asm_max[i]).filled(1.0).flatten())
    opn_max_data=np.log10(ma.masked_where(org_max[i]<=0.0,opn_max[i]).filled(1.0).flatten())
    ax1.plot(org_max_data,asm_max_data,marker="o",fillstyle="none",markersize=3.0,linestyle="none",linewidth=0,color=colorsA[i],zorder=112)
    #===============
    slope, intercept, r, p, se = stats.linregress(org_max_data, asm_max_data)
    datamax=max(np.max(org_max_data),np.max(asm_max_data))
    datamin=min(np.min(org_max_data),np.min(asm_max_data))
    ax1.plot(np.array([datamin,datamax]),intercept + slope*np.array([datamin,datamax]), color=colorsA[i], label='fitted line',linewidth=0.3,linestyle="-")
    r_sqr="$r^2 (assim)$: %3.2f"%(r)
    r_asm=np.corrcoef(org_max_data, asm_max_data)[0,1]
    s_asm=slope
    print r_asm, s_asm
    # ax1.text(0.8,0.2,r_sqr,ha="left",va="center",transform=ax1.transAxes,fontsize=6)

    # Open
    ax1.plot(org_max_data,opn_max_data,marker="s",fillstyle="none",markersize=3.0,linestyle="none",linewidth=0,color=colorsO[i],zorder=111)
    #===============
    slope, intercept, r, p, se = stats.linregress(org_max_data, opn_max_data)
    datamax=max(np.max(org_max_data),np.max(asm_max_data))
    datamin=min(np.min(org_max_data),np.min(asm_max_data))
    ax1.plot(np.array([datamin,datamax]),intercept + slope*np.array([datamin,datamax]), color=colorsO[i], label='fitted line',linewidth=0.3,linestyle="--")
    r_sqr="$r^2 (open)$: %3.2f"%(r)
    r_opn=np.corrcoef(org_max_data, opn_max_data)[0,1]
    s_opn=slope
    # ax1.text(0.8,0.1,r_sqr,ha="left",va="center",transform=ax1.transAxes,fontsize=6)
    #====
    print "%5s%7.2f%7.2f%7.2f%7.2f"%(labels[i],r_asm,r_opn,s_asm,s_opn)

ax1.plot([0.0,6.0],[0.0,6.0],linestyle="--",linewidth=1.2,color="grey",zorder=110)
# ax1.plot(ma.masked_equal(org_max*flag,0.0),ma.masked_equal(asm_max*flag,0.0),marker="o",markersize=1.0,linestyle="none",linewidth=0,color="#ff8021")
# ax1.plot(ma.masked_equal(org_max*flag,0.0),ma.masked_equal(opn_max*flag,0.0),marker="o",markersize=1.0,linestyle="none",linewidth=0,color="#4dc7ec")


ax1.set_ylabel('Simulated / Assimilated', color='k',fontsize=8)
ax1.set_xlabel('GRDC', color='k',fontsize=8)
ax1.tick_params(labelsize=5)
ax1.set_xticklabels([r"$10^{%d}$"%(i) for i in np.arange(0,6+1,1)])
ax1.set_yticklabels([r"$10^{%d}$"%(i) for i in np.arange(0,6+1,1)])
ax1.set_ylim(ymin=0.0,ymax=6.0) #math.log10(260000.0)) #max(np.max(org_max),np.max(asm_max)))
ax1.set_xlim(xmin=0.0,xmax=6.0) #math.log10(260000.0)) #max(np.max(org_max),np.max(asm_max)))
ax1.set_title("Annual Maximum", loc="center")

#==== write in table ====

ax2 = fig.add_subplot(G[0,1])
print "Annual Minimum"
print "%10s%10s%10s%10s%10s"%("Experiment", "r^2(assim)", "r^2(open)","slope(assim)","slope(open)")
for i in np.arange(lexp):
    # Assim
    # print ma.masked_less_equal(org_min[i],0.0).flatten(), colors[i]
    org_min_data=np.log10(ma.masked_less_equal(org_min[i],0.0).filled(1.0).flatten())
    asm_min_data=np.log10(ma.masked_where(org_min[i]<=0.0,asm_min[i]).filled(1.0).flatten())
    opn_min_data=np.log10(ma.masked_where(org_min[i]<=0.0,opn_min[i]).filled(1.0).flatten())
    ax2.plot(org_min_data,asm_min_data,marker="o",fillstyle="none",markersize=3.0,linestyle="none",linewidth=0,color=colorsA[i],zorder=112)
    #===============
    slope, intercept, r, p, se = stats.linregress(org_min_data, asm_min_data)
    datamax=max(np.max(org_min_data),np.max(asm_min_data))
    datamin=min(np.min(org_min_data),np.min(asm_min_data))
    ax2.plot(np.array([datamin,datamax]),intercept + slope*np.array([datamin,datamax]), color=colorsA[i], label='fitted line',linewidth=0.3,linestyle="-")
    r_sqr="$r^2 (assim)$: %3.2f"%(r)
    r_asm=np.corrcoef(org_min_data, asm_min_data)[0,1]
    s_asm=slope
    # ax2.text(0.8,0.2,r_sqr,ha="left",va="center",transform=ax2.transAxes,fontsize=6)
        
    # Open
    ax2.plot(org_min_data,opn_min_data,marker="s",fillstyle="none",markersize=3.0,linestyle="none",linewidth=0,color=colorsO[i],zorder=111)
    #===============
    slope, intercept, r, p, se = stats.linregress(org_min_data, opn_min_data)
    datamax=max(np.max(org_min_data),np.max(asm_min_data))
    datamin=min(np.min(org_min_data),np.min(asm_min_data))
    ax2.plot(np.array([datamin,datamax]),intercept + slope*np.array([datamin,datamax]), color=colorsO[i], label='fitted line',linewidth=0.3,linestyle="--")
    r_sqr="$r^2 (open)$: %3.2f"%(r)
    r_opn=np.corrcoef(org_min_data, opn_min_data)[0,1]
    s_opn=slope
    # ax2.text(0.8,0.1,r_sqr,ha="left",va="center",transform=ax2.transAxes,fontsize=6)

    #====
    print "%5s%7.2f%7.2f%7.2f%7.2f"%(labels[i],r_asm,r_opn,s_asm,s_opn)

ax2.plot([0.0,6.0],[0.0,6.0],linestyle="--",linewidth=1.2,color="grey",zorder=110)

# ax2.set_ylabel('Simulated/Assimilated', color='k',fontsize=5)
ax2.set_xlabel('GRDC', color='k',fontsize=8)
ax2.tick_params(labelsize=5)
ax2.set_xticklabels([r"$10^{%d}$"%(i) for i in np.arange(0,6+1,1)])
ax2.set_yticklabels([r"$10^{%d}$"%(i) for i in np.arange(0,6+1,1)])
ax2.set_ylim(ymin=0.0,ymax=6.0) #math.log10(130000.0)) #max(np.max(org_min),np.max(asm_min)))
ax2.set_xlim(xmin=0.0,xmax=6.0) #math.log10(130000.0)) #max(np.max(org_min),np.max(asm_min)))
ax2.set_title("Annual Minimum", loc="center")
# corr_val, r_squared =measures(org_min,asm_min)
# print ("min: ",r_squared)
# r_sqr="$r^2 (assim)$: %3.2f"%(r_squared)
# ax2.text(0.8,0.2,r_sqr,ha="left",va="center",transform=ax2.transAxes,fontsize=6)
# corr_val, r_squared =measures(org_min,opn_min)
# print ("min: ",r_squared)
# r_sqr="$r^2 (open)$: %3.2f"%(r_squared)
# ax2.text(0.8,0.1,r_sqr,ha="left",va="center",transform=ax2.transAxes,fontsize=6)

# print (ma.masked_greater(asm_min,0.0).compress)
# print (ma.masked_greater(opn_min.flatten(),0.0).compress)
# figname=expname+"_"+"dis_annual_max_min"
patches=[]
colors=["xkcd:lime green","xkcd:emerald","xkcd:salmon","xkcd:burnt umber"]
labels=[labels[0]+"(Open-loop)", labels[0]+"(Assimilated)",labels[1]+"(Open-loop)", labels[1]+"(Assimilated)"]
for i,exp in enumerate(labels):
    patches.append(mpatches.Patch(color=colors[i], label=labels[i]))
legend1=plt.legend(handles=patches,bbox_to_anchor=(0.1,-0.1,0.6,0.1), loc="lower left",
           bbox_transform=fig.transFigure, ncol=2,  borderaxespad=0.0, frameon=False)#
#=========================================
# legend 
feature1=mlines.Line2D([], [], color='grey', marker='o',
                          markersize=5, label='Assimilation',linewidth=0.0)
feature2=mlines.Line2D([], [], color='grey', marker='s',
                          markersize=5, label='Open-loop',linewidth=0.0)
# plt.axis('off')
# legend=plt.legend(handles=[feature1,feature2],bbox_to_anchor=(0.8, 0.001, 0.15, 0.1), loc=3,
#            ncol=1,  borderaxespad=0.0)#,fancybox=False) mode="expand",
legend2=plt.legend(handles=[feature1,feature2],bbox_to_anchor=(0.8,-0.1), loc="lower right",
           bbox_transform=fig.transFigure, ncol=1,  borderaxespad=0.0, frameon=False)#
# plt.gca().add_artist(legend2)
plt.gca().add_artist(legend1)
#--
print ("./figures/"+figname+".pdf")
# print ("./figures/"+figname+".png")
plt.savefig("./figures/"+figname+".jpg",dpi=800,bbox_inches="tight", pad_inches=0.01)
plt.savefig("./figures/"+figname+".pdf",dpi=800,bbox_inches="tight", pad_inches=0.01)
plt.savefig("./figures/"+figname+".png",dpi=800,bbox_inches="tight", pad_inches=0.01)