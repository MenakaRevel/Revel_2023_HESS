#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import datetime
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
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
#========================================
mk_dir("./figures")
mk_dir("./figures/"+expname)
#----
fname=CaMa_dir+"/map/"+mapname+"/params.txt"
f=open(fname,"r")
lines=f.readlines()
f.close()
#-------
nx     = int(filter(None, re.split(" ",lines[0]))[0])
ny     = int(filter(None, re.split(" ",lines[1]))[0])
gsize  = float(filter(None, re.split(" ",lines[3]))[0])
print mapname, nx, ny
#----
syear,smonth,sdate=2009,1,1
eyear,emonth,edate=2015,1,1 #2015
#month=1
#date=1
start_dt=datetime.date(syear,smonth,sdate)
end_dt=datetime.date(eyear,emonth,edate)
size=60

start=0
last=(end_dt-start_dt).days
#last=365#int(argvs[1])
#if calendar.isleap(year):
#    last=366
#else:
#    last=365

#last=89
N=int(last)

green2="greenyellow"
green ="green"
#colors = pc.color_assimd()
staid=[]
pname=[]
xlist=[]
ylist=[]
river=[]
#--
#rivernames  = ["LENA","NIGER","CONGO","OB","MISSISSIPPI","MEKONG","AMAZON","MEKONG","IRRAWADDY","VOLGA", "NIGER","YUKON","DANUBE"] #,"INDUS"] #["AMAZONAS"]#["CONGO"]#
# rivernames  = ["AMAZON"]
# rivernames = grdc.grdc_river_name_v396(CaMa_dir,mapname)
# rivernames  = ["AMAZON","NEGRO, RIO","PURUS, RIO","MADEIRA, RIO","JURUA, RIO","TAPAJOS, RIO","XINGU, RIO","CURUA, RIO","JAPURA, RIO"]
# rivernames  = ["NEGRO, RIO"]
rivernames  = ["AMAZON","NEGRO, RIO","PURUS, RIO","MADEIRA, RIO","JURUA, RIO"
              ,"TAPAJOS, RIO","XINGU, RIO","CURUA, RIO","JAPURA, RIO","BRANCO, RIO"
              ,"JAVARI, RIO","IRIRI, RIO","JURUENA, RIO","ACRE, RIO","BENI, RIO"
              ,"MAMORE, RIO","GUAPORE, RIO","ARINOS, RIO","TROMBETAS, RIO"]
for rivername in rivernames:
  grdc_id,station_loc,x_list,y_list = grdc.get_grdc_loc_v396(rivername,CaMa_dir,mapname)
#   print rivername, grdc_id,station_loc
  river.append([rivername]*len(station_loc))
  staid.append(grdc_id)
  pname.append(station_loc)
  xlist.append(x_list)
  ylist.append(y_list)
#--
river=([flatten for inner in river for flatten in inner])
staid=([flatten for inner in staid for flatten in inner])
pname=([flatten for inner in pname for flatten in inner])
print len(pname), len(xlist)
xlist=([flatten for inner in xlist for flatten in inner])
ylist=([flatten for inner in ylist for flatten in inner])


pnum=len(pname)
NN=eyear-syear

# opn_max=np.ctypeslib.as_ctypes(np.zeros([NN,pnum],np.float32))
# opn_min=np.ctypeslib.as_ctypes(np.zeros([NN,pnum],np.float32))
# shared_array_opn_max  = sharedctypes.RawArray(opn_max._type_, opn_max)
# shared_array_opn_min  = sharedctypes.RawArray(opn_min._type_, opn_min)

# asm_max=np.ctypeslib.as_ctypes(np.zeros([NN,pnum],np.float32))
# asm_min=np.ctypeslib.as_ctypes(np.zeros([NN,pnum],np.float32))
# shared_array_asm_max  = sharedctypes.RawArray(asm_max._type_, asm_max)
# shared_array_asm_min  = sharedctypes.RawArray(asm_min._type_, asm_min)

# for parallel calcualtion
inputlist=[]
for year in np.arange(syear,eyear):
    inputlist.append(year)

def read_data(year):
    year=int(year)
    #--
    tmp_opn_max  = np.ctypeslib.as_array(shared_array_opn_max)
    tmp_opn_min  = np.ctypeslib.as_array(shared_array_opn_min)
    tmp_asm_max  = np.ctypeslib.as_array(shared_array_asm_max)
    tmp_asm_min  = np.ctypeslib.as_array(shared_array_asm_min)

    start_date=datetime.date(year,1,1)
    days=365
    if calendar.isleap(year):
        days=366
    tmp_opn=np.zeros([days,ens_mem,pnum],np.float32)
    tmp_asm=np.zeros([days,ens_mem,pnum],np.float32)
    for day in np.arange(0,days):
        target_dt=start_date+datetime.timedelta(days=day)
        yyyy='%04d' % (target_dt.year)
        mm='%02d' % (target_dt.month)
        dd='%02d' % (target_dt.day)
        for num in np.arange(1,ens_mem+1):
            numch="%03d"%(num)
            # corrpted
            fname=DA_dir+"/out/"+expname+"/assim_out/outflw/open/outflw"+yyyy+mm+dd+"_"+numch+".bin"
            opnfile=np.fromfile(fname,np.float32).reshape([ny,nx])
            # assimilated
            fname=DA_dir+"/out/"+expname+"/assim_out/outflw/assim/outflw"+yyyy+mm+dd+"_"+numch+".bin"
            asmfile=np.fromfile(fname,np.float32).reshape([ny,nx])
            for point in np.arange(pnum):
                ix1,iy1,ix2,iy2=grdc.get_grdc_station_v396(pname[point],CaMa_dir,mapname)
                if ix2 == -9999 or iy2 == -9999:
                    tmp_asm[day,num-1,point]=asmfile[iy1,ix1]
                    tmp_opn[day,num-1,point]=opnfile[iy1,ix1]
                else:
                    tmp_opn[day,num-1,point]=opnfile[iy1,ix1]+opnfile[iy2,ix2]
                    tmp_asm[day,num-1,point]=asmfile[iy1,ix1]+asmfile[iy2,ix2]
    #===
    tmp_opn_max[year-syear,:]=np.max(np.mean(tmp_opn,axis=1),axis=0)
    tmp_opn_min[year-syear,:]=np.min(np.mean(tmp_opn,axis=1),axis=0)
    tmp_asm_max[year-syear,:]=np.max(np.mean(tmp_asm,axis=1),axis=0)
    tmp_asm_min[year-syear,:]=np.min(np.mean(tmp_asm,axis=1),axis=0)
# #--------
# para_flag=0
# if para_flag==1:
#     p       = Pool(20)
#     res     = p.map(read_data, inputlist)
#     opn_max = np.ctypeslib.as_array(shared_array_opn_max)
#     opn_min = np.ctypeslib.as_array(shared_array_opn_min)
#     asm_max = np.ctypeslib.as_array(shared_array_asm_max)
#     asm_min = np.ctypeslib.as_array(shared_array_asm_min)
#     p.terminate()
# else:
#     for inpi in np.arange(eyear-syear):
#         res = read_data(inputlist[inpi])
#         opn_max = np.ctypeslib.as_array(shared_array_opn_max)
#         opn_min = np.ctypeslib.as_array(shared_array_opn_min)
#         asm_max = np.ctypeslib.as_array(shared_array_asm_max)
#         asm_min = np.ctypeslib.as_array(shared_array_asm_min)
# print pnum
# pnum=pnum-1
# print pnum
org_max=np.ones([NN,pnum])*-9999.0
org_min=np.ones([NN,pnum])*-9999.0
opn_max=np.ones([NN,pnum])*-9999.0
opn_min=np.ones([NN,pnum])*-9999.0
asm_max=np.ones([NN,pnum])*-9999.0
asm_min=np.ones([NN,pnum])*-9999.0
for point in np.arange(pnum):
    # if pname[point]=="SERRINHA":
    #     continue
    # print (pname[point])
    for year in np.arange(syear,eyear):
        asm,opn=read_data_txt(staid[point],expname,year,year)
        # asm,opn=read_dis(expname,staid[point])
        # opn_max[year-syear,point]=np.max(opn)
        # opn_min[year-syear,point]=np.min(opn)
        # asm_max[year-syear,point]=np.max(asm)
        # asm_min[year-syear,point]=np.min(asm)
        org=grdc.grdc_dis(staid[point],year,year)
        if np.sum((org!=-9999.0)*1.0) >= 365:
            # print (np.sum((org!=-9999.0)*1.0))
            org_max[year-syear,point]=np.max(org)
            org_min[year-syear,point]=np.min(org) #ma.masked_equal(org,-9999.0).filled(1.e20))
            opn_max[year-syear,point]=np.max(opn)
            opn_min[year-syear,point]=np.min(opn)
            asm_max[year-syear,point]=np.max(asm)
            asm_min[year-syear,point]=np.min(asm)
            # print (np.max(asm), np.min(asm), np.max(org), np.min(ma.masked_equal(org,-9999.0).filled(1.e20)))
        else:
            # print (np.sum((org!=-9999.0)*1.0))
            org_max[year-syear,point]=-9999.0
            org_min[year-syear,point]=-9999.0
            opn_max[year-syear,point]=-9999.0
            opn_min[year-syear,point]=-9999.0
            asm_max[year-syear,point]=-9999.0
            asm_min[year-syear,point]=-9999.0
        if pname[point]=="SERRINHA":
            # print ("remove: ",pname[point])
            org_max[year-syear,point]=-9999.0
            org_min[year-syear,point]=-9999.0
            opn_max[year-syear,point]=-9999.0
            opn_min[year-syear,point]=-9999.0
            asm_max[year-syear,point]=-9999.0
            asm_min[year-syear,point]=-9999.0
        # print (pname[point], year, org_max[year-syear,point],org_min[year-syear,point])
        # print (pname[point], year, org_max[year-syear,point], asm_max[year-syear,point], opn_max[year-syear,point])
# make figure
va_margin= 0.0#1.38#inch 
ho_margin= 0.0#1.18#inch
hgt=(11.69 - 2*va_margin)*(1.0/3.0)
wdt=(8.27 - 2*ho_margin)*(1.0)
#fig=plt.figure(figsize=(8.27,11.69))
fig=plt.figure(figsize=(wdt,hgt))
#fig.suptitle("auto-correalated area")#"maximum autocorrelation length")
#G = gridspec.GridSpec(2,1)
G = gridspec.GridSpec(1,2)
#--boxplot
#ax = fig.add_subplot(G[1,0])
ax1 = fig.add_subplot(G[0,0])
flag=(org_max!=-9999.0)*(asm_max!=opn_max)*1.0
# ax1.plot(org_max,asm_max,marker="o",linestyle="none",linewidth=0,color="#ff8021")
ax1.plot(ma.masked_less_equal(org_max,0.0).flatten(),ma.masked_where(org_max<=0.0,asm_max).flatten(),marker="o",fillstyle="none",markersize=5.0,linestyle="none",linewidth=0,color="#ff8021",zorder=112)
ax1.plot(ma.masked_less_equal(org_max,0.0).flatten(),ma.masked_where(org_max<=0.0,opn_max).flatten(),marker="o",fillstyle="none",markersize=5.0,linestyle="none",linewidth=0,color="#4dc7ec",zorder=111)
ax1.plot([0.0,260000.0],[0.0,260000.0],linestyle="--",linewidth=0.5,color="k",zorder=110)
# ax1.plot(ma.masked_equal(org_max*flag,0.0),ma.masked_equal(asm_max*flag,0.0),marker="o",markersize=1.0,linestyle="none",linewidth=0,color="#ff8021")
# ax1.plot(ma.masked_equal(org_max*flag,0.0),ma.masked_equal(opn_max*flag,0.0),marker="o",markersize=1.0,linestyle="none",linewidth=0,color="#4dc7ec")

# Assim
slope, intercept, r, p, se = stats.linregress(ma.masked_less_equal(org_max,0.0).flatten(), ma.masked_where(org_max<=0.0,asm_max).flatten())
datamax=max(np.max(ma.masked_less_equal(org_max,0.0).flatten()),np.max(ma.masked_where(org_max<=0.0,asm_max).flatten()))
datamin=min(np.min(ma.masked_less_equal(org_max,0.0).flatten()),np.min(ma.masked_where(org_max<=0.0,asm_max).flatten()))
ax1.plot(np.array([datamin,datamax]),intercept + slope*np.array([datamin,datamax]), color="#ff8021", label='fitted line')
r_sqr="$r^2 (assim)$: %3.2f"%(r)
# ax1.text(0.8,0.2,r_sqr,ha="left",va="center",transform=ax1.transAxes,fontsize=6)

# Open
slope, intercept, r, p, se = stats.linregress(ma.masked_less_equal(org_max,0.0).flatten(), ma.masked_where(org_max<=0.0,opn_max).flatten())
datamax=max(np.max(ma.masked_less_equal(org_max,0.0).flatten()),np.max(ma.masked_where(org_max<=0.0,asm_max).flatten()))
datamin=min(np.min(ma.masked_less_equal(org_max,0.0).flatten()),np.min(ma.masked_where(org_max<=0.0,asm_max).flatten()))
ax1.plot(np.array([datamin,datamax]),intercept + slope*np.array([datamin,datamax]), color="#4dc7ec", label='fitted line')
r_sqr="$r^2 (open)$: %3.2f"%(r)
# ax1.text(0.8,0.1,r_sqr,ha="left",va="center",transform=ax1.transAxes,fontsize=6)


ax1.set_ylabel('Simulated / Assimilated', color='k',fontsize=8)
ax1.set_xlabel('GRDC', color='k',fontsize=8)
ax1.tick_params(labelsize=5)
ax1.set_ylim(ymin=0.0,ymax=260000.0) #max(np.max(org_max),np.max(asm_max)))
ax1.set_xlim(xmin=0.0,xmax=260000.0) #max(np.max(org_max),np.max(asm_max)))
ax1.set_title("Annual Maximum", loc="center")

ax2 = fig.add_subplot(G[0,1])
ax2.plot(ma.masked_less_equal(org_min,0.0).flatten(),ma.masked_where(org_min<=0.0,asm_min).flatten(),marker="o",fillstyle="none",markersize=5.0,linestyle="none",linewidth=0,color="#ff8021",zorder=112)
ax2.plot(ma.masked_less_equal(org_min,0.0).flatten(),ma.masked_where(org_min<=0.0,opn_min).flatten(),marker="o",fillstyle="none",markersize=5.0,linestyle="none",linewidth=0,color="#4dc7ec",zorder=111)
ax2.plot([0.0,130000.0],[0.0,130000.0],linestyle="--",linewidth=0.5,color="k")

# Assim
slope, intercept, r, p, se = stats.linregress(ma.masked_less_equal(org_min,0.0).flatten(), ma.masked_where(org_min<=0.0,asm_min).flatten())
datamax=max(np.max(ma.masked_less_equal(org_min,0.0).flatten()),np.max(ma.masked_where(org_min<=0.0,asm_min).flatten()))
datamin=min(np.min(ma.masked_less_equal(org_min,0.0).flatten()),np.min(ma.masked_where(org_min<=0.0,asm_min).flatten()))
ax2.plot(np.array([datamin,datamax]),intercept + slope*np.array([datamin,datamax]), color="#ff8021", label='fitted line')
r_sqr="$r^2 (assim)$: %3.2f"%(r)
# ax2.text(0.8,0.2,r_sqr,ha="left",va="center",transform=ax2.transAxes,fontsize=6)

# Open
slope, intercept, r, p, se = stats.linregress(ma.masked_less_equal(org_min,0.0).flatten(), ma.masked_where(org_min<=0.0,opn_min).flatten())
datamax=max(np.max(ma.masked_less_equal(org_min,0.0).flatten()),np.max(ma.masked_where(org_min<=0.0,asm_min).flatten()))
datamin=min(np.min(ma.masked_less_equal(org_min,0.0).flatten()),np.min(ma.masked_where(org_min<=0.0,asm_min).flatten()))
ax2.plot(np.array([datamin,datamax]),intercept + slope*np.array([datamin,datamax]), color="#4dc7ec", label='fitted line')
r_sqr="$r^2 (open)$: %3.2f"%(r)
# ax2.text(0.8,0.1,r_sqr,ha="left",va="center",transform=ax2.transAxes,fontsize=6)


# ax2.set_ylabel('Simulated/Assimilated', color='k',fontsize=5)
ax2.set_xlabel('GRDC', color='k',fontsize=8)
ax2.tick_params(labelsize=5)
ax2.set_ylim(ymin=0.0,ymax=130000.0) #max(np.max(org_min),np.max(asm_min)))
ax2.set_xlim(xmin=0.0,xmax=130000.0) #max(np.max(org_min),np.max(asm_min)))
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
figname=expname+"_"+"dis_annual_max_min"
#--
print ("./figures/"+figname+".pdf")
plt.savefig("./figures/"+figname+".pdf",dpi=800,bbox_inches="tight", pad_inches=0.05)
print ("./figures/"+figname+".png")
plt.savefig("./figures/"+figname+".png",dpi=800,bbox_inches="tight", pad_inches=0.05)