#!/opt/local/bin/python
# -*- coding: utf-8 -*-

from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
# import datetime
from matplotlib.colors import LogNorm,Normalize,ListedColormap,BoundaryNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
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
import warnings
from adjustText import adjust_text
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression

from statistics import *
import read_grdc as grdc
import my_colorbar as mbar
import read_hydroweb as hweb
from river_patch import patchms
#====================================================================
#=== functions ===
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
def ret_metric(exp,metric,thr=0.0):
    i=1
    # GRDC
    obslist="/cluster/data6/menaka/HydroDA/dat/grdc_"+mapname+".txt"
    with open(obslist,"r") as f:
        lines=f.readlines()
    for line in lines[1::]:
        line    = re.split(";",line)
        line    = list(filter(None, line))
        # print (line)
        num     = line[0].strip()
        basin   = line[1].split(",")[0].strip()
        stream  = line[2].strip()
        ix1     = int(line[3])
        iy1     = int(line[4])
        ix2     = int(line[5])
        iy2     = int(line[6])
        staid   = int(num)
        # print (num, basin ,stream)
        #-------------------------
        # if basin not in basins:
        #     if np.sum((org!=-9999.0)*1.0) >= 365*1:
        #         nolist.append(basin)
        #     continue
        #-------------------------
        if uparea[iy1-1,ix1-1] < 1.0e9:
            continue
        #-------------------------
        if rivermap[iy1-1,ix1-1] !=1.0:
            continue
        #-------------------------
        if satcov[iy1-1,ix1-1] !=1.0:
            continue
        #--------------
        org=grdc.grdc_dis(num,syear,eyear)
        org=np.array(org)
        #--------------
        if np.sum((org!=-9999.0)*1.0) < 365*1:
            # print ("no obs: ",stream, np.sum((org!=-9999.0)*1.0), "< 365 days")
            continue
        #--------------
        asm, opn = read_dis(exp,num)
        if metric=="NSEAI":
            # NSEAI
            NSEasm=NS(asm,org)
            NSEopn=NS(opn,org)
            NAI=(NSEasm-NSEopn)/(1.0-NSEopn+1.0e-20)
            # if NAI < thr:
            #     print basin, stream, NAI, NAI<0.0
            if NAI < thr: # thr:
                #-----------------
                # print stream
                VSstation,flag,egm08,egm96,dist=near_VS(ix1-1,iy1-1,20,nextxy,uparea,nxtdst,nx,ny)
                # print stream, flag, dist
                if flag==0:
                    # print stream, flag
                    continue
                if dist > 100.0: #25:
                    continue
                # print stream, "-", VSstation, dist
                # print i,stream,ix1-1, iy1-1, ix2-1,iy2-1, VSstation,flag
                try:
                    asm1, opn1 = read_wse(exp,VSstation,"simulation")
                except:
                    asm1, opn1 = read_wse(exp,VSstation,"validation")
                #--------------
                org1=hweb.HydroWeb_continous_WSE(VSstation,syear,1,1,eyear,12,31,egm08,egm96)
                org1=np.array(org1)#+np.array(EGM08)-np.array(EGM96)
                rmsea=RMSE(asm1,org1)
                rmseo=RMSE(opn1,org1)
                biasa=BIAS(asm1,org1)
                biaso=BIAS(opn1,org1)
                delAa=dAmplitude(asm1,org1)
                delAo=dAmplitude(opn1,org1)
                dpeaka=dpeak_time(asm1,org1)
                dpeako=dpeak_time(opn1,org1)
                print basin,";", stream,";", VSstation,";", dist,";", NAI,";", rmseo,";", biaso,";", delAo,";", dpeako #np.mean(opn1),np.mean(np.ma.masked_less(org1,0.0)) #num, basin, VSstation
        if metric=="delta_r":
            # delat r
            CORasm=correlation(asm,org)
            CORopn=correlation(opn,org)
            DCOR=CORasm-CORopn
            if DCOR < thr:
                print num, basin, stream, VSstation, DCOR, CORasm, CORopn
        if metric=="rISS":
            # rISS
            ua, la, uo, lo = read_ul_all(exp,num)
            ISSasm=ISS(ua,la,org,0.05)
            ISSopn=ISS(uo,lo,org,0.05)
            rISS  =(ISSasm-ISSopn)/(ISSopn)
            if rISS > thr:
                print num, basin, stream, rISS, ISSasm, ISSopn
    return 0
#====================================================================
def patchMS(ix,iy,ngrids,nextxy,uparea,nxtdst,nx,ny):
    outdir='/cluster/data6/menaka/Empirical_LocalPatch'
    nextX=nextxy[0]
    nextY=nextxy[1]
    # print "**********",ix,iy,nx,ny #outdir,
    xlist,ylist,dlist=patchms(ix+1,iy+1,ngrids,nextX.T,nextY.T,uparea.T,nxtdst.T) #,outdir)
    xlist=xlist[xlist>0]
    ylist=ylist[ylist>0]
    dlist=dlist[xlist>0]
    # print xlist-1, ylist-1, dlist
    return xlist-1, ylist-1, dlist
#====================================================================
def near_VS(ix0,iy0,ngrids,nextxy,uparea,nxtdst,nx,ny):
    # print ix0,iy0
    xlist0,ylist0,dlist0=patchMS(ix0,iy0,ngrids,nextxy,uparea,nxtdst,nx,ny)
    # print xlist
    # print ylist
    # print dlist
    indexs=np.argsort(dlist0)
    xlist=xlist0[indexs]
    ylist=ylist0[indexs]
    dlist=dlist0[indexs]
    flag=0
    dist0=1.0e20
    station0="none"
    eg08=0.0
    eg96=0.0
    for iXX,iYY,iDD in zip(xlist,ylist,dlist):
        # print "***",iXX, iYY, iDD
        obslist="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+mapname+"_QC0_simulation.txt"
        with open(obslist,"r") as f:
            lines=f.readlines()
        # metric_frag=[]
        for line in lines[1::]:
            line    = re.split(" ",line)
            line    = list(filter(None, line))
            # print line
            num     = line[0].strip()
            station = line[1].strip()
            lon     = float(line[2])
            lat     = float(line[3])
            ix      = int(line[4]) - 1
            iy      = int(line[5]) - 1
            ele     = float(line[6])
            ele_dif = float(line[7])
            EGM08   = float(line[8])
            EGM96   = float(line[9])
            sat     = line[10].strip()
            #-------------------------
            if ix == iXX and iy == iYY:
                flag=1
                station0=station
                eg08=EGM08
                eg96=EGM96
                dist0=iDD
                # print "L390",iXX, iYY, station, iDD
                # break
        #=============================
        obslist="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+mapname+"_QC0_validation.txt"
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
            ix      = int(line[4]) - 1
            iy      = int(line[5]) - 1
            ele     = float(line[6])
            ele_dif = float(line[7])
            EGM08   = float(line[8])
            EGM96   = float(line[9])
            sat     = line[10].strip()
            #-------------------------
            if ix == iXX and iy == iYY:
                if flag==0:
                    flag=1
                    station0=station
                    eg08=EGM08
                    eg96=EGM96
                    dist0=iDD
                    # print "L419", iXX, iYY, station, iDD
                    # break
                elif flag==1 and iDD<dist0:
                    flag=1
                    station0=station
                    eg08=EGM08
                    eg96=EGM96
                    dist0=iDD
                    # print "L426",iXX, iYY, station, iDD
                    # break
        if flag==1:
            break
    return station0, flag, eg08, eg96, dist0
#===================
# CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"
#map="glb_15min"
# map="glb_06min"
mapname="amz_06min"
# lexp=["DIR_WSE_E2O_HWEB_001","ANO_WSE_E2O_HWEB_004","NOM_WSE_E2O_HWEB_004"]
# figname="fig11-all_relative_metrics"
#=============================
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
nxtdst = CaMa_dir+"/map/"+mapname+"/nxtdst.bin"
nextxy = np.fromfile(nextxy,np.int32).reshape(2,ny,nx)
rivwth = np.fromfile(rivwth,np.float32).reshape(ny,nx)
rivhgt = np.fromfile(rivhgt,np.float32).reshape(ny,nx)
rivlen = np.fromfile(rivlen,np.float32).reshape(ny,nx)
elevtn = np.fromfile(elevtn,np.float32).reshape(ny,nx)
lonlat = np.fromfile(lonlat,np.float32).reshape(2,ny,nx)
uparea = np.fromfile(uparea,np.float32).reshape(ny,nx)
nxtdst = np.fromfile(nxtdst,np.float32).reshape(ny,nx)
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
syear=2009
eyear=2014
if __name__ == "__main__":
    # expname="NOM_WSE_E2O_HWEB_004"
    # expname="DIR_WSE_E2O_HWEB_001"
    expname="DIR_WSE_E2O_HWEB_002"
    # expname="ANO_WSE_E2O_HWEB_004"
    print expname
    ret_metric(expname,"NSEAI",0.0)