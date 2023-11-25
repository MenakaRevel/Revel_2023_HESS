#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import datetime
from matplotlib.colors import LogNorm,Normalize,ListedColormap,BoundaryNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import sys
import os
import math
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
import seaborn as sns
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages

import read_grdc as grdc
import my_colorbar as mbar
import read_hydroweb as hweb
# from read_sfcelv import read_sfcelv, read_sfcelv_multi
from river_patch import upstream
from statistics import *
from get_DA_data import *
#====================================================================
def upgrid(ix,iy,nextxy,uparea,rivseq,ups=1):
    nextX=nextxy[0]
    nextY=nextxy[1]
    uXX = ix
    uYY = iy
    tmp_XX = ix 
    tmp_YY = iy
    for _ in np.arange(ups):
        uXX, uYY = upstream(tmp_XX,tmp_YY,nextX.T,nextY.T,uparea.T)
        tmp_XX, tmp_YY = uXX, uYY
        if uXX < 0:
            break
        if rivseq[uYY-1,uXX-1] == 1:
            break
    return uXX, uYY
#====================================================================
def dwgrid(ix,iy,nextxy,dws=1):
    dXX = ix
    dYY = iy
    tmp_XX = ix 
    tmp_YY = iy 
    for _ in np.arange(dws):
        dXX=nextxy[0,tmp_YY-1,tmp_XX-1]
        dYY=nextxy[1,tmp_YY-1,tmp_XX-1]
        tmp_XX, tmp_YY = dXX, dYY
        if nextxy[0,tmp_YY-1,tmp_XX-1] < 0:
            break
    return dXX, dYY
#====================================================================
def power_law(x, a, b):
    return a * np.power(x, b)
#====================================================================
def mk_fig(ldis,lwse,syear,eyear,metric,chtitle):
    """
    make figures
    1. H-Q curve
    2. WSE graph
    3. Q graph
    """
    colors=["k","#D81B60","#FFC107","#004D40","grey"]
    # # #==
    # # com_data=(lwse[0]!=-9999.0)*(ldis[0]>0.0)*1.0
    # # locw=np.where(lwse[0]!=-9999.0)[0]
    # # orgw=ma.masked_equal(lwse[0],-9999.0).compressed()
    # # wse=ma.masked_where(com_data==0.0,lwse[0]).compressed()
    # # # discharge observations
    # # # locd=np.where(disobs!=-9999.0)[0]
    # # locd=np.arange(len(ldis[0]))#np.where(disobs!=-9999.0)[0]
    # # orgd=ma.masked_less(ldis[0],-100.0)
    # # dis=ma.masked_where(com_data==0.0,ldis[0]).compressed()
    # Figure
    hgt=11.69
    wdt=8.27
    plt.clf()
    fig=plt.figure(figsize=(wdt, hgt))
    #plt.title(pname[point][0],fontsize=12)
    G   = gridspec.GridSpec(3,2)
    ax0 = fig.add_subplot(G[0,0])
    ax1 = fig.add_subplot(G[0,1])
    ax2 = fig.add_subplot(G[1,:])
    ax3 = fig.add_subplot(G[2,:])
    #----
    #ax0.set_title(re.split("_",wselist[ixiy][0])[1],fontsize=12)
    # chtitle=dislist[ixiy][ii].strip()+"-"+wselist[ixiy][0].strip()
    # ax0.set_title(chtitle,fontsize=12)
    ax0.text(1.0,1.2,chtitle,horizontalalignment='center',
        verticalalignment='center',transform=ax0.transAxes,fontsize=14)
    ax0.axis('off')
    tab=metric
    tab_2 = [['%.2f' % j for j in i] for i in tab]
    ax0.table(cellText=tab_2,rowLabels=["NSEAI (Q)","pBias (Q)","$r$Amp (Q)","rRMSE (WSE)","Bias (WSE)","$\Delta$Amp (WSE)"],colLabels=["DIR","ANO","NOM"],loc='center')
    #----
    # ax1.plot(dis,wse,color=colors[0],linestyle='None',linewidth=0,marker="o",fillstyle="none",markersize=5)
    # # plot WSE
    # ax2.plot(locw,orgw,color=colors[0],linestyle='None',linewidth=0,marker="o",fillstyle="none",markersize=5)
    # ax3.plot(locd,orgd,color=colors[0])
    # for assimilated/simulated
    # print ("len", len(lwse), len(ldis))
    for i in np.arange(0,len(ldis)):
        com_data=(lwse[i]!=-9999.0)*(ldis[i]>0.0)*1.0
        # WSE
        locw=np.where(lwse[i]!=-9999.0)[0]
        dwse=ma.masked_equal(lwse[i],-9999.0).compressed()
        wse=ma.masked_where(com_data==0.0,lwse[i]).compressed()
        # Q
        locd=np.arange(len(ldis[i]))#np.where(disobs!=-9999.0)[0]
        ddis=ma.masked_less(ldis[i],-100.0)
        dis=ma.masked_where(com_data==0.0,ldis[i]).compressed()
        #========
        # plot 
        #========
        # H-Q curve
        # print ("maximum values",max(dis), max(wse))
        # print ("minimum values",min(dis), min(wse))
        ax1.plot(dis,wse,color=colors[i],linestyle='None',linewidth=0,marker="o",fillstyle="none",markersize=0.5,zorder=99-i,alpha=0.1)
        # fit a power_law
        try:
            popt, pcov = curve_fit(power_law,dis,wse)
            index=np.linspace(min(dis),max(dis),1000)
            ax1.plot(index, power_law(index, *popt),color=colors[i],linewidth=1.0,zorder=115+i) #, label='power law')
            print ("power law: ", popt) #, power_law(dis, *popt))
        except:
            print ("no optimal fit found......")
        #========
        # WSE 
        if i==0:
            ax2.plot(locw,dwse,color=colors[i],linestyle='None',linewidth=0,marker="o",fillstyle="none",markersize=5)
        else:    
            ax2.plot(locw,dwse,color=colors[i],linewidth=0.3)
        #========
        # Q
        if i==0:
            ax3.plot(locd,ddis,color=colors[i],linewidth=1.0)
        else:
            ax3.plot(locd,ddis,color=colors[i],linewidth=0.3)
    ##########
    # print (np.min(ma.masked_less(np.array(lwse),0.0)))
    ymin=np.min(ma.masked_less(np.array(lwse),0.0))
    ymax=np.max(np.array(lwse))
    xmin=np.min(ma.masked_less(np.array(ldis),0.0))
    xmax=np.max(np.array(ldis))
    #-----------
    ax1.set_ylim(ymin=ymin,ymax=ymax)
    ax1.set_xlim(xmin=xmin,xmax=xmax)
    ax1.set_xlabel("discharge $(m^3/s)$", color='k',fontsize=10)
    ax1.set_ylabel("WSE $(m)$", color='k',fontsize=10)
    ##########
    N=len(lwse[0])
    ax2.set_xlim(xmin=0,xmax=N+1)
    xxlab=np.arange(syear,eyear+2,1)
    dt=int(math.ceil(((eyear-syear)+2)/1.0))
    xxlist=np.linspace(0,N,dt,endpoint=True)
    ax2.set_ylabel('WSE $(m)$', color='k',fontsize=10)
    ax2.set_xticks(xxlist)
    ax2.set_xticklabels([])
    # ax2.set_xticklabels(xxlab,fontsize=8)
    ##########
    N=len(ldis[0])
    ax3.set_xlim(xmin=0,xmax=N+1)
    xxlab=np.arange(syear,eyear+2,1)
    dt=int(math.ceil(((eyear-syear)+2)/1.0))
    xxlist=np.linspace(0,N,dt,endpoint=True)
    ax3.set_ylabel('discharge $(m^3/s)$', color='k',fontsize=10)
    ax3.set_xlabel('Years', color='k',fontsize=10)
    ax3.set_xticks(xxlist)
    ax3.set_xticklabels(xxlab,fontsize=8)
    ##########
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()
#====================================================================
# argv=sys.argv
# syear=int(argv[1])
# eyear=int(argv[2])
# namea=argv[3]
# nameb=argv[4]
# CaMa_dir=argv[5]
# mapname=argv[6]
# figname=argv[7]
# indirs=argv[8]
# numvar=2
syear=2009
eyear=2014
mapname="amz_06min"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"
lexp=["DIR_WSE_E2O_HWEB_001","ANO_WSE_E2O_HWEB_004","NOM_WSE_E2O_HWEB_004"]
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
#=============================
nextxy = CaMa_dir+"/map/"+mapname+"/nextxy.bin"
rivwth = CaMa_dir+"/map/"+mapname+"/rivwth.bin"
rivhgt = CaMa_dir+"/map/"+mapname+"/rivhgt.bin"
rivlen = CaMa_dir+"/map/"+mapname+"/rivlen.bin"
elevtn = CaMa_dir+"/map/"+mapname+"/elevtn.bin"
lonlat = CaMa_dir+"/map/"+mapname+"/lonlat.bin"
uparea = CaMa_dir+"/map/"+mapname+"/uparea.bin"
nxtdst = CaMa_dir+"/map/"+mapname+"/nxtdst.bin"
rivseq = CaMa_dir+"/map/"+mapname+"/rivseq.bin"
nextxy = np.fromfile(nextxy,np.int32).reshape(2,ny,nx)
rivwth = np.fromfile(rivwth,np.float32).reshape(ny,nx)
rivhgt = np.fromfile(rivhgt,np.float32).reshape(ny,nx)
rivlen = np.fromfile(rivlen,np.float32).reshape(ny,nx)
elevtn = np.fromfile(elevtn,np.float32).reshape(ny,nx)
lonlat = np.fromfile(lonlat,np.float32).reshape(2,ny,nx)
uparea = np.fromfile(uparea,np.float32).reshape(ny,nx)
nxtdst = np.fromfile(nxtdst,np.float32).reshape(ny,nx)
rivseq = np.fromfile(rivseq,np.int32).reshape(ny,nx)
#=============================
# mean discharge
meanQ = "/work/a06/menaka/Empirical_LocalPatch/CaMa_out/amz_06min_S14FD/outflw_mean_1958-2013.bin"
meanQ = np.fromfile(meanQ,np.float32).reshape(ny,nx)
#******************************
#==============================
# wse
#==============================
wsexx=[]
wseyy=[]
wselist={}
wsenum={}
lEGM08={}
lEGM96={}
#==============================
wsetxt="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+mapname+"_QC0.txt"
#==============================
with open(wsetxt,"r") as fwse:
    linewse=fwse.readlines()
for point,line in enumerate(linewse[1::]):
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
    #----------------------
    ixiy="%05d%05d"%(ix,iy)
    if ixiy in wselist.keys():
        wselist[ixiy].append(station)
        wsenum[ixiy].append(point)
    else:
        wselist[ixiy]=[station]
        wsenum[ixiy]=[point]
    #-----------------------
    lEGM08[station]=EGM08
    lEGM96[station]=EGM96
#******************************
#==============================
# discharge
#==============================
disxx=[]
disyy=[]
dislist={}
disnum={}
#==============================
distxt="/cluster/data6/menaka/HydroDA/dat/grdc_"+mapname+".txt"
with open(distxt,"r") as fdis:
    linedis=fdis.readlines()
for point,line in enumerate(linedis[1::]):
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
    #----------------------
    ixiy="%05d%05d"%(ix1,iy1)
    if ixiy in dislist.keys():
        dislist[ixiy].append(num)
        disnum[ixiy].append(point)
    else:
        dislist[ixiy]=[num]
        disnum[ixiy]=[point]
#==============================
pdfname="./pdfdoc/dis_wse_relation.pdf"
#==============================
thr=0.05 # threshold for discharge difference
print (len(dislist.keys()), len(wselist.keys()))
with PdfPages(pdfname) as pdf:
    for ixiy in wselist.keys():
        #print dislist[ixiy]
        wsekey=ixiy
        diskey=ixiy
        if ixiy not in dislist.keys():
            # find next river grid
            # print (ixiy)
            ix=int(ixiy[0:5])
            iy=int(ixiy[5::])
            # print (ix,iy)
            #----------------------
            # upstream
            for gds in np.arange(10):
                dx, dy = dwgrid(ix,iy,nextxy,dws=gds+1)
                Q_diff = abs(meanQ[iy-1,ix-1]-meanQ[dy-1,dx-1])/(meanQ[iy-1,ix-1]+1e-20)
                if Q_diff <= thr:
                    new_ixiy="%05d%05d"%(dx-1,dy-1)
                    if new_ixiy in dislist.keys():
                        print ("downstream GRDC found")
                        diskey=new_ixiy
                        break
                ux, uy = upgrid(ix,iy,nextxy,uparea,rivseq,ups=gds+1)
                Q_diff = abs(meanQ[iy-1,ix-1]-meanQ[uy-1,ux-1])/(meanQ[iy-1,ix-1]+1e-20)
                if Q_diff <= thr:
                    new_ixiy="%05d%05d"%(ux-1,uy-1)
                    if new_ixiy in dislist.keys():
                        print ("upstream GRDC found")
                        diskey=new_ixiy
                        break
            if diskey not in dislist.keys():
                continue
            # for _ in np.arange(4):
            #     ix1, iy1 = tmp_ix, tmp_iy
            #     new_ixiy="%05d%05d"%(ix1,iy1)
            #     if new_ixiy in dislist.keys():
            #         print ("downstream GRDC found")
            #         diskey=new_ixiy
            #         break
            #     ix1=nextxy[0,tmp_iy,tmp_ix]
            #     iy1=nextxy[1,tmp_iy,tmp_ix]
            #     new_ixiy="%05d%05d"%(ix1,iy1)
            #     if new_ixiy not in dislist.keys():
            #         ix=ix1
            #         iy=iy1
            #         ix1=nextxy[0,iy-1,ix-1]
            #         iy1=nextxy[1,iy-1,ix-1]
            #         new_ixiy="%05d%05d"%(ix1,iy1)
            #         if new_ixiy not in dislist.keys():
            #             continue
            #         else:
            #             print ("downstream GRDC found")
            #             diskey=new_ixiy
            #     else:
            #         print ("downstream GRDC found")
            #         diskey=new_ixiy
            # print ("no WSE data")
            # continue
        # print ("==============================================================")
        # print (dislist[ixiy], wselist[ixiy])
        # print dislist[ixiy][0], wselist[ixiy][0]
        # wse observations
        # for wse in wsenum[ixiy]:
        #     print wse
        # wsepoint=wsenum[ixiy][0]
        # dispoint=disnum[ixiy][0]
        # grdc_id=grdc.get_id(dislist[ixiy][0],mapname)
        #==================================
        iQ=0
        lQ=0.0
        dnum=len(dislist[diskey])
        for i in np.arange(dnum):
            GRDCID=dislist[diskey][i]
            disobs=grdc.grdc_dis(GRDCID,syear=syear,eyear=eyear)
            if np.sum((disobs>0.0)*1.0) > lQ:
                iQ=i
                lQ=np.sum((disobs>0.0)*1.0)
            # print ("Q",iQ, lQ)
        #==================================
        iW=0
        lW=0.0
        snum=len(wselist[wsekey])
        for i in np.arange(snum):
            VS=wselist[wsekey][i]
            sfcobs=hweb.HydroWeb_continous_WSE(VS,syear=syear,eyear=eyear,egm08=lEGM08[VS],egm96=lEGM96[VS])
            if np.sum((sfcobs!=-9999.0)*1.0) > lW:
                iW=i
                lW=np.sum((sfcobs>0.0)*1.0)
            # print ("WSE",iW, lW)
            # print np.sum((sfcobs!=-9999.0)*1.0), np.sum((disobs>0.0)*1.0)
        # print (len(sfcobs),len(disobs))
        # # com_data=(sfcobs!=-9999.0)*(disobs>0.0)*1.0
        # # if np.sum(np.nan_to_num(com_data)) == 0.0:
        # #     print ("==============================================================")
        # #     print ("no common data",dislist[ixiy], wselist[ixiy])
        # #     continue
        flag=1
        # print (sfcobs, disobs)
        # observed values
        Q_list=[disobs]
        W_list=[sfcobs]
        # read assimilated and open-loop
        orgd=disobs
        orgw=sfcobs
        lNSEAI=[]
        lRMSE=[]
        lbiasQ=[]
        lampdQ=[]
        lbiasW=[]
        lampdW=[]
        for exp in lexp:
            #===================
            # discharge
            asm, opn = read_dis(exp,GRDCID)
            #===============
            NSEasm=NS(asm,orgd)
            NSEopn=NS(opn,orgd)
            NSEAI=(NSEasm-NSEopn)/(1.0-NSEopn+1.0e-20)
            # print ("NSEAI", NSEAI)
            lNSEAI.append(NSEAI)
            lbiasQ.append(pBIAS(opn,orgd))
            lampdQ.append(rAmplitude(opn,orgd))
            #
            Q_list.append(asm)
            if exp == lexp[-1]:
                Q_list.append(opn)
            #===================
            # WSE
            try:
                asm, opn = read_wse(exp,VS,"simulation")
            except:
                try:
                    asm, opn = read_wse(exp,VS,"validation")
                except:
                    print ("no file", exp, VS)
                    flag=-9
                    # continue
            W_list.append(asm)
            if exp == lexp[-1]:
                W_list.append(opn)
            #===============
            RMSEasm=RMSE(asm,orgw)
            RMSEopn=RMSE(opn,orgw)
            rRMSE=(RMSEasm-RMSEopn)/(RMSEopn+1e-20)
            # print ("rRMSE", rRMSE)
            lRMSE.append(rRMSE)
            lbiasW.append(BIAS(opn,orgw))
            lampdW.append(dAmplitude(opn,orgw))
        if flag==-9:
            continue
        print ("==============================================================")
        print (dislist[diskey], wselist[wsekey])
        # make figure
        chtitle=grdc.get_station(int(dislist[diskey][iQ].strip()))+"("+dislist[diskey][iQ]+") | "+wselist[ixiy][iW].strip()
        # chtitle=(dislist[diskey][iQ].strip())+" | "+wselist[wsekey][iW].strip()
        # print (np.shape(Q_list))
        metric=[lNSEAI,lbiasQ,lampdQ,lRMSE,lbiasW,lampdW]
        mk_fig(Q_list,W_list,syear,eyear,metric,chtitle)
    # set the file's metadata via the PdfPages object:
    d = pdf.infodict()
    d['Title'] = 'Q - WSE comparison'
    d['Author'] = 'Menaka Revel'
    d['Subject'] = 'Q - WSE comparison'
    d['Keywords'] = 'Q - WSE comparison'
    d['CreationDate'] = datetime.datetime.today()
    d['ModDate'] = datetime.datetime.today()