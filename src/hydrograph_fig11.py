#!/opt/local/bin/python
# -*- coding: utf-8 -*-
'''
making hydrographs for give GRDC ID for figure 11 in manuscript
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
from matplotlib.ticker import MaxNLocator,FormatStrFormatter
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
# import cal_stat as stat
#import plot_colors as pc
#========================================
#====  functions for making figures  ====
#========================================
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
argv=sys.argv
# syear=int(argv[1])
# eyear=int(argv[2])
# CaMa_dir=argv[3]
# mapname=argv[4]
# expname=argv[5]
# ncpus=int(sys.argv[6])
#====================================================================
# station=argv[1]
# expname=argv[2]
station="3620000"
expname="DIR_WSE_E2O_HWEB_001"
syear=2009
eyear=2009
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
mapname="amz_06min"
# expname="NOM_WSE_E2O_HWEB_006"
# expname="DIR_WSE_E2O_HWEB_001"
# expname="DIR_WSE_E2O_HWEB_002"
#====================================================================
os.system("mkdir -p ./pptimage")
start_dt=datetime.date(syear,1,1)
end_dt=datetime.date(eyear,12,31)
size=60

start=0
last=(end_dt-start_dt).days + 1
N=int(last)
#====================================================================
# obslist="/cluster/data6/menaka/HydroDA/dat/grdc_"+mapname+".txt"
# with open(obslist,"r") as f:
#     lines=f.readlines()
# pname=[]
# for line in lines[1::]:
#     line    = re.split(";",line)
#     line    = list(filter(None, line))
#     # print (line)
#     num     = line[0].strip()
#     basin   = line[1].strip()
#     stream  = line[2].strip()
#     ix1     = int(line[3])
#     iy1     = int(line[4])
#     ix2     = int(line[5])
#     iy2     = int(line[6])
#     staid   = int(num)
# pname.append(num)
#====================================================================
def make_fig(station,exp1,exp2,exp3,ax=None,ylab='discharge (m$^3$/s)',xlab='year'):
    ax=ax or plt.gca()
    # read data
    ua, la, uo, lo = read_ul_all(expname,station)
    asm1, opn1 = read_dis(exp1,station)
    asm2, opn2 = read_dis(exp2,station)
    asm3, opn3 = read_dis(exp3,station)
    org=grdc.grdc_dis(station,syear,eyear)
    org=np.array(org)
    #=============
    # make figure
    #=============
    labels=["GRDC","open-loop","Direct","Anomaly","Normalized"]
    # colors=["#34495e","grey","#004488","#ddaa33","#ba5566"]
    colors=["#34495e","grey","#D81B60","#FFC107","#004D40"]
    lines=[ax.plot(np.arange(start,last),ma.masked_less(org,0.0),label=labels[0],color=colors[0],linewidth=3.0,zorder=101)[0]] #,marker = "o",markevery=swt[point])
    # draw mean of ensembles
    lines.append(ax.plot(np.arange(start,last),opn1,label=labels[1],color=colors[1],linewidth=0.8,linestyle="--",alpha=1,zorder=104)[0])
    lines.append(ax.plot(np.arange(start,last),asm1,label=labels[2],color=colors[2],linewidth=0.3,alpha=1,zorder=106)[0])
    lines.append(ax.plot(np.arange(start,last),asm2,label=labels[3],color=colors[3],linewidth=0.3,alpha=1,zorder=106)[0])
    lines.append(ax.plot(np.arange(start,last),asm3,label=labels[4],color=colors[4],linewidth=0.3,alpha=1,zorder=106)[0])
    # ax.fill_between(np.arange(start,last),lo,uo,color="#4dc7ec",alpha=0.2,zorder=102)
    # ax.fill_between(np.arange(start,last),la,ua,color="#ff8021",alpha=0.2,zorder=103)
    # print ua
    # print asm
    # print la
    #    plt.ylim(ymin=)
    #=======================================
    #xxlist=np.linspace(0,N,(eyear-syear)+1)
    #xlab=np.arange(syear,eyear+1,1)
    #xxlab=[calendar.month_name[i][:3] for i in range(1,13)]
    if eyear-syear <= 1:
        dtt=1
        dt=1 #int(math.ceil(((eyear-syear)+1)/dtt))
        # xxlist=np.linspace(0,N,14,endpoint=True)
        xxlist=np.cumsum(np.concatenate(([0],[calendar.monthrange(syear, i)[1] for i in range(1,13)])))
        xxlab=[calendar.month_name[i][:3] for i in range(1,13)]
    else:
        dtt=1 
        dt=(eyear-syear)+1
        xxlist=np.linspace(0,N,dt,endpoint=True)
        xxlab=np.arange(syear,eyear+1,dtt)
    #=======================================
    # Make the y-axis label, ticks and tick labels match the line color.
    if ylab:
        ax.set_ylabel(ylab, color='k',fontsize=10)
    #====
    if xlab:
        ax.set_xlabel(xlab, color='k',fontsize=10)
        ax.set_xticks(xxlist)
        ax.set_xticklabels(xxlab,fontsize=6,rotation=45)
    else:
        ax.set_xticks(xxlist)
        ax.set_xticklabels([],fontsize=6)

    ax.set_xlim(xmin=0,xmax=last+1)
    ax.set_ylim(ymin=0,ymax=3.0e5)
    ax.tick_params('y', colors='k')
    # scentific notaion
    ax.yaxis.set_major_locator(MaxNLocator(nbins=1,integer=True))
    # ax.yaxis.set_major_formatter(FormatStrFormatter('%1.0e'))
    ax.ticklabel_format(style="sci",axis="y",scilimits=(0,0))
    ax.yaxis.major.formatter._useMathText=True 
    ax.yaxis.offsetText.set_fontsize(10) 
    # ax.yaxis.label.set_size(10)
    ax.yaxis.get_offset_text().set_x(-0.05)
    #=====
    return 0
#============
if __name__ == '__main__':
    # figname="fe11-hydrograph_bathy_runoff"
    figname="figs7-hydrograph_bathy_runoff"
    fig=plt.figure() #figsize=(wdt,hgt))
    #fig.suptitle("auto-correalated area")#"maximum autocorrelation length")
    #G = gridspec.GridSpec(2,1)
    G = gridspec.GridSpec(ncols=2,nrows=2)
    ax1 = fig.add_subplot(G[0,0])
    ax2 = fig.add_subplot(G[0,1])
    ax3 = fig.add_subplot(G[1,0])
    ax4 = fig.add_subplot(G[1,1])
    #making panels
    # station="3629001"
    station="3626000"
    make_fig(station,"DIR_WSE_ECMWF_HWEB_011","ANO_WSE_ECMWF_HWEB_011","NOM_WSE_ECMWF_HWEB_011",ax=ax1,xlab=False)
    make_fig(station,"DIR_WSE_ECMWF_HWEB_012","ANO_WSE_ECMWF_HWEB_012","NOM_WSE_ECMWF_HWEB_012",ax=ax2,ylab=False,xlab=False)
    make_fig(station,"DIR_WSE_ECMWF_HWEB_013","ANO_WSE_ECMWF_HWEB_013","NOM_WSE_ECMWF_HWEB_013",ax=ax3)
    make_fig(station,"DIR_WSE_ECMWF_HWEB_014","ANO_WSE_ECMWF_HWEB_014","NOM_WSE_ECMWF_HWEB_014",ax=ax4,ylab=False)
    #
    ax1.text(-0.20,0.50,"without runoff bias",rotation=90,ha="right",va="center",transform=ax1.transAxes,fontsize=10)
    ax1.text(0.5,1.15,"without bathymetry error",rotation=0,ha="center",va="bottom",transform=ax1.transAxes,fontsize=10)
    ax2.text(0.5,1.15,"with bathymetry error",rotation=0,ha="center",va="bottom",transform=ax2.transAxes,fontsize=10)
    ax3.text(-0.20,0.50,"with runoff bias",rotation=90,ha="right",va="center",transform=ax3.transAxes,fontsize=10)
    #==
    features=[]
    labels=["observation","open-loop","Direct","Anomaly","Normalized"]
    # colors=["#34495e","grey","#004488","#ddaa33","#ba5566"]
    colors=["#34495e","grey","#D81B60","#FFC107","#004D40"]
    for i in np.arange(len(labels)):
        label=labels[i]
        features.append(mlines.Line2D([], [], color=colors[i], marker=None,
                            label=labels[i],linewidth=2.5)) #edgecolors="k",
    legend=plt.legend(handles=features,bbox_to_anchor=(0.5,0.0), loc="upper center",
           bbox_transform=fig.transFigure, ncol=len(labels),  borderaxespad=0.0, frameon=False)#
    #==
    plt.savefig("./figures/"+figname+".pdf",dpi=800,bbox_inches="tight", pad_inches=0.01)
    plt.savefig("./figures/"+figname+".png",dpi=800,bbox_inches="tight", pad_inches=0.01)
    plt.savefig("./figures/"+figname+".jpg",dpi=800,bbox_inches="tight", pad_inches=0.01)