#!/opt/local/bin/python
# -*- coding: utf-8 -*-
'''
making hydrographs for give GRDC ID with confidence intervals
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
eyear=2014
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
mapname="amz_06min"
# expname="NOM_WSE_E2O_HWEB_006"
# expname="DIR_WSE_E2O_HWEB_001"
# expname="DIR_WSE_E2O_HWEB_002"
#====================================================================

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
def make_fig(station,expname):
    # read data
    ua, la, uo, lo = read_ul_all(expname,station)
    asm, opn = read_dis(expname,station)
    org=grdc.grdc_dis(station,syear,eyear)
    org=np.array(org)
    #=============
    # make figure
    #=============
    plt.close()
    #labels=["GRDC","corrupted","assimilated"]
    labels=["GRDC","simulated","assimilated"]
    fig, ax1 = plt.subplots()
    lines=[ax1.plot(np.arange(start,last),ma.masked_less(org,0.0),label="GRDC",color="#34495e",linewidth=3.0,zorder=101)[0]] #,marker = "o",markevery=swt[point])
    # draw mean of ensembles
    lines.append(ax1.plot(np.arange(start,last),opn,label="corrupted",color="#4dc7ec",linewidth=1.0,alpha=1,zorder=104)[0])
    lines.append(ax1.plot(np.arange(start,last),asm,label="assimilated",color="#ff8021",linewidth=1.0,alpha=1,zorder=106)[0])
    ax1.fill_between(np.arange(start,last),lo,uo,color="#4dc7ec",alpha=0.2,zorder=102)
    ax1.fill_between(np.arange(start,last),la,ua,color="#ff8021",alpha=0.2,zorder=103)
    # print ua
    # print asm
    # print la
    #    plt.ylim(ymin=)
    # Make the y-axis label, ticks and tick labels match the line color.
    ax1.set_ylabel('discharge (m$^3$/s)', color='k',fontsize=10)
    ax1.set_xlim(xmin=0,xmax=last+1)
    ax1.tick_params('y', colors='k')
    # scentific notaion
    ax1.yaxis.set_major_locator(MaxNLocator(nbins=1,integer=True))
    # ax1.yaxis.set_major_formatter(FormatStrFormatter('%1.0e'))
    ax1.ticklabel_format(style="sci",axis="y",scilimits=(0,0))
    ax1.yaxis.major.formatter._useMathText=True 
    ax1.yaxis.offsetText.set_fontsize(10) 
    # ax1.yaxis.label.set_size(10)
    ax1.yaxis.get_offset_text().set_x(-0.05)
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
    ax1.set_xticks(xxlist)
    ax1.set_xticklabels(xxlab,fontsize=10)
    plt.legend(lines,labels,ncol=1,loc='upper right') #, bbox_to_anchor=(1.0, 1.0),transform=ax1.transAxes)
    # station_loc_list=pname[point].split("/")
    # station_name="-".join(station_loc_list) 
    # print ('save',river[point] , station_name)
    # plt.savefig(assim_out+"/figures/disgraph/"+river[point]+"-"+station_name+".png",dpi=500)
    plt.show()
    return 0
#============
if __name__ == '__main__':
    # make_fig(station,expname)
    make_fig("3621400","NOM_WSE_E2O_HWEB_004")