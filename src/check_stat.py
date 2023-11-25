#!/opt/local/bin/python
# -*- coding: utf-8 -*-
'''
check statistic values
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

#==========
ny,nx=250, 350
ens_mem=49
odir="/cluster/data6/menaka/HydroDA/dat"
means=["sfcelv_49_ECMWF_amz_06min_2000-2014", # original statistic
"sfcelv_courrpt_ECMWF_amz_06min_2000-2014", #  courrpted bathymetry statistic
"sfcelv_bias_ECMWF_amz_06min_2000-2014", # biased runoff statistic
"sfcelv_bias_corrupt_ECMWF_amz_06min_2000-2014"] # biased courrpted bathymetry statistic
station="3629001"
ix=245
iy=70
valss=[]
for mean in means:
    vals=[]
    for ens in np.arange(1,ens_mem+1):
        fname=odir+"/mean_"+mean+"_%03d.bin"%(ens)
        sfcelv=np.fromfile(fname,np.float32).reshape(ny,nx)
        val=sfcelv[iy-1,ix-1]
        vals.append(val)
    valss.append(vals)

valss=np.array(valss)
valss=np.transpose(valss)
# print (valss)

# plot boxplot using seaborn
colors=["w","#004488","#ddaa33","#ba5566"]
# cmap   = matplotlib.colors.ListedColormap(matplotlib.cm.get_cmap("viridis_r").colors[:len(basins)])
va_margin= 0.0#1.38#inch 
ho_margin= 0.0#1.18#inch
hgt=(11.69 - 2*va_margin)*(8.27/11.69)*(2.0/3.0)
wdt=(8.27 - 2*ho_margin)*(8.27/8.27)
#fig=plt.figure(figsize=(8.27,11.69))
fig=plt.figure(figsize=(wdt,hgt))
#fig.suptitle("auto-correalated area")#"maximum autocorrelation length")
#G = gridspec.GridSpec(2,1)
G = gridspec.GridSpec(ncols=1,nrows=1)
ax = fig.add_subplot(G[0,0])
#
flierprops = dict(marker='o', markerfacecolor='none', markersize=8,linestyle='none', markeredgecolor='k')
boxprops = dict(color='grey')#facecolor='none'
whiskerprops = dict(color='grey',linestyle="--")
capprops = dict(color='grey')
medianprops = dict(color='grey',linestyle="-",linewidth=1.0)
meanprops = dict(marker='D', markeredgecolor='black',markerfacecolor='green',markersize=8)
#==========
box=sns.boxplot(ax=ax,data=valss, fliersize=0.0, palette=colors, whis=1.5\
    ,meanline=True, width=0.8, linewidth=0.3, dodge=True\
    ,meanprops=meanprops,capprops=capprops,medianprops=medianprops)
# median line
ax.axhline(y=0.0, linestyle="--", color="k",linewidth=0.5)
ax.set_xticklabels(["Original","Corrupt", "Bias", "Bias-Corrupt"],rotation=90)
plt.show()