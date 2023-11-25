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
import warnings;warnings.filterwarnings('ignore')
from adjustText import adjust_text
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression

from statistics import *
import read_grdc as grdc
import my_colorbar as mbar
import read_hydroweb as hweb
from river_patch import patchms
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
syear=2009
eyear=2014
# CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"
#map="glb_15min"
# map="glb_06min"
mapname="amz_06min"
# expname="NOM_WSE_E2O_HWEB_004"
lexp=["DIR_WSE_E2O_HWEB_001","ANO_WSE_E2O_HWEB_004","NOM_WSE_E2O_HWEB_004",]
exp="NOM_WSE_E2O_HWEB_004"
figname="figs13-deltaA_relative_metrics"
#====================================================================
obslist="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_"+mapname+"_QC0_simulation.txt"
with open(obslist,"r") as f:
    lines=f.readlines()
metric_frag=[]
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
    #=========
    try:
        asm1, opn1 = read_wse(exp,station,"simulation")
    except:
        asm1, opn1 = read_wse(exp,station,"validation")

    org1=hweb.HydroWeb_continous_WSE(station,syear,1,1,eyear,12,31,EGM08,EGM96)
    org1=np.array(org)#+np.array(EGM08)-np.array(EGM96)
    delAa=rAmplitude(asm1,org1)
    delAo=rAmplitude(opn1,org1)
    print delAa, delAo