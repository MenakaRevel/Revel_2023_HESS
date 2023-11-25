#!/opt/local/bin/python
# -*- coding: utf-8 -*-

from configparser import Interpolation
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import datetime
from matplotlib.colors import LogNorm,Normalize,ListedColormap,LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import matplotlib.colors as colors
import sys
import os
import calendar
from multiprocessing import Pool
from multiprocessing import Process
from multiprocessing import sharedctypes
from numpy import ma
import re
import my_colorbar as mbar
import cartopy.crs as ccrs
import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
import os

import read_grdc as grdc
#########################
#====================================================================
def get_grdc_list(mapname):
    obslist="/cluster/data6/menaka/HydroDA/dat/grdc_"+mapname+".txt"
    with open(obslist,"r") as f:
        lines=f.readlines()
    pname=[]
    IX1=[]
    IY1=[]
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
        pname.append(num)
        IX1.append(ix1)
        IY1.append(iy1)
    return pname,IX1,IY1
#==================================
def vec_par(LEVEL,ax=None):
    ax=ax or plt.gca()
    txt="tmp_%02d.txt"%(LEVEL)
    os.system("./bin/print_rivvec tmp1.txt 1 "+str(LEVEL)+" > "+txt)
    width=(float(LEVEL)**sup)*w
    #print LEVEL, width#, lon1,lat1,lon2-lon1,lat2-lat1#x1[0],y1[0],x1[1]-x1[0],y1[1]-y1[0]
    # open tmp2.txt
    f = open(txt,"r")
    lines = f.readlines()
    f.close()
    #print LEVEL, width, lines, txt
    #---
    for line in lines:
        line    = filter(None, re.split(" ",line))
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
        colorVal="b"#"w" 
        #print (lon1,lon2,lat1,lat2,width)
        plot_ax(lon1,lon2,lat1,lat2,width,colorVal,ax)
#==================================
def plot_ax(lon1,lon2,lat1,lat2,width,colorVal,ax=None):
    ax=ax or plt.gca()
    return ax.plot([lon1,lon2],[lat1,lat2],color=colorVal,linewidth=width,zorder=105,alpha=alpha)
#==================================
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap 
#==================================
def createColormapFromGPF(file_, get_dict=False):
    print (file_)
    # data = sp.loadtxt(file_)
    data = np.loadtxt(file_)
    cdict = {'red': np.take(data, (0, 1, 1), axis=1),
             'green': np.take(data, (0, 2, 2), axis=1),
             'blue': np.take(data, (0, 3, 3), axis=1)}
    name = os.path.splitext(os.path.basename(file_))[0]
    if get_dict:
        return cdict
    else:
        return colors.LinearSegmentedColormap(name, cdict)
#==================================
argv=sys.argv
syear=int(argv[1])
eyear=int(argv[2])
CaMa_dir=argv[3]
mapname=argv[4]
figname=argv[5]
ncpus=int(sys.argv[6])
ens_mem=21
seaborn_map=True
#==================================
# syear=2009
# eyear=2014
# DA_dir="/cluster/data6/menaka/HydroDA"
# CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"
# mapname="amz_06min"
obslist="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_amz_06min_QC0.txt"
ens_mem=49
lexp=6
upthr=1e9 #m2
#----
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
#--------------------------------------------
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
# elevtn=elevtn*((rivermap==1.0)*1.0)
sea=(elevtn>0.0)*1.0
#===================
satcov="./dat/satellite_coverage.bin"
satcov=np.fromfile(satcov,np.float32).reshape(ny,nx)
#--------------------------------------------
# land="#C0C0C0"
# water="#FFFFFF"
land="w"
water="b"

londiff=(east-west)*4
latdiff=(north-south)*4

npix=(90-north)*4
spix=(90-south)*4
wpix=(180+west)*4
epix=(180+east)*4

#cmap=make_colormap(colors_list)
cmap=mbar.colormap("H01")
cmap.set_under("w",alpha=0)
cmapL=cmap #cm.get_cmap("rainbow_r")
vmin=0.0
vmax=1.0
norm=Normalize(vmin=vmin,vmax=vmax)
#==
# print (createColormapFromGPF("./dat/DEM-poster.gpf"))
cmap = createColormapFromGPF("./dat/DEM-poster.gpf")
# cmap = createColormapFromGPF("./dat/DEM-print.gpf")
# cmap = createColormapFromGPF("./dat/wiki-schwarzwald-d020.gpf")
# cmap = plt.get_cmap('terrain')
# cmap = truncate_colormap(cmap, 0.25, 0.9)
cmap.set_under("b")
#
cmapL=ListedColormap(["b"])
#------
# river width
sup=2
w=0.02
alpha=1
width=0.5
#------
print (west,east,south,north)
resol=1
# figure in A4 size
va_margin= 0.0#1.38#inch 
ho_margin= 0.0#1.18#inch
hgt=(11.69 - 2*va_margin)*(1.0/3.0)
wdt=(8.27 - 2*ho_margin)*(1.0/2.0)
fig=plt.figure(figsize=(wdt,hgt))
G = gridspec.GridSpec(1,1)
ax=fig.add_subplot(G[0,0])
m = Basemap(projection='cyl',llcrnrlat=south-2,urcrnrlat=north+2,llcrnrlon=west,urcrnrlon=east, lat_ts=0,resolution='c',ax=ax)
#m.drawcoastlines( linewidth=0.1, color='k' )
# m.fillcontinents(color=land,lake_color=water,zorder=99)
# im=plt.scatter([],[],c=[],cmap=cmap,s=0.1,vmin=vmin,vmax=vmax,norm=norm,zorder=101)
# im.set_visible(False)
# m.drawparallels(np.arange(south,north+0.1,5), labels = [1,0,0,0], fontsize=10,linewidth=0,zorder=102)
# m.drawmeridians(np.arange(west,east+0.1,5), labels = [0,0,0,1], fontsize=10,linewidth=0,zorder=102)
# plot elevation
lat = np.arange(north, south+gsize/2.0, -gsize)
lon = np.arange(west, east+gsize/2.0, gsize)
lon, lat = np.meshgrid(lon, lat)
im=m.pcolormesh(lon, lat, ma.masked_less_equal(elevtn,0.0),latlon=True, cmap=cmap,zorder=101)
m.pcolormesh(lon, lat, ma.masked_greater(sea,0.0),latlon=True, cmap=cmapL,zorder=100)
# im=m.imshow(ma.masked_less_equal(elevtn,0.0),origin="upper",interpolation="nearest",cmap=cmap)#,norm=colors.LogNorm())
# amazon shapefile
m.readshapefile('./etc/Amazon_basin/Amazon_boundry', 'Amazon_boundry', drawbounds = False)
for info, shape in zip(m.Amazon_boundry, m.Amazon_boundry):
    # if info['nombre'] == 'Amazon_boundry':
    x, y = zip(*shape) 
    m.plot(x, y, marker=None,color='k',linewidth=1.0,zorder=102)
#--
box="%f %f %f %f"%(west,east,north,south) 
os.system("./bin/txt_vector "+box+" "+CaMa_dir+" "+mapname+" > tmp1.txt") 
# #map(vec_par,np.arange(1,10+1,1))
map(vec_par,np.arange(2,10+1,1))
# map(vec_par,np.arange(5,10+1,1))

# #######
# lname, xlist, ylist = get_grdc_list(mapname) 
# pnum = len(lname)
# #--
# for point in np.arange(pnum):
#     ix =xlist[point]
#     iy =ylist[point]
#     num=lname[point]
#     upa=uparea[iy-1,ix-1]
#     # if upa < upthr:
#     #     continue
#     #-------------------------
#     if uparea[iy-1,ix-1] < 1.0e9:
#         continue
#     #-------------------------
#     if rivermap[iy-1,ix-1] !=1.0:
#         continue
#     # #-------------------------
#     # if satcov[iy-1,ix-1] !=1.0:
#     #     continue
#     #--------------
#     org=grdc.grdc_dis(num,syear,eyear)
#     if np.sum((org!=-9999.0)*1.0) < 365*1:
#         continue
#     lon=lonlat[0,iy-1,ix-1]
#     lat=lonlat[1,iy-1,ix-1]
#     if satcov[iy-1,ix-1]==1.0:
#         ax.scatter(lon,lat,s=10,marker="o",edgecolors="k", linewidth=0.5, facecolors="xkcd:blue green",zorder=106)
#     else:
#         ax.scatter(lon,lat,s=10,marker="s",edgecolors="k", linewidth=0.5, facecolors="xkcd:grapefruit",zorder=106)
# #=======
# # legend 
# feature1=mlines.Line2D([], [], color='xkcd:blue green', marker='s',
#                           markersize=5, label='Inside Satellite Observations',linewidth=0.0)
# feature2=mlines.Line2D([], [], color='xkcd:grapefruit', marker='o',
#                           markersize=5, label='Outside Satellite Observations',linewidth=0.0)

# legend=plt.legend(handles=[feature1,feature2],bbox_to_anchor=(0.5,0.01), loc="lower center",
#            bbox_transform=fig.transFigure, ncol=1,  borderaxespad=0.0, frameon=False)#
#======================
# colorbar
cb=m.colorbar(im,size="2%")
cb.ax.tick_params(labelsize=8)
#======================
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
#======================
plt.savefig("./figures/"+figname+".pdf",dpi=800,bbox_inches="tight", pad_inches=0.05)
plt.savefig("./figures/"+figname+".jpg",dpi=800,bbox_inches="tight", pad_inches=0.05)
plt.savefig("./figures/"+figname+".png",dpi=800,bbox_inches="tight", pad_inches=0.05)
os.system("rm -r tmp*.txt")