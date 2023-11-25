#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import datetime
from matplotlib.colors import LogNorm,Normalize,ListedColormap
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
import my_colorbar as mbar
import cartopy.crs as ccrs
import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
import os
#########################
def get_HydroWeb(fname):
	#---
    lname=[]
    l_lon=[]
    l_lat=[]
    xlist=[]
    ylist=[]
    llsat=[]
    # fname=DA_dir+"/dat/HydroWeb_alloc_"+mapname+".txt"
    # fname=DA_dir+"/dat/HydroWeb_alloc_"+mapname+"_new.txt"
    # fname=DA_dir+"/dat/HydroWeb_alloc_"+mapname+"_QC0.txt"
    with open(fname,"r") as f:
        lines=f.readlines()
    for line in lines[1::]:
        line    = re.split(" ",line)
        line    = list(filter(None, line))
        # print (line)
        #================
        num     = line[0]
        station = line[1]
        riv     = re.split("_",station)[1]
        lon     = float(line[2])
        lat     = float(line[3])
        ix      = int(line[4]) - 1
        iy      = int(line[5]) - 1
        sat     = line[10].rstrip()
        #------
        lname.append(station)
        l_lon.append(lon)
        l_lat.append(lat)
        xlist.append(ix)
        ylist.append(iy)
        llsat.append(sat)
    return lname, l_lon, l_lat, xlist, ylist, llsat
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

        if rivermap[iy,ix] == 0:
            continue

        if lon1-lon2 > 180.0:
            # print (lon1, lon2)
            lon2=180.0
        elif lon2-lon1> 180.0:
            # print (lon1,lon2)
            lon2=-180.0
        #--------
        colorVal="grey"#"w" 
        #print (lon1,lon2,lat1,lat2,width)
        plot_ax(lon1,lon2,lat1,lat2,width,colorVal,ax)
#==================================
def plot_ax(lon1,lon2,lat1,lat2,width,colorVal,ax=None):
    ax=ax or plt.gca()
    return ax.plot([lon1,lon2],[lat1,lat2],color=colorVal,linewidth=width,zorder=105,alpha=alpha)
#==================================
argv=sys.argv
syear=int(argv[1])
eyear=int(argv[2])
CaMa_dir=argv[3]
mapname=argv[4]
obslist1=argv[5]
obslist2=argv[6]
figname=argv[7]
ncpus=int(sys.argv[8])
#==================
# DA_dir="/cluster/data6/menaka/HydroDA"
# CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"
# mapname="amz_06min"
# obslist="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_amz_06min_QC0.txt"
ens_mem=49
lexp=6
upz=1e09 #m2
#----
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
#===================
satcov="./dat/satellite_coverage.bin"
satcov=np.fromfile(satcov,np.float32).reshape(ny,nx)
#===================
land="#C0C0C0"
water="#FFFFFF"

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
# # im=plt.scatter([],[],c=[],cmap=cmap,s=0.1,vmin=vmin,vmax=vmax,norm=norm,zorder=101)
# # im.set_visible(False)
# m.drawparallels(np.arange(south,north+0.1,5), labels = [1,0,0,0], fontsize=10,linewidth=0,zorder=102)
# m.drawmeridians(np.arange(west,east+0.1,5), labels = [0,0,0,1], fontsize=10,linewidth=0,zorder=102)
# amazon shapefile
m.readshapefile('./etc/Amazon_basin/Amazon_boundry', 'Amazon_boundry', drawbounds = False)
for info, shape in zip(m.Amazon_boundry, m.Amazon_boundry):
    # if info['nombre'] == 'Amazon_boundry':
    x, y = zip(*shape) 
    m.plot(x, y, marker=None,color='grey',linewidth=1.0)
#--
box="%f %f %f %f"%(west,east,north,south) 
os.system("./bin/txt_vector "+box+" "+CaMa_dir+" "+mapname+" > tmp1.txt") 
# #map(vec_par,np.arange(1,10+1,1))
map(vec_par,np.arange(2,10+1,1))
#######
for obs in ["sim","val"]:
    c="g"
    m="o"
    if obs=="sim":
        obslist=obslist1
        marker="o"
    else:
        obslist=obslist2
        marker="s"
    #===================
    lname, l_lon, l_lat, xlist, ylist, llsat = get_HydroWeb(fname=obslist)
    pnum = len(lname)
    #--
    for point in np.arange(pnum):
        lon=l_lon[point]
        lat=l_lat[point]
        ix =xlist[point]
        iy =ylist[point]
        sat=llsat[point]
        upa=uparea[iy,ix]
        # print sat
        # c=cmapL(norm(NSEAI))
        # if upa >= upz:
        #     m="D"
        # else:
        #     m="o"
        if sat == "ENVISAT":
            # c="g"
            # c="xkcd:electric blue"
            # c="xkcd:steel blue"
            # c="#95d0fc"
            # c="#0165fc"
            c="#069af3"
        elif sat == "JASON2":
            # c="r"
            # c="xkcd:goldenrod"
            # c="xkcd:pinkish"
            # c="#e50000" #"#cb416b" #
            # c="#ff000d"
            c="#ff000d"
        elif sat == "JASON3":
            # c="b"
            # c="#e50000" #"#cb416b" #
            c="#ff000d"
        elif sat == "SENTINEL3A" or sat == "SENTINEL3B":
            # c="k"
            c="#cb416b"
        ax.scatter(lon,lat,s=10,marker=marker,edgecolors="k",linewidth=0.3, facecolors=c,zorder=106)
#=======
labels=["ENVISAT", "Jason 2/3"]#, "Jason 3"]#, "SENTINEL3A/B"]
# colors=["g","r"]#,"r"]#,"k"]
# colors=["xkcd:electric blue","xkcd:goldenrod"]
# colors=["xkcd:steel blue","xkcd:pinkish"]
# colors=["#95d0fc","#cb416b"]
# colors=["#95d0fc","#e50000"]
# colors=["#0165fc","#ff000d"]
colors=["#069af3","#ff000d"]
patches=[]
for point in np.arange(len(labels)):
    patches.append(mpatches.Patch(color=colors[point],label=labels[point]))
#--
legend1=plt.legend(handles=patches,bbox_to_anchor=(0.9,0.13), loc="lower right",
           bbox_transform=fig.transFigure, ncol=1,  borderaxespad=0.0, frameon=False)#

# Add the legend manually to the Axes.
plt.gca().add_artist(legend1)

# cbar=m.colorbar(im,"right",size="2%",ticks=np.arange(vmin,vmax+0.001,0.2))
#plt.title(stitle)
# legend 
feature1=mlines.Line2D([], [], color='grey', marker='s',
                          markersize=5, label='Assimilation',linewidth=0.0)
feature2=mlines.Line2D([], [], color='grey', marker='o',
                          markersize=5, label='Validation',linewidth=0.0)
# plt.axis('off')
# legend=plt.legend(handles=[feature1,feature2],bbox_to_anchor=(0.8, 0.001, 0.15, 0.1), loc=3,
#            ncol=1,  borderaxespad=0.0)#,fancybox=False) mode="expand",
legend2=plt.legend(handles=[feature1,feature2],bbox_to_anchor=(0.1,0.13), loc="lower left",
           bbox_transform=fig.transFigure, ncol=1,  borderaxespad=0.0, frameon=False)#
#======================
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
#======================
# plt.savefig("./figures/AltiVS_scatter_QC0.pdf",dpi=500,bbox_inches="tight", pad_inches=0.05)
# plt.savefig("./figures/AltiVS_scatter_QC0.jpg",dpi=500,bbox_inches="tight", pad_inches=0.05)
plt.savefig("./figures/"+figname+".jpg",dpi=500,bbox_inches="tight", pad_inches=0.0)
plt.savefig("./figures/"+figname+".png",dpi=500,bbox_inches="tight", pad_inches=0.0)
plt.savefig("./figures/"+figname+".pdf",dpi=800,bbox_inches="tight", pad_inches=0.0)
os.system("rm -r tmp*.txt")