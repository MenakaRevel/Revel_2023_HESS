#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os
from multiprocessing import Pool
import re
from numpy import ma
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm,Normalize,ListedColormap,BoundaryNorm
from mpl_toolkits.basemap import Basemap
#====================================================================
def vec_par(LEVEL,ax=None):
    ax=ax or plt.gca()
    txt=re.split("-",figname)[0]+"_%02d.txt"%(LEVEL)
    os.system("./bin/print_rivvec "+tmp0+" 1 "+str(LEVEL)+" > "+txt)
    width=(float(LEVEL)**sup)*w
    #print LEVEL, width#, lon1,lat1,lon2-lon1,lat2-lat1#x1[0],y1[0],x1[1]-x1[0],y1[1]-y1[0]
    # open tmp2.txt
    with open(txt,"r") as f:
        lines = f.readlines()

    #print LEVEL, width, lines, txt
    #---
    for line in lines:
        line = filter(None, re.split(" ",line))
        lon1 = float(line[0])
        lat1 = float(line[1])
        lon2 = float(line[3])
        lat2 = float(line[4])

        #- higher resolution data
        ixx = int((lon1  - west)*60.0)
        iyy = int((-lat1 + north)*60.0)

        #----
        ix =catmxy[0,iyy,ixx] 
        iy =catmxy[1,iyy,ixx]

        #==
        # ix = int((lon1 - west)*(1/gsize))
        # iy = int((-lat1 + north)*(1/gsize))

        if rivermap[iy-1,ix-1] == 0:
            continue

        if lon1-lon2 > 180.0:
            # print (lon1, lon2)
            lon2=180.0
        elif lon2-lon1> 180.0:
            # print (lon1,lon2)
            lon2=-180.0
        #--------
        # colorVal="grey"#"k" 
        colorVal=cmap(norm(data[iy-1,ix-1]))
        #print (lon1,lon2,lat1,lat2,width)
        plot_ax(lon1,lon2,lat1,lat2,width,colorVal,ax)
#====================================================================
def plot_ax(lon1,lon2,lat1,lat2,width,colorVal,ax=None):
    ax=ax or plt.gca()
    return ax.plot([lon1,lon2],[lat1,lat2],color=colorVal,linewidth=width,zorder=105,alpha=alpha)
#====================================================================
#
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"
mapname="amz_06min"
coeff=0.25
figname="figs5-diff_corrupt_rivhgt"
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
#higher resolution data
catmxy = CaMa_dir+"/map/"+mapname+"/1min/1min.catmxy.bin"
catmxy = np.fromfile(catmxy,np.int16).reshape(2,ny*6,nx*6)
#===================
rivnum="/cluster/data6/menaka/HydroDA/dat/rivnum_"+mapname+".bin"
rivnum=np.fromfile(rivnum,np.int32).reshape(ny,nx)
rivermap=((nextxy[0]>0)*(rivnum==1))*1.0
#----
rivgt_corrupt=(rivhgt>=5.0)*(rivhgt+rivhgt*coeff)+(rivhgt<5.0)*rivhgt
rivgt_corrupt.tofile(CaMa_dir+"/map/"+mapname+"/rivhgt_corrupt.bin")
#====================================
#---------- making fig --------------
#====================================
#files for making river network
tmp0=re.split("-",figname)[0]+".txt"

land="#C0C0C0"
water="#FFFFFF"
#==========================
# river width
sup=3 #2
w=0.002 #0.02
alpha=1
width=0.5
#------
west0=west
east0=east
south0=south-2
north0=north+2
# figure in A4 size
# norm=norm=Normalize(vmin=vmin,vmax=vmax)
vmin=-10.0
vmax=10.0
# cmap=plt.cm.get_cmap("bwr_r")
cmap=plt.cm.get_cmap("coolwarm_r")
norm=BoundaryNorm(np.arange(vmin,vmax+0.1,1),cmap.N)
# cmap   = matplotlib.colors.ListedColormap(matplotlib.cm.get_cmap("viridis_r").colors[:len(basins)])
va_margin= 0.0#1.38#inch 
ho_margin= 0.0#1.18#inch
hgt=(11.69 - 2*va_margin)*(8.27/11.69)*(1.0/3.0)
wdt=(8.27 - 2*ho_margin)*(8.27/8.27)
#fig=plt.figure(figsize=(8.27,11.69))
fig=plt.figure(figsize=(wdt,hgt))
#fig.suptitle("auto-correalated area")#"maximum autocorrelation length")
#G = gridspec.GridSpec(2,1)
G = gridspec.GridSpec(ncols=1,nrows=1)
ax= fig.add_subplot(G[0,0])
m = Basemap(projection='cyl',llcrnrlat=south0,urcrnrlat=north0,llcrnrlon=west0,urcrnrlon=east0, lat_ts=0,resolution='c',ax=ax)
# amazon shapefile
m.readshapefile('./etc/Amazon_basin/Amazon_boundry', 'Amazon_boundry', drawbounds = False)
for info, shape in zip(m.Amazon_boundry, m.Amazon_boundry):
    # if info['nombre'] == 'Amazon_boundry':
    x, y = zip(*shape) 
    m.plot(x, y, marker=None,color='grey',linewidth=1.0)
box="%f %f %f %f"%(west,east,north,south) 
os.system("./bin/txt_vector "+box+" "+CaMa_dir+" "+mapname+" > "+tmp0) 
data=ma.masked_where(rivhgt<0.0,rivgt_corrupt-rivhgt).filled(0.0)
#map(vec_par,np.arange(1,10+1,1))
map(vec_par,np.arange(2,10+1,1))
# map(vec_par,np.arange(5,10+1,1))
# im=m.imshow(ma.masked_where(rivhgt<0.0,rivgt_corrupt-rivhgt),origin="upper",interpolation="nearest",cmap=cmap,vmin=vmin,vmax=vmax)
#======================
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
#======================
# plt.colorbar(im)
# im=m.imshow([],[],origin="upper",interpolation="nearest",cmap=cmap,vmin=vmin,vmax=vmax)
im=plt.scatter([],[],c=[],cmap=cmap,s=0.1,vmin=vmin,vmax=vmax,norm=norm,zorder=101)
cbar=m.colorbar(im,"right",size="2%",extend="both",ticks=np.arange(vmin,vmax+0.001,2.0))
cbar.set_label("Difference of river channel depth\n (courrpted - original)")
# plt.show()
#===
plt.savefig("./figures/"+figname+".pdf",dpi=800,bbox_inches="tight", pad_inches=0.05)
plt.savefig("./figures/"+figname+".png",dpi=800,bbox_inches="tight", pad_inches=0.05)
plt.savefig("./figures/"+figname+".jpg",dpi=800,bbox_inches="tight", pad_inches=0.05)
os.system("rm -r "+re.split("-",figname)[0]+"*.txt")