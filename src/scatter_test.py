#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
# import datetime
from matplotlib.colors import LogNorm,Normalize,ListedColormap,BoundaryNorm
# from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import sys
import os
# import calendar
# from multiprocessing import Pool
# from multiprocessing import Process
# from multiprocessing import sharedctypes
from numpy import ma
import re
# import my_colorbar as mbar
# import cartopy.crs as ccrs
# import cartopy
# from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
# import cartopy.feature as cfeature
import os
import string
# figure in A4 size
va_margin= 0.0#1.38#inch 
ho_margin= 0.0#1.18#inch
hgt=(11.69 - 2*va_margin)*(1.0/4.0)*2.0
wdt=(8.27 - 2*ho_margin)
fig=plt.figure(figsize=(wdt,hgt))
G = gridspec.GridSpec(2,2)
#=========================
west0=-80
east0=-45
south0=-22
north0=7
cmap=cm.get_cmap("PRGn")
vmin=-1.0
vmax=1.0
norm=Normalize(vmin=vmin,vmax=vmax)
land="#C0C0C0"
water="#FFFFFF"
#=========================
ax=[]
for i,exp in enumerate(range(2)):
    ax.append(fig.add_subplot(G[i,0]))
    m = Basemap(projection='cyl',llcrnrlat=south0,urcrnrlat=north0,llcrnrlon=west0,urcrnrlon=east0, lat_ts=0,resolution='c',ax=ax[2*i])
    # m.drawcoastlines( linewidth=0.1, color='k' )
    # m.fillcontinents(color=land,lake_color=water,zorder=99)
    # m.drawmapboundary(fill_color=water,zorder=100)
    im=plt.scatter([],[],c=[],cmap=cmap,s=0.1,vmin=vmin,vmax=vmax,norm=norm,zorder=101)
    im.set_visible(False)
    # m.drawparallels(np.arange(south,north+0.1,5), labels = [1,0,0,0], fontsize=10,linewidth=0,zorder=102)
    # m.drawmeridians(np.arange(west,east+0.1,5), labels = [0,0,0,1], fontsize=10,linewidth=0,zorder=102)
    #--
    m.readshapefile('./etc/Amazon_basin/Amazon_boundry', 'Amazon_boundry', drawbounds = False)
    for info, shape in zip(m.Amazon_boundry, m.Amazon_boundry):
        # if info['nombre'] == 'Amazon_boundry':
        x, y = zip(*shape) 
        m.plot(x, y, marker=None,color='grey',linewidth=1.0)
    #===
    ax[2*i].spines['top'].set_visible(False)
    ax[2*i].spines['right'].set_visible(False)
    ax[2*i].spines['bottom'].set_visible(False)
    ax[2*i].spines['left'].set_visible(False)
    if i==0:
        ax[2*i].text(0.5,1.05,"$Q$",ha="center",va="center",transform=ax[2*i].transAxes,fontsize=10)
    ax[2*i].text(-0.05,1.05,"%s)"%(string.ascii_lowercase[2*i]),ha="left",va="center",transform=ax[2*i].transAxes,fontsize=10)
    
    ax.append(fig.add_subplot(G[i,1]))
    m = Basemap(projection='cyl',llcrnrlat=south0,urcrnrlat=north0,llcrnrlon=west0,urcrnrlon=east0, lat_ts=0,resolution='c',ax=ax[2*i+1])
    # m.drawcoastlines( linewidth=0.1, color='k' )
    # m.fillcontinents(color=land,lake_color=water,zorder=99)
    # m.drawmapboundary(fill_color=water,zorder=100)
    im=plt.scatter([],[],c=[],cmap=cmap,s=0.1,vmin=vmin,vmax=vmax,norm=norm,zorder=101)
    im.set_visible(False)
    # m.drawparallels(np.arange(south,north+0.1,5), labels = [1,0,0,0], fontsize=10,linewidth=0,zorder=102)
    # m.drawmeridians(np.arange(west,east+0.1,5), labels = [0,0,0,1], fontsize=10,linewidth=0,zorder=102)
    #--
    m.readshapefile('./etc/Amazon_basin/Amazon_boundry', 'Amazon_boundry', drawbounds = False)
    for info, shape in zip(m.Amazon_boundry, m.Amazon_boundry):
        # if info['nombre'] == 'Amazon_boundry':
        x, y = zip(*shape) 
        m.plot(x, y, marker=None,color='grey',linewidth=1.0)
    #===
    ax[2*i+1].spines['top'].set_visible(False)
    ax[2*i+1].spines['right'].set_visible(False)
    ax[2*i+1].spines['bottom'].set_visible(False)
    ax[2*i+1].spines['left'].set_visible(False)
    if i==0:
        ax[2*i+1].text(0.5,1.05,"$WSE$",ha="center",va="center",transform=ax[2*i+1].transAxes,fontsize=10)
    ax[2*i+1].text(-0.05,1.05,"%s)"%(string.ascii_lowercase[2*i+1]),ha="left",va="center",transform=ax[2*i+1].transAxes,fontsize=10)

#--
caxQ=fig.add_axes([0.1,0.05,0.35,0.01])
cbarQ=plt.colorbar(im,cax=caxQ,ticks=[-0.4,-0.2,0.0,0.2,0.4],orientation='horizontal',extend="both")#np.arange(-0.5,0.5+0.001,0.25)
cbarQ.set_label("$\Delta$$r$") #$(r_{assim} - r_{open})$") correlation
cbarQ.ax.tick_params(labelsize=6.0)
# cbar.set_label("$r $ $correlation$") #$(r_{assim} - r_{open})$")
# caxW=fig.add_axes([0.55,0.05,0.35,0.01])
# cbarW=plt.colorbar(imW,cax=caxW,ticks=[-0.8,-0.4,0.0,0.4,0.8],orientation='horizontal',extend="both") #np.arange(-1.0,1.0+0.001,0.2)
# cbarW.set_label("$rRMSE$") #$(r_{assim} - r_{open})$")
# cbarW.ax.tick_params(labelsize=6.0)

caxW=fig.add_axes([0.55,0.05,0.20,0.01])
cbarW=plt.colorbar(im,cax=caxW,ticks=[-0.8,-0.4,0.0,0.4,0.8],orientation='horizontal',extend="both") #np.arange(-1.0,1.0+0.001,0.2)
cbarW.set_label("$rRMSE$") #$(r_{assim} - r_{open})$")
cbarW.ax.tick_params(labelsize=6.0)

# legend 
feature1=mlines.Line2D([], [], color='grey', marker='s',
                          markersize=5, label='Assimilation',linewidth=0.0)
feature2=mlines.Line2D([], [], color='grey', marker='o',
                          markersize=5, label='Validation',linewidth=0.0)
# plt.axis('off')
# legend=plt.legend(handles=[feature1,feature2],bbox_to_anchor=(0.8, 0.001, 0.15, 0.1), loc=3,
#            ncol=1,  borderaxespad=0.0)#,fancybox=False) mode="expand",
legend=plt.legend(handles=[feature1,feature2],bbox_to_anchor=(0.9,0.02), loc="lower right",
           bbox_transform=fig.transFigure, ncol=1,  borderaxespad=0.0, frameon=False)#
legend.get_frame().set_alpha(None)
# legend.get_frame().set_facecolor((0, 0, 1, 0.1))
# legend.set_facecolor('white')
# legend.set_edgecolor('white')

#plt.title(stitle)
plt.savefig("figures/amazon_test.png",dpi=800,bbox_inches="tight", pad_inches=0.0)
plt.savefig("figures/amazon_test.png",dpi=800,bbox_inches="tight", pad_inches=0.0)
# plt.show()