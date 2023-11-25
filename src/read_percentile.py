#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import datetime
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
import sys
import os
import calendar
from multiprocessing import Pool
from multiprocessing import Process
from multiprocessing import sharedctypes
from numpy import ma
import re
import math

import read_grdc as grdc
#====================================================================
def ulbound(data,alpha=0.05):
    a1=(1-(alpha/2.0))*100.0
    a2=(alpha/2.0)*100.0
    upb=np.percentile(data,a1)
    lwb=np.percentile(data,a2)
    return upb, lwb
#====================================================================
argv=sys.argv
syear=int(argv[1])
eyear=int(argv[2])
CaMa_dir=argv[3]
mapname=argv[4]
expname=sys.argv[5]
ncpus=int(sys.argv[6])

ens_mem=49
indir="/cluster/data7/menaka/HydroDA/out"
#==================================
# syear=2019#int(argv[1])
# eyear=2019#int(argv[2])
# CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
# mapname="amz_06min"
# expname="NOM_WSE_E2O_HWEB_006"
# ncpus=5
#====================================================================
fname=CaMa_dir+"/map/"+mapname+"/params.txt"
with open(fname,"r") as f:
    lines=f.readlines()
#-------
nx     = int(filter(None, re.split(" ",lines[0]))[0])
ny     = int(filter(None, re.split(" ",lines[1]))[0])
gsize  = float(filter(None, re.split(" ",lines[3]))[0])
#----
start_dt=datetime.date(syear,1,1)
end_dt=datetime.date(eyear,12,31)
size=60

start=0
last=(end_dt-start_dt).days + 1
#====================================================================
obslist="/cluster/data6/menaka/HydroDA/dat/grdc_"+mapname+".txt"
with open(obslist,"r") as f:
    lines=f.readlines()
pname=[]
for line in lines[1::]:
    line    = re.split(";",line)
    line    = list(filter(None, line))
    # print (line)
    num     = line[0].strip()
    # basin   = line[1].strip()
    # stream  = line[2].strip()
    # ix1     = int(line[3])
    # iy1     = int(line[4])
    # ix2     = int(line[5])
    # iy2     = int(line[6])
    # staid   = int(num)
    pname.append(num)
#====================================================================
org=[]
opn=[]
asm=[]

pnum=len(pname)
N=int(last)

# multiprocessing array
opn=np.ctypeslib.as_ctypes(np.zeros([N,ens_mem,pnum],np.float32))
shared_array_opn  = sharedctypes.RawArray(opn._type_, opn)
asm=np.ctypeslib.as_ctypes(np.zeros([N,ens_mem,pnum],np.float32))
shared_array_asm  = sharedctypes.RawArray(asm._type_, asm)

# for parallel calcualtion
inputlist=[]
for day in np.arange(start,last):
    target_dt=start_dt+datetime.timedelta(days=day)
    yyyy='%04d' % (target_dt.year)
    mm='%02d' % (target_dt.month)
    dd='%02d' % (target_dt.day)
    for num in np.arange(1,ens_mem+1):
        numch='%03d'%num
        inputlist.append([yyyy,mm,dd,numch])
        #print (yyyy,mm,dd,numch)

def read_data(inputlist):
    yyyy = inputlist[0]
    mm   = inputlist[1]
    dd   = inputlist[2]
    numch= inputlist[3]
    # print (yyyy,mm,dd,numch)
    #--
    tmp_opn  = np.ctypeslib.as_array(shared_array_opn)
    tmp_asm  = np.ctypeslib.as_array(shared_array_asm)

    # year, mon, day
    year=int(yyyy)
    mon=int(mm)
    day=int(dd)
    num=int(numch)-1
    #--
    target_dt=datetime.date(year,mon,day)
    dt=(target_dt-start_dt).days
    # corrpted
    fname=indir+"/"+expname+"/assim_out/outflw/open/outflw"+yyyy+mm+dd+"_"+numch+".bin"
    # print (fname)
    #fname=assim_out+"/assim_out/rivout/open/rivout"+yyyy+mm+dd+"_"+numch+".bin"
    opnfile=np.fromfile(fname,np.float32).reshape([ny,nx])
    # assimilated
    fname=indir+"/"+expname+"/assim_out/outflw/assim/outflw"+yyyy+mm+dd+"_"+numch+".bin"
    #fname=assim_out+"/assim_out/rivout/assim/rivout"+yyyy+mm+dd+"_"+numch+".bin"
    asmfile=np.fromfile(fname,np.float32).reshape([ny,nx])
    #-------------
    for point in np.arange(pnum):
        ix1,iy1,ix2,iy2=grdc.get_grdc_station_v396(pname[point],fname=obslist)
        if ix2 == -9999 or iy2 == -9999:
            tmp_opn[dt,num,point]=opnfile[iy1,ix1]
            tmp_asm[dt,num,point]=asmfile[iy1,ix1]
        else:
            tmp_opn[dt,num,point]=opnfile[iy1,ix1]+opnfile[iy2,ix2]
            tmp_asm[dt,num,point]=asmfile[iy1,ix1]+asmfile[iy2,ix2]
#--------
p   = Pool(ncpus)
res = p.map(read_data, inputlist)
opn = np.ctypeslib.as_array(shared_array_opn)
asm = np.ctypeslib.as_array(shared_array_asm)
p.terminate()
#====================================================================
def write_txt(inputlist):
    point=int(inputlist[0])
    syear=int(inputlist[1])
    eyear=int(inputlist[2])
    station=pname[point]
    st_dt=datetime.date(syear,1,1)
    ed_dt=datetime.date(eyear,12,31)
    start=0
    last=(ed_dt-st_dt).days + 1
    fname="./txt/"+expname+"/Q_percentile/"+station+".txt"
    with open(fname,"w") as f:
        for day in np.arange(start,last):
            t_dt=st_dt+datetime.timedelta(days=1)
            upba, lwba=ulbound(asm[day,:,point])
            upbo, lwbo=ulbound(opn[day,:,point])
            line="%04d%02d%02d%10.2f%10.2f%10.2f%10.2f\n"%(t_dt.year,t_dt.month,t_dt.day,upba,lwba,upbo,lwbo)
            print line
            f.write(line)
    return 0

# for parallel 
inputlist=[]
for point in np.arange(pnum):
    inputlist.append(["%02d"%(point),"%04d"%(syear),"%04d"%(eyear)])

p = Pool(ncpus)
p.map(write_txt, inputlist)
p.terminate()