# usr/lib/python 

import numpy as np
import re
import shutil
import os
import datetime

#====================================================================
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
def read_wse(experiment,station,type0):
    asm=[]
    opn=[]
    fname="./txt/"+experiment+"/wse."+type0+"/"+station+".txt"
    with open(fname,"r") as f:
        lines=f.readlines()
    for line in lines:
        line    = re.split(" ",line)
        line    = list(filter(None, line))
        asm.append(float(line[1]))
        opn.append(float(line[2]))
    return np.array(asm), np.array(opn)
#====================================================================
def get_GRDClist(fname):
    # obslist="/cluster/data6/menaka/HydroDA/dat/grdc_"+mapname+".txt"
    with open(fname,"r") as f:
        lines=f.readlines()
    stationlist=[]
    for line in lines[1::]:
        line    = re.split(";",line)
        line    = list(filter(None, line))
        # print (line)
        num     = line[0].strip()
        # basin   = line[1].strip()
        # stream  = line[2].strip()
        ix1     = int(line[3])
        iy1     = int(line[4])
        # ix2     = int(line[5])
        # iy2     = int(line[6])
        # staid   = int(num)
        #-------------------------
        if rivermap[iy1-1,ix1-1] !=1.0:
            continue
        stationlist.append(num)
    return np.array(stationlist)
#====================================================================