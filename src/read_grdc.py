# usr/lib/python 

import numpy as np
import re
import shutil
import os
import datetime
#--
#os.system("ln -sf ../params.py params.py")
#shutil.copy("../gosh/params.py","params.py")
# import params as pm
#--
def get_grdc_loc(name,info = "a",CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514",mapname="glb_15min"):
  #--ask the river name and a or b
  # a - most downsream loc
  # b - all locations  
  grdc = CaMa_dir + "/map/"+mapname+"/grdc_loc.txt"
  #grdc = "../data/grdc_loc.txt"
  f = open(grdc,"r")
  lines = f.readlines()
  f.close()
  
  station_loc = []
  x_list      = []
  y_list      = [] 
  #---
  for line in lines:
    line    = list(filter(None, re.split(" ",line)))
    grdc_id = line[0]
    river   = line[1]
    station = line[2]
    d_info  = line[3]
    #--
    lon     = float(line[4])
    lat     = float(line[5])
    #--
    ix      = int(line[6])-1
    iy      = int(line[7])-1
    if name == river:
        if info == "a":
            if d_info=="a":
                station_loc.append(station)
                x_list.append(ix)
                y_list.append(iy)
        else:
            station_loc.append(station)
            x_list.append(ix)
            y_list.append(iy)
  
  return station_loc,x_list,y_list
#------------
def get_grdc_loc_dic(name,info = "a",CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514",mapname="glb_15min"):
  #--aske the river name and a or b
  # a - most downsream loc
  # b - all locations  
  grdc = CaMa_dir + "/map/"+mapname+"/grdc_loc.txt"
  #grdc = CaMa_dir + "/map/glb_15min/grdc_loc.txt"
  #grdc = "../data/grdc_loc.txt"
  f = open(grdc,"r")
  lines = f.readlines()
  f.close()
  
  station_loc = {}
#  print station_loc 
  for line in lines:
    line    = list(filter(None, re.split(" ",line)))
    grdc_id = line[0]
    river   = line[1]
    station = line[2]
    d_info  = line[3]
    station_loc[station] = ()
    #--
    lon     = float(line[4])
    lat     = float(line[5])
    #--
    ix      = int(line[6])-1
    iy      = int(line[7])-1
    if info == d_info and info == "a":
      station_loc[station] = (ix,iy)
    else:
      station_loc[station] = (ix,iy)
  
#  print station_loc   
  return station_loc   

#----
def grdc_river_name(CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514",mapname="glb_15min"):
    # get river names
    grdc = CaMa_dir + "/map/"+mapname+"/grdc_loc.txt"
    #grdc = CaMa_dir + "/map/glb_15min/grdc_loc.txt"
    #grdc = "../data/grdc_loc.txt"
    f = open(grdc,"r")
    lines = f.readlines()
    f.close()

    rivername =  []

    for line in lines:
      line    = list(filter(None, re.split(" ",line)))
      grdc_id = line[0]
      river   = line[1]
      d_info  = line[3]
      if d_info == "a":
        rivername.append(river)

    return rivername
#--
def get_grdc_loc_latlon(name,info = "a",CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514",mapname="glb_15min"):
  #--aske the river name and a or b
  # a - most downsream loc
  # b - all locations  
  #grdc = CaMa_dir + "/map/glb_15min/grdc_loc.txt"
  #grdc = "../data/grdc_loc.txt"
  grdc = CaMa_dir + "/map/"+mapname+"/grdc_loc.txt"
  f = open(grdc,"r")
  lines = f.readlines()
  f.close()
  
  station_loc = []
  x_list      = []
  y_list      = [] 
  #---
  for line in lines:
    line    = list(filter(None, re.split(" ",line)))
    grdc_id = line[0]
    river   = line[1]
    station = line[2]
    d_info  = line[3]
    #--
    lon     = float(line[4])
    lat     = float(line[5])
    #--
    ix      = int(line[6])-1
    iy      = int(line[7])-1
    if info == d_info and info == "a":
      station_loc.append(station)
      x_list.append(lon)
      y_list.append(lat)
    else:
      station_loc.append(station)
      x_list.append(lon)
      y_list.append(lat)
  
  return station_loc,x_list,y_list
#---
def get_id(name,CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514",mapname="glb_15min"):
  #--get GRDC id
  #grdc = CaMa_dir + "/map/glb_15min/grdc_loc.txt"
  #grdc = "../data/grdc_loc.txt"
  grdc = CaMa_dir + "/map/"+mapname+"/grdc_loc.txt"
  
  f = open(grdc,"r")
  lines = f.readlines()
  f.close()
  
  gid=-9999
  #---
  for line in lines:
    line    = list(filter(None, re.split(" ",line)))
    grdc_id = int(line[0])
    river   = line[1]
    station = line[2]
    d_info  = line[3]
    #--
    
    if name == station:
      gid=grdc_id
    
  return gid
#--
def get_station(gid,CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4",mapname="glb_15min"):
  #--get GRDC id
  #grdc = CaMa_dir + "/map/glb_15min/grdc_loc.txt"
  #grdc = "../data/grdc_loc.txt"
  # grdc = CaMa_dir + "/map/"+mapname+"/grdc_loc.txt"
  grdc = "/cluster/data6/menaka/HydroDA/dat/grdc_amz_06min.txt"
  f = open(grdc,"r")
  lines = f.readlines()
  f.close()
  
  gid=int(gid)
  # gid=-9999
  #---
  for line in lines[1::]:
    line    = list(filter(None, re.split(";",line)))
    grdc_id = int(line[0])
    river   = line[1]
    station = line[2]
    d_info  = line[3]
    #--
    
    if grdc_id == gid:
      sname=station
    
  return sname
#--
def get_loc_v394(gid,CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514",mapname="glb_15min"):
  #--ask the river name and a or b
  # a - most downsream loc
  # b - all locations  
  #grdc = CaMa_dir + "/map/glb_15min/grdc_loc.txt"
  #grdc = "../data/GRDC_alloc.txt"
  # grdc = CaMa_dir + "/map/"+mapname+"/grdc_loc.txt"
  # grdc = "/cluster/data6/menaka/HydroDA/dat/grdc_"+mapname+".txt"
  #----
  f = open(grdc,"r")
  lines = f.readlines()
  f.close()
  ix=0
  iy=0
  gid=int(gid)
  #---
  for line in lines[1::]:
    line    = list(filter(None, re.split(" ",line)))
    grdc_id = int(line[0])
    u_info  = line[7]
    
    #--
    if gid==grdc_id:
      #print "get_loc_v394", grdc_id 
      ix      = int(line[8])-1
      iy      = int(line[9])-1
  
  return ix,iy
#--
def get_grdc_loc_v396(name,CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514",mapname="glb_15min"):
  #--ask the river name and a or b
  #  all locations
  #grdc = CaMa_dir + "/map/glb_15min/grdc_loc.txt"
  grdc = CaMa_dir + "/map/"+mapname+"/grdc_loc.txt"
  f = open(grdc,"r")
  lines = f.readlines()
  f.close()
  #--
  id_list     = []
  station_loc = []
  x_list      = []
  y_list      = []
  #---
  for line in lines[1::]:
    #print line
    line    = list(filter(None, re.split(";",line)))
    #print line
    grdc_id = line[0]
    river   = line[1].strip()
    station = line[2].strip()
    ix      = int(line[3])-1
    iy      = int(line[4])-1
    ix2     = int(line[5])
    #print river
    # if ix2 != -9999:
    #     continue
    if river == name:
        # print (river, station)
        id_list.append(grdc_id)
        station_loc.append(station)
        x_list.append(ix)
        y_list.append(iy)

  return id_list,station_loc,x_list,y_list
#------------
#--
def get_grdc_latlon_v396(name,CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514",mapname="glb_15min"):
  #--ask the river name and a or b
  #  all locations
  #grdc = CaMa_dir + "/map/glb_15min/grdc_loc.txt"
  grdc = CaMa_dir + "/map/"+mapname+"/grdc_loc.txt"
  f = open(grdc,"r")
  lines = f.readlines()
  f.close()
  #--
  id_list     = []
  station_loc = []
  x_list      = []
  y_list      = []
  lons        = []
  lats        = []
  #---
  for line in lines[1::]:
    #print line
    line    = list(filter(None, re.split(";",line)))
    #print line
    grdc_id = line[0]
    river   = line[1].strip()
    station = line[2].strip()
    ix      = int(line[3])-1
    iy      = int(line[4])-1
    ix2     = int(line[5])
    #print river
    if ix2 != -9999:
        continue
    if river == name:
        # print (river, name)
        id_list.append(grdc_id)
        station_loc.append(station)
        x_list.append(ix)
        y_list.append(iy)
        lons.append()

  return id_list,station_loc,x_list,y_list
#------------
def grdc_river_name_v396(CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514",mapname="glb_15min"):
  #grdc = CaMa_dir + "/map/glb_15min/grdc_loc.txt"
  grdc = CaMa_dir + "/map/"+mapname+"/grdc_loc.txt"
  f = open(grdc,"r")
  lines = f.readlines()
  f.close()
  
  rivername =  []
  
  for line in lines[1::]:
    line    = list(filter(None, re.split(";",line)))
    grdc_id = line[0]
    river   = line[1].strip()
    #d_info  = line[3]
    #if d_info == "a":
    if river not in rivername:
        rivername.append(river)

  return rivername
#--
def get_grdc_station_v396(name,fname="/cluster/data6/menaka/CaMa-Flood_v396a_20200514/map/glb_15min/grdc_loc.txt"):
#CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514",mapname="glb_15min"):
  #--ask the river name and a or b
  #  all locations
  #grdc = CaMa_dir + "/map/glb_15min/grdc_loc.txt"
  # grdc = CaMa_dir + "/map/"+mapname+"/grdc_loc.txt"
  f = open(fname,"r")
  lines = f.readlines()
  f.close()

  station_loc = []
  x_list      = []
  y_list      = []
  #---
  for line in lines[1::]:
    #print line
    line    = list(filter(None, re.split(";",line)))
    #print line
    grdc_id = line[0]
    river   = line[1].strip()
    station = line[2].strip()
    ix1     = int(line[3])-1
    iy1     = int(line[4])-1
    ix2     = int(line[5])
    iy2     = int(line[6])
    #print river
    if ix2 != -9999:
        ix2 = ix2-1
    if iy2 != -9999:
        iy2 = iy2-1
    # if station == name:
    if grdc_id == name:
        return ix1,iy1,ix2,iy2
#------------
def grdc_dis(grdc_id,syear=2000,eyear=2020,smon=1,emon=12,sday=1,eday=31):
  start_dt=datetime.date(syear,smon,sday)
  last_dt=datetime.date(eyear,emon,eday)
  #--
  start=0
  last=(last_dt-start_dt).days + 1
  #iid=grdc_id()[name]
  #iid="%d"%(grdc_id)
  iid=grdc_id
  # read grdc q
  # grdc ="/cluster/data6/menaka/GRDC_2019/"+iid+"_Q_Day.Cmd.txt"
  grdc ="/cluster/data7/menaka/GRDC_2021/"+iid+"_Q_Day.Cmd.txt"
  if not os.path.exists(grdc):
      return np.ones([last],np.float32)*-9999.0
  else:
      with open(grdc,"r") as f: #,"r"
        lines = f.readlines()
      dis = {}
      for line in lines[37::]:
        line     = list(filter(None, re.split(";",line)))
        yyyymmdd = list(filter(None, re.split("-",line[0])))
        #print yyyymmdd
        yyyy     = int(yyyymmdd[0])
        mm       = int(yyyymmdd[1])
        dd       = int(yyyymmdd[2])
        #---
        #print start_dt.year,start_dt.month,start_dt.day
        if start_dt <= datetime.date(yyyy,mm,dd) and last_dt  >= datetime.date(yyyy,mm,dd):
          # if float(line[2])==-999.000:
          if float(line[2])<=0.0: # updated
            dis[yyyy,mm,dd]=-9999.0
          else:
            dis[yyyy,mm,dd]=float(line[2])
          #print float(line[2])
        elif last_dt  < datetime.date(yyyy,mm,dd):
          break
      #---
      start=0
      last=(last_dt-start_dt).days + 1
      Q=[]
      for day in np.arange(start,last):
        target_dt=start_dt+datetime.timedelta(days=day)
        if (target_dt.year,target_dt.month,target_dt.day) in dis.keys():
          Q.append(dis[target_dt.year,target_dt.month,target_dt.day])
        else:
          Q.append(-9999.0)
      return np.array(Q)
#--
