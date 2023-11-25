#!/opt/local/bin/python
# -*- coding: utf-8 -*-
'''
Statistics for streamflow caclculations
'''
import numpy as np
from numpy import ma
import datetime
import calendar
import warnings;warnings.filterwarnings('ignore')

#========================================
#====  functions calculating stats  ====
#========================================
def filter_nan(s,o):
    """
    this functions removed the data  from simulated and observed data
    where ever the observed data contains nan
    """
    data = np.array([s.flatten(),o.flatten()])
    data = np.transpose(data)
    data = data[~np.isnan(data).any(1)]

    return data[:,0],data[:,1]
#====================================================================
def NS(s,o):
    """
    Nash Sutcliffe efficiency coefficient
    input:
        s: simulated
        o: observed
    output:
        ns: Nash Sutcliffe efficient coefficient
    """
    s,o = filter_nan(s,o)
    o=ma.masked_where(o<=0.0,o).filled(0.0)
    s=ma.masked_where(o<=0.0,s).filled(0.0)
    o=np.compress(o>0.0,o)
    s=np.compress(o>0.0,s) 
    return 1 - sum((s-o)**2)/(sum((o-np.mean(o))**2)+1e-20)
#====================================================================
def correlation(s,o):
    """
    correlation coefficient
    input:
        s: simulated
        o: observed
    output:
        correlation: correlation coefficient
    """
    s,o = filter_nan(s,o)
    o=ma.masked_where(o<=0.0,o).filled(0.0)
    s=ma.masked_where(o<=0.0,s).filled(0.0)
    o=np.compress(o>0.0,o)
    s=np.compress(o>0.0,s)
    if s.size == 0:
        corr = 0.0 #np.NaN
    else:
        corr = np.corrcoef(o, s)[0,1]
        
    return corr
#========================================
def KGE(s,o):
    """
	Kling Gupta Efficiency (Kling et al., 2012, http://dx.doi.org/10.1016/j.jhydrol.2012.01.011)
	input:
        s: simulated
        o: observed
    output:
        KGE: Kling Gupta Efficiency
    """
    o=ma.masked_where(o<=0.0,o).filled(0.0)
    s=ma.masked_where(o<=0.0,s).filled(0.0)
    o=np.compress(o>0.0,o)
    s=np.compress(o>0.0,s)
    s,o = filter_nan(s,o)
    B = np.mean(s) / np.mean(o)
    y = (np.std(s) / np.mean(s)) / (np.std(o) / np.mean(o))
    r = np.corrcoef(o, s)[0,1]
    return 1 - np.sqrt((r - 1) ** 2 + (B - 1) ** 2 + (y - 1) ** 2)
#==========================================================
def RMSE(s,o):
    """
    Root Mean Squre Error
    input:
        s: simulated
        o: observed
    output:
        RMSE: Root Mean Squre Error
    """
    o=ma.masked_where(o==-9999.0,o).filled(0.0)
    s=ma.masked_where(o==-9999.0,s).filled(0.0)
    o=np.compress(o>0.0,o)
    s=np.compress(o>0.0,s)
    s,o = filter_nan(s,o)
    # return np.sqrt(np.mean((s-o)**2))
    return np.sqrt(np.ma.mean(np.ma.masked_where(o<=0.0,(s-o)**2)))
#==========================================================
def NRMSE(s,o):
    """
    Normalised Root Mean Squre Error
    input:
        s: simulated
        o: observed
    output:
        NRMSE: Normalised Root Mean Squre Error
    """
    o=ma.masked_where(o==-9999.0,o).filled(0.0)
    s=ma.masked_where(o==-9999.0,s).filled(0.0)
    o=np.compress(o>0.0,o)
    s=np.compress(o>0.0,s)
    s,o = filter_nan(s,o)
    # return np.sqrt(np.mean((s-o)**2))
    return np.sqrt(np.ma.mean(np.ma.masked_where(o<=0.0,(s-o)**2)))/np.mean(o)
#==========================================================
def pBIAS(s,o):
    """
    Precentage Bias
    input:
        s: simulated
        o: observed
    output:
        pBias: Precentage Bias
    """
    s,o = filter_nan(s,o)
    o=ma.masked_where(o==-9999.0,o).filled(0.0)
    s=ma.masked_where(o==-9999.0,s).filled(0.0)
    o=np.compress(o>0.0,o)
    s=np.compress(o>0.0,s)
    return abs((np.mean(s)-np.mean(o))/(np.mean(o)+1e-20))
#==========================================================
def BIAS(s,o):
    """
    Bias
    input:
        s: simulated
        o: observed
    output:
        Bias: Bias
    """
    s,o = filter_nan(s,o)
    o=ma.masked_where(o==-9999.0,o).filled(0.0)
    s=ma.masked_where(o==-9999.0,s).filled(0.0)
    o=np.compress(o>0.0,o)
    s=np.compress(o>0.0,s)
    # o=np.compress(o==-9999.0,o)
    # s=np.compress(o==-9999.0,s)
    return abs((np.mean(s)-np.mean(o)))
#====================================================================
def rAmplitude(s,o,syear=2009,eyear=2014):
    """
    Realtive Amplitude Difference
    input:
        s: simulated
        o: observed
    output:
        rAmplitude: Realtive Amplitude Difference
    """ 
    s,o = filter_nan(s,o)
    # o=ma.masked_where(o==-9999.0,o).filled(0.0)
    # s=ma.masked_where(o==-9999.0,s).filled(0.0)
    # o=np.compress(o>0.0,o)
    # s=np.compress(o>0.0,s)
    alti_diff=[]
    for year in np.arange(syear,eyear+1):
        date1 = datetime.date(year, 1, 1) # month and day are 1-base
        date2 = datetime.date(year, 12, 31)
        days  = (date2-date1).days + 1
        if year == syear:
            st_dt = 0
            ed_dt = days
        elif year == eyear:
            st_dt = ed_dt
            ed_dt = -1
        else:
            st_dt = ed_dt
            ed_dt = ed_dt + days
        # print st_dt, ed_dt
        maxloc = np.argmax(s[st_dt:ed_dt])
        minloc = np.argmin(s[st_dt:ed_dt])
        smax   = np.amax(s[st_dt:ed_dt])
        smin   = np.amin(s[st_dt:ed_dt])
        maxloc1=max(maxloc-15,0)
        maxloc2=min(maxloc+15,len(s)-1)
        minloc1=max(minloc-15,0)
        minloc2=min(minloc+15,len(s)-1)
        maxarry=ma.masked_equal(o[st_dt:ed_dt],-9999.0).filled(-9999.0)
        minarry=ma.masked_equal(o[st_dt:ed_dt],-9999.0).filled(9999.0)
        omax   = np.amax(maxarry[maxloc1:maxloc2])
        omin   = np.amin(minarry[minloc1:minloc2])
        if omax == -9999.0 or omin == 9999.0:
            continue
        alti_diff.append(((smax-smin)-(omax-omin))/((omax-omin)+1e-20))
    return np.mean(alti_diff)
#====================================================================
def dAmplitude(s,o,syear=2009,eyear=2014):
    """
    Amplitude Difference
    input:
        s: simulated
        o: observed
    output:
        dAmplitude: Amplitude Difference
    """ 
    s,o = filter_nan(s,o)
    sim_amp=[]
    obs_amp=[]
    for year in np.arange(syear,eyear+1):
        date1 = datetime.date(year, 1, 1) # month and day are 1-base
        date2 = datetime.date(year, 12, 31)
        days  = (date2-date1).days + 1
        if year == syear:
            st_dt = 0
            ed_dt = days
        elif year == eyear:
            st_dt = ed_dt
            ed_dt = -1
        else:
            st_dt = ed_dt
            ed_dt = ed_dt + days
        # max min simuulated
        smax   = np.amax(s[st_dt:ed_dt])
        smin   = np.amin(s[st_dt:ed_dt])
        sim_amp.append(smax-smin)
        # max min observed
        maxarry=ma.masked_equal(o[st_dt:ed_dt],-9999.0).filled(-9999.0)
        minarry=ma.masked_equal(o[st_dt:ed_dt],-9999.0).filled(9999.0)
        #===
        omax   = np.amax(maxarry)
        omin   = np.amin(minarry)
        if omax == -9999.0 or omin == 9999.0:
            continue
        maxloc = np.argmax(maxarry)
        minloc = np.argmin(minarry)
        obs_amp.append(omax-omin)
    # return abs(np.mean(np.array(sim_amp))-np.mean(np.array(obs_amp)))
    return np.mean(np.array(sim_amp))-np.mean(np.array(obs_amp))
#========================================
def dpeak_time(s,o,syear=2009,eyear=2014):
    """
    Peak Timing Difference
    input:
        s: simulated
        o: observed
    output:
        dpeak_time: Peak Timing Difference
    """ 
    s,o = filter_nan(s,o)
    sim_peak=[]
    obs_peak=[]
    for year in np.arange(syear,eyear+1):
        date1 = datetime.date(year, 1, 1) # month and day are 1-base
        date2 = datetime.date(year, 12, 31)
        days  = (date2-date1).days + 1
        if year == syear:
            st_dt = 0
            ed_dt = days
        elif year == eyear:
            st_dt = ed_dt
            ed_dt = -1
        else:
            st_dt = ed_dt
            ed_dt = ed_dt + days
        smaxloc= np.argmax(s[st_dt:ed_dt])
        sim_peak.append(smaxloc)
        maxarry=ma.masked_equal(o[st_dt:ed_dt],-9999.0).filled(-9999.0)
        omaxloc= np.argmax(maxarry)
        if omaxloc == None:
            continue
        obs_peak.append(omaxloc)
    # return abs(np.mean(np.array(sim_peak))-np.mean(np.array(obs_peak)))
    return np.mean(np.array(sim_peak))-np.mean(np.array(obs_peak))
#========================================
def ISS(l,u,o,alpha=0.05):
    """
    Interval Skill Score
    input:
        l: lower bound
        u: upper bound
    output:
        ISS: Interval Skill Score 
    """
    u,l = filter_nan(u,l)
    o=ma.masked_where(o<=0.0,o).filled(0.0)
    u=ma.masked_where(o<=0.0,u).filled(0.0)
    l=ma.masked_where(o<=0.0,l).filled(0.0)
    o=np.compress(o>0.0,o)
    u=np.compress(o>0.0,u)
    l=np.compress(o>0.0,l)
    pnum=len(o)
    issalpha=[]
    for i in np.arange(pnum):
        if l[i] <= o[i] <= u[i]:
            issalpha.append(u[i]-l[i])
        elif o[i] < l[i]:
            issalpha.append((u[i]-l[i])+(2.0/alpha)*(l[i]-o[i]))
        elif u[i] < o[i]:
            issalpha.append((u[i]-l[i])+(2.0/alpha)*(o[i]-u[i]))
        else:
            issalpha.append(0.0)
    issalpha=np.array(issalpha)
    return np.sum(issalpha)
#====================================================================
def sharpness(l,u,o):
    """
    Sharpness
    input:
        l: lower bound
        u: upper bound
    output:
        Sharpness 
    """
    u,l = filter_nan(u,l)
    o=ma.masked_where(o<=0.0,o).filled(0.0)
    u=ma.masked_where(o<=0.0,u).filled(0.0)
    l=ma.masked_where(o<=0.0,l).filled(0.0)
    o=np.compress(o>0.0,o)
    u=np.compress(o>0.0,u)
    l=np.compress(o>0.0,l)
    return np.sum(u-l)
#====================================================================
def reliability(l,u,o):
    """
    Reliability
    input:
        l: lower bound
        u: upper bound
    output:
        Reliability 
    """
    u,l = filter_nan(u,l)
    o=ma.masked_where(o<=0.0,o).filled(0.0)
    u=ma.masked_where(o<=0.0,u).filled(0.0)
    l=ma.masked_where(o<=0.0,l).filled(0.0)
    o=np.compress(o>0.0,o)
    u=np.compress(o>0.0,u)
    l=np.compress(o>0.0,l)
    pnum=len(o)
    rcount=0
    for i in np.arange(pnum):
        if l[i] <= o[i] <= u[i]:
            rcount=rcount+1
    return float(rcount)/float(pnum)
#====================================================================
# End
#====================================================================