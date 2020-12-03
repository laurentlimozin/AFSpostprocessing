#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 15:40:32 2020

@author: laurent
""""""

Created on Sun Sep 29 19:33:24 2019
@author: Laurent
#Open one tdms file created by Labview AFS Acquisition software\\
#Export trajectories and force vs time in readable csv format (for Charlie)\\
#Measure acquisition frequency, period and force for square signal\\
#Perform reference on a set of beads and anchor point\\
#Perform anchoring point on low or zero power periods\\
#Measure end-to-end distance (Length) with Claire's formula\\
#Detect jumps in square power AFS experiments on JDNA\\
#Draw and fit survival curves at constant force\\
#Compute maps and histogram of x, y, z at low and high force\\
#Displays all XY traj with number on the same graph\\
#Display 3D scatter plot\\
#Trace plot X,Y,Z vs T by number SingleTrace(n) or coordinates\\

20200322 v3.6   Tests different smoothing values NZlavg for anchor point
                Separation of histograms Open , Close, Low states and gaussian fits
                Power Spectrum in Open states and fit
20200402 v3.61  Automatic save of graphics and fit results with standardized name
                Faxen correction for z
20200617 v3.7   Extended information saved in csv file
                New representation of spectrum by binning frequencies
                Implemented stop frame
                Anchor point can be calculated on zero force interval; use AnchorPointState variable
                Faxen correction for x,y
20200620 v3.8   Use of step_detect gaussian1dfilter (copy step_detect.py in main folder)
         v3.81  Step visualization in FindJumps2() and graph saving
         v3.83  PullAngle2 calculated between high force and zero force state
20200703 v3.85  Correct bug for start / stop
         v3.86  Records absolute coordinates; calculate fit error on force
20200716 v3.87  Corection bugs for start/stop
20200826 v3.88  Correction to exclude Nan for Median XY position
20201009 v3.9   Permits to select the first long step for the spectrum calculation (set AnchorPointStateList=['0'] )
                Take into account the NoOpen cases in the survival curve
20201025 v3.91  Revised NoOpen/NoClose cases with introduction of a threshold MiddleOpenClose=900nm
                Some discrepancy for step_detectON=0 or 1 (which may depend on windowsize_detect and NZavg)
20201203 v4.0   Full revision of code
                Adapt to last version of nptdms package

"""""

import numpy as np
from nptdms import TdmsFile #from nptdms import tdms  # pip install nptdms
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit
from scipy import signal
from matplotlib.pyplot import cm
import pandas as pd
import os
import step_detect_2 as sd
from scipy.ndimage.filters import gaussian_filter1d

np.random.seed(444); plt.close('all')
import warnings; warnings.filterwarnings("ignore")
def pfl(f): return '{:.2E}'.format(f)

#path1="C:/DATA/Experiences Diverses/AFS/Data/Debug20201002/"      #☺ data folder
path1="/home/laurent/DATA/Experiences Diverses/AFS/Data/" #☺ data folder  #path1="D:/Jim/JDNA/2020-08-24/FOV3/"  
tdmsfilename="20200828-121251.tdms"          #tdmsfilename="20200903-144604.tdms"    tdmsfilename="20200826-120242.tdms"
OutputFolder='Outputtest/'
if not os.path.exists(path1+OutputFolder): os.makedirs(path1+OutputFolder)

wbref=[6,28,49]      # reference beads  wbref=[17,29,43]               # reference beads
b=[23]               # test beads  b=[6, 8, 21, 41, 46, 58, 73, 91]     #b=[58, 73, 91] b=[73, 91] b=[5]

#RangeNZlavg=[20,40,60, 80, 100, 120, 140, 160, 200, 240, 280, 320, 400, 480, 560, 640, 720, 800, 900, 1000]
RangeNZlavg=[400]       # size of moving average for smoothing z before anchor point determination

# initialiation for coordinate reading
nbead=len(b); fname=[""]*nbead; start=np.zeros(nbead, dtype=int); stop=np.zeros(nbead, dtype=int); load=np.zeros(nbead, dtype=int)
#for nf in range(nbead): fname[nf] = tdmsfilename; start[nf]=15000; stop[nf]=1600000; load[nf]=nf==0  # stop[nf]=1600000
for nf in range(nbead): fname[nf] = tdmsfilename; start[nf]=0.*60*1000/18; stop[nf]=175*60*1000/18; load[nf]=nf==0  # stop[nf]=1600000
range_anchorpoint=(0,100)
range_spectrum=(1000, 2000)

##========================= ========================= =========================
loaddata=True  # to load the tdms file
RefBeadZcorrection=True # correction with reference bead
RefBeadXYcorrection=True # correction with reference bead
AnchorPointStateList=['0'] # correction with AnchorPoint as determined by anchor point  'low' or '0' force ; or 'None': no anchor point 
export=False # to produce csv files of raw data (for Charlie)
GraphTdms=False  # graph of frame rate and force data (duration and amplitude of 2 steps)
GraphDhisto=True; GraphTupdown=False; GraphTZ_Power=False
HistoZ=True; GraphXYZ=False # characteristics of trajectory of the selected beads
pprint=True; # printing detailed infos in console window
GraphXYall=False; iminXY=0; imaxXY=0 # display xy trajectories with number on a single graph
SaveGraph=True; OutFormat=".jpg" ; CloseAfterSave=False # alternative ".jpg" or ".pdf" 
DisplaySpectrum=True
DisplayCloseOpen=True
DisplayGaussianFits=False
DisplayJumpGraphs=True # graphs of jump analysis
DrawSuperposal=False # graph of superposition of jump traces
CalculateSurvival=False  # calculate and draw survival curve
tshiftsurvival=0; 
DisplaySurvival=True
Display3D=False
DisplayplotvariableNZlavg=True
step_detectON=True; plot_detect=False; nbplot_detect=15; thres_detect=0.5; windowsize_detect=20 # parameters for gaussian step detection
NZavg=100; NZlavg=400; dzjump=80.; NZrefavg=5000; tol=0.01; sample=1; sampleGraph=100 # parameters for direct step detection
MiddleOpenClose=900
BeadRad=790. # Bead Radius in nm
fmin=0.1; fmax=15.   # frequency range for spectrum fit (Hz)
p1Lo0=1000; p1Lc0=800; p1Ll0=500; p1Zl0=500  # guess for histogram fit

##========================= ========================= =========================

def track(d0, i):    # retrieves coordinates / data of bead i
    X = d0['ROI{:04d}'.format(i)+' X (nm)'][:]
    Y = d0['ROI{:04d}'.format(i)+' Y (nm)'][:]
    Z = d0['ROI{:04d}'.format(i)+' Z (nm)'][:]
    MinLUT = d0['ROI{:04d}'.format(i)+' MinLUT'][:]
    SSIM = d0['ROI{:04d}'.format(i)+' SSIM'][:]
    return (X,Y,Z, MinLUT, SSIM)

def events(d2):  # retrieve events from tdms file Time/Name/Data
    T=d2['Time (ms)'][:]
    N=d2['Name'][:]
    D=d2['Data'][:]
    return (T,N,D)

def BuildP(T, TE, NE, DE):      # builds power time trace form event times
    I=(NE=='SigGen.Power (%)')
    TP=TE[I]; DP=DE[I]
    print('BuildP:', len(TP), 'events') #   print(TP)
    P=np.zeros(len(T)); j1=0
    for i in range(1,len(TP)-1):
   #     j=int(np.where(T==TP[i])[0])
        j=int(np.where(T==TP[i])[0][0])
        P[j1:j]=DP[i-1]; j1=j
    return P

def time(d0, Graph):    # retrieves time trace from tdms file and measures time step
    T = d0['Time (ms)'][:]
    print('Total points:', len(T))
    dT=np.gradient(T)
    if Graph:
        figname=refname+'dTrawhisto'
        plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
        ax=sns.distplot( dT, kde=False, bins=np.linspace(0.,20.,num=200) )
        ax.set_xlabel("Time step (ms)");  ax.set_ylabel("Number"); ax.set_yscale('log')
        plt.ylim(0.1, 1.e7)
        if SaveGraph: plt.savefig(path1+OutputFolder+figname+OutFormat)
#    mostfreq=np.argmax(np.bincount(dT))
    print('dTrawhisto (frame time ms)', np.amin(dT), np.amax(dT), '  frequency (Hz)', 1000/np.amax(dT))   
    return T, 1000/np.amax(dT)

def shortentrack(X0, Y0, Z0, T0, P0, MinLUT0, i0,j0):   # cut trajectories between 0 and i0  # TO OPTIMIZE !
    X=np.delete(X0, range(i0));  X=np.delete(X, range(j0-i0,len(X0)-i0))
    Y=np.delete(Y0, range(i0));  Y=np.delete(Y, range(j0-i0,len(Y0)-i0))
    Z=np.delete(Z0, range(i0));  Z=np.delete(Z, range(j0-i0,len(Z0)-i0))
    T=np.delete(T0, range(i0));  T=np.delete(T, range(j0-i0,len(T0)-i0))
    P=np.delete(P0, range(i0));  P=np.delete(P, range(j0-i0,len(P0)-i0))
    MinLUT=np.delete(MinLUT0, range(i0));  MinLUT=np.delete(MinLUT, range(j0-i0,len(MinLUT0)-i0))
    return (X,Y,Z,T,P, MinLUT)

def period(TE, Graph, refname):  # retrieves low power and high power duration in case of binary power
    dTE=np.gradient(TE)
    dTE=-TE[:-1] + TE[1:]
    if Graph:
        figname=refname+'dTEventsHisto'
        plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
        ax=sns.distplot( dTE, kde=False, bins=np.linspace(0.,200000.,num=200) )
        ax.set_xlabel("Event Time step (ms)");  ax.set_ylabel("Number"); ax.set_yscale('log')
        plt.ylim(0.1, 1.e3)
        if SaveGraph: plt.savefig(path1+OutputFolder+figname+OutFormat)
    bins=np.linspace(0, 200000 ,num=1001)
    hist, bin_edges = np.histogram(dTE,bins=bins)
    hist[np.argmax(hist)]=0
    TE1=bin_edges[np.argmax(hist)]
    hist[np.argmax(hist)]=0
    TE2=bin_edges[np.argmax(hist)]
    TEhigh=max(TE1,TE2)
    TElow=min(TE1,TE2)
    print('TEhigh (ms):', TEhigh, "TElow (ms):", TElow)
    return TEhigh, TElow, dTE

def MakeHistoP(P, Graph):       # histogram of power trace and retrieval of low/high values for binary power
    bins=np.linspace(0, 40 ,num=40001)
    if Graph:
        figname=refname+'Phisto'
        plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
        ax=sns.distplot( P, kde=False, bins=bins )
        ax.set_xlabel("Power (%)");  ax.set_ylabel("Number")
        if SaveGraph: plt.savefig(path1+OutputFolder+figname+OutFormat)
    hist, bin_edges = np.histogram(P,bins=bins)    #   np.bincount(P)
    Phigh=bin_edges[np.argmax(hist)]
    hist[np.argmax(hist)]=0
    Plow=bin_edges[np.argmax(hist)]
    if Plow==0:
        hist[np.argmax(hist)]=0
        Plow=bin_edges[np.argmax(hist)]
    print('Phigh:', Phigh, "Plow:", Plow)
    return Phigh, Plow, hist, bin_edges

def Rollavg(Z, N, mode):      # rolling average on N time steps
    Zavg=np.convolve(Z, np.ones((N,))/N, mode=mode)
    n=int(N/2.)
    Zavg[range(n)]=Zavg[n+1]        # np.nan
    Zavg[range(-1,-n-1,-1)]= Zavg[-n]    #   np.nan
    return Zavg

def MakeGraphXYall(d0, imax, n, refrms, name, corr=False, Xr0avg=None, Yr0avg=None):   
# Displays all XY traj with number and SD on the same graph
    if corr: name=name+" Corr"
    figname='XYall_'+name
    plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
    xmin=0; xmax=350000; ymin=0; ymax=250000
    for i in range(1,imax):
        (X,Y,Z,MinLUT,SSIM255)=track(d0, i)
        if corr: X=X-Xr0avg; Y=Y-Yr0avg
        xx=X[::n]; yy=Y[::n]
        ds=np.sqrt((X-X.mean())**2+(Y-Y.mean())**2)
        dsmean=ds.mean()
        if np.isnan(dsmean): lab=str(i)+'/'+'Nan'
        else: lab=str(i)+'/'+str(int(dsmean))
        if dsmean<refrms: lab=lab+'ref'
        print('plot xy traj:'+lab)
        ax.scatter(xx, yy, marker='.', alpha=0.5, label=str(i))
        plt.text(xx[-1], yy[-1], lab, fontsize=6)
    ax.axis([xmin, xmax, ymin, ymax])
    if corr: ax.axis([xmin-Xr0avg.mean(), xmax-Xr0avg.mean(), ymin-Yr0avg.mean(), ymax-Yr0avg.mean()])
    if SaveGraph: plt.savefig(path1+OutputFolder+figname+OutFormat)

def exporttrack(X,Y,Z,T,P, datasave):    # exports raw traces as csv file
    file = open(datasave, "w")
    file.write('X_nm'+','+'Y_nm'+','+'Z_nm'+','+'T_ms'+','+'P_%'+'\n')
    for i in range(len(X)):
        file.write(str(X[i])+','+str(Y[i])+ ','+str(Z[i])+ ','+str(T[i])+ ','+str(P[i])+"\n")
    file.close()
    
def Plot3D(x,y,z,n, refname):    # 3D plot
    fig = plt.figure(refname, figsize=(7,7), dpi=100)
    xx=x[::n]; yy=y[::n]; zz=z[::n]
    ax = fig.add_subplot(111, projection='3d')
    # For each set of style and range settings, plot n random points in the box
    # defined by x in [23, 32], y in [0, 100], z in [zlow, zhigh].
    for m, zlow, zhigh in [('.', -50, -25), ('.', -30, -5)]:
        ax.scatter(xx, yy, zz, marker=m)
    if SaveGraph: plt.savefig(path1+OutputFolder+refname+OutFormat)
    
def MakeHistoZ(Zh, Zl, labelZ1, labelZ2, zmin, zmax, refname):  # histogram of 2 Z traces for comparison
    figname=refname+'Zhisto'
    plt.figure(figname, figsize=(4,4), dpi=100); ax = plt.gca()
    if zmin==0 and zmax==0:
        zmin=min(np.amin(Zh), np.amin(Zl))
        zmax=max(np.amax(Zh), np.amax(Zl))
    bins=np.linspace(zmin,zmax,num=200)
    ax=sns.distplot( Zh, kde=False, bins=bins, label=labelZ1 )
    ax=sns.distplot( Zl, kde=False, bins=bins, label=labelZ2 )
    ax.set_xlabel("Z (nm)");  ax.set_ylabel("Number")
    plt.xlim(zmin,zmax); ax.legend(); ax.set_yscale('log')
    print(labelZ1, len(Zh), labelZ2, len(Zl))
    if SaveGraph: plt.savefig(path1+OutputFolder+figname+OutFormat)
    
def MakeGraphXY(X1, Y1, X2, Y2, xmin, xmax, ymin, ymax, label1, label2, refname, line1=False, line2=False, log=False):
    figname=refname+'_XY'
    plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
    if line1: 
        ax.plot(X1, Y1, alpha=0.5, label=label1)
    else:
        ax.scatter(X1, Y1, marker='.', alpha=0.1, label=label1)
#    ax=sns.kdeplot(X1,Y1, cmap="Reds", shade=True, bw=.15, legend=False)
    if line2:
        ax.plot(X2, Y2, alpha=0.5, label=label2)
    else:
        ax.scatter(X2, Y2, marker='.', alpha=0.1, label=label2)
    if log: ax.set_yscale('log')
    ax.legend(fontsize=6); ax.axis([xmin, xmax, ymin, ymax])
    if SaveGraph: plt.savefig(path1+OutputFolder+figname+OutFormat)

def SingleTrace(d0, n,m, X,Y,Z,T, MinLUT):    # SingleTrace(30,10, 0,0,0,0,0)
# Traces time traces and histograms for bead n (step m) or coordinates if n<=0
    if n>0:
        (X,Y,Z,MinLUT,SSIM)=track(d0, n); lab=str(n)        # test bead
        T = tdms_file.object('Tracking data', 'Time (ms)').data
    T_s=(1/1000.)*T
    if n==0: lab='Corr by ref'
    if n==-1: lab='ref'
    if n==-2: lab='After anchor point'
    figname=refname+'SingleTrace'+lab
    fig=plt.figure(figname, figsize=(9,6), dpi=100)
    spec = fig.add_gridspec(ncols=2, nrows=2, width_ratios=[3, 1], height_ratios=[1, 1])
    print('Averages: X=', X.mean(), 'Y=', Y.mean(), 'Z=', Z.mean())
    ax = fig.add_subplot(spec[0, 0])#    fig.add_subplot(2,1,1); ax = plt.gca()
    ax.scatter(T_s[::m], (X-X.mean())[::m], marker='.', alpha=0.5, label='X')
    ax.scatter(T_s[::m], (Y-Y.mean())[::m], marker='.', alpha=0.5, label='Y')
    ax.legend(fontsize=6); ax.axis([np.amin(T_s), np.amax(T_s), -1000, 1000])
    ax.set_xlabel("time (s)"); ax.set_ylabel("X,Y (nm)")
    if np.isnan(X.mean()) or np.isnan(Y.mean()): labx='Nan'; laby='Nan'
    else: labx=str(int(X.mean())); laby=str(int(Y.mean()))
    ax.annotate('Xavg='+labx+'/Yavg='+laby , (0.1, 0.5), xycoords='axes fraction', va='center')
    ax = fig.add_subplot(spec[0, 1])
    ax=sns.distplot( X-X.mean(), kde=False, bins=np.linspace(-1000,1000,num=200) )
    ax=sns.distplot( Y-Y.mean(), kde=False, bins=np.linspace(-1000,1000,num=200) )
    ax.set_xlabel("X,Y (nm)"); ax.set_yscale('log'); plt.ylim(1, 1.e5)#;  ax.set_ylabel("Number")
    ax = fig.add_subplot(spec[1, 0]); ax2=ax.twinx()
    ax.scatter(T_s[::m], (Z-Z.mean())[::m], marker='.', alpha=0.5, label='Z')
    ax2.scatter(T_s[::m], MinLUT[::m], marker='.', color='k', alpha=0.5, label='MinLUT')
    ax.legend(fontsize=6, loc='upper left'); ax.axis([np.amin(T_s), np.amax(T_s), -1000, 1000])
    ax2.set_ylim(0., 20.); ax2.legend(fontsize=6, loc='upper right'); ax2.set_ylabel('MinLUT')
    ax.set_xlabel("time (s)"); ax.set_ylabel("Z-Zavg (nm)") 
    if np.isnan(Z.mean()): labz='Nan'
    else: labz=str(int(Z.mean()))
    ax.annotate('Zavg='+labz , (0.1, 0.5), xycoords='axes fraction', va='center')
    ax = fig.add_subplot(spec[1, 1])
    ax=sns.distplot( Z, kde=False, bins=np.linspace(Z.mean()-1000,Z.mean()+1000,num=200) )
    ax.set_xlabel("Z (nm)");  ax.set_ylabel("Number"); ax.set_yscale('log'); plt.ylim(1, 1.e5)
    if SaveGraph: plt.savefig(path1+OutputFolder+figname+OutFormat)
    
##========================= ========================= =========================
##=========================   # MAIN LOOP  # ==================================
listresultsTotal=[]
color=iter(cm.rainbow(np.linspace(0, 1, 4*nbead*len(RangeNZlavg))))

if nbead>1: figAllSpec, axAllSpec = plt.subplots(nbead, 3, num='figAllSpec', figsize=(6, 2*nbead), dpi=100)

for ibead in range(nbead):        ## Loop on test beads
    print('=============================================================')
    print('NZavg=',NZavg, 'dzjump=', dzjump, 'tol=',tol, 'NZrefavg=',NZrefavg)
    if loaddata and load[ibead]: print('Loading data', fname[ibead]); tdms_file = TdmsFile.open(path1+fname[ibead])  # alternative TdmsFile.read(path1+fname[ibead])
    if loaddata and not load[ibead]: print('Data already loaded', fname[ibead])
    tdms_groups = tdms_file.groups(); # if pprint: print(tdms_groups)
    tdms_track = tdms_groups[0]; tdms_temp = tdms_groups[1]; tdms_event = tdms_groups[2]
    dmax=1000000.; imax=(len(tdms_track)-1)//5; dX=1000; dY=1000    
    print('Nb of traj=', imax)
    refname=fname[ibead]+'_Bead'+str(b[ibead])+'_'; print(refname)
    if GraphXYall and ibead==0: MakeGraphXYall(tdms_track, imax, 100, 100, fname[ibead])

    (X0,Y0,Z0,MinLUT0,SSIM0)=track(tdms_track, b[ibead])        # test bead
    X0med=np.median(X0[~np.isnan(X0)]); Y0med=np.median(Y0[~np.isnan(Y0)])
    if load[ibead]:
        print("Calculating pool of reference bead(s)", wbref )
        if len(wbref)>1:
            (Xr0,Yr0,Zr0,MinLUTr0,SSIMr0)=track(tdms_track, wbref[0])
            for i in range(1,len(wbref)):
                (Xr0_,Yr0_,Zr0_,MinLUTr0_,SSIMr0_)=track(tdms_track, wbref[i])
                Xr0=Xr0+Xr0_; Yr0=Yr0+Yr0_; Zr0=Zr0+Zr0_
            Xr0=Xr0/len(wbref); Yr0=Yr0/len(wbref); Zr0=Zr0/len(wbref)
        else:
            (Xr0,Yr0,Zr0,MinLUTr0,SSIMr0)=track(tdms_track, wbref[0])   # reference bead
    else: print("Pool of reference bead(s) already calculated") 
    (TE0,NE0,DE0)=events(tdms_event)
    T0, fs = time(tdms_track, GraphTdms)
    P0=BuildP(T0, TE0, NE0, DE0)
    if GraphTdms:
        plt.figure('Power', figsize=(6,6), dpi=100); ax = plt.gca()
        ax.scatter(T0, P0, marker='.', c='b', alpha=0.5)
 
    # make sure that chosen time interval does not start or stop at high power        
    print('Power at start', start[ibead],':', P0[start[ibead]])
    print('Power at stop:', stop[ibead],':', P0[stop[ibead]])
    if P0[start[ibead]]>0.1:
        wstart=10*np.arange(1000)+start[ibead]; wPstart=P0[wstart]
        i0=np.argmin(wPstart); start[ibead]= wstart[i0+5]
    if P0[stop[ibead]]>0.1:
        wstop=10*np.arange(1000)+stop[ibead]; wPstop=P0[wstop]
        i0=np.argmin(wPstop); stop[ibead]= wstop[i0+5]
    print('Power at new start', start[ibead],':', P0[start[ibead]])
    print('Power at new stop:', stop[ibead],':', P0[stop[ibead]])

    print('Tracks shortening from frame', start[ibead], 'to ',stop[ibead])
    # shorten test bead
    (X,Y,Zu,T,P, MinLUT)= shortentrack(X0, Y0, Z0, T0, P0, MinLUT0, start[ibead], stop[ibead])
    # shorten reference bead
    (Xr,Yr,Zr,Tr,Pr, MinLUTr)= shortentrack(Xr0, Yr0, Zr0, T0, P0, MinLUTr0, start[ibead], stop[ibead])
     
    mask=(TE0>=T0[start[ibead]])&(TE0<T0[stop[ibead]])
    TE=TE0[mask]; NE=NE0[mask]; DE=DE0[mask];
    
    TEhigh, TElow, dTE=period(TE, GraphTdms, refname)        
        
    Phigh, Plow ,hist,bin_edges= MakeHistoP(P, GraphTdms)
    # https://stackoverflow.com/questions/13728392/moving-average-or-running-mean
    mode='same'     # edge modes = ['full', 'same', 'valid']
    Zravg=Rollavg(Zr, NZrefavg, mode)        # rolling average of shorten reference bead
    Zr0avg=Rollavg(Zr0, NZrefavg, mode)       # rolling average of reference bead
    if RefBeadZcorrection:
        print('Zcorrection by reference bead(s): Zravg=', Zravg[0] )#int(bref[ibead]) )
        Zc=Zu-Zravg; Z0c=Zr0-Zr0avg; Zec=Z0-Zr0avg
    else: print('No Zcorrection'); Zc=Zu; Z0c=Zr0; Zec=Z0
    Xr0avg=Rollavg(Xr0, NZrefavg, mode); Yr0avg=Rollavg(Yr0, NZrefavg, mode) 
    Xravg=Rollavg(Xr, NZrefavg, mode); Yravg=Rollavg(Yr, NZrefavg, mode) 
    if RefBeadXYcorrection:
         print('XYcorrection by reference bead(s): Xravg=', Xravg[0],' Yravg=', Yravg[0] )
         Xc=X-Xravg; Yc=Y-Yravg; Xec=X0-Xr0avg; Yec=Y0-Yr0avg
    else: print('No XYcorrection'); Xc=X; Yc=Y; Xec=X0; Yec=Y0

    if GraphXYall and ibead==0: MakeGraphXYall(tdms_track, imax, 100, 100, fname[ibead], True, Xr0avg, Yr0avg)

    indnonan=(~np.isnan(Xc))*(~np.isnan(Yc))
    Xc=Xc[indnonan]; Yc=Yc[indnonan]; Zc=Zc[indnonan]
    MinLUT=MinLUT[indnonan]; T=T[indnonan]; P=P[indnonan]
    listresults=[]
    if export:    # exports raw data for third party analysis
        exporttrack(X,Y,Zu,T,P, path1+'Export_'+refname+'.csv')
        exporttrack(Xr,Yr,Zr,Tr,P, path1+'Export_'+refname+'_r'+str(wbref)+'.csv')        
       
    indhigh=(P>=Phigh*(1-tol))*(P<=Phigh*(1+tol))
    indlow=(P>=Plow*(1-tol))*(P<=Plow*(1+tol)); #indlow=(P<Phigh*0.9)+(P>Pmax)
    ind0force=(P==0)
    indeforce0=(P0==0)
    maskcustom=(TE0>=T0[range_anchorpoint[0]])&(TE0<T0[range_anchorpoint[1]])
        
    p1Lo=p1Lo0; p1Lc=p1Lc0; p1Ll=p1Ll0; p1Zl=p1Zl0  # initialization of gaussian center for fit guess
    for NZlavg in RangeNZlavg:      # test loop for various choices of NZlavg
        print('******* NZlavg=', NZlavg, '*******')
        Zcavg=np.convolve(Zc, np.ones((NZlavg,))/NZlavg, mode=mode); dZ=Zc-Zcavg
        Zecavg=np.convolve(Zec, np.ones((NZlavg,))/NZlavg, mode=mode); dZe=Zec-Zecavg
         
        x = np.random.randn(len(P))/10
        Zch=Zc[indhigh]; Xch=Xc[indhigh]; Ych=Yc[indhigh]; dZh=dZ[indhigh]; Ph=P[indhigh]+x[indhigh]
        Zcl=Zcavg[indlow]; Xcl=Xc[indlow]; Ycl=Yc[indlow]; dZl=dZ[indlow]; Pl=P[indlow]+x[indlow]
        Zc0=Zcavg[ind0force]; Xc0=Xc[ind0force]; Yc0=Yc[ind0force]; dZ0=dZ[ind0force]
        Zec0=Zecavg[indeforce0]; Xec0=Xec[indeforce0]; Yec0=Yec[indeforce0]; dZe0=dZe[indeforce0]
        Zec0custom=Zecavg[maskcustom]; Xec0custom=Xec[maskcustom]; Yec0custom=Yec[maskcustom]; dZe0custom=dZe[maskcustom]
              
        indstart=(T0<start[ibead]); Z0start=Z0c[indstart]
        Nl=len(Zc[indlow]); Nh=len(Zc[indhigh]); N0=len(Zc[ind0force]); Ntot=len(Zc)
        print('States Populations: low=', Nl, ' high=', Nh, ' 0force=', N0, ' total=', Ntot)        
        print('States Fractions: low=', Nl/Ntot, ' high=', Nh/Ntot, ' 0force=', N0/Ntot, 'total=', (Nl+Nh+N0)/Ntot)        
        for AnchorPointState in AnchorPointStateList:      # ['low','0','custom']     
            print('Anchor Point Reference Power:', AnchorPointState)
            if AnchorPointState=='None':
                print('No anchoring point'); Zs=Zc; Xs=Xc; Ys=Yc
            else:
                if AnchorPointState=='low':
                    Zc_=Zcl; Xc_=Xcl; Yc_=Ycl
                elif AnchorPointState=='0':
                    NZlavg+=10
                    Zc_=Zc0; Xc_=Xc0; Yc_=Yc0   
                    Zc_=Zec0; Xc_=Xec0; Yc_=Yec0
                elif AnchorPointState=='custom':
                    Zc_=Zec0custom; Xc_=Xec0custom; Yc_=Yec0custom
                Zc_sort=np.sort(Zc_)       #     cleanedZclsort = [x for x in Zclsort if str(x) != 'nan']
                if len(Zc_sort)>0: Za=np.median(Zc_sort[np.arange(sample)])  # Za=np.amin(Zcl) #-BeadRad
                else:
                    print('Zc_ null wave:  Bad Anchor Point Reference Choice')
                    break
#                print('Anchoring point at', AnchorPointState, 'Force. Length calculation', 'min=', pfl(np.min(Zc_)), 'min_',sample, pfl(Za))
                Xa=np.median(Xc_); Ya=np.median(Yc_)     
                CA=np.sqrt((Xc-Xa)**2+(Yc-Ya)**2+(Zc-Za)**2)
                Zs=Zc-Za; Xs=Xc-Xa; Ys=Yc-Ya
            Xsl=Xs[indlow]; Ysl=Ys[indlow]; Zsl=Zs[indlow]; Tl=T[indlow]; MinLUTl=MinLUT[indlow]
            Xsh=Xs[indhigh]; Ysh=Ys[indhigh] ; Zsh=Zs[indhigh]; Th=T[indhigh]; MinLUTh=MinLUT[indhigh]
            Xs0=Ys[ind0force]; Ys0=Ys[ind0force]; Zs0=Zs[ind0force]       
            Length=np.sqrt(Xs**2+Ys**2+(Zs+BeadRad)**2)-BeadRad; Lengthl=Length[indlow]; Lengthh=Length[indhigh]; Length0=Length[ind0force]
            PullAngle=(180/np.pi)*np.arcsin( np.sqrt( (np.median(Xsh)-np.median(Xsl))**2 + (np.median(Ysh)-np.median(Ysl))**2  )/np.median(Lengthh))
            PullAngle2=(180/np.pi)*np.arcsin( np.sqrt( (np.median(Xsh)-np.median(Xs0))**2 + (np.median(Ysh)-np.median(Ys0))**2  )/np.median(Lengthh))            
            print('High length=', pfl(np.median(Lengthh)), 'Low length=', pfl(np.median(Lengthl)), ' Pulling Angle=', pfl(PullAngle), ' Pulling Angle2=', pfl(PullAngle2))
            if Display3D: Plot3D(Xs,Ys,Zs, 30, refname+'XsYsZs')
            if HistoZ:
        #        MakeHistoZ(Zcl, Z0start, 'Zcl', 'Z0start', 0,0, refname+'Z0startZcl')
        #        MakeHistoZ(Zcl, Zch, 'Zcl', 'Zch', 0,0, refname+'ZchZcl')
        #        MakeHistoZ(dZl, dZh, 'dZl', 'dZh', -200, 200, refname+'dZhdZl')
        #        MakeHistoZ(Zc, Zu, 'Zc', 'Zu', 0, 0, refname+'ZcZu')
        #        MakeHistoZ(Zu, Zr, 'Zu', 'Zr', 0,0, refname+'ZuZr')
             #   MakeHistoZ(Zc, Zs, 'Zc', 'Zs', -100,1200, refname+'ZcZs')
                MakeHistoZ(Lengthh, Lengthl, 'Lengthh', 'Lengthl', -100,1200, refname+'LengthhLengthl')
        #        MakeHistoZ(Length2h, Length2l, 'Lengthh', 'Length2l', 0,0, refname+'LengthhLengthl')
            if GraphXYZ:
                MakeGraphXY(Xsl, Ysl, Xsh, Ysh, -dX, dX, -dY, dY, 'XYl', 'XYh', refname+'XYslXYsh')
        #        MakeGraphXY(Pl, Zcl, Ph, Zch, -1, 5,-500, 2000, 'PZcl', 'PZch', refname+'PZclPZch')
        #        MakeGraphXY(Pl, Zsl, Ph, Zsh, -1, 5,-500, 2000, 'PZsl', 'PZsh', refname+'PZslPZsh')
        #        MakeGraphXY(T, Zc, T, Zs, 0, 6500000, -500, 2000, 'Zc', 'Zs', refname+'TZcZs')
                MakeGraphXY(Tl, Zsl, Th, Zsh, 0, T[-1], -100, 1200, 'Xs', 'Ys', refname+'TZslZsh')
        #        MakeGraphXY(Zsl, MinLUTl, Zsh, MinLUTh, -100, 1200, 0, 12, 'MinLUTl', 'MinLUTh', refname+'ZslLUTlZshLUTh')
        #        MakeGraphXY(dZl, MinLUTl, dZh, MinLUTh, -300, 300, 0, 12, 'MinLUTl', 'MinLUTh', refname+'dZslLUTldZshLUTh')
        #        MakeGraphXY(Lengthl, MinLUTl, Lengthh, MinLUTh, -100, 1000, 0, 12, 'MinLUTl', 'MinLUTh', refname+'lengthlLUTllengthhLUTh')
        #        MakeGraphXY(Length2l, MinLUTl, Length2h, MinLUTh, -100, 1000, 0, 12, 'MinLUTl', 'MinLUTh', refname+'lengthlLUTllengthhLUTh') 
        #        SingleTrace(tdms_track,0, sampleGraph, Xc, Yc, Zc, T, MinLUT )
        #        SingleTrace(tdms_track,-1,sampleGraph, Xr, Yr, Zr, T, MinLUT )
                SingleTrace(tdms_track, -2,sampleGraph, Xs, Ys, Zs, T, MinLUT )
                MakeGraphXY(T, Length, T, Zs, 0, T[-1], 0, 1200, 'Length', 'Zs', refname+'Length')
 
