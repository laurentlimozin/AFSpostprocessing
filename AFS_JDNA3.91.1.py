# -*- coding: utf-8 -*-
"""
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
"""
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
import warnings
warnings.filterwarnings("ignore")
def pfl(f): return '{:.2E}'.format(f)

path1="C:/Data/"
path1="C:/DATA/Experiences Diverses/AFS/Data/Debug20201002/"      #☺ data folder
path1="/home/laurent/DATA/Experiences Diverses/AFS/Data/"#"Debug20201002/"      #☺ data folder
#path1="D:/Jim/JDNA/2020-08-24/FOV3/"      #☺ data folder
tdmsfilename="20200903-144604.tdms"
tdmsfilename="20200826-120242.tdms"
tdmsfilename="20200828-121251.tdms"
#tdmsfilename="20200825-070046.tdms"
#path1="E:/Jim/AFS/20200625/"      #☺ data folder
#tdmsfilename="20200625-204554.tdms"
OutputFolder='Outputtest/'   # please create it first in he data folder
#OutputFolder=path1+'Outputtest20200701/'   # please create it first in he data folder
if not os.path.exists(path1+OutputFolder): os.makedirs(path1+OutputFolder)

wbref=[17,29,43]               # reference beads
wbref=[6,28,49]               # reference beads
b=[6, 8, 21, 41, 46, 58, 73, 91]    # test beads
#b=[58, 73, 91]
b=[73, 91]
b=[5]
b=[23]

RangeNZlavg=[20,40,60, 80, 100, 120, 140, 160, 200, 240, 280, 320, 400, 480, 560, 640, 720, 800, 900, 1000]
RangeNZlavg=[400]       # size of moving average for smoothing z before anchor point determination

# initialiation for coordibnate reading
nbead=len(b); fname=[""]*nbead; start=np.zeros(nbead, dtype=int); stop=np.zeros(nbead, dtype=int); load=np.zeros(nbead, dtype=int)
#PP=np.zeros(nbead, dtype=float); FpN=np.zeros(nbead, dtype=float)
#for nf in range(nbead): fname[nf] = tdmsfilename; start[nf]=15000; stop[nf]=1600000; load[nf]=nf==0  # stop[nf]=1600000
for nf in range(nbead): fname[nf] = tdmsfilename; start[nf]=20.1*60*1000/18; stop[nf]=175*60*1000/18; load[nf]=nf==0  # stop[nf]=1600000
for nf in range(nbead): fname[nf] = tdmsfilename; start[nf]=0.*60*1000/18; stop[nf]=175*60*1000/18; load[nf]=nf==0  # stop[nf]=1600000

#ForceCalibration=0.39881  # pN per %power
#for i in range(nbead): FpN[i]=ForceCalibration*PP[i]

##========================= ========================= =========================
loaddata=True  # to load the tdms file
oldtdsmspackage=False
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
CalculateSurvival=True  # calculate and draw survival curve
tshiftsurvival=0; 
DisplaySurvival=True
Display3D=False
DisplayplotvariableNZlavg=True
step_detectON=True; plot_detect=True; nbplot_detect=15; thres_detect=0.5; windowsize_detect=20 # parameters for gaussian step detection
NZavg=100; NZlavg=400; dzjump=80.; NZrefavg=5000; tol=0.01; sample=1; sampleGraph=100 # parameters for direct step detection
MiddleOpenClose=900
BeadRad=790. # Bead Radius in nm
fmin=0.1; fmax=15.   # frequency range for spectrum fit (Hz)
p1Lo0=1000; p1Lc0=800; p1Ll0=500; p1Zl0=500  # guess for histogram fit

##========================= ========================= =========================

def SingleTrace(n,m, X,Y,Z,T, MinLUT):    # SingleTrace(30,10, 0,0,0,0,0)
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

def printinfo(tdms_file):   # displays info about tdms file in console
    root_object = tdms_file.object()
 # Iterate over all items in the properties dictionary and print them
    for name, value in root_object.properties.items():
        print("{0}: {1}".format(name, value))
    print(type(d0), len(d0))
    print(d0[0], d0[1], d0[2], d0[3], d0[4], d0[5])
    print(type(d1), len(d1))
    print(d1[0], d1[1])
    print(type(d2), len(d2))
    print(d2[0], d2[1], d2[2])

def track(d0, i):    # retrieves coordinates / data of bead i
    X = d0['ROI{:04d}'.format(i)+' X (nm)'][:]
    Y = d0['ROI{:04d}'.format(i)+' Y (nm)'][:]
    Z = d0['ROI{:04d}'.format(i)+' Z (nm)'][:]
    MinLUT = d0['ROI{:04d}'.format(i)+' MinLUT'][:]
    SSIM = d0['ROI{:04d}'.format(i)+' SSIM'][:]
    return (X,Y,Z, MinLUT, SSIM)

def shortentrack(X0, Y0, Z0, T0, P0, MinLUT0, i0,j0):   # cut trajectories between 0 and i0
    X=np.delete(X0, range(i0));  X=np.delete(X, range(j0-i0,len(X0)-i0))
    Y=np.delete(Y0, range(i0));  Y=np.delete(Y, range(j0-i0,len(Y0)-i0))
    Z=np.delete(Z0, range(i0));  Z=np.delete(Z, range(j0-i0,len(Z0)-i0))
    T=np.delete(T0, range(i0));  T=np.delete(T, range(j0-i0,len(T0)-i0))
    P=np.delete(P0, range(i0));  P=np.delete(P, range(j0-i0,len(P0)-i0))
    MinLUT=np.delete(MinLUT0, range(i0));  MinLUT=np.delete(MinLUT, range(j0-i0,len(MinLUT0)-i0))
#    print('before',len(X0), 'after', len(X))
    return (X,Y,Z,T,P, MinLUT)

def events(d2):  # retrieve events from tdms file Time/Name/Data
    T=d2['Time (ms)'][:]
    N=d2['Name'][:]
    D=d2['Data'][:]
    # N=np.array(tdms_file.object('Events', 'Name').data)
    # D=tdms_file.object('Events', 'Data').data
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

def MakeHistoP(P, Graph):       # histogram of power trace and retrieval of low/high values for binary power
    bins=np.linspace(0, 40 ,num=40001)
    if Graph:
        figname=refname+'Phisto'
        plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
  #      ax=sns.distplot( P, kde=False, bins=np.linspace(0, 10 ,num=1000) )
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

def time(Graph):    # retrieves time trace from tdms file and measures time step
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

def period(TE, Graph, refname):  # retrieves low power and high power duration in case of binary power
    dTE=np.gradient(TE)
  #  np.count_nonzero(x[:-1] < x[1:])
    dTE=-TE[:-1] + TE[1:]
    if Graph:
        figname=refname+'dTEventsHisto'
        plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
        ax=sns.distplot( dTE, kde=False, bins=np.linspace(0.,200000.,num=200) )
        ax.set_xlabel("Event Time step (ms)");  ax.set_ylabel("Number"); ax.set_yscale('log')
        plt.ylim(0.1, 1.e3)
        if SaveGraph: plt.savefig(path1+OutputFolder+figname+OutFormat)
    bins=np.linspace(0, 200000 ,num=1001)
    hist, bin_edges = np.histogram(dTE,bins=bins)    #   np.bincount(P)
    #TE0=bin_edges[np.argmax(hist)]
    hist[np.argmax(hist)]=0
    TE1=bin_edges[np.argmax(hist)]
    hist[np.argmax(hist)]=0
    TE2=bin_edges[np.argmax(hist)]
    TEhigh=max(TE1,TE2)
    TElow=min(TE1,TE2)
    print('TEhigh (ms):', TEhigh, "TElow (ms):", TElow)
#    print('periodhisto', np.amin(dTE), np.amax(dTE))
    return TEhigh, TElow, dTE

def exporttrack(X,Y,Z,T,P, datasave):    # exports raw traces as csv file
    file = open(datasave, "w")
    file.write('X_nm'+','+'Y_nm'+','+'Z_nm'+','+'T_ms'+','+'P_%'+'\n')
    for i in range(len(X)):
        file.write(str(X[i])+','+str(Y[i])+ ','+str(Z[i])+ ','+str(T[i])+ ','+str(P[i])+"\n")
    file.close()

def MakeHistodZ(dZ, NZavg, refname):    # histogram of z fluctuations with rolling average on NZAvg time steps
    figname=refname+'Zhisto'
    plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
    ax=sns.distplot( dZ, kde=False, bins=np.linspace(-100,100,num=200) )
    ax.set_xlabel("Z fluctuation around "+ str(NZavg)+" steps time average (nm)");  ax.set_ylabel("Number")
    plt.xlim(-100, 100.)
    print('dZrawhisto', np.amin(dZ), np.amax(dZ))
    if SaveGraph: plt.savefig(path1+OutputFolder+figname+OutFormat)

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
#☺    return (sns.distplot( Zl, kde=False, bins=bins, label=labelZ2 )).get_lines()[0].get_data()

#   FindJumps2(T, Zs, Zsavg, NZavg, TE, DE, NE, Phigh, Plow, TEhigh, tol, dzjump, GraphDhisto, refname, pprint)
def FindJumps2(T, Zs, Zsavg, NZavg, TE, DE, NE, Phigh, Plow, TEhigh, tol, dzjump, GraphDhisto, refname, pprint):
    ISup=(NE=='SigGen.Power (%)')*(Phigh*(1-tol)<DE)*(DE<Phigh*(1+tol))
#    ISdown=(NE=='SigGen.Power (%)')*(DE<=Plow*(1+tol))
    ISdown=(NE=='SigGen.Power (%)')*(Plow*(1-tol)<DE)*(DE<Plow*(1+tol))
    ISdown=(NE=='SigGen.Power (%)')*(DE<=Plow*(1+tol))
   # ISdown=(NE=='SigGen.Power (%)')*(Plow*(1-tol*5)<DE)*(DE<Plow*(1+tol*5))
    Tup=(TE[ISup]).astype(float); Dup=DE[ISup]; nup=len(Tup); print('Nb steps UP:', nup )
    Tdown=(TE[ISdown]).astype(float); Ddown=DE[ISdown]; ndown=len(Tdown); print('Nb steps DOWN:', ndown )
    if ndown==nup+1 and Tdown[0]<Tup[0]:
        Tdown=Tdown[1:ndown]; Ddown=Ddown[1:ndown]; ndown=len(Tdown); print('remove first steps DOWN:', ndown  )
    if nup>ndown: Tup=Tup[0:nup-1]; Dup=Dup[0:nup-1]; nup=len(Tup); print('reduce Nb steps UP:', nup  )
    if pprint:
        print('Tup', len(Tup), Tup); print('Tdown', len(Tdown), Tdown)
    if GraphDhisto:
        plt.figure('Duphisto', figsize=(6,6), dpi=100); sns.distplot( Dup, kde=False, bins=20 )
        plt.figure('Ddownhisto', figsize=(6,6), dpi=100); sns.distplot( Ddown, kde=False, bins=20 )    
    if DrawSuperposal: plt.figure('Superposal'+refname, figsize=(6,6), dpi=100); ax1 = plt.gca()
    Tupdown=np.zeros(nup); Tuphigh=np.zeros(nup); Zupzero=np.zeros(nup)
    Zuplow=np.zeros(nup); Zuphigh=np.zeros(nup); Zupmid=np.zeros(nup)
    Tupjump=np.zeros(nup); Zupjump=np.zeros(nup)
    Tupjump_=np.zeros(nup); Zupjump_=np.zeros(nup)
    countOpen=0; countNoClose=0; countNoOpen=0
    countOpen_=0; countNoClose_=0; countNoOpen_=0

    for i in range(nup):  # miss the last step to avoid pb when stoped at high force
        for j in range(ndown):
            if Tup[i]<Tdown[j]:
                Tupdown[i]=Tdown[j]     #       print(i,Tup[i],j,Tdown[j])
                break
        u0=np.where(T==Tup[i])[0]; u1=np.where(T==Tupdown[i])[0]
        if np.size(u0)>0 and np.size(u1)>0:
            j0=int(u0); j1=int(u1)     #    print(i,j, j0, j1)
        else:
            print('    next step'); break
        Zupzero[i]=np.amin(Zs[range(j0-1-5*NZavg,j0-1)])
        Zuplow[i]=np.mean(Zs[range(j0-1-NZavg,j0-1)])    # Zs[j0-1]
        Zupmid[i]=np.mean(Zs[range(j0,j0+NZavg)])   #Zs[j0]
        Zuphigh[i]=Zsavg[j0+NZavg]; Tuphigh[i]=T[j0+NZavg]      
        if DrawSuperposal and i>5:
            j0m=j0-500; j1p=j1+100
      #      Zs1=Zs[range(j0m,j1p)]-np.amin(Zs[range(j0m,j0)])
            Zs1=Zs[range(j0m,j1p)]-np.mean(Zs[range(j0m,j0)])
            T1=T[range(j0m,j1p)]-T[j0]
            Zs1avg=np.convolve(Zs1, np.ones((NZavg,))/NZavg, mode=mode)
            ax1.scatter(T1, Zs1, c='r', alpha=0.01)
            ax1.plot(T1, Zs1avg, c='k', alpha=0.5)
 #       if not step_detectON:  ################################################"
        print('-----------------------------------------------') 
        for j in range(j0+1,j1):
            if Zsavg[j]>Zuphigh[i]+dzjump:
                Tupjump[i]=T[j]; Zupjump[i]=np.amin(Zs[range(j+NZavg,j+2*NZavg)]); countOpen+=1; State='Open '; Count=countOpen
                break
        if j==j0+1 or j==j1-1:
            if Zuphigh[i]>MiddleOpenClose: Tupjump[i]=Tup[i]; Zupjump[i]=Zuphigh[i]+200; countNoClose+=1; State='NoClo'; Count=countNoClose            
            if Zuphigh[i]<MiddleOpenClose: Tupjump[i]=TEhigh+Tup[i]; Zupjump[i]=Zuphigh[i]; countNoOpen+=1; State='NoOpen'; Count=countNoOpen
        if pprint: print(i,j, int(Tupjump[i]), int(Zupjump[i]),'Direct', State, Count, ' ', 
                         int((Tupjump[i]-Tup[i])/1000), int(Zupjump[i]-Zupmid[i]))
#        if step_detectON:  ################################################"
        shift=windowsize_detect
        data=Zs[j0+shift:j1]; tdata=range(j0+shift,j1)
        dg1=gaussian_filter1d(data, windowsize_detect, order=1); dg1 /= np.abs(dg1).max()
        steps_dg1 = sd.find_steps(dg1, thres_detect)  #  print(type(steps_dg1),len(steps_dg1), steps_dg1[0])
        step_sizes, step_error = sd.get_step_sizes(data, steps_dg1, window=shift)  #      print(step_sizes, step_error)
#        step=np.mean(step_sizes); error=np.mean(step_error)
#        base=np.mean([data[step] for step in steps_dg1]); dbase=np.sqrt(np.var([data[step] for step in steps_dg1]))
        if len(steps_dg1)==1:# and step>dzjump:
            j=steps_dg1[0]+j0+shift
            Tupjump_[i]=T[j]; Zupjump_[i]=step_sizes[0]+Zupmid[i]; countOpen_+=1; Count_=countOpen_; State_='Open '
        else:
            if Zuphigh[i]>MiddleOpenClose: Tupjump_[i]=Tup[i]; Zupjump_[i]=Zuphigh[i]+200; countNoClose_+=1; State_='NoClo'; Count_=countNoClose            
            if Zuphigh[i]<MiddleOpenClose: Tupjump_[i]=TEhigh+Tup[i]; Zupjump_[i]=Zuphigh[i]; countNoOpen_+=1; State_='NoOpen'; Count_=countNoOpen
        if pprint: print(i,j, int(Tupjump_[i]), int(Zupjump_[i]),'Gaus1d', State_, Count_, len(steps_dg1),
                         int((Tupjump_[i]-Tup[i])/1000), int(Zupjump_[i]-Zupmid[i]))
        if step_detectON and plot_detect:
            ii=i%nbplot_detect                
            if ii==0: 
                figname=refname+'step_detect'+str(i)
                fig, axii= plt.subplots(nbplot_detect, 1, num=refname+'step_detect'+str(i), figsize=(nbplot_detect*1,8), dpi=100)
            ax=axii[ii]; ax2 = ax.twinx()
            ax.plot(tdata, dg1, c='b', alpha=0.7, label='dg1'); ax.legend(fontsize=6)
            ax2.plot(tdata, data, c='k', alpha=0.5, label='data'+str(i)); ax2.legend(loc=0,fontsize=6) 
            ax.set(ylim=(-1., 1.)); ax2.set(ylim=(400., 1000.))
            for k in steps_dg1:
      #          ax.plot([T[k+j0+shift],T[k+j0+shift]], [-1., 1.], c='r', alpha=0.5, label='step')
                ax.plot([k+j0+shift,k+j0+shift], [-1., 1.], c='r', alpha=0.5, label='step')
            if ii==nbplot_detect-1 and SaveGraph:
                plt.savefig(path1+OutputFolder+figname+OutFormat)
                if CloseAfterSave: plt.close()

    if step_detectON and plot_detect:
        plt.figure(refname+'Direct_vs_Gaus1d_T', figsize=(6,6), dpi=100); ax3 = plt.gca()
        ax3.scatter((Tupjump-Tup)/1000, (Tupjump_-Tup)/1000, c='r', alpha=0.5)
        plt.figure(refname+'Direct_vs_Gaus1d_Z', figsize=(6,6), dpi=100); ax3 = plt.gca()
        ax3.scatter(Zupjump_-Zupmid, Zupjump-Zupmid, c='b', alpha=0.5)
    I0=Tupjump==0; Tupjump[I0]=np.nan; Zupjump[I0]=np.nan#; Zupmid[I0]=np.nan
    
    print('countOpen=', countOpen, 'countNoClose=', countNoClose, 'countNoOpen=', countNoOpen)
    print('Total=', countOpen+countNoClose+countNoOpen)
    if not step_detectON:
        return (Tup, Dup, Tdown, Ddown, Tupdown, Tupjump, Tuphigh, Zupjump, Zuplow, 
                Zuphigh, Zupmid, Zupzero, countOpen, countNoClose, countNoOpen)
    if step_detectON:
        print('countOpen_=', countOpen_, 'countNoClose_=', countNoClose_, 'countNoOpen_=', countNoOpen_)
        print('Total_=', countOpen_+countNoClose_+countNoOpen_)
        return (Tup, Dup, Tdown, Ddown, Tupdown, Tupjump_, Tuphigh, Zupjump_, Zuplow, 
                Zuphigh, Zupmid, Zupzero, countOpen_, countNoClose_, countNoOpen_)

def Rollavg(Z, N):      # rolling average on N time steps
    Zavg=np.convolve(Z, np.ones((N,))/N, mode='same')
    n=int(N/2.)
    Zavg[range(n)]=Zavg[n+1]        # np.nan
    Zavg[range(-1,-n-1,-1)]= Zavg[-n]    #   np.nan
    return Zavg
 
def FitExp(x, *p): return p[0]+(p[1]-p[0])*np.exp(x/p[2])
def FitExp1(x, *p): return p[1]+(1-p[1])*np.exp(-x/p[0])
def FitExp0(x, *p): return np.exp(-x/p[0])
def FitLin(x, *p): return p[1]+x*p[0]
def FitLin1(x, *p): return x*p[0]

def survival(Tupjump_Tuphigh, TEhigh, countNoOpen, shift, refname, color):
# calculates, fits and plots survival curve until TEhigh/1000.-shift
    print('Survival:', len(Tupjump_Tuphigh),'countNoOpen=', countNoOpen )
    dur2s=np.sort(Tupjump_Tuphigh[Tupjump_Tuphigh>0])
#    dur2s=np.sort(np.append(Tupjump_Tuphigh, np.ones(countNoOpen)*TEhigh))
    print('nb events:', len(dur2s), 'without Nan:', len(dur2s[~np.isnan(dur2s)]) )
    dur2s = dur2s[~np.isnan(dur2s)]
    dur2s_y=np.arange(len(dur2s))
    dur2s_y=1-dur2s_y/len(dur2s)
    maxtfit=1000#?TEhigh/1000.-shift
    x2 = dur2s[dur2s<maxtfit]
    y2 = dur2s_y[dur2s<maxtfit]   #    y2log = dur2s_logy[dur2s<maxtfit]
#    dur2s_logy=np.log10(dur2s_y)
  #  plt.figure('Survival_'+refname, figsize=(6,6), dpi=100); ax = plt.gca()
    if DisplaySurvival:
        figname=refname+'Survival'
        plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
        ax.scatter(dur2s, dur2s_y, c=color, marker='o', alpha=0.5, label=refname)
        ax.set_xlabel("Time (s)");  ax.set_ylabel("Survival Fraction")
        ax.axis([0., TEhigh/1000., 0.01, 1.]); ax.set_yscale('log')
    maxtfit=TEhigh/1000.-shift
    x2 = dur2s[dur2s<maxtfit]
    y2 = dur2s_y[dur2s<maxtfit]   #    y2log = dur2s_logy[dur2s<maxtfit]
    while True:
        try:
            initialguessEq=[1.]
            pEq0, pcovEq0=curve_fit(FitExp0, x2, y2, p0=initialguessEq)
            break
        except RuntimeError:
            print("No convergence"); pEq0=[0.]
            break
    y2_fit0=FitExp0(x2, pEq0[0])
    while True:
        try:
            initialguessEq=[1., 0.]
            pEq, pcovEq=curve_fit(FitExp1, x2, y2, p0=initialguessEq)
            break
        except RuntimeError:
            print("No convergence"); pEq=[0.,0.]
            break
    y2_fit=FitExp1(x2, pEq[0], pEq[1])
    if DisplaySurvival:
        ax.plot(x2, y2_fit, c='r',alpha=0.5, label='kd 2param fit='+str(1/pEq[0]))
        ax.plot(x2, y2_fit0, c='k',alpha=0.5, label='kd 1param fit='+str(1/pEq[0]))
        ax.legend(fontsize=6)
        if SaveGraph: plt.savefig(path1+OutputFolder+figname+OutFormat)
    print('Survival 1param fit 0-', maxtfit, 's: offrate=', '{:.2E}'.format(1/pEq[0]), 'tau (s)=', '{:.2E}'.format(pEq[0]))
    print('Survival 2param fit 0-', maxtfit, 's: offrate=', '{:.2E}'.format(1/pEq[0]), 'tau (s)=', '{:.2E}'.format(pEq[0]), 'P0 (s)=', '{:.2E}'.format(pEq[1]))
    return(pEq[0], pEq[1], x2, y2, y2_fit)
    
def Plot3D(x,y,z,n, refname):    # 3D plot
    fig = plt.figure(refname, figsize=(7,7), dpi=100)
    xx=x[::n]; yy=y[::n]; zz=z[::n]
    ax = fig.add_subplot(111, projection='3d')
    # For each set of style and range settings, plot n random points in the box
    # defined by x in [23, 32], y in [0, 100], z in [zlow, zhigh].
    for m, zlow, zhigh in [('.', -50, -25), ('.', -30, -5)]:
        ax.scatter(xx, yy, zz, marker=m)
    if SaveGraph: plt.savefig(path1+OutputFolder+refname+OutFormat)

def JumpGraphs(refname, Zs, Zsavg, Zupjump_Zupmid, Tupjump_Tuphigh, Tup, Zuplow,Tupjump, Zupjump, 
               Zupmid, Tuphigh, Zuphigh, Zupzero, Dup, Pmax, T, P):
    figname=refname+'dZjumphisto_'
    plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
    ax=sns.distplot( Zupjump_Zupmid, kde=False, bins=np.linspace(-100.,500.,num=20) )
    ax.set_xlabel("Jump height (nm)");  ax.set_ylabel("Number"); plt.xlim(0, 500.);
    if SaveGraph: plt.savefig(path1+OutputFolder+figname+OutFormat)
    
    figname=refname+'JumpHigh dZ vs dT'
    plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
    ax.scatter(Tupjump_Tuphigh, Zupjump_Zupmid, c='b', marker='o', alpha=0.5, label='dZ vs dT')
    ax.set_xlabel("dt jump (s)");  ax.set_ylabel("dZ jump (nm)")
    plt.figure('JumpHigh dZ vs Zzero_'+refname, figsize=(6,6), dpi=100); ax = plt.gca()
    ax.scatter(Zupzero, Zupjump_Zupmid, c='b', marker='o', alpha=0.5, label='dZ vs Zzero')
    ax.set_xlabel("Z zero (nm)");  ax.set_ylabel("dZ jump (nm)")    
    if SaveGraph: plt.savefig(path1+OutputFolder+figname+OutFormat)

    figname=refname+'TZ'
    plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
    #ax.scatter(T, Z, c='c', marker='o', alpha=0.5, label='Z')
#    ax.plot(T, Zs, c='c', alpha=0.5, label='Zs')
    ax.scatter(Tclosed, Zsclosed, marker='.', alpha=0.2, label='Zsclosed')
    ax.scatter(Topen, Zsopen, marker='.', alpha=0.2, label='Zsopen')
    ax.scatter(Tl, Zsl, marker='.', alpha=0.2, label='Zsl')
    ax.plot(T, Zsavg, c='k', alpha=0.5, label='Zsavg')
    ax.scatter(Tup, Zuplow, c='y', marker='o', alpha=0.8, label='Zuplow')
    ax.scatter(Tupjump, Zupjump, c='r', marker='o', alpha=0.8, label='Zupjump')
    ax.scatter(Tuphigh, Zuphigh, c='g', marker='o', alpha=0.8, label='Zuphigh')
    ax.scatter(Tuphigh, Zupmid, c='b', marker='o', alpha=0.8, label='Zupmid')
    ax.legend(loc=0, fontsize=6)
    zmin=np.amin(Zs); zmax=np.amax(Zs)  #   zmin=-500; zmax=1500
    ax.set_ylim(zmin, zmax)
    ax.set_xlabel("Time (ms)");  ax.set_ylabel("Z")
    ax2 = ax.twinx(); ax2.set_ylim(-0.25, Pmax); plt.grid()
    ax2.scatter(Tup, Dup, c='b', marker='.', alpha=0.5, label='Dup')
    ax2.plot(T, P, c='b', alpha=0.5, label='P'); ax2.legend(fontsize=6)
    ax2.set_ylabel("Power (%)")
    if SaveGraph: plt.savefig(path1+OutputFolder+figname+OutFormat)
    if CloseAfterSave: plt.close()

def MakeGraphXYall(d0, imax, n, refrms, name, corr, Xr0avg, Yr0avg):   
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

def MakeGraphTZ_Power():
    figname=refname+'TZ_Power_'
    plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
    ax.scatter(Tup, Zuplow, c='k', marker='o', alpha=0.5, label='Zuplow')
    ax.scatter(Tup, Zupmid, c='m', marker='o', alpha=0.5, label='Zupmid')
    ax.scatter(Tup, Zuphigh, c='g', marker='o', alpha=0.5, label='Zuphigh')
    zmin=np.amin(Zuplow); zmax=np.amax(Zuphigh)
    ax.set_ylim(zmin, zmax); ax.legend(fontsize=6)
#    ax2 = ax.twinx(); ax2.set_ylim(-0.25, PP[ibead]*1.5); plt.grid()
    ax2 = ax.twinx(); ax2.set_ylim(-0.25, 10); plt.grid()
    ax2.scatter(Tup, Dup, c='b', marker='.', alpha=0.5, label='Dup')
    ax2.plot(T, P, c='b', alpha=0.5, label='P'); ax2.legend(fontsize=6)
    if SaveGraph: plt.savefig(path1+OutputFolder+figname+OutFormat)

def MakeGraphTupdown():
    figname=refname+'Tupdown vs Tup_'
    plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
    ax.scatter(Tup, Tupdown, c='b', marker='o', alpha=0.5, label='Tupdown')
    ax.scatter(Tup, Tupjump, c='r', marker='o', alpha=0.5, label='Tupjump')
    ax.plot(Tup, Tup, c='k', alpha=0.5, label='Tup')
    ax.legend(fontsize=6)
    plt.figure('Tupdown_Tup histo_'+refname, figsize=(6,6), dpi=100); ax = plt.gca()
    ax=sns.distplot( Tupdown_Tup, kde=False, bins=np.linspace(-10.,200.,num=20) )
    plt.figure('Tup-Tupdown_'+refname, figsize=(6,6), dpi=100); ax = plt.gca()
    ax.scatter(np.arange(nup), Tup, c='r', marker='o', alpha=0.5, label='Tup')
    ax.scatter(np.arange(nup), Tupdown, c='b', marker='o', alpha=0.5, label='Tupdown')
    ax.scatter(np.arange(nup), Tupjump, c='k', marker='o', alpha=0.5, label='Tupjump')
    ax2 = ax.twinx(); plt.grid(); ax.legend(fontsize=6)
    ax2.scatter(np.arange(nup), Tupdown_Tup, c='g', marker='.', alpha=0.5, label='Tupdown_Tup')
    ax2.legend(fontsize=6)
    ax.set_xlabel("Event Number");  ax.set_ylabel("Time (ms)")
    if SaveGraph: plt.savefig(path1+OutputFolder+figname+OutFormat)
    #ax2.set_ylim(0.,0.75)

def Spectrum(xx, axis, p1Zo, fs, label, display, fmin, fmax, axtable=None):
    friction0 = 6*np.pi*1.e-9*BeadRad       # units pN.s/nm
    if axis=='XY':
        friction = friction0 / ( 1 - (9/16)*(BeadRad/(p1Zo+BeadRad)) + (1/8)*(BeadRad/(p1Zo+BeadRad))**3 )
    elif axis=='Z':
        friction = friction0 / ( 1 - (9/8)*(BeadRad/(p1Zo+BeadRad)) + (1/2)*(BeadRad/(p1Zo+BeadRad))**3 )
    def FitSpectrum(x, *p):
        return p[1]/(2*np.pi**2)/( x**2 + (p[0]/(2*np.pi*friction))**2 )
    f, Pxx_spec = signal.periodogram(xx, fs, scaling='density') 
    Pxxt_spec= Pxx_spec[(f>fmin)&(f<fmax)]; ft=f[(f>fmin)&(f<fmax)]
#    Pxxt_spec2= Pxx_spec[(Pxx_spec>0)&(f>0)]; ft2=f[(Pxx_spec>0)&(f>0)]
#    dfP = pd.concat([pd.DataFrame(np.log10(ft2[::10])), pd.DataFrame(np.log10(Pxxt_spec2[::10]))], axis=1, keys=['d', 'P'])
#    ax = sns.jointplot(x='d', y='P', data=dfP, kind="kde")
    nbins=101; fbins=np.logspace(-2,2,nbins); Pbins=np.ones(nbins); dPbins=np.zeros(nbins)
    for m in range(nbins-1):
        u=Pxx_spec[(f>=fbins[m])&(f<fbins[m+1])]
        Pbins[m]=np.mean(u)
        dPbins[m]=np.std(u)/np.sqrt(len((u[~np.isnan(u)])))
#    ax = sns.dist(Pxxt_spec)
#    ax=sns.distplot( Pxxt_spec, kde=False, axlabel=False, bins=my_bins, label="test" )
#    plt.ylim(1e-3, 1e4); plt.xlim(1e-2, 1e2)# plt.set_xscale('log'); plt.set_yscale('log') 
#    plt.ylim(-3., 4.); plt.xlim(-2., 2.)# plt.set_xscale('log'); plt.set_yscale('log') 
    if display:        
        figname=refname+'PowerSpectrum'+label
        plt.figure(figname+'bis', figsize=(6,6), dpi=100); ax = plt.gca()
        if axtable!=None:
            wax=[ax, axtable]
        else:
            wax=[ax]   
        for axbis in wax:
            axbis.errorbar(fbins, Pbins, dPbins, marker='.', c='b', alpha=0.2)
            axbis.set_ylim(1e1, 1e4); axbis.set_xscale('log'); axbis.set_xlim(1e-2, 1e2); axbis.set_yscale('log')    
            axbis.set_xlabel('frequency [Hz]'); ax.set_ylabel('PSD [nm²/Hz] '+label);    #     ax.set_ylabel('spectrum [nm²] '+label); 
            if axbis==axtable: axbis.axes.get_xaxis().set_visible(False); axbis.axes.get_yaxis().set_visible(False)
#        plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
#        ax.scatter(f[::10], Pxx_spec[::10], marker='.', c='b', alpha=0.1)
#        ax.set_ylim(1e-3, 1e4); ax.set_xscale('log'); ax.set_xlim(1e-2, 1e2); ax.set_yscale('log')        
#        ax.set_xlabel('frequency [Hz]'); ax.set_ylabel('PSD [nm²/Hz] '+label);    #     ax.set_ylabel('spectrum [nm²] '+label); 
    while True:
        try:
            pEq, pcovEq=curve_fit(FitSpectrum, ft, Pxxt_spec, p0=[1.e-3, 1.e5])
            eEq=np.sqrt(np.diag(pcovEq))
            break
        except RuntimeError:
            print("No convergence"); pEq=[np.nan, np.nan]; eEq=[np.nan, np.nan]; break 
    pEq[0]=np.abs(pEq[0])
    FitPxxt_spec=FitSpectrum(ft, pEq[0], pEq[1])
    if display:
#        ax.plot(ft, FitPxxt_spec, c='r', alpha=0.3)
        for axbis in wax: axbis.plot(ft, FitPxxt_spec, c='r', alpha=0.8)
        if SaveGraph: plt.savefig(path1+OutputFolder+figname+OutFormat)
        if CloseAfterSave: plt.close()
    print('Spectrum'+label, ' k (pN/nm)=',pfl(pEq[0]),' D (nm²/s)=',pfl(pEq[1]))
    return pEq[0], pEq[1], eEq[0], eEq[1], friction

def gauss_function(x, *p): return p[0]*np.exp(-(x-p[1])**2/(2*p[2]**2))

def FitHistoGaussian(Length, label, display, delta):
#    bins=np.linspace(-500, 1500 ,num=2001)
    bins=np.linspace(-1500, 2500 ,num=2001)
#    print(sum(np.isinf(Length)==True), sum(np.isinf(Length)==True))
    Hy, Hx = np.histogram(Length, bins=bins, density=True)
    Hymax=np.argmax(Hy); xHymax=bins[Hymax]   #    print(Hymax, xHymax)
    indxfit=np.abs(Hx-xHymax)<delta
    Hxt=Hx[indxfit]; Hyt=Hy[indxfit[1:]]
#    print(len(Hxt),len(Hyt))
    if len(Hxt)>len(Hyt): Hxt=Hxt[1:]
#    print(sum(np.isinf(Hxt)==True), sum(np.isinf(Hyt)==True))
#    print(sum(np.isnan(Hxt)==True), sum(np.isnan(Hyt)==True))
    pEq, pcovEq=curve_fit(gauss_function, Hxt, Hyt, p0=[1., xHymax,20])
    pEq[2]=np.abs(pEq[2])
    FitHy=gauss_function(Hxt, pEq[0], pEq[1], pEq[2])
    print(label+'Mod',xHymax, '  Gaussian fit: Amp=', pfl(pEq[0]), 'Avg=', pfl(pEq[1]), 'SD=', pfl(abs(pEq[2])))
    if display:
#        MakeGraphXY(Hxt, Hyt, Hxt, FitHy, -500, 1500, 0.00001, 0.02, label, 'density', refname+'Hist'+label, line2=True, log=True)
        MakeGraphXY(Hxt, Hyt, Hxt, FitHy, -1500, 2500, 0.00001, 0.02, label, 'density', refname+'Hist'+label, line2=True, log=True)
    return pEq[0], pEq[1], pEq[2], xHymax

def plotresults(df):
    figname=refname+'Variation_NZlavg'
    plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
    for y in ['Lo_nm', 'Lc_nm', 'Ll_nm', 'Zl_nm']:
        ax.errorbar(df['NZlavg'], df['p1'+y], df['p2'+y], label=y, marker='o', ms=10)
    ax.legend(fontsize=6); ax.set_xlabel("NZlavg");  ax.set_ylabel("Peak length or Z_nm")
    ax.axis([0, 1000, -100, 1100])
    if SaveGraph: plt.savefig(path1+OutputFolder+figname+OutFormat)


##========================= ========================= =========================
##=========================   # MAIN LOOP  # ==================================
listresultsTotal=[]
color=iter(cm.rainbow(np.linspace(0, 1, 4*nbead*len(RangeNZlavg))))

if nbead>1: figAllSpec, axAllSpec = plt.subplots(nbead, 3, num='figAllSpec', figsize=(6, 2*nbead), dpi=100)

for ibead in range(nbead):        ## Loop on test beads
    print('=============================================================')
    print('NZavg=',NZavg, 'dzjump=', dzjump, 'tol=',tol, 'NZrefavg=',NZrefavg)
#    if loaddata and load[ibead]: print('Loading data', fname[ibead]); tdms_file = TdmsFile.read(path1+fname[ibead])
    if loaddata and load[ibead]: print('Loading data', fname[ibead]); tdms_file = TdmsFile.open(path1+fname[ibead])
    if loaddata and not load[ibead]: print('Data already loaded', fname[ibead])
    tdms_groups = tdms_file.groups()
#    d0 = tdms_file.group_channels(tdms_groups[0]); d1 = tdms_file.group_channels(tdms_groups[1]); d2 = tdms_file.group_channels(tdms_groups[2])
    d0 = tdms_groups[0]; d1 = tdms_groups[1]; d2 = tdms_groups[2]
    if pprint: print(tdms_groups)
    dmax=1000000.; imax=(len(d0)-1)//5; dX=1000; dY=1000    
    print('Nb of traj=', imax)
    refname=fname[ibead]+'_Bead'+str(b[ibead])+'_'; print(refname)

    (X0,Y0,Z0,MinLUT0,SSIM0)=track(d0, b[ibead])        # test bead
#    X0med=np.median(X0); Y0med=np.median(Y0)
    X0med=np.median(X0[~np.isnan(X0)]); Y0med=np.median(Y0[~np.isnan(Y0)])
    if load[ibead]:
        print("Calculating pool of reference bead(s)", wbref )
        if len(wbref)>1:
            (Xr0,Yr0,Zr0,MinLUTr0,SSIMr0)=track(d0, wbref[0])
            for i in range(1,len(wbref)):
                (Xr0_,Yr0_,Zr0_,MinLUTr0_,SSIMr0_)=track(d0, wbref[i])
                Xr0=Xr0+Xr0_; Yr0=Yr0+Yr0_; Zr0=Zr0+Zr0_
            Xr0=Xr0/len(wbref); Yr0=Yr0/len(wbref); Zr0=Zr0/len(wbref)
        else:
            (Xr0,Yr0,Zr0,MinLUTr0,SSIMr0)=track(d0, wbref[0])   # reference bead
    else: print("Pool of reference bead(s) already calculated") 
    (TE0,NE0,DE0)=events(d2)
    T0, fs = time(GraphTdms)
    P0=BuildP(T0, TE0, NE0, DE0)
    if GraphTdms:
        plt.figure('Power', figsize=(6,6), dpi=100); ax = plt.gca()
        ax.scatter(T0, P0, marker='.', c='b', alpha=0.5)

    print('Power at start', start[ibead],':', P0[start[ibead]])
    print('Power at stop:', stop[ibead],':', P0[stop[ibead]])
    if P0[start[ibead]]>0.1:
        wstart=10*np.arange(1000)+start[ibead]; wPstart=P0[wstart]
        #print(wstart); print(wPstart)
        i0=np.argmin(wPstart); start[ibead]= wstart[i0+5]
   #     print(i0)
    if P0[stop[ibead]]>0.1:
        wstop=10*np.arange(1000)+stop[ibead]; wPstop=P0[wstop]
        #print(wstop); print(wPstop)
        i0=np.argmin(wPstop); stop[ibead]= wstop[i0+5]
   #     print(i0)
    print('Power at new start', start[ibead],':', P0[start[ibead]])
    print('Power at new stop:', stop[ibead],':', P0[stop[ibead]])

#    if P0[stop[ibead]]>0.1:
    print('Tracks shortening from frame', start[ibead], 'to ',stop[ibead])
    # shorten test bead
    (X,Y,Zu,T,P, MinLUT)= shortentrack(X0, Y0, Z0, T0, P0, MinLUT0, start[ibead], stop[ibead])
    # shorten reference bead
    (Xr,Yr,Zr,Tr,Pr, MinLUTr)= shortentrack(Xr0, Yr0, Zr0, T0, P0, MinLUTr0, start[ibead], stop[ibead])
     
    cutE=range(np.argmax(TE0>=T0[start[ibead]]))
    mask=(TE0>=T0[start[ibead]])&(TE0<T0[stop[ibead]])
    TE=TE0[mask]; NE=NE0[mask]; DE=DE0[mask];
#    TE=np.delete(TE0,cutE); NE=np.delete(NE0,cutE); DE=np.delete(DE0,cutE)
    
    TEhigh, TElow, dTE=period(TE, GraphTdms, refname)

    Phigh, Plow ,hist,bin_edges= MakeHistoP(P, GraphTdms)
    # https://stackoverflow.com/questions/13728392/moving-average-or-running-mean
    mode='same'   # edge modes = ['full', 'same', 'valid']  
    Zravg=Rollavg(Zr, NZrefavg )        # rolling average of shorten reference bead
    Zr0avg=Rollavg(Zr0, NZrefavg )       # rolling average of reference bead
    if RefBeadZcorrection:
        print('Zcorrection by reference bead(s): Zravg=', Zravg[0] )#int(bref[ibead]) )
        Zc=Zu-Zravg; Z0c=Zr0-Zr0avg; Zec=Z0-Zr0avg
   #     Zc=Zu-Zravg+BeadRad; Z0c=Zr0-Zr0avg+BeadRad
    else: print('No Zcorrection'); Zc=Zu; Z0c=Zr0; Zec=Z0
    Xr0avg=Rollavg(Xr0, NZrefavg ); Yr0avg=Rollavg(Yr0, NZrefavg ) 
    Xravg=Rollavg(Xr, NZrefavg ); Yravg=Rollavg(Yr, NZrefavg ) 
    if RefBeadXYcorrection:
         print('XYcorrection by reference bead(s): Xravg=', Xravg[0],' Yravg=', Yravg[0] )
         Xc=X-Xravg; Yc=Y-Yravg; Xec=X0-Xr0avg; Yec=Y0-Yr0avg
    else: print('No XYcorrection'); Xc=X; Yc=Y; Xec=X0; Yec=Y0

#    if GraphXYall: MakeGraphXYall(tdms_file, imax, 100, 100, fname[ibead], False, Xr0avg, Yr0avg)
    if GraphXYall: MakeGraphXYall(tdms_file, imax, 100, 100, fname[ibead], True, Xr0avg, Yr0avg)

    indnonan=(~np.isnan(Xc))*(~np.isnan(Yc))
    Xc=Xc[indnonan]; Yc=Yc[indnonan]; Zc=Zc[indnonan]
    MinLUT=MinLUT[indnonan]; T=T[indnonan]; P=P[indnonan]
    listresults=[]
    if export:    # exports raw data for third party analysis
        exporttrack(X,Y,Zu,T,P, path1+'Export_'+refname+'.csv')
        exporttrack(Xr,Yr,Zr,Tr,P, path1+'Export_'+refname+'_r'+str(wbref)+'.csv')
            
    p1Lo=p1Lo0; p1Lc=p1Lc0; p1Ll=p1Ll0; p1Zl=p1Zl0  # initialization of gaussian center for fit guess
    for NZlavg in RangeNZlavg:      # test loop for various choices of NZlavg
        print('******* NZlavg=', NZlavg, '*******')
        Zcavg=np.convolve(Zc, np.ones((NZlavg,))/NZlavg, mode=mode); dZ=Zc-Zcavg
        
        Zecavg=np.convolve(Zec, np.ones((NZlavg,))/NZlavg, mode=mode); dZe=Zec-Zecavg
        
        indhigh=(P>=Phigh*(1-tol))*(P<=Phigh*(1+tol))
        indlow=(P>=Plow*(1-tol))*(P<=Plow*(1+tol)); #indlow=(P<Phigh*0.9)+(P>Pmax)
        ind0force=(P==0)
        indeforce0=(P0==0)
        x = np.random.randn(len(P))/10
        Zch=Zc[indhigh]; Xch=Xc[indhigh]; Ych=Yc[indhigh]; dZh=dZ[indhigh]; Ph=P[indhigh]+x[indhigh]
        Zcl=Zcavg[indlow]; Xcl=Xc[indlow]; Ycl=Yc[indlow]; dZl=dZ[indlow]; Pl=P[indlow]+x[indlow]
        Zc0=Zcavg[ind0force]; Xc0=Xc[ind0force]; Yc0=Yc[ind0force]; dZ0=dZ[ind0force]
        
        Zec0=Zecavg[indeforce0]; Xec0=Xec[indeforce0]; Yec0=Yec[indeforce0]; dZe0=dZe[indeforce0]
        
        indstart=(T0<start[ibead]); Z0start=Z0c[indstart]
        Nl=len(Zc[indlow]); Nh=len(Zc[indhigh]); N0=len(Zc[ind0force]); Ntot=len(Zc)
        print('States Populations: low=', Nl, ' high=', Nh, ' 0force=', N0, ' total=', Ntot)        
        print('States Fractions: low=', Nl/Ntot, ' high=', Nh/Ntot, ' 0force=', N0/Ntot, 'total=', (Nl+Nh+N0)/Ntot)        
        for AnchorPointState in AnchorPointStateList:      # ['low','0']     
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
        #        SingleTrace(0, sampleGraph, Xc, Yc, Zc, T, MinLUT )
        #        SingleTrace(-1,sampleGraph, Xr, Yr, Zr, T, MinLUT )
                SingleTrace(-2,sampleGraph, Xs, Ys, Zs, T, MinLUT )
                MakeGraphXY(T, Length, T, Zs, 0, T[-1], 0, 1200, 'Length', 'Zs', refname+'Length')
    
            print('Opening analysis********')
            Zsavg=np.convolve(Zs, np.ones((NZavg,))/NZavg, mode=mode)
            DeltaCOFit=0
            for i in [0]:         #  [0,1] first pass dzjump, 2nd pass DeltaCOFit
                dzjump_tmp = dzjump*(i==0) + DeltaCOFit*0.7*(i==1)
                print('Round ',i,' /0  dzjump=',dzjump_tmp)
                (Tup, Dup, Tdown, Ddown, Tupdown, Tupjump, Tuphigh, Zupjump, Zuplow, Zuphigh, Zupmid,Zupzero, countOpen, countNoClose,
                     countNoOpen)=FindJumps2(T, Zs, Zsavg, NZavg, TE, DE, NE, Phigh, Plow, TEhigh, tol, dzjump_tmp, GraphDhisto, refname, pprint)
     #           nup=len(Tup); ndown=len(Tdown)
                Tupdown_Tup=(Tupdown-Tup)/1000.; Tupjump_Tup=(Tupjump-Tup)/1000.
                Zupjump_Zuphigh=Zupjump-Zuphigh; Tupjump_Tuphigh=(Tupjump-Tuphigh)/1000.
                Zupjump_Zupmid=Zupjump-Zupmid
        
                indclosed=P<-10; indopen=P<-10      # labeling of tracking data with open and closed states
                for nj in range(len(Tupjump)):
                    indclosed=indclosed+(T<Tupjump[nj])*(T>Tuphigh[nj])*indhigh*(Tupjump[nj]>Tup[nj])
           #         indopen=indopen+(T>Tupjump[nj])*(T<Tdown[nj+1])*indhigh
       #             indopen=indopen+(T>Tupjump[nj])*(T<Tdown[nj])*indhigh
                    indopen=indopen+(T>Tupjump[nj])*(T<Tupdown[nj])*indhigh
                Zsclosed=Zs[indclosed]; Tclosed=T[indclosed]; Xsclosed=Xs[indclosed]; Ysclosed=Ys[indclosed]; Lengthclosed=Length[indclosed]
                Zsopen=Zs[indopen]; Topen=T[indopen]; Xsopen=Xs[indopen]; Ysopen=Ys[indopen]; Lengthopen=Length[indopen]
                Nc=len(Zsclosed); No=len(Zsopen)
                print('States Population: close=', Nc, ' open=', No, ' high=', Nh)
                print('States Fractions: close=', Nc/Nh, ' open=', No/Nh, ' (Nc+No)/Nh=', (Nc+No)/Nh)
    
                print('Gaussian fits ********')
                if i==0:
                    p0Lo, p1Lo, p2Lo, pLo = FitHistoGaussian(Lengthopen, 'Lengthopen nm ', DisplayGaussianFits , 100 )
                    p0Lc, p1Lc, p2Lc, pLc = FitHistoGaussian(Lengthclosed, 'Lengthclosed nm ', DisplayGaussianFits, 100 )
                    p0Ll, p1Ll, p2Ll, pLl = FitHistoGaussian(Lengthl, 'Lengthlow nm ', False, 300 )
                    if len(Length0)>0: p0L0, p1L0, p2L0, pL0 = FitHistoGaussian(Length0, 'Length0 nm ', False, 300 )
                    p0Zl, p1Zl, p2Zl, pZl = FitHistoGaussian(Zsl, 'Zsl nm ', False, 300 )
                    p0Xo, p1Xo, p2Xo, pXo = FitHistoGaussian(Xsopen, 'Xsopen nm ', False , 300 )
                    p0Yo, p1Yo, p2Yo, pYo = FitHistoGaussian(Ysopen, 'Ysopen nm ', False , 300 )
                p0Zo, p1Zo, p2Zo, pZo = FitHistoGaussian(Zsopen, 'Zsopen nm ', False , 150 )
                p0Zc, p1Zc, p2Zc, pZc = FitHistoGaussian(Zsclosed, 'Zsclosed nm ', False , 150 )
    #            p0X0, p1X0, p2X0, pX0 = FitHistoGaussian(Xs0, 'Xs0force nm ', False , 300 )
    #            p0Y0, p1Y0, p2Y0, pY0 = FitHistoGaussian(Ys0, 'Ys0force nm ', False , 300 )
    #            p0Z0, p1Z0, p2Z0, pZ0 = FitHistoGaussian(Zs0, 'Zs0force nm ', False , 150 )
                DeltaCOFit=p1Zo-p1Zc   

            if DisplayCloseOpen:
                figname=refname+'XY_'+'TZclosedopen'
                plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca(); ax2=ax.twinx()
                ax.scatter(Tclosed, Zsclosed, marker='.', alpha=0.2, label='Zsclosed')
                ax.scatter(Topen, Zsopen, marker='.', alpha=0.2, label='Zsopen')
                ax.scatter(Tl, Zsl, marker='.', alpha=0.2, label='Zsl')
                ax.plot(T, Zsavg, c='m', alpha=0.5, label='Zsavg')
                ax2.plot(T, P, c='k', alpha=0.3, label='P'), ax2.legend(fontsize=6); 
                ax.legend(fontsize=6); ax.axis([0, T[-1], -100, 1200])
                if SaveGraph: plt.savefig(path1+OutputFolder+figname+OutFormat)
                if CloseAfterSave: plt.close()
            MakeHistoZ(Zsclosed, Zsopen, 'Zsclosed', 'Zsopen', -100, 1200, refname+'ZsclosedZsopen')
            
            print('Power Spectra fits********')
            if nbead>1: axtable0=axAllSpec[ibead,0]; axtable1=axAllSpec[ibead,1]; axtable2=axAllSpec[ibead,2]
            else: axtable0=None; axtable1=None; axtable2=None
            kx, Dx, dkx, dDx, frictionXY = Spectrum(Xsopen, 'XY', p1Zo, fs,'XsopenRef'+AnchorPointState, DisplaySpectrum, fmin, fmax, axtable=axtable0) 
            ky, Dy, dky, dDy, frictionXY = Spectrum(Ysopen, 'XY', p1Zo, fs,'YsopenRef'+AnchorPointState, DisplaySpectrum, fmin, fmax, axtable=axtable1)
            kz, Dz, dkz, dDz, frictionZ = Spectrum(Zsopen, 'Z', p1Zo, fs,'ZsopenRef'+AnchorPointState, DisplaySpectrum, fmin, fmax, axtable=axtable2)
            Fx=kx*(BeadRad+p1Lo);  Fy=ky*(BeadRad+p1Lo); Fz=kz*(BeadRad+p1Lo)
            dFx=dkx*(BeadRad+p1Lo);  dFy=dky*(BeadRad+p1Lo); dFz=dkz*(BeadRad+p1Lo)
            dFy=np.abs(dkx*(BeadRad+p1Lo))+np.abs(kx*p2Lo)
            dFx=np.abs(dky*(BeadRad+p1Lo))+np.abs(ky*p2Lo)
            dFz=np.abs(dkz*(BeadRad+p1Lo))+np.abs(kz*p2Lo)
            if nbead>1: 
                axAllSpec[ibead,0].set_title(str(b[ibead])+' Fx='+str(pfl(Fx)), fontsize=6)
                axAllSpec[ibead,1].set_title(str(b[ibead])+' Fy='+str(pfl(Fy)), fontsize=6)
                axAllSpec[ibead,2].set_title(str(b[ibead])+' Fz='+str(pfl(Fz)), fontsize=6)
            print("Corner frequencies (Hz) x,y,z: ", pfl(kx*Dx/4), pfl(ky*Dy/4), pfl(kz*Dz/4))
            print("Forces (pN) x,y,z: ", pfl(Fx), pfl(Fy), pfl(Fz), "Forces Errors (pN) x,y,z: ", pfl(dFx), pfl(dFy), pfl(dFz))
            Dxytheo=4/frictionXY; Dztheo=4/frictionZ
            print("D_by_Dtheo x, y,z: ", pfl(Dx/Dxytheo), pfl(Dy/Dxytheo), pfl(Dz/Dztheo))
            FSDXFit=4*(p1Lo+BeadRad)/p2Xo**2; FSDXMod=4*(pLo+BeadRad)/p2Xo**2
            FSDYFit=4*(p1Lo+BeadRad)/p2Yo**2; FSDYMod=4*(pLo+BeadRad)/p2Yo**2    
 
            if GraphTZ_Power: MakeGraphTZ_Power()
            if GraphTupdown: MakeGraphTupdown()
            if CalculateSurvival: (pE0, pE1, x2, y2, y2_fit)=survival(Tupjump_Tuphigh, TEhigh, countNoOpen, tshiftsurvival, refname, next(color))
     #       if CalculateSurvival: (pE0, pE1, x2, y2, y2_fit)=survival(Tupjump_Tuphigh, TEhigh, 2., '', next(color))   # to plot all beads on th same figure
            if DisplayJumpGraphs:
                JumpGraphs(refname, Zs, Zsavg, Zupjump_Zupmid, Tupjump_Tuphigh, Tup, Zuplow, Tupjump, Zupjump, Zupmid,
                           Tuphigh, Zuphigh, Zupzero, Dup, Phigh*(1+tol), T, P)
        
            results=(tdmsfilename, b[ibead], wbref, BeadRad, X0med, Y0med, NZavg, NZlavg, dzjump, NZrefavg, tol, sample, start[ibead], stop[ibead], fs, TEhigh, TElow, Phigh, Plow,
                     PullAngle, PullAngle2, p1Lo, p1Lc, p1Ll, p1Zl,p2Lo, p2Lc, p2Ll, p2Zl, pLo, pLc, pLl, pZl, 
                     p0Xo, p1Xo, p2Xo, pXo, p0Yo, p1Yo, p2Yo, pYo, p0Zo, p1Zo, p2Zo, pZo, fmin, fmax,
                     kx, ky, kz, Dx, Dy, Dz, Dxytheo, Dztheo, Fx, Fy, Fz, dFx, dFy, dFz, FSDXFit, FSDXMod, FSDYFit, FSDYMod,
                     step_detectON, thres_detect, windowsize_detect, countOpen, countNoClose, countNoOpen, Nc, No, Nh, 1/pE0, pE1, Zupjump_Zupmid, Tupjump_Tuphigh)  
            listresults.append(results)
            listresultsTotal.append(results)
            columns=['tdmsfilename', 'bead', 'refbeads', 'BeadRad_nm', 'X0med', 'Y0med', 'NZavg', 'NZlavg', 'dzjump_nm', 'NZrefavg', 'tol', 'sample', 'start_frame', 'stop_frame', 'fs_Hz', 'TEhigh_ms', 'TElow_ms', 'Phigh_%', 'Plow_%',
                     'PullAngle_deg', 'PullAngle2_deg', 'p1Lo_nm' ,'p1Lc_nm','p1Ll_nm','p1Zl_nm', 'p2Lo_nm' ,'p2Lc_nm','p2Ll_nm','p2Zl_nm','pLo_nm', 'pLc_nm', 'pLl_nm', 'pZl_nm',
                     'p0Xo_nm', 'p1Xo_nm', 'p2Xo_nm', 'pXo_nm', 'p0Yo_nm', 'p1Yo_nm', 'p2Yo_nm', 'pYo_nm', 'p0Zo_nm', 'p1Zo_nm', 'p2Zo_nm', 'pZo_nm', 'fmin_Hz', 'fmax_Hz',
                     'kx_pN_per_nm', 'ky_pN_per_nm', 'kz_pN_per_nm', 'Dx_nm2_per_s', 'Dy_nm2_per_s', 'Dz_nm2_per_s', 'Dxytheo_nm2_per_s',
                     'Dztheo', 'FSpecx_pN', 'FSpecy_pN', 'FSpecz_pN', 'dFSpecx_pN', 'dFSpecy_pN', 'dFSpecz_pN', 'FSDXFit_pN', 'FSDXMod_pN', 'FSDYFit_pN', 'FSDYMod_pN',
                     'step_detectON', 'thres_detect', 'windowsize_detect', 'countOpen', 'countNoClose', 'countNoOpen', 'Nc', 'No', 'Nh', 'Offrate_1_per_s', 'Offrate_P0', 'Zupjump_Zupmid', 'Tupjump_Tuphigh']
    df = pd.DataFrame(listresults, columns = columns)
    plotresults(df)
if nbead>1 and SaveGraph: plt.figure('figAllSpec'); plt.savefig(path1+OutputFolder+fname[0]+'figAllSpec'+OutFormat)
dfTotal = pd.DataFrame(listresultsTotal, columns =columns)
dfTotal.to_csv(path1+OutputFolder+fname[0]+'FitsResults.csv')

#dfTotal=pd.read_csv('E:/Laurent/AFS/Output20200319/dfTotal.csv', index_col=False, sep=',', skiprows=0)
listFig1=['Peak Length Open', 'Delta CloseOpen', 'Delta vs Peak Open', 'SD Open', 'ForceSDX', 'D_by_Dtheo'] 
listFig2=['ForceSpectreXY','ForceX SD vs PSD', 'ForceSDXY']
nFig1=len(listFig1); ax1=[]; figname1=[]; color=iter(cm.rainbow(np.linspace(0, 1, nbead)))
nFig2=len(listFig2); ax2=[]; figname2=[]
if DisplayplotvariableNZlavg:
    for i in range(nFig1): figname1.append( plt.figure(listFig1[i], figsize=(6,6), dpi=100) ); ax1.append(plt.gca())
for i in range(nFig2): figname2.append( plt.figure(listFig2[i], figsize=(6,6), dpi=100) ); ax2.append(plt.gca())

for ib in b:
    c=[next(color),]
    dfT=dfTotal.where(dfTotal['bead']==ib)
    DeltaCOFit=dfT['p1Lo_nm']-dfT['p1Lc_nm']; DeltaCOMod=dfT['pLo_nm']-dfT['pLc_nm']

#    ForceSDXFit=4*(dfT['p1Lo']+BeadRad)/dfT['p2Xo']**2; ForceSDXMod=4*(dfT['pLo']+BeadRad)/dfT['p2Xo']**2
#    ForceSDYFit=4*(dfT['p1Lo']+BeadRad)/dfT['p2Yo']**2; ForceSDYMod=4*(dfT['pLo']+BeadRad)/dfT['p2Yo']**2   
 #   ForceSDXFit=dfT['FSDXFit']; ForceSDXMod=dfT['FSDXMod']
#    ForceSDYFit=dfT['FSDYFit']; ForceSDYMod=dfT['FSDYMod'] 
    #    ForceSpecX=(dfT['p1Lo']+BeadRad)*np.abs(dfT['kx']); ForceSpecY=(dfT['p1Lo']+BeadRad)*np.abs(dfT['ky'])
 #   ForceSpecX=dfT['Fx']; ForceSpecY=dfT['Fy']

    if DisplayplotvariableNZlavg:
        ax1[0].scatter(dfT['NZlavg'], dfT['p1Lo_nm'], marker='.', c=c, alpha=0.5, label=str(ib)+' Fit')
        ax1[0].scatter(dfT['NZlavg'], dfT['pLo_nm'], marker='+', c=c, alpha=0.5, label=str(ib)+' Mod')
        ax1[0].set_xlabel("NZlavg");  ax1[0].set_ylabel("Peak Length Open_nm")    
        ax1[1].scatter(dfT['NZlavg'], DeltaCOFit, marker='.', c=c, alpha=0.5, label=str(ib)+' Fit')
        ax1[1].scatter(dfT['NZlavg'], DeltaCOMod, marker='+', c=c, alpha=0.5, label=str(ib)+' Mod')
        ax1[1].set_xlabel("NZlavg");  ax1[1].set_ylabel("Delta CloseOpen_nm")
        ax1[4].scatter(dfT['NZlavg'], dfT['FSDXFit_pN'], marker='.', c=c, alpha=0.5, label=str(ib)+' Fit')
        ax1[4].scatter(dfT['NZlavg'], dfT['FSDXMod_pN'], marker='+', c=c, alpha=0.5, label=str(ib)+' Mod')
        ax1[4].set_xlabel("NZlavg");  ax1[4].set_ylabel("Force SDX _pN")
        ax1[5].scatter(dfT['NZlavg'], np.abs(dfT['Dx_nm2_per_s'])/dfT['Dxytheo_nm2_per_s'], marker='.', c=c, alpha=0.5, label=str(ib)+' Dx')
        ax1[5].scatter(dfT['NZlavg'], np.abs(dfT['Dy_nm2_per_s'])/dfT['Dxytheo_nm2_per_s'], marker='.', c=c, alpha=0.5, label=str(ib)+' Dy')
        ax1[5].set_xlabel("NZlavg");  ax1[5].set_ylabel("D/Dtheo")
    
        ax1[2].scatter(dfT['p1Lo_nm'], DeltaCOFit, marker='.', c=c, alpha=0.5, label=str(ib)+' Fit')
        ax1[2].scatter(dfT['pLo_nm'], DeltaCOMod, marker='+', c=c, alpha=0.5, label=str(ib)+' Mod')
        ax1[2].set_xlabel("Peak Length Open_nm");  ax1[2].set_ylabel("Delta CloseOpen_nm")
        ax1[3].scatter(dfT['p1Lo_nm'], dfT['p2Lo_nm'], marker='.', c=c, alpha=0.5, label=str(ib)+' Fit')
        ax1[3].scatter(dfT['pLo_nm'], dfT['p2Lo_nm'], marker='+', c=c, alpha=0.5, label=str(ib)+' Fit')
        ax1[3].set_xlabel("Peak Length Open_nm");  ax1[3].set_ylabel("SD Open_nm")

    ax2[2].scatter(dfT['FSDXFit_pN'], dfT['FSDYFit_pN'], marker='.', c=c, alpha=0.5, label=str(ib)+' Fit')
    ax2[2].scatter(dfT['FSDXMod_pN'], dfT['FSDYMod_pN'] , marker='+', c=c, alpha=0.5, label=str(ib)+' Mod')
    ax2[2].set_xlabel("Force SDX_pN");  ax2[2].set_ylabel("Force SDY_pN")
#    ax2[0].scatter(dfT['FSpecx_pN'], dfT['FSpecy_pN'], marker='.', c=c, alpha=0.5, label=str(ib)+' Angle='+str(dfT['PullAngle_deg'][dfT.idxmin()[0]]))
    ax2[0].scatter(dfT['FSpecx_pN'], dfT['FSpecy_pN'], marker='.', c=c, alpha=0.5, label=str(ib)+' Angle='+str(dfT['PullAngle_deg'][0]))
    ax2[0].set_xlabel("ForceSpecX_pN");  ax2[0].set_ylabel("ForceSpecY_pN")
    ax2[1].scatter(dfT['FSpecx_pN'], dfT['FSDXFit_pN'], marker='.', c=c, alpha=0.5, label=str(ib))
    ax2[1].set_xlabel("ForceSpec_pN");  ax2[1].set_ylabel("ForceSD_pN")

for axi in [ax2[2], ax2[0], ax2[1]]: axi.plot(np.linspace(0, 2, 10), np.linspace(0, 2, 10), c='k', alpha=0.5, label=' y=x')    
if DisplayplotvariableNZlavg:
    for i in range(nFig1):
        ax1[i].legend(fontsize=6)
        if SaveGraph: plt.figure(listFig1[i]); plt.savefig(path1+OutputFolder+fname[0]+listFig1[i]+OutFormat)
for i in range(nFig2):
    ax2[i].legend(fontsize=6)
    if SaveGraph: plt.figure(listFig2[i]); plt.savefig(path1+OutputFolder+fname[0]+listFig2[i]+OutFormat)
