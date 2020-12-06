#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""""
Created on Sun Sep 29 19:33:24 2019
@author: Laurent Limozin - Laboratoire Adgesion Inflammation
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
         v3.81  Step visualization in FindJumps() and graph saving
         v3.83  PullAngle2 calculated between high force and zero force state
20200703 v3.85  Correct bug for start / stop
         v3.86  Records absolute coordinates; calculate fit error on force
20200716 v3.87  Corection bugs for start/stop
20200826 v3.88  Correction to exclude Nan for Median XY position
20201009 v3.9   Permits to select the first long step for the spectrum calculation (set AnchorPointStateList=['0 force'] )
                Take into account the NoOpen cases in the survival curve
20201025 v3.91  Revised NoOpen/NoClose cases with introduction of a threshold MiddleOpenClose=900nm
                Some discrepancy for step_detectON=0 or 1 (which may depend on windowsize_detect and NZavg)
20201206 v4.0   Full revision of code; All graphic functions shifted to AFS_JDNA_Graphs
                Adapt to last version of nptdms package
                Option of custom interval for anchor point and power spectrum (cf graph 'Various intervals')
                    range_anchorpoint=(start,stop) range_spectrum=(start,stop) and boolean Select_range_spectrum

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
import AFS_JDNA_Graphs as AFSG
from scipy.ndimage.filters import gaussian_filter1d
import warnings
warnings.filterwarnings("ignore")
np.random.seed(444); plt.close('all')
#def pfl(f): return '{:.2E}'.format(f)

#path1="C:/DATA/Experiences Diverses/AFS/Data/Debug20201002/"      #☺ data folder
path1="/home/laurent/DATA/Experiences Diverses/AFS/Data/" #☺ data folder  #path1="D:/Jim/JDNA/2020-08-24/FOV3/"  
tdmsfilename="20200828-121251.tdms"          #tdmsfilename="20200903-144604.tdms"    tdmsfilename="20200826-120242.tdms"
OutputFolder='Outputtest/'
if not os.path.exists(path1+OutputFolder): os.makedirs(path1+OutputFolder)

wbref=[6,28,49]      # reference beads  wbref=[17,29,43]               # reference beads
b=[23]               # test beads  b=[6, 8, 21, 41, 46, 58, 73, 91]     #b=[58, 73, 91] b=[73, 91] b=[5]

#RangeNZlavg=[20,40,60, 80, 100, 120, 140, 160, 200, 240, 280, 320, 400, 480, 560, 640, 720, 800, 900, 1000]
RangeNZlavg=[400]       # size of moving average for smoothing z before anchor point determination

# initialiation for various intervals
nbead=len(b); fname=[""]*nbead; start=np.zeros(nbead, dtype=int); stop=np.zeros(nbead, dtype=int); load=np.zeros(nbead, dtype=int)
#for nf in range(nbead): fname[nf] = tdmsfilename; start[nf]=15000; stop[nf]=1600000; load[nf]=nf==0  # stop[nf]=1600000
for nf in range(nbead): fname[nf] = tdmsfilename; start[nf]=30*60*1000/18; stop[nf]=175*60*1000/18; load[nf]=nf==0  # stop[nf]=1600000
range_anchorpoint=(0,int(180000/18)); range_spectrum=(int(0.9e6/18), int(1.1e6/18)); Select_range_spectrum=False

##========================= ========================= =========================
# General settings
loaddata=True  # to load the tdms file
RefBeadZcorrection=True # correction with reference bead
RefBeadXYcorrection=True # correction with reference bead
AnchorPointStateList=['0 force'] # correction with AnchorPoint as determined by anchor point  'low force' or '0 force' or 'custom' ; or 'None': no anchor point 
export=False # to produce csv files of raw data (for Charlie)
CalculateSurvival=False  # calculate and draw survival curve
tshiftsurvival=0; 
step_detectON=True; plot_detect=False; nbplot_detect=15; thres_detect=0.5; windowsize_detect=20 # parameters for gaussian step detection
NZavg=100; NZlavg=400; dzjump=80.; NZrefavg=5000; tol=0.01; sample=1; sampleGraph=100 # parameters for direct step detection
MiddleOpenClose=900
BeadRad=790. # Bead Radius in nm
fmin=0.1; fmax=15.   # frequency range for spectrum fit (Hz)
p1Lo0=1000; p1Lc0=800; p1Ll0=500; p1Zl0=500  # guess for histogram fit

# Print and display options
pprint=False; # printing detailed infos in console window
SaveGraph=True; OutFormat=".jpg" ; CloseAfterSave=False # alternative ".jpg" or ".pdf" 
GraphTdms=False  # graph of frame rate and force data (duration and amplitude of 2 steps)
GraphDhisto=False; GraphTupdown=False; GraphTZ_Power=True
HistoZ=True; GraphXYZ=False # characteristics of trajectory of the selected beads
GraphXYall=False; iminXY=0; imaxXY=0 # display xy trajectories with number on a single graph
DisplaySpectrum=True
DisplayCloseOpen=True
DisplayGaussianFits=False
DisplayJumpGraphs=True # graphs of jump analysis
DrawSuperposal=False # graph of superposition of jump traces
DisplaySurvival=True
Display3D=False
DisplayplotvariableNZlavg=False

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
    print('dTrawhisto (frame time ms)', np.amin(dT), np.amax(dT), '  frequency (Hz)', "%.3f" % (1000/np.amax(dT)))   
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

def exporttrack(X,Y,Z,T,P, datasave):    # exports raw traces as csv file
    file = open(datasave, "w")
    file.write('X_nm'+','+'Y_nm'+','+'Z_nm'+','+'T_ms'+','+'P_%'+'\n')
    for i in range(len(X)):
        file.write(str(X[i])+','+str(Y[i])+ ','+str(Z[i])+ ','+str(T[i])+ ','+str(P[i])+"\n")
    file.close()
    
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
    
def FindJumps(T, Zs, Zsavg, NZavg, TE, DE, NE, Phigh, Plow, TEhigh, tol, dzjump, GraphDhisto, refname, pprint):
    ISup=(NE=='SigGen.Power (%)')*(Phigh*(1-tol)<DE)*(DE<Phigh*(1+tol))
    ISdown=(NE=='SigGen.Power (%)')*(Plow*(1-tol)<DE)*(DE<Plow*(1+tol))
#    ISdown=(NE=='SigGen.Power (%)')*(DE<=Plow*(1+tol))
    Tup=(TE[ISup]).astype(float); Dup=DE[ISup]; nup=len(Tup); print('Nb steps UP:', nup, end='  ' )
    Tdown=(TE[ISdown]).astype(float); Ddown=DE[ISdown]; ndown=len(Tdown); print('Nb steps DOWN:', ndown )
    if ndown==nup+1 and Tdown[0]<Tup[0]:
        Tdown=Tdown[1:ndown]; Ddown=Ddown[1:ndown]; ndown=len(Tdown); print('remove first steps DOWN:', ndown  )
    if nup>ndown: Tup=Tup[0:nup-1]; Dup=Dup[0:nup-1]; nup=len(Tup); print('reduce Nb steps UP:', nup  )
    if pprint: print('Tup', len(Tup), Tup); print('Tdown', len(Tdown), Tdown)
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
            Zs1=Zs[range(j0m,j1p)]-np.mean(Zs[range(j0m,j0)])
            T1=T[range(j0m,j1p)]-T[j0]
            Zs1avg=np.convolve(Zs1, np.ones((NZavg,))/NZavg, mode=mode)
            ax1.scatter(T1, Zs1, c='r', alpha=0.01)
            ax1.plot(T1, Zs1avg, c='k', alpha=0.5)
        if pprint: print('-----------------------------------------------') 
        for j in range(j0+1,j1):
            if Zsavg[j]>Zuphigh[i]+dzjump:
                Tupjump[i]=T[j]; Zupjump[i]=np.amin(Zs[range(j+NZavg,j+2*NZavg)]); countOpen+=1; State='Open '; Count=countOpen
                break
        if j==j0+1 or j==j1-1:
            if Zuphigh[i]>MiddleOpenClose: Tupjump[i]=Tup[i]; Zupjump[i]=Zuphigh[i]+200; countNoClose+=1; State='NoClo'; Count=countNoClose            
            if Zuphigh[i]<MiddleOpenClose: Tupjump[i]=TEhigh+Tup[i]; Zupjump[i]=Zuphigh[i]; countNoOpen+=1; State='NoOpen'; Count=countNoOpen
        if pprint: print(i,j, int(Tupjump[i]), int(Zupjump[i]),'Direct', State, Count, ' ', 
                         int((Tupjump[i]-Tup[i])/1000), int(Zupjump[i]-Zupmid[i]))
        shift=windowsize_detect
        data=Zs[j0+shift:j1]; tdata=range(j0+shift,j1)
        dg1=gaussian_filter1d(data, windowsize_detect, order=1); dg1 /= np.abs(dg1).max()
        steps_dg1 = sd.find_steps(dg1, thres_detect)  #  print(type(steps_dg1),len(steps_dg1), steps_dg1[0])
        step_sizes, step_error = sd.get_step_sizes(data, steps_dg1, window=shift)  #      print(step_sizes, step_error)

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
            for k in steps_dg1: ax.plot([k+j0+shift,k+j0+shift], [-1., 1.], c='r', alpha=0.5, label='step')
            if ii==nbplot_detect-1 and SaveGraph:
                plt.savefig(path1+OutputFolder+figname+OutFormat)
                if CloseAfterSave: plt.close()

    if step_detectON and plot_detect:
        plt.figure(refname+'Direct_vs_Gaus1d_T', figsize=(6,6), dpi=100); ax3 = plt.gca()
        ax3.scatter((Tupjump-Tup)/1000, (Tupjump_-Tup)/1000, c='r', alpha=0.5)
        plt.figure(refname+'Direct_vs_Gaus1d_Z', figsize=(6,6), dpi=100); ax3 = plt.gca()
        ax3.scatter(Zupjump_-Zupmid, Zupjump-Zupmid, c='b', alpha=0.5)
    I0=Tupjump==0; Tupjump[I0]=np.nan; Zupjump[I0]=np.nan#; Zupmid[I0]=np.nan
    
    print('Direct countOpen=', countOpen, 'countNoClose=', countNoClose, 'countNoOpen=', countNoOpen, end='')
    print(' Total=', countOpen+countNoClose+countNoOpen)
    if not step_detectON:
        return (Tup, Dup, Tdown, Ddown, Tupdown, Tupjump, Tuphigh, Zupjump, Zuplow, 
                Zuphigh, Zupmid, Zupzero, countOpen, countNoClose, countNoOpen)
    if step_detectON:
        print('step_detect countOpen_=', countOpen_, 'countNoClose_=', countNoClose_, 'countNoOpen_=', countNoOpen_, end='')
        print(' Total_=', countOpen_+countNoClose_+countNoOpen_)
        return (Tup, Dup, Tdown, Ddown, Tupdown, Tupjump_, Tuphigh, Zupjump_, Zuplow, 
                Zuphigh, Zupmid, Zupzero, countOpen_, countNoClose_, countNoOpen_)

def FitHistoGaussian(Length, label, display, delta):

    def gauss_function(x, *p): return p[0]*np.exp(-(x-p[1])**2/(2*p[2]**2))

    bins=np.linspace(-1500, 2500 ,num=2001)
    Hy, Hx = np.histogram(Length, bins=bins, density=True)
    Hymax=np.argmax(Hy); xHymax=bins[Hymax]   #    print(Hymax, xHymax)
    indxfit=np.abs(Hx-xHymax)<delta
    Hxt=Hx[indxfit]; Hyt=Hy[indxfit[1:]]
    if len(Hxt)>len(Hyt): Hxt=Hxt[1:]
    pEq, pcovEq=curve_fit(gauss_function, Hxt, Hyt, p0=[1., xHymax,20])
    pEq[2]=np.abs(pEq[2])
    FitHy=gauss_function(Hxt, pEq[0], pEq[1], pEq[2])
    print(label+'Mod',xHymax, '  Gaussian fit: Amp=', "%.3f" % pEq[0], 'Avg=', "%.3f" % pEq[1] , 'SD=', "%.3f" % abs(pEq[2]))
    if display: AFSG.AFSG.MakeGraphXY(Hxt, Hyt, Hxt, FitHy, -1500, 2500, 0.00001, 0.02, label, 'density',
                            refname+'Hist'+label, SaveGraph, path1+OutputFolder, OutFormat, line2=True, log=True)
    return pEq[0], pEq[1], pEq[2], xHymax

def Spectrum(xx, axis, p1Zo, fs, label, display, fmin, fmax, axtable=None):
    friction0 = 6*np.pi*1.e-9*BeadRad       # units pN.s/nm
    if axis=='XY':
        friction = friction0 / ( 1 - (9/16)*(BeadRad/(p1Zo+BeadRad)) + (1/8)*(BeadRad/(p1Zo+BeadRad))**3 )
    elif axis=='Z':
        friction = friction0 / ( 1 - (9/8)*(BeadRad/(p1Zo+BeadRad)) + (1/2)*(BeadRad/(p1Zo+BeadRad))**3 )

    def FitSpectrum(x, *p): return p[1]/(2*np.pi**2)/( x**2 + (p[0]/(2*np.pi*friction))**2 )

    f, Pxx_spec = signal.periodogram(xx, fs, scaling='density') 
    Pxxt_spec= Pxx_spec[(f>fmin)&(f<fmax)]; ft=f[(f>fmin)&(f<fmax)]
    nbins=101; fbins=np.logspace(-2,2,nbins); Pbins=np.ones(nbins); dPbins=np.zeros(nbins)
    for m in range(nbins-1):
        u=Pxx_spec[(f>=fbins[m])&(f<fbins[m+1])]
        Pbins[m]=np.mean(u)
        dPbins[m]=np.std(u)/np.sqrt(len((u[~np.isnan(u)])))
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
        for axbis in wax: axbis.plot(ft, FitPxxt_spec, c='r', alpha=0.8)
        if SaveGraph: plt.savefig(path1+OutputFolder+figname+OutFormat)
        if CloseAfterSave: plt.close()
    print('Spectrum'+label, ' k (pN/nm)=',"%.5f" %  (pEq[0]),' D (nm²/s)=',"%.3f" %  (pEq[1]))
    return pEq[0], pEq[1], eEq[0], eEq[1], friction

def survival(Tupjump_Tuphigh, TEhigh, countNoOpen, shift, refname, color):
# calculates, fits and plots survival curve until TEhigh/1000.-shift

    def FitExp1(x, *p): return p[1]+(1-p[1])*np.exp(-x/p[0])
    def FitExp0(x, *p): return np.exp(-x/p[0])

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

##========================= ========================= =========================
##=========================   # MAIN LOOP  # ==================================
listresultsTotal=[]
color=iter(cm.rainbow(np.linspace(0, 1, 4*nbead*len(RangeNZlavg))))

if nbead>1: figAllSpec, axAllSpec = plt.subplots(nbead, 3, num='figAllSpec', figsize=(6, 2*nbead), dpi=100)

for ibead in range(nbead):        ## Loop on test beads
    print('======================BEAD ', b[ibead] ,' ===========================')
    print('NZavg=',NZavg, 'dzjump=', dzjump, 'tol=',tol, 'NZrefavg=',NZrefavg)
    if loaddata and load[ibead]: print('Loading data', fname[ibead]); tdms_file = TdmsFile.open(path1+fname[ibead])  # alternative TdmsFile.read(path1+fname[ibead])
    if loaddata and not load[ibead]: print('Data already loaded', fname[ibead])
    tdms_groups = tdms_file.groups(); # if pprint: print(tdms_groups)
    tdms_track = tdms_groups[0]; tdms_temp = tdms_groups[1]; tdms_event = tdms_groups[2]
    dmax=1000000.; imax=(len(tdms_track)-1)//5; dX=1000; dY=1000    
    print('Nb of traj=', imax)
    refname=fname[ibead]+'_Bead'+str(b[ibead])+'_'; print(refname)
    if GraphXYall and ibead==0:
        AFSG.MakeGraphXYall(tdms_track, imax, 100, 100, fname[ibead], path1+OutputFolder, OutFormat, SaveGraph)

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
    print('Power at    start', start[ibead],':', P0[start[ibead]], end=' ')
    print('Power at    stop:', stop[ibead],':', P0[stop[ibead]])
    if P0[start[ibead]]>0.1:
        wstart=10*np.arange(1000)+start[ibead]; wPstart=P0[wstart]
        i0=np.argmin(wPstart); start[ibead]= wstart[i0+5]
    if P0[stop[ibead]]>0.1:
        wstop=10*np.arange(1000)+stop[ibead]; wPstop=P0[wstop]
        i0=np.argmin(wPstop); stop[ibead]= wstop[i0+5]
    print('Power at new start', start[ibead],':', P0[start[ibead]], end=' ')
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
        print('Zcorrection by reference bead(s): Zravg=', "%.3f" % Zravg[0] )#int(bref[ibead]) )
        Zc=Zu-Zravg; Z0c=Zr0-Zr0avg; Zec=Z0-Zr0avg
    else: print('No Zcorrection'); Zc=Zu; Z0c=Zr0; Zec=Z0
    Xr0avg=Rollavg(Xr0, NZrefavg, mode); Yr0avg=Rollavg(Yr0, NZrefavg, mode) 
    Xravg=Rollavg(Xr, NZrefavg, mode); Yravg=Rollavg(Yr, NZrefavg, mode) 
    if RefBeadXYcorrection:
         print('XYcorrection by reference bead(s): Xravg=', "%.3f" % Xravg[0],' Yravg=', "%.3f" % Yravg[0] )
         Xc=X-Xravg; Yc=Y-Yravg; Xec=X0-Xr0avg; Yec=Y0-Yr0avg
    else: print('No XYcorrection'); Xc=X; Yc=Y; Xec=X0; Yec=Y0

    if GraphXYall and ibead==0:
        AFSG.MakeGraphXYall(tdms_track, imax, 100, 100, fname[ibead], path1+OutputFolder,
                            SaveGraph, OutFormat, SaveGraph, True, Xr0avg, Yr0avg)

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
    mask_anchor=(T0>=T0[range_anchorpoint[0]])&(T0<T0[range_anchorpoint[1]])
    mask_spectrum=(T0>=T0[range_spectrum[0]])&(T0<T0[range_spectrum[1]])
        
    p1Lo=p1Lo0; p1Lc=p1Lc0; p1Ll=p1Ll0; p1Zl=p1Zl0  # initialization of gaussian center for fit guess
    for NZlavg in RangeNZlavg:      # test loop for various choices of NZlavg
        if len(RangeNZlavg)>1: print('=============== NZlavg=', NZlavg, '===============')
        Zcavg=np.convolve(Zc, np.ones((NZlavg,))/NZlavg, mode=mode); dZ=Zc-Zcavg
        Zecavg=np.convolve(Zec, np.ones((NZlavg,))/NZlavg, mode=mode); dZe=Zec-Zecavg
         
        x = np.random.randn(len(P))/10
        Zch=Zc[indhigh]; Xch=Xc[indhigh]; Ych=Yc[indhigh]; dZh=dZ[indhigh]; Ph=P[indhigh]+x[indhigh]
        Zcl=Zcavg[indlow]; Xcl=Xc[indlow]; Ycl=Yc[indlow]; dZl=dZ[indlow]; Pl=P[indlow]+x[indlow]
        Zc0=Zcavg[ind0force]; Xc0=Xc[ind0force]; Yc0=Yc[ind0force]; dZ0=dZ[ind0force]
        Zec0=Zecavg[indeforce0]; Xec0=Xec[indeforce0]; Yec0=Yec[indeforce0]; dZe0=dZe[indeforce0]
        Zec0custom=Zecavg[mask_anchor]; Xec0custom=Xec[mask_anchor]; Yec0custom=Yec[mask_anchor]; dZe0custom=dZe[mask_anchor]

        Zcavg=np.convolve(Zc, np.ones((NZavg,))/NZavg, mode=mode)
        plt.figure('Various intervals', figsize=(6,6), dpi=100); ax = plt.gca()
        ax.plot(T, Zcavg, c='k', alpha=0.1, label='Zcavg')
        ax.scatter(T[indhigh], np.median(Zch)*np.ones(len(T))[indhigh], marker='.', s=3, alpha=0.5, label='high power')
        ax.scatter(T[indlow], np.median(Zcl)*np.ones(len(T))[indlow], marker='.', s=3, alpha=0.5, label='low power')
        ax.scatter(T[ind0force], np.median(Zc0)*np.ones(len(T))[ind0force], marker='.', s=3, alpha=0.5, label='0 force')
        ax.scatter(T0[indeforce0], np.median(Zec0)*np.ones(len(T0))[indeforce0], marker='.', s=3, alpha=0.5, label='0 force full range')
        ax.scatter(T0[mask_anchor], np.median(Zec0custom)*np.ones(len(T0))[mask_anchor], marker='.', s=8, alpha=0.5, label='mask_anchor')
        if Select_range_spectrum: ax.scatter(T0[mask_spectrum], np.zeros(len(T0))[mask_spectrum], marker='.', s=8, alpha=0.5, label='mask_spectrum')
        ax.legend(fontsize=10, markerscale=3)

        indstart=(T0<start[ibead]); Z0start=Z0c[indstart]
        Nl=len(Zc[indlow]); Nh=len(Zc[indhigh]); N0=len(Zc[ind0force]); Ntot=len(Zc)
        print('States Populations: low=', Nl, ' high=', Nh, ' 0force=', N0, ' total=', Ntot)        
        print('States Fractions: low=',"%.3f" %  (Nl/Ntot), ' high=',"%.3f" %  (Nh/Ntot), ' 0force=',"%.3f" %  (N0/Ntot), 'total=',"%.3f" %  ((Nl+Nh+N0)/Ntot))        
        for AnchorPointState in AnchorPointStateList:      # ['low force','0 force','custom']     
            print('=============== Anchor Point Reference:', AnchorPointState, end=' ')
            if AnchorPointState=='None':
                print('No anchoring point'); Zs=Zc; Xs=Xc; Ys=Yc
            else:
                if AnchorPointState=='low force':
                    Zc_=Zcl; Xc_=Xcl; Yc_=Ycl; print("") 
                elif AnchorPointState=='0 force':
                    NZlavg+=10; print("") 
                    Zc_=Zc0; Xc_=Xc0; Yc_=Yc0   
                    Zc_=Zec0; Xc_=Xec0; Yc_=Yec0
                elif AnchorPointState=='custom':
                    Zc_=Zec0custom; Xc_=Xec0custom; Yc_=Yec0custom
                    print(range_anchorpoint[0], range_anchorpoint[1])
                Zc_sort=np.sort(Zc_)       #     cleanedZclsort = [x for x in Zclsort if str(x) != 'nan']
                if len(Zc_sort)>0: Za=np.median(Zc_sort[np.arange(sample)])  # Za=np.amin(Zcl) #-BeadRad
                else:
                    print('Zc_ null wave:  Bad Anchor Point Reference Choice')
                    break
                print('Anchoring point at', AnchorPointState, '. Length calculation', 'min=', "%.3f" % np.min(Zc_), 'min_',sample, "%.3f" %Za)
                Xa=np.median(Xc_); Ya=np.median(Yc_)     
                CA=np.sqrt((Xc-Xa)**2+(Yc-Ya)**2+(Zc-Za)**2)
                Zs=Zc-Za; Xs=Xc-Xa; Ys=Yc-Ya
                Zs00=Zec-Za; Xs00=Xec-Xa; Ys00=Yec-Ya
            Xsl=Xs[indlow]; Ysl=Ys[indlow]; Zsl=Zs[indlow]; Tl=T[indlow]; MinLUTl=MinLUT[indlow]
            Xsh=Xs[indhigh]; Ysh=Ys[indhigh] ; Zsh=Zs[indhigh]; Th=T[indhigh]; MinLUTh=MinLUT[indhigh]
            Xs0=Ys[ind0force]; Ys0=Ys[ind0force]; Zs0=Zs[ind0force]       
            Length=np.sqrt(Xs**2+Ys**2+(Zs+BeadRad)**2)-BeadRad; Lengthl=Length[indlow]; Lengthh=Length[indhigh]; Length0=Length[ind0force]
            PullAngle=(180/np.pi)*np.arcsin( np.sqrt( (np.median(Xsh)-np.median(Xsl))**2 + (np.median(Ysh)-np.median(Ysl))**2  )/np.median(Lengthh))
            PullAngle2=(180/np.pi)*np.arcsin( np.sqrt( (np.median(Xsh)-np.median(Xs0))**2 + (np.median(Ysh)-np.median(Ys0))**2  )/np.median(Lengthh))            
            print('High length=', "%.3f" % np.median(Lengthh), 'Low length=', "%.3f" % np.median(Lengthl), ' Pulling Angle=', "%.3f" % PullAngle, ' Pulling Angle2=', "%.3f" % PullAngle2)

            if Display3D: AFSG.Plot3D(Xs,Ys,Zs, 30, refname+'XsYsZs', path1+OutputFolder, OutFormat, SaveGraph)
            if HistoZ: MakeHistoZ(Lengthh, Lengthl, 'Lengthh', 'Lengthl', -100,1200, refname+'LengthhLengthl')
            if GraphXYZ:
                AFSG.MakeGraphXY(Xsl, Ysl, Xsh, Ysh, -dX, dX, -dY, dY, 'XYl', 'XYh', refname+'XYslXYsh', SaveGraph, path1+OutputFolder, OutFormat)
                AFSG.MakeGraphXY(Tl, Zsl, Th, Zsh, 0, T[-1], -100, 1200, 'Xs', 'Ys', refname+'TZslZsh', SaveGraph, path1+OutputFolder, OutFormat)
                AFSG.SingleTrace(refname, tdms_track, -2, sampleGraph, Xs, Ys, Zs, T, MinLUT, SaveGraph, path1+OutputFolder, OutFormat )
                AFSG.MakeGraphXY(T, Length, T, Zs, 0, T[-1], 0, 1200, 'Length', 'Zs', refname+'Length', SaveGraph, path1+OutputFolder, OutFormat)
              
            print('=============== Opening analysis ===============')       
            Zsavg=np.convolve(Zs, np.ones((NZavg,))/NZavg, mode=mode)
            DeltaCOFit=0
            for i in [0]:         #  [0,1] first pass dzjump, 2nd pass DeltaCOFit
                dzjump_tmp = dzjump*(i==0) + DeltaCOFit*0.7*(i==1)
                print('Round ',i,' /0  dzjump=',dzjump_tmp)
                Pdown=0.01
                (Tup, Dup, Tdown, Ddown, Tupdown, Tupjump, Tuphigh, Zupjump, Zuplow, Zuphigh, Zupmid,Zupzero, countOpen, countNoClose,
                     countNoOpen)=FindJumps(T, Zs, Zsavg, NZavg, TE, DE, NE, Phigh, Pdown, TEhigh, tol, dzjump_tmp, GraphDhisto, refname, pprint)
                Tupdown_Tup=(Tupdown-Tup)/1000.; Tupjump_Tup=(Tupjump-Tup)/1000.
                Zupjump_Zuphigh=Zupjump-Zuphigh; Tupjump_Tuphigh=(Tupjump-Tuphigh)/1000.
                Zupjump_Zupmid=Zupjump-Zupmid
        
                indclosed=P<-10; indopen=P<-10      # labeling of tracking data with open and closed states
                for nj in range(len(Tupjump)):
                    indclosed=indclosed+(T<Tupjump[nj])*(T>Tuphigh[nj])*indhigh*(Tupjump[nj]>Tup[nj])
                    indopen=indopen+(T>Tupjump[nj])*(T<Tupdown[nj])*indhigh
                Zsclosed=Zs[indclosed]; Tclosed=T[indclosed]; Xsclosed=Xs[indclosed]; Ysclosed=Ys[indclosed]; Lengthclosed=Length[indclosed]
                Zsopen=Zs[indopen]; Topen=T[indopen]; Xsopen=Xs[indopen]; Ysopen=Ys[indopen]; Lengthopen=Length[indopen]
                Nc=len(Zsclosed); No=len(Zsopen)
                print('States Population: close=', Nc, ' open=', No, ' high=', Nh)
                print('States Fractions: close=',"%.3f" %  (Nc/Nh), ' open=', "%.3f" %  (No/Nh), ' (Nc+No)/Nh=', "%.3f" %  ((Nc+No)/Nh))
    
                print('=============== Gaussian fits ===============')
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
                DeltaCOFit=p1Zo-p1Zc   

            if DisplayCloseOpen: AFSG.MakeDisplayCloseOpen(refname, Tclosed, Zsclosed, Topen, Zsopen, Tl, Zsl, T, Zsavg, P, SaveGraph, CloseAfterSave, path1+OutputFolder, OutFormat)
            
            MakeHistoZ(Zsclosed, Zsopen, 'Zsclosed', 'Zsopen', -100, 1200, refname+'ZsclosedZsopen')

            print('=============== Power Spectra fits ===============')
            if Select_range_spectrum:
                Zsopen_=Zs00[mask_spectrum]; Xsopen_=Xs00[mask_spectrum];  Ysopen_=Ys00[mask_spectrum]
                print('Select_range_spectrum', range_spectrum[0], range_spectrum[1])
            else:
                Zsopen_=Zs[indopen]; Xsopen_=Xsopen;  Ysopen_=Ysopen
            if nbead>1: axtable0=axAllSpec[ibead,0]; axtable1=axAllSpec[ibead,1]; axtable2=axAllSpec[ibead,2]
            else: axtable0=None; axtable1=None; axtable2=None
            kx, Dx, dkx, dDx, frictionXY = Spectrum(Xsopen_, 'XY', p1Zo, fs,'XsopenRef'+AnchorPointState, DisplaySpectrum, fmin, fmax, axtable=axtable0) 
            ky, Dy, dky, dDy, frictionXY = Spectrum(Ysopen_, 'XY', p1Zo, fs,'YsopenRef'+AnchorPointState, DisplaySpectrum, fmin, fmax, axtable=axtable1)
            kz, Dz, dkz, dDz, frictionZ = Spectrum(Zsopen_, 'Z', p1Zo, fs,'ZsopenRef'+AnchorPointState, DisplaySpectrum, fmin, fmax, axtable=axtable2)
            Fx=kx*(BeadRad+p1Lo); dFx=np.abs(dkx*(BeadRad+p1Lo))+np.abs(kx*p2Lo)
            Fy=ky*(BeadRad+p1Lo); dFy=np.abs(dky*(BeadRad+p1Lo))+np.abs(ky*p2Lo)
            Fz=kz*(BeadRad+p1Lo); dFz=np.abs(dkz*(BeadRad+p1Lo))+np.abs(kz*p2Lo)
            if nbead>1: 
                axAllSpec[ibead,0].set_title(str(b[ibead])+' Fx='+"%.3f" % Fx, fontsize=6)
                axAllSpec[ibead,1].set_title(str(b[ibead])+' Fy='+"%.3f" % Fy, fontsize=6)
                axAllSpec[ibead,2].set_title(str(b[ibead])+' Fz='+"%.3f" % Fz, fontsize=6)
            print("Corner frequencies (Hz) x,y,z: ", "%.3f" %  (kx*Dx/4), "%.3f" %  (ky*Dy/4), "%.3f" %  (kz*Dz/4))
            print("Forces (pN) x,y,z: ", "%.3f" % Fx, "%.3f" % Fy, "%.3f" % Fz, 
                  "Forces Errors (pN) x,y,z: ", "%.3f" % dFx, "%.3f" % dFy, "%.3f" % dFz)
            Dxytheo=4/frictionXY; Dztheo=4/frictionZ
            print("D_by_Dtheo x, y,z: ", "%.3f" %  (Dx/Dxytheo), "%.3f" %  (Dy/Dxytheo), "%.3f" %  (Dz/Dztheo))
            FSDXFit=4*(p1Lo+BeadRad)/p2Xo**2; FSDXMod=4*(pLo+BeadRad)/p2Xo**2
            FSDYFit=4*(p1Lo+BeadRad)/p2Yo**2; FSDYMod=4*(pLo+BeadRad)/p2Yo**2    
 
            if GraphTZ_Power: AFSG.MakeGraphTZ_Power(refname, Tup, Zuplow, Zupmid, Zuphigh, Dup, T, P, SaveGraph, path1+OutputFolder, OutFormat)
            if GraphTupdown: AFSG.MakeGraphTupdown(refname, Tup, Tupdown, Tupjump, Tupdown_Tup, SaveGraph, path1+OutputFolder, OutFormat)
            if CalculateSurvival:
                (pE0, pE1, x2, y2, y2_fit)=survival(Tupjump_Tuphigh, TEhigh, countNoOpen, tshiftsurvival, refname, next(color))
            else:
                pE0=1.; pE1=0
     #       if CalculateSurvival: (pE0, pE1, x2, y2, y2_fit)=survival(Tupjump_Tuphigh, TEhigh, 2., '', next(color))   # to plot all beads on th same figure
            if DisplayJumpGraphs:
                AFSG.JumpGraphs(refname, Zs, Zsavg, Zupjump_Zupmid, Tupjump_Tuphigh, Tup, Zuplow, Tupjump, Zupjump, Zupmid,
                           Tuphigh, Zuphigh, Zupzero, Dup, Phigh*(1+tol), T, P, Tclosed, Zsclosed, Topen, Zsopen,
                           Tl, Zsl, SaveGraph, CloseAfterSave, path1+OutputFolder, OutFormat)
        
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
    AFSG.plotresults(refname, df, SaveGraph, path1+OutputFolder, OutFormat)
    
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
    
    AFSG.MakeDisplayplotvariableNZlavg(DisplayplotvariableNZlavg, ax1, ax2, dfT, c, ib, DeltaCOFit, DeltaCOMod)
    
    for axi in [ax2[2], ax2[0], ax2[1]]: axi.plot(np.linspace(0, 2, 10), np.linspace(0, 2, 10), c='k', alpha=0.5, label=' y=x')    

    for i in range(nFig1):
        ax1[i].legend(fontsize=6)
        if SaveGraph: plt.figure(listFig1[i]); plt.savefig(path1+OutputFolder+fname[0]+listFig1[i]+OutFormat)
    for i in range(nFig2):
        ax2[i].legend(fontsize=6)
        if SaveGraph: plt.figure(listFig2[i]); plt.savefig(path1+OutputFolder+fname[0]+listFig2[i]+OutFormat)
