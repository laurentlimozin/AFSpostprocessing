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
20201009 v3.9   Permits to select the first long step for the spectrum calculation (set AnchorPointStateList=['0_force'] )
                Take into account the NoOpen cases in the survival curve
20201025 v3.91  Revised NoOpen/NoClose cases with introduction of a threshold MiddleOpenClose=900nm
                Some discrepancy for step_detectON=0 or 1 (which may depend on windowsize_detect and NZavg)
20201206 v4.0   Full revision of code; All graphic functions shifted to AFS_JDNA_Graphs
                Adapt to last version of nptdms package
                Option of custom interval for anchor point and power spectrum (cf graph 'Various intervals')
                    range_anchorpoint=(start,stop) range_spectrum=(start,stop) and boolean Select_range_spectrum
20201217 v4.1   Option of reference trace calculated on low_force cycle positions: set boolean selfreference=True
20201221 v4.11  Reads file x_x_param.csv for paramterization of analysis: set loadparameter=True
                Save FitsResults.csv after each bead
         v4.12  Tuning of step_detect and 2 iterations on FindJumps to discriminate low and high states 20210105
20210105 v4.13  Integrates Jim's add-ons to compute XY fluctuations
20210108 v4.14  Contour plot for fluctuations (uses updates AFS_JDNA_Graphs)
20210109 v4.15  Add (mean) center on contour plots; correct bug to use fmin_Hz, fmax_Hz from param.csv
20210111 v4.16  Figures Transparent bgd. Superimposes centers of fluctuations contour plot;
                list of contour data saved as contourdata for each figure , containing the contour number, its length and center of mass
20210112 v4.17  Jim implements the total extension and the pulling angle based on the anchor point from the contour.
                improved MakeHistoP() function using signal.find_peaks
20210113 v4.18  Jim typo corrected - set nb of contours in XY maps = levels=20 
20210114 v4.19  Labelling of cycles and possibility to eliminate some selected in the param file 
20210117 v4.20  Refinement of step detection with weightsecondpass_step = 0.7 (if ==0 step guess zjump imposed )
                nstepmax=5  if more steps are detectred this is considered as noclose or noopen
20210118 v4.21  Modification BuildZinterpol() for selfreference
20210119 v4.22  Export of spectrum and fit as csv file. Correct bugs name and BeadRadius in Resultsfile.csv
20210120 v4.23  several test beads in parameter file works now  / standardisation of file names for figures and data export
                export of close state duration and cycle number / center of mass and center of outer plot on XY fluctuation graph
20210120 v4.24  Improved baseline for self reference BuildZinterpol()             
20210121 v4.25  Modified detection of high force steps to handle multiple power traces used for force calibration : set boolean forcecalib=True
                correct bug in BuiltP (last cycle missed) and introduces selfreference for X,Y
20210123 v4.26  Correct print / calculation of corner frequency
                Correct X,Y baseline for selfreference case, using position at high force steps
                manual anchor point: default manual_anchor_point=False , coordinates_anchor_point=(0,0)
                anchor point from contour center calculated in interval range_anchorpoint_min; set AnchorPointStateList=custom
                forcecalib=True to analyze force calibration tdms file
                    (each power has to be measured separately by setting range_anchorpoint_min and range_spectrum_min)
20210125 v4.27  graph of Y selfreference; modify Buildinterpol() for backward or forward filling of first stretch
                adjust origin of baseline corrected data and raw data
20210126 v4.28  save results.csv in mode forcecalib=True
20210127 v4.29  modify origin of baseline corrected data and raw data x,y,z; Spectrum fit function=Sitter2015
                correction line 689 Zec=Z0-Z10avg   # corr LL 27012021
                correction 
20210129 v4.30  radial fluctuations r; spectrum 1D (Schäffer) or 2D (Sitter)      
20210201 v4.31  modify Buildinterpol() keeping first/last cycle (Z, Outer=True) or interpolating them (X,Y Outer=False) 
20210203 v4.32  top graph option in MakeDisplayCloseOpen; print jump infos in history without rounding
20210210 v4.33  spectrum with r fluctuations ; layout with x,y,z,r spectra for all beads.
20210210 v4.34  Three command lines (640-642) added to exclude the segment containing spikes at the end of a trace
                safety interval for gradient calculation in time()
20210407 v4.35  covariance matrix for AnchorPoint XY
                reference bead added to selfreference graph for comparison (if selfreference=True)
                phi pulling angle in horizontal plane polar coordinates (PullPhi and PullPhiC)
                spectrum on closed state: fitting parameters kx2, ky2, Dx2, Dy2; Forces Fx2, Fy2, Fx2C, Fy2C
                    estimate of force F3 based on difference open/closed state
                    to do:  check line 697 and work on self reference with file 20200824-201201
20210409 v4.36 SymmetryFactor = sqrt(lambdamax/lambdamin) with lambdas eigenvalues of covariance matrix
20210413 v4.37 Correction signe F3
20210514 v4.38 Fix bug duration countNoOpen, by extending maximal step duration to 1000 s; adjusted guess for tau survival
20210614 v4.39 calculates Allan deviation: https://allantools.readthedocs.io/en/latest/readme_copy.html
                select AllenComputation = True and a typical range nAllen(nstart, nend)
               bypass calibration in closed state for forcecalib case
20211021 v4.41 Accounts for cases where Plow=0
20211023 v4.42 Improves detection of Plow/Phig and Plow/Phigh
20211025 v4.43 Corrects bug in tdmsfilename type, in output paramaters if CalculateSurvival==False
20211106 v4.44 force calibration runs with power as low as 0.2%
20211108 v4.45 Variable RefBeadZcorrection in param.csv 
20211109 v4.46 Select intersection of spectrummask and highpower range for measuring length in force calibration mode

20220314 v4.55 Exports raw spectra for calibration mode
20220314 v4.56 Exports raw spectra for calibration mode and cyclic mode
20220408 v4.57 Error management in fits. Git push
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
import cv2
import allantools
import AFS_JDNA_Graphs as AFSG
from scipy.ndimage.filters import gaussian_filter1d, uniform_filter1d
import warnings
import sys
import ast
#import dataframe_image as dfi ; dfi.export(df, 'dataframe.png')

fontsizelarge = 14
plt.rc('font', size=fontsizelarge)

warnings.filterwarnings("ignore")
np.random.seed(444); plt.close('all')
#def pfl(f): return '{:.2E}'.format(f)

#path1="C:/DATA/Experiences Diverses/AFS/Data/Debug20201002/"      #☺ data folder
path1="/home/laurent/DATA/Experiences Diverses/AFS/Data/" #☺ data folder  #path1="D:/Jim/JDNA/2020-08-24/FOV3/"  
#path1="D:/Jim/JDNA/Rapamycin/2019-11-21-25/FOV4/4%/"
#path1="D:/Jim/JDNA/Nef19/2020-12-28_jDNA_3um silica beads/"
#path1="D:/Jim/JDNA/Rapamycin/2020-09-01/"
#path1="D:/Jim/JDNA/Rapamycin/2020-08-24/FOV1/"
#path1="D:/Jim/JDNA/Rapamycin/2020-12-08_1.5m silica beads/FOV2/"
#path1="E:/Jim/2020-12-08_1.5m silica beads/FOV2/"
#tdmsfilename="20200901-194315"          #tdmsfilename="20200903-144604.tdms"    tdmsfilename="20200826-120242.tdms"
#tdmsfilename="20191124-151720"          #tdmsfilename="20200903-144604.tdms"    tdmsfilename="20200826-120242.tdms"
#tdmsfilename="20210102-1013125" #tdmsfilename="20201231-204406" #tdmsfilename="20201203-194431"
#tdmsfilename="20200828-154633"
#tdmsfilename="20201208-170730"
#tdmsfilename="20201208-210756" #tdmsfilename="20201210-1002859"
#tdmsfilename="20201210-151130"
#tdmsfilename="20201230-121830"
#tdmsfilename="20201210-170702"
#tdmsfilename="20201210-123017"
#tdmsfilename="20201211-131102"
#tdmsfilename="20201214-1074703"
#tdmsfilename="20200825-1003446"
#tdmsfilename="20200825-183857"
tdmsfilename="20200825-114338"
#tdmsfilename="20210423-220614"
#tdmsfilename="20210424-1071541"
#tdmsfilename="20210423-152842"
#tdmsfilename="20210928-153622"
#tdmsfilename="20210928-200240"
#tdmsfilename="20211020-182424"
#tdmsfilename="20210929-225744"
#tdmsfilename="20210929-191511"
#tdmsfilename="20211008-205604"
#tdmsfilename="20211015-154115"
#tdmsfilename="20211015-225305"
#tdmsfilename="20211015-150329"
#tdmsfilename="20211007-123758"
#tdmsfilename="20211009-162627"
#tdmsfilename="20211001-120231"
tdmsfilename="20210929-005129"
#tdmsfilename="20200903-231501"
tdmsfilename="20211008-205604"
tdmsfilename="20211001-184602"
#tdmsfilename="20211007-130743"
tdmsfilename="20200824-201201"

OutputFolder='Outputtest/'
if not os.path.exists(path1+OutputFolder): os.makedirs(path1+OutputFolder)

bref=[6,28,49]      # reference beads  bref=[17,29,43]               # reference beads
b=[17, 78, 79, 1, 1, 1]               # test beads  b=[6, 8, 21, 41, 46, 58, 73, 91]     #b=[58, 73, 91] b=[73, 91] b=[5]

bad_cycles=[]
#RangeNZlavg=[20,40,60, 80, 100, 120, 140, 160, 200, 240, 280, 320, 400, 480, 560, 640, 720, 800, 900, 1000]
RangeNZlavg=[400]       # size of moving average for smoothing z before anchor point determination

# initialiation for various intervals
nbead=len(b); fname=[""]*nbead; start=np.zeros(nbead, dtype=int); stop=np.zeros(nbead, dtype=int); load=np.zeros(nbead, dtype=int)
start_min=np.zeros(nbead, dtype=int); stop_min=np.zeros(nbead, dtype=int)
#for nf in range(nbead): fname[nf] = tdmsfilename; start[nf]=15000; stop[nf]=1600000; load[nf]=nf==0  # stop[nf]=1600000
for nf in range(nbead): fname[nf] = tdmsfilename; start[nf]=30*60*1000/18; stop[nf]=175*60*1000/18; load[nf]=nf==0  # stop[nf]=1600000
range_anchorpoint_min=(0,100); range_spectrum_min=(0, 100); Select_range_spectrum=False
#for nf in range(nbead): fname[nf] = tdmsfilename+".tdms"; start[nf]=6*60*1000/20; stop[nf]=306*60*1000/20; load[nf]=nf==0  # stop[nf]=1600000
for nf in range(nbead): fname[nf] = tdmsfilename+".tdms"; start_min[nf]=6; stop_min[nf]=306*60*1000/20; load[nf]=nf==0  # stop[nf]=1600000

##========================= ========================= =========================
# General settings
loaddata=True  # to load the tdms file
loadparameter=True  # to load the analysis parameter file
forcecalib=False   # to analyze series of single steps with variable power
CleanOutliers=True
RefBeadZcorrection=True # correction with reference bead or selfreference
RefBeadXYcorrection=False # correction with reference bead or selfreference
AnchorPointStateList=['low_force'] # correction with AnchorPoint as determined by anchor point
        #   'low_force' or '0_force' or 'custom' ; or 'None': no anchor point 
manual_anchor_point=False; coordinates_anchor_point=(0,0)
export=False # to produce csv files of raw data (for Charlie)
CalculateSurvival=True  # calculate and draw survival curve
tshiftsurvival=0; 
step_detectON=True; plot_detect=True; nbplot_detect=15; thres_detect=0.5; windowsize_detect=25 # parameters for gaussian step detection
dzjump=120.; weightsecondpass_step = 0.7; nstepmax=5
NZavg=100; NZlavg=400; NZrefavg=5000; tol=0.01; sample=1; sampleGraph=100 # parameters for direct step detection
MiddleOpenClose=900  #  900
MaxZHistoPlot=2500
BeadRad=1500. # Bead Radius in nm
Temperature_C = 25; kBT_pN_nm= 1.38e-23*(Temperature_C+273)*1.e12*1.e9
fmin_Hz=0.1; fmax_Hz=15.   # frequency range for spectrum fit (Hz)
p1Lo0=1000; p1Lc0=800; p1Ll0=500; p1Zl0=500  # guess for histogram fit
selfreference=False
AllenComputation = False; nAllen=(1,400000); tauAllen = np.logspace(0, 2, 50); stepsAllen=[1,2,3,4]
mode='reflect' #'nearest'  # mode for rolling average   # edge modes = ['full', 'same', 'valid', 'nearst']
# see details on https://stackoverflow.com/questions/13728392/moving-average-or-running-mean

# Print and display options
pprint=False; # printing detailed infos in console window
pprintstep=True
SaveGraph=True; OutFormat=".png" ; CloseAfterSave=False # alternative ".jpg" or ".pdf" 
GraphTdms=False  # graph of frame rate and force data (duration and amplitude of 2 steps)
GraphDhisto=False; GraphTupdown=False; GraphTZ_Power=True
HistoZ=True; GraphXYZ=True # characteristics of trajectory of the selected beads
GraphXYall=False; iminXY=0; imaxXY=0 # display xy trajectories with number on a single graph
DisplaySpectrum=True
DisplayCloseOpen=True
DisplayGaussianFits=False
DisplayJumpGraphs=False # graphs of jump analysis
DrawSuperposal=False # graph of superposition of jump traces
DisplaySurvival=True
Display3D=False
DisplayplotvariableNZlavg=False

param_name=tdmsfilename+"_param.csv"
#param_name=tdmsfilename+"_param_force_calibration.csv"
#param_name=tdmsfilename+"_param_cyclic.csv"

if loadparameter and os.path.isfile(path1+param_name):
    print('load '+param_name)
    dfparam = pd.read_csv(path1+param_name, encoding = "ISO-8859-1")
    print(pd.DataFrame(dfparam[['parameter','value','type']]))

    for i, row in dfparam.iterrows():
        p = row['parameter']; v = row['value']; t = row['type']
  #      print(i, p, type(p), v, type(v), t, type(t))
   #     if t == 's': exec("%s = %s" % (p,v))  # str added on 20211025
        if t == 's': exec("%s = '%s'" % (p,v))  # str added on 20211025
        if t == 'd': exec("%s = int(%s)" % (p,v))
        if t == 'f': exec("%s = float(%s)" % (p,v))
        if t == 'l': exec("%s = '%s'.split()" % (p,v))
        if t == 't': exec("%s = tuple(%s)" % (p,v))
     #   exec("print(%s , type(%s))" %(p,p))
Select_range_spectrum = Select_range_spectrum == 'True'
CalculateSurvival = CalculateSurvival == 'True'   # 20211025 transforms strings in booleans
selfreference = selfreference == 'True'
forcecalib = forcecalib == 'True'
manual_anchor_point = manual_anchor_point == 'True'
RefBeadXYcorrection = RefBeadXYcorrection == 'True'
RefBeadZcorrection = RefBeadZcorrection == 'True'
tdmsfilename = tdmsfilename.replace('-', '_')
for i, bi in enumerate(b): b[i]=int(bi) 
for i, bi in enumerate(bref): bref[i]=int(bi)
for i, bi in enumerate(bad_cycles): bad_cycles[i]=int(bi)
#print(tdmsfilename, type(tdmsfilename))
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
    print('Build Power trace:', len(TP), 'events') #   print(TP)
    P=np.zeros(len(T)); j1=0
    for i in range(1,len(TP)):
        j=int(np.where(T==TP[i])[0][0])
        P[j1:j]=DP[i-1]; j1=j
    return P

def Buildinterpol(Z, T, TE, NE, DE, Outer=True, MinDDP=0.1):      # builds power time trace form event times
    I=(NE=='SigGen.Power (%)')
    DP=DE[I]; TP=TE[I]
    DDP=np.diff(DP); DDP = np.append(DDP,MinDDP)
    IDP = np.abs(DDP)>=MinDDP
    backstep = 2000
#    print(len(DP), len(DDP))
    TP=TP[IDP]
    print('Buildinterpol:', len(TP), 'events') #   print(TP)
    Z1=np.zeros(len(T)); j1=0
    for i in range(1,len(TP)):
        j = int(np.where(T==TP[i])[0][0])
        if i<len(TP)-1: jmax = int(np.where(T==TP[i+1])[0][0])
        if pprint: print(i, T[j1], T[j], TP[i], DP[i])
        if i==1:
            if Outer: Z1[j1:j]=Z[j1:j]
            if not Outer: Z1[j1:j]=np.mean(Z[j+1: jmax])
        else:
            Z1[j1:j]=np.mean(Z[max(0, j1-backstep): j1-1])
  #          Z1[j1:j]=np.mean(Z[jmin: j1-1])
            if i==len(TP)-1:
                if Outer: Z1[j:len(Z1)]=Z[j:len(Z)]
                if not Outer: Z1[j:len(Z1)]=np.mean(Z[max(0, j-backstep): j-1])
        j1=j
    return Z1

def time(d0, Graph):    # retrieves time trace from tdms file and measures time step
    T = d0['Time (ms)'][:]
    print('Total points:', len(T))
    dT=np.gradient(T)[(np.gradient(T)>0) & (np.gradient(T)<1000)]
    print('Total points gradient:', len(dT))
    avfreq = 1000/np.mean(dT)
    print('dTrawhisto (frame time ms)', np.amin(dT), np.amax(dT), ' average frequency (Hz)', "%.3f" % avfreq)   
    if Graph:
        figname=refname+'dTrawhisto'
        plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
        ax=sns.distplot( dT, kde=False, bins=np.linspace(0.,20.,num=200) )
        ax.set_xlabel("Time step (ms)");  ax.set_ylabel("Number"); ax.set_yscale('log')
        plt.ylim(0.1, 1.e7)
        if SaveGraph: plt.savefig(path1+OutputFolder+figname+OutFormat, transparent=True)
    return T, avfreq

def shortentrack(X0, Y0, Z0, T0, P0, MinLUT0, i0,j0):   # cut trajectories between 0 and i0  # TO OPTIMIZE !
    X=np.delete(X0, range(i0));  X=np.delete(X, range(j0-i0,len(X0)-i0))
    Y=np.delete(Y0, range(i0));  Y=np.delete(Y, range(j0-i0,len(Y0)-i0))
    Z=np.delete(Z0, range(i0));  Z=np.delete(Z, range(j0-i0,len(Z0)-i0))
    T=np.delete(T0, range(i0));  T=np.delete(T, range(j0-i0,len(T0)-i0))
    P=np.delete(P0, range(i0));  P=np.delete(P, range(j0-i0,len(P0)-i0))
    MinLUT=np.delete(MinLUT0, range(i0));  MinLUT=np.delete(MinLUT, range(j0-i0,len(MinLUT0)-i0))
    return (X,Y,Z,T,P, MinLUT)

def period(TE, NE, DE, Graph, refname, height):  # retrieves low power and high power duration in case of binary power
    TE_ = TE[(NE=='SigGen.Power (%)')]
    DE_ = DE[(NE=='SigGen.Power (%)')]
    dTE = -TE_[:-1] + TE_[1:]  # dTE=np.gradient(TE)
    dDE = -DE_[:-1] + DE_[1:]  # dTE=np.gradient(TE)
    dTEUp = dTE[dDE>0]
    dTEDo = dTE[dDE<0]
    bins=np.linspace(0, 10000000 , num=10001)
    histUp, bin_edges = np.histogram(dTEUp,bins=bins)
    peaksUp, propsUp  = signal.find_peaks(histUp, height=height)
    histDo, bin_edges = np.histogram(dTEDo,bins=bins)
    peaksDo, propsDo  = signal.find_peaks(histDo, height=height)
    # print(histUp, histDo)
    # print(np.unique(histUp), np.unique(histDo))
    # print(peaksDo, peaksUp)
    if Graph:
        figname=refname+'dTEventsHisto'
        plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
        ax=sns.distplot( dTE, kde=False, bins=np.linspace(0.,1000000.,num=1001) )
        ax.set_xlabel("Event Time step (ms)");  ax.set_ylabel("Number"); ax.set_yscale('log')
        plt.ylim(0.1, 1.e3)
        if SaveGraph: plt.savefig(path1+OutputFolder+figname+OutFormat, transparent=True)
    TEhigh = bin_edges[peaksDo[0]]; TElow = bin_edges[peaksUp[0]]
    print('period detection with height=', height)
    print('TEhigh (ms):', TEhigh, "TElow (ms):", TElow)
    return TEhigh, TElow

def MakeHistoP(P, NE, DE, Graph):       # histogram of power trace and retrieval of low/high values for binary power
    bins=np.linspace(-0.001, 100 , num=100002) # bins=np.linspace(-1, 100 ,num=101001)
    if Graph:
        figname=refname+'Phisto'
        plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
        ax=sns.distplot( P, kde=False, bins=bins )
        ax.set_xlabel("Power (%)");  ax.set_ylabel("Number")
        if SaveGraph: plt.savefig(path1+OutputFolder+figname+OutFormat, transparent=True)
    DE_ = DE[(NE=='SigGen.Power (%)')]
    hist, bin_edges = np.histogram(DE_, bins=bins)    #   np.bincount(P)       
    peaks , _ = signal.find_peaks(hist, height=4)
    hpeaks = bins[peaks]; npeaks = hist[peaks]
    print(hpeaks)
 #   u = np.sort(np.unique(P)); Plow = u[-2]; Phigh = u[-1]
    Plow = np.min(bins[peaks]); Phigh = np.max(bins[peaks])
    print('Phigh:', Phigh, "Plow:", Plow)
    return Phigh, Plow

def Allen(y, x, nAllen, name, steps):
    print('Allen computation for '+name)   
    figname=refname+"Allen_"+name; plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
    for s in steps:
        if nAllen[0]+(s+1)*nAllen[1]<len(y):
            ya = y[nAllen[0]+s*nAllen[1]:nAllen[0]+(s+1)*nAllen[1]]
            (t2, ad, ade, adn) = allantools.oadev(ya, rate=fs, data_type="freq", taus=x)  # Compute the overlapping ADEV
 #       (t2, ad, ade, adn) = allantools.oadev(ya, rate=r)  # Compute the overlapping ADEV
#    (t2, ad, ade, adn) = allantools.oadev(ya, rate=fs, data_type="freq", taus=ta)  # Compute the overlapping ADEV
            ax.loglog(t2/fs, ad)
    ax.set_ylabel(name+" axis Allan Deviation (nm)");  ax.set_xlabel("tau (s)"); plt.tight_layout()
    if SaveGraph: plt.savefig(path1+OutputFolder+figname+OutFormat)

def Rollavg(Z, N, mode):      # rolling average on N time steps
 #   Zavg=np.convolve(Z, np.ones((N,))/N, mode=mode)
    Zavg = uniform_filter1d(Z, size=N, mode=mode)
 #   n=int(N/2.); Zavg[range(n)]=Zavg[n+1]; Zavg[range(-1,-n-1,-1)]= Zavg[-n]
 #   if mode=='valid': Z=Z[0: len(Z)-N+1];# Z=Z[0: len(Z)-N+1]
    return Z, Zavg

def exporttrack(X,Y,Z,T,P, datasave):    # exports raw traces as csv file
    file = open(datasave, "w")
    file.write('X_nm'+','+'Y_nm'+','+'Z_nm'+','+'T_ms'+','+'P_%'+'\n')
    for i in range(len(X)):
        file.write(str(X[i])+','+str(Y[i])+ ','+str(Z[i])+ ','+str(T[i])+ ','+str(P[i])+"\n")
    file.close()

def exportsurvival(T, datasave):    # exports raw traces as csv file
    file = open(datasave, "w")
    file.write('i_cycle'+','+'duration_s'+'\n')
    for i in range(len(T)):
        file.write(str(i)+','+str(T[i])+"\n")
    file.close()

def exportspectrum(ft, FitPxxt_spec, Xfbins, Pbins, dPbins, datasave, colname = ''):    # exports raw traces as csv file

    file = open(datasave, "w")
    if colname == '': colname='f_fit_Hz,P_fit,Xfbins_Hz,Pbins,dPbins'
    file.write(colname+'\n')
    for i in range(len(ft)):
        line = str(ft[i])+ ','+str(FitPxxt_spec[i])
        if i<len(Xfbins): 
            line += ',' +str(Xfbins[i])+',' + str(Pbins[i])+ ','+str(dPbins[i])+"\n"
        else:
            line += "\n"
        file.write(line)
    file.close()
    
def from_np_array(array_string):
    array_string = ','.join(array_string.replace('[ ', '[').split())
    return np.array(ast.literal_eval(array_string))

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
    if SaveGraph: plt.savefig(path1+OutputFolder+figname+OutFormat, transparent=True)
    
def FindJumps(T, Zs, Zsavg, NZavg, TE, DE, NE, Phigh, Plow, TEhigh, tol, dzjump, GraphDhisto, refname, pprint):
    ISup=(NE=='SigGen.Power (%)')*(Phigh*(1-tol)<DE)*(DE<Phigh*(1+tol))
    if forcecalib: ISup=(NE=='SigGen.Power (%)')*(Phigh*(1-tol)<DE)
    ISdown=(NE=='SigGen.Power (%)')*(Plow*(1-tol)<=DE)*(DE<=Plow*(1+tol))
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
       #     Zs1avg=np.convolve(Zs1, np.ones((NZavg,))/NZavg, mode='same')
            Zs1avg=uniform_filter1d(Zs1, size=NZavg, mode=mode)
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

 #       Zsavg[j0+NZavg]
 #       if len(steps_dg1)==1:# and step>dzjump:
        if len(steps_dg1)<nstepmax and len(steps_dg1)>0:
            idx = np.argmax(step_sizes)
            j=steps_dg1[idx]+j0+shift
            Tupjump_[i]=T[j]; Zupjump_[i]=step_sizes[idx]+Zupmid[i]; countOpen_+=1; Count_=countOpen_; State_='Open '
        else:
            Zuphigh[i]=np.mean(Zs[range(j0,j0+int(TEhigh/1000*fs))])
            if pprint: print(i, 'multiple steps. z=', Zuphigh[i],'of size:', step_sizes)
       #     if Zuphigh[i]>=MiddleOpenClose: Tupjump_[i]=Tup[i]-10000; Zupjump_[i]=Zuphigh[i]+200; countNoClose_+=1; State_='NoClo'; Count_=countNoClose            
            if Zuphigh[i]>=MiddleOpenClose: Tupjump_[i]=Tup[i]; Zupjump_[i]=Zuphigh[i]+200; countNoClose_+=1; State_='NoClo'; Count_=countNoClose            
            if Zuphigh[i]<MiddleOpenClose: Tupjump_[i]=TEhigh+Tup[i]; Zupjump_[i]=Zuphigh[i]; countNoOpen_+=1; State_='NoOpen'; Count_=countNoOpen
#        if pprintstep: print(i,j, int(Tupjump_[i]), int(Zupjump_[i]),'Gaus1d', State_, Count_, len(steps_dg1),
 #                        int((Tupjump_[i]-Tup[i])/1000), int(Zupjump_[i]-Zupmid[i]))
        if pprintstep: print(i,j, Tupjump_[i], Zupjump_[i],'Gaus1d', State_, Count_, len(steps_dg1),
                         (Tupjump_[i]-Tup[i])/1000, Zupjump_[i]-Zupmid[i])
        if step_detectON and plot_detect:
            ii=i%nbplot_detect                
            if ii==0: 
                figname=refname+'step_detect'+str(i)
                fig, axii= plt.subplots(nbplot_detect, 1, num=refname+'step_detect'+str(i), figsize=(nbplot_detect*1,8), dpi=100)
            ax=axii[ii]; ax2 = ax.twinx()
            ax.plot(tdata, dg1, c='b', alpha=0.7, label=str(len(steps_dg1))); ax.legend(loc="upper right", fontsize=6)
            ax2.plot(tdata, data, c='k', alpha=0.5, label='data'+str(i)); ax2.legend(loc="upper left",fontsize=6) 
            ax.set(ylim=(-1., 1.)); ax2.set(ylim=(400., 1000.))
            for k in steps_dg1: ax.plot([k+j0+shift,k+j0+shift], [-1., 1.], c='r', alpha=0.5, label='step')
            if ii==nbplot_detect-1 and SaveGraph:
                plt.savefig(path1+OutputFolder+figname+OutFormat, transparent=True)
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

    bins=np.linspace(-15000, 25000 ,num=20001)
    Hy, Hx = np.histogram(Length, bins=bins, density=True)
    Hymax=np.argmax(Hy); xHymax=bins[Hymax]   #    print(Hymax, xHymax)
    indxfit=np.abs(Hx-xHymax)<delta
    Hxt=Hx[indxfit]; Hyt=Hy[indxfit[1:]]
    if len(Hxt)>len(Hyt): Hxt=Hxt[1:]
    while True:
         try:
             pEq, pcovEq=curve_fit(gauss_function, Hxt, Hyt, p0=[1., xHymax,20])
             break
         except RuntimeError:
             print("No convergence"); pEq=[0., 0., 1.]
             break
         except TypeError:
             print("Improper input"); pEq=[0., 0., 1.]
             break
         except ValueError:
             print("Data contains Nan"); pEq=[0., 0., 1.]
             break            
    pEq[2]=np.abs(pEq[2])
    FitHy=gauss_function(Hxt, pEq[0], pEq[1], pEq[2])
    print(label+'Mod',xHymax, '  Gaussian fit: Amp=', "%.3f" % pEq[0], 'Avg=', "%.3f" % pEq[1] , 'SD=', "%.3f" % abs(pEq[2]))
    if display: AFSG.MakeGraphXY(Hxt, Hyt, Hxt, FitHy, -15000, 25000, 0.00001, 0.02, label, 'density',
                            refname+'Hist'+label, SaveGraph, path1+OutputFolder, OutFormat, line2=True, log=True)
    return pEq[0], pEq[1], pEq[2], xHymax

def Spectrum(xx, axis, p1Zo, fs, label, display, fmin_Hz, fmax_Hz, exportname=None, axtable=None):
    friction0 = 6*np.pi*1.e-9*BeadRad       # units pN.s/nm
    if axis=='XY' or axis=='R':
        friction = friction0 / ( 1 - (9/16)*(BeadRad/(p1Zo+BeadRad)) + (1/8)*(BeadRad/(p1Zo+BeadRad))**3 )
    elif axis=='Z':
        friction = friction0 / ( 1 - (9/8)*(BeadRad/(p1Zo+BeadRad)) + (1/2)*(BeadRad/(p1Zo+BeadRad))**3 )

#    if axis=='XY' or axis=='Z':
      #  def FitSpectrum(x, *p): return p[1]/(np.pi**2)/( x**2 + (p[0]/(2*np.pi*kBT_pN_nm/p[1]))**2 )   # corrected 28/01/2022 # Sitter 2015 one sided version
    def FitSpectrum(x, *p): return 4*kBT_pN_nm**2/(p[1]*p[0]*p[0]) * 1/( 1+ (x*2*np.pi*kBT_pN_nm/(p[1]*p[0]))**2 )  # corrected 28/01/2022 # Daldrop 20215 one sided
#    if axis=='R':
     #   def FitSpectrum(x, *p): return p[1]/(np.pi**2)/( x**2 + (p[0]/(2*np.pi*kBT_pN_nm/p[1]))**2 )   # corrected 28/01/2022 # Sitter 2015 one sided version
#        def FitSpectrum(x, *p): return 4*kBT_pN_nm**2/(p[1]*p[0]*p[0]) * 1/( 1+ (x*2*np.pi*kBT_pN_nm/(p[1]*p[0]))**2 )  

    f, Pxx_spec = signal.periodogram(xx, fs, scaling='density') 
#    Pxxt_spec= Pxx_spec[(f>fmin_Hz)&(f<fmax_Hz)]; ft=f[(f>fmin_Hz)&(f<fmax_Hz)]
    Pxxt_spec= Pxx_spec[(f>=fmin_Hz)&(f<=fmax_Hz)]; ft=f[(f>=fmin_Hz)&(f<=fmax_Hz)]  # correction 28/01/2022
    nbins=101; fbins=np.logspace(-2,3,nbins); Pbins=np.ones(nbins); dPbins=np.zeros(nbins)
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
            axbis.set_ylim(1e-1, 1e4); axbis.set_xscale('log'); axbis.set_xlim(1e-2, 1e3); axbis.set_yscale('log')    
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
    fc = pEq[0]*pEq[1]/(2*np.pi*kBT_pN_nm)   #  modification 26/01/2022    fc=pEq[0]/(2*np.pi*friction)
    print('friction / friction_bulk =', friction/friction0, ' Corner frequency fc=', fc)
    FitPxxt_spec=FitSpectrum(ft, pEq[0], pEq[1])

    if display:
        for axbis in wax: axbis.plot(ft, FitPxxt_spec, c='r', alpha=0.8)
        if SaveGraph: plt.savefig(path1+OutputFolder+figname+OutFormat, transparent=True)
        if CloseAfterSave: plt.close()
    if exportname!=None: exportspectrum(ft, FitPxxt_spec, fbins, Pbins, dPbins, exportname+label+'.csv')
    print('Spectrum'+label, ' k (pN/nm)=',"%.5f" %  (pEq[0]),' D (µm²/s)=',"%.3f" %  (pEq[1]*1.e-6), 'D/Dtheo=', pEq[1]*friction/kBT_pN_nm), 'fc (Hz)= ',"%.3f" %  fc
    return pEq[0], pEq[1], eEq[0], eEq[1], fc, friction, ft, Pxxt_spec, f, Pxx_spec

def survival(Tupjump_Tuphigh, TEhigh, countNoOpen, shift, refname, color, exportname=None):
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
    if DisplaySurvival:
        figname=refname+'Survival'
        plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
        ax.scatter(dur2s, dur2s_y, c=color, marker='o', alpha=0.5, label=refname)
        ax.set_xlabel("Duration of closed state $\Delta$t (s)");  ax.set_ylabel("Survival Fraction")
        ax.axis([0., TEhigh/1000., 0.01, 1.]); ax.set_yscale('log')
    maxtfit=0.95*TEhigh/1000-shift
    x2 = dur2s[dur2s<maxtfit]
    y2 = dur2s_y[dur2s<maxtfit]   #    y2log = dur2s_logy[dur2s<maxtfit]
    while True:
        try:
            initialguessEq=[20]
            pEq0, pcovEq0=curve_fit(FitExp0, x2, y2, p0=initialguessEq)
            print(pEq0[0], pcovEq0[0][0])
            break
        except RuntimeError:
            print("No convergence"); pEq0=[1.]; pcovEq0=[[0]]
            break
        except ValueError:
            print("Empty data"); pEq0=[1.]; pcovEq0=[[0]]
            break
    y2_fit0=FitExp0(x2, pEq0[0])
    while True:
        try:
            initialguessEq=[20, 0.]
            pEq, pcovEq=curve_fit(FitExp1, x2, y2, p0=initialguessEq)
            break
        except RuntimeError:
            print("No convergence"); pEq=[1.,0.]; pcovEq=[[0,0],[0,0]]
            break
        except ValueError:
            print("Empty data"); pEq=[1.,0.]; pcovEq=[[0,0],[0,0]]
            break
    y2_fit=FitExp1(x2, pEq[0], pEq[1])
    if DisplaySurvival:
   #     ax.plot(x2, y2_fit, c='r',alpha=0.5, label='kd 2param fit='+"%.6f" % (1/pEq[0]))   
        ax.plot(x2, y2_fit0, c='k',alpha=0.5, label='kd 1param fit='+"%.6f" % (1/pEq0[0]))
   #     ax.legend(fontsize=6)
        plt.tight_layout()
        if SaveGraph: plt.savefig(path1+OutputFolder+figname+OutFormat, transparent=False)
    if exportname!=None: exportsurvival(Tupjump_Tuphigh, exportname+'.csv')

    print('Survival 1param fit 0-', maxtfit, 's: offrate=', '{:.2E}'.format(1/pEq0[0]), '+/-', '{:.2E}'.format(np.sqrt(pcovEq0[0][0])/pEq0[0]**2) , 'tau (s)=', '{:.2E}'.format(pEq0[0]))
    print('Survival 2param fit 0-', maxtfit, 's: offrate=', '{:.2E}'.format(1/pEq[0]), 'tau (s)=', '{:.2E}'.format(pEq[0]), 'P0 (s)=', '{:.2E}'.format(pEq[1]))
    return(pEq0[0], np.sqrt(pcovEq0[0][0]), pEq[0], pEq[1], x2, y2, y2_fit0, y2_fit)

##========================= ========================= =========================
##=========================   # MAIN LOOP  # ==================================
listresultsTotal=[]
#listresultsTotalSpec=[]
nbead=len(b)
color=iter(cm.rainbow(np.linspace(0, 1, 4*nbead*len(RangeNZlavg))))

if nbead>1: figAllSpec, axAllSpec = plt.subplots(nbead, 4, num='figAllSpec', figsize=(6, 2*nbead), dpi=100)

for ibead in range(nbead):        ## Loop on test beads
    print('======================BEAD ', b[ibead] ,' ===========================')
    print('NZavg=',NZavg, 'dzjump=', dzjump, 'tol=',tol, 'NZrefavg=',NZrefavg)
    if loaddata and load[ibead]: print('Loading data', fname[ibead]); tdms_file = TdmsFile.open(path1+fname[ibead])  # alternative TdmsFile.read(path1+fname[ibead])
    if loaddata and not load[ibead]: print('Data already loaded', fname[ibead])
    tdms_groups = tdms_file.groups(); # if pprint: print(tdms_groups)
    tdms_track = tdms_groups[0]; tdms_temp = tdms_groups[1]; tdms_event = tdms_groups[2]
    dmax=1000000.; imax=(len(tdms_track)-1)//5; dX=1500; dY=1500    
    print('Nb of traj=', imax)
    refname=fname[ibead][:-5]+'_Bead'+str(b[ibead])+'_'; print(refname)
    if GraphXYall and ibead==0:
 #       AFSG.MakeGraphXYall(tdms_track, imax, 100, 100, fname[ibead], path1+OutputFolder, OutFormat, SaveGraph)
        AFSG.MakeGraphXYall(tdms_track, imax, 100, 100, fname[ibead][:-5], path1+OutputFolder, OutFormat, SaveGraph, pprint=pprint)

    (X0,Y0,Z0,MinLUT0,SSIM0)=track(tdms_track, b[ibead])        # test bead
    X0med=np.median(X0[~np.isnan(X0)]); Y0med=np.median(Y0[~np.isnan(Y0)])
    if load[ibead]:
        print("Calculating pool of reference bead(s)", bref )
        if len(bref)>1:
            (Xr0,Yr0,Zr0,MinLUTr0,SSIMr0)=track(tdms_track, bref[0])
            for i in range(1,len(bref)):
                (Xr0_,Yr0_,Zr0_,MinLUTr0_,SSIMr0_)=track(tdms_track, bref[i])
                Xr0=Xr0+Xr0_; Yr0=Yr0+Yr0_; Zr0=Zr0+Zr0_
            Xr0=Xr0/len(bref); Yr0=Yr0/len(bref); Zr0=Zr0/len(bref)
        else:
            (Xr0,Yr0,Zr0,MinLUTr0,SSIMr0)=track(tdms_track, bref[0])   # reference bead
    else: print("Pool of reference bead(s) already calculated") 
    (TE0,NE0,DE0) = events(tdms_event)
    T0, fs = time(tdms_track, GraphTdms)
    start[ibead]=start_min*60*fs; stop[ibead]=stop_min*60*fs
    range_anchorpoint=(int(range_anchorpoint_min[0]*60*fs),int(range_anchorpoint_min[1]*60*fs))
    range_spectrum=(int(range_spectrum_min[0]*60*fs),int(range_spectrum_min[1]*60*fs))
    
    # Clean outliers
    if CleanOutliers:
        mm=0.5; mmabs=2000; print('Clean outliers nan or +/-', mm)
        i0 = np.isnan(X0) + np.isnan(Y0) + np.isnan(Z0) 
        mX0 = np.median(X0[~i0]);  X0[i0] = mX0
        mY0 = np.median(Y0[~i0]);  Y0[i0] = mY0
        mZ0 = np.median(Z0[~i0]);  Z0[i0] = mZ0
        if pprint: print('median X0=', mX0, 'median Y0=', mY0, 'median Z0=', mZ0 )
    #    j0 = (X0>(1+mm)*mX0) + (X0<(1-mm)*mX0) + (Y0>(1+mm)*mY0) + (Y0<(1-mm)*mY0) + (Z0>(1+mm)*mZ0) + (Z0<(1-mm)*mZ0)
        j0 = (X0>mX0+mmabs) + (X0<mX0-mmabs) + (Y0>mY0+mmabs) + (Y0<mY0-mmabs) + (Z0>mZ0+mmabs) + (Z0<mZ0-mmabs)
        X0[j0] = mX0; Y0[j0] = mY0; Z0[j0] = mZ0; 
    #    m0 = np.median(X0[~i0]);  X0[i0] = m0*(1+.01*np.random.rand(len(X0[i0]))); X0[ (X0>(1+mm)*m0) | (X0<(1-mm)*m0) ] = m0
        m0 = np.median(MinLUT0[~i0]);  MinLUT0[i0] = m0 ; MinLUT0[ (MinLUT0>(1+mm)*m0) | (MinLUT0<(1-mm)*m0) ] = m0
        m0 = np.median(SSIM0[~i0]);  SSIM0[i0] = m0 ; SSIM0[ (SSIM0>(1+mm)*m0) | (SSIM0<(1-mm)*m0) ] = m0
    
    P0 = BuildP(T0, TE0, NE0, DE0)
    # three command lines added to remove the segment containing spikes at the end of a trace
    # Tr0 = np.copy(T0); Pr0 = np.copy(P0); MinLUTr0 = np.copy(MinLUT0)
    # (X0,Y0,Z0,T0,P0, MinLUT0) = shortentrack(np.copy(X0), np.copy(Y0), np.copy(Z0), np.copy(T0), np.copy(P0), np.copy(MinLUT0), 0, 1020000)
    # (Xr0,Yr0,Zr0,T0,P0, MinLUT0) = shortentrack(np.copy(Xr0), np.copy(Yr0), np.copy(Zr0), Tr0, Pr0, MinLUTr0, 0, 1020000)

    if GraphTdms:
        plt.figure('Power', figsize=(6,6), dpi=100); ax = plt.gca()
        ax.scatter(T0, P0, marker='.', c='b', alpha=0.5)
 
    # make sure that chosen time interval does not start or stop at high power       
    if not forcecalib:
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

    # interpolated traces for selfreference
    Z10 = Buildinterpol(Z0, T0, TE0, NE0, DE0)
    X10 = Buildinterpol(X0, T0, TE0, NE0, DE0, Outer=False)
    Y10 = Buildinterpol(Y0, T0, TE0, NE0, DE0, Outer=False)

    print('Tracks shortening from frame', start[ibead], 'to ',stop[ibead])
    (X1,Y1,Z1,T1,P1, MinLUT1)= shortentrack(X10, Y10, Z10, np.copy(T0), np.copy(P0), np.copy(MinLUT0), start[ibead], stop[ibead]) # shorten selfreference trace
    (X,Y,Zu,T,P, MinLUT)= shortentrack(X0, Y0, Z0, T0, P0, MinLUT0, start[ibead], stop[ibead])  # shorten test bead    
    (Xr,Yr,Zr,Tr,Pr, MinLUTr)= shortentrack(Xr0, Yr0, Zr0, T0, P0, MinLUTr0, start[ibead], stop[ibead])   # shorten reference bead
     
    mask=(TE0>=T0[start[ibead]])&(TE0<T0[stop[ibead]])
    TE=TE0[mask]; NE=NE0[mask]; DE=DE0[mask]
    
    height=2
    if forcecalib:
        height=1
        (TEhigh, TElow) = (0,0)
    else:
        TEhigh, TElow = period(TE, NE, DE, GraphTdms, refname, height)        
        
    if not forcecalib:
        Phigh, Plow = MakeHistoP(P, NE, DE, GraphTdms)
    elif forcecalib:
        Phigh=0.2; Plow=0
        print('force calibration mode, enforced values:  Phigh=', Phigh, ' Plow=', Plow)
  
    Zr, Zravg = Rollavg(Zr, NZrefavg, mode)        # rolling average of shorten reference bead
    Zr0, Zr0avg = Rollavg(Zr0, NZrefavg, mode)       # rolling average of reference bead
    Xr, Xravg = Rollavg(Xr, NZrefavg, mode); Yr, Yravg = Rollavg(Yr, NZrefavg, mode) 
    Xr0avg, Xr0 = Rollavg(Xr0, NZrefavg, mode); Yr0avg,Yr0=Rollavg(Yr0, NZrefavg, mode) 

    # selfreference calculation based on low/high force states
    indlow=(P>=Plow*(1-tol))*(P<=Plow*(1+tol))
    indhigh=(P>=Phigh*(1-tol))*(P<=Phigh*(1+tol))
    if forcecalib: indhigh=(P>=Phigh*(1-tol))
    NZrefavgLarge = 10*4*(1+9*int(forcecalib==True))*NZrefavg  # prefactor 10*
    Z1[indlow] = Zu[indlow]; Z1, Z1avg = Rollavg(Z1, NZrefavgLarge, mode)    # keep low force state for selfreference Z
    X1[indhigh] = X[indhigh]; X1, X1avg = Rollavg(X1, NZrefavgLarge, mode)   # keep high force state for selfreference X
    Y1[indhigh] = Y[indhigh]; Y1, Y1avg = Rollavg(Y1, NZrefavgLarge, mode)   # keep high force state for selfreference Y
    # full length selfreference
    ind0low=(P0>=Plow*(1-tol))*(P0<=Plow*(1+tol))
    ind0high=(P0>=Phigh*(1-tol))*(P0<=Phigh*(1+tol))
    if forcecalib: ind0high=(P0>=Phigh*(1-tol))
#    Z10[start[ibead]:stop[ibead]]=Z1; 
    Z10[ind0low] = Z0[ind0low]; Z10, Z10avg = Rollavg(Z10, NZrefavgLarge, mode)
#    X10[start[ibead]:stop[ibead]]=X1
    X10[ind0high] = X0[ind0high] ; X10, X10avg = Rollavg(X10, NZrefavgLarge, mode)
#    Y10[start[ibead]:stop[ibead]]=Y1;
    Y10[ind0high] = Y0[ind0high] ; Y10, Y10avg = Rollavg(Y10, NZrefavgLarge, mode)
   
    ind0force=(P==0)
    indeforce0=(P0==0)
    mask_anchor0=(T0>=T0[range_anchorpoint[0]])&(T0<T0[range_anchorpoint[1]])
    mask_spectrum0=(T0>=T0[range_spectrum[0]])&(T0<T0[range_spectrum[1]])
    
    if RefBeadZcorrection:
        if selfreference==True:
            print('Baseline Zcorrection by low_force self reference' )          
            Zc = Zu-Z1avg+Z1avg[0]; Zec = Z0-Z10avg+Z10avg[0]   # corr LL 27012021   # Z0c=Zr0; 
            Zref=np.array([]); Zref0=np.array([])
        else:
            print('Baseline Zcorrection by reference bead(s): Zravg=', "%.3f" % Zravg[0] )#int(bref[ibead]) )
            Zc = Zu-Zravg+Zravg[0]; Zec = Z0-Zr0avg+Zr0avg[0]   # ; Z0c=Zr0-Zr0avg
            Zref=Zravg-Zravg[0]+Z1avg[0]; Zref0=Zr0avg-Zr0avg[0]+Z10avg[0]
        AFSG.GraphSelfReference(T, Zu, Z1, Z1avg, Zref, 'Z', refname+'_self reference_', path1+OutputFolder,  OutFormat, SaveGraph)
        AFSG.GraphSelfReference(T0, Z0, Z10, Z10avg, Zref0, 'Z0', refname+'_self reference_', path1+OutputFolder,  OutFormat, SaveGraph)
                
        if AllenComputation:
            Allen(Z0, tauAllen, nAllen, 'Z0', stepsAllen)
            Allen(X0, tauAllen, nAllen, 'X0', stepsAllen)
            Allen(Y0, tauAllen, nAllen, 'Y0', stepsAllen)
        
    else: print('No Baseline Zcorrection'); Zc = Zu; Zec = Z0      # ; Z0c =Zr0

    if RefBeadXYcorrection:
        if selfreference:
            print('Baseline XYcorrection by high force self reference: Xravg=', "%.3f" % X1avg[0],' Yravg=', "%.3f" % Y1avg[0] )
            Xbl = X1avg-X1avg[0]; Ybl = Y1avg-Y1avg[0]; X0bl = X10avg-X10avg[0]; Y0bl = Y10avg-Y10avg[0]
            Xref=np.array([]); Xref0=np.array([]); Yref=np.array([]); Yref0=np.array([])
        else:
            print('Baseline XYcorrection by reference bead(s): Xravg=', "%.3f" % Xravg[0],' Yravg=', "%.3f" % Yravg[0] )
            Xbl = Xravg-Xravg[0]; Ybl = Yravg-Yravg[0]; X0bl = Xr0avg-Xr0avg[0]; Y0bl = Yr0avg-Yr0avg[0]
            Xref=Xravg-Xravg[0]+X1avg[0]; Xref0=Xr0avg-Xr0avg[0]+X10avg[0]
            Yref=Yravg-Yravg[0]+Y1avg[0]; Yref0=Yr0avg-Yr0avg[0]+Y10avg[0]
        AFSG.GraphSelfReference(T, X, X1, X1avg, Xref, 'X', refname+'_self reference_', path1+OutputFolder,  OutFormat, SaveGraph)
        AFSG.GraphSelfReference(T, Y, Y1, Y1avg, Yref, 'Y', refname+'_self reference_', path1+OutputFolder,  OutFormat, SaveGraph)
        AFSG.GraphSelfReference(T0, X0, X10, X10avg, Xref0, 'X0', refname+'_self reference_', path1+OutputFolder,  OutFormat, SaveGraph)
        AFSG.GraphSelfReference(T0, Y0, Y10, Y10avg, Yref0, 'Y0', refname+'_self reference_', path1+OutputFolder,  OutFormat, SaveGraph)
    
    else: print('No Baseline XYcorrection'); Xbl = 0; Ybl = 0; X0bl = 0; Y0bl = 0
    
    Xc = X-Xbl; Yc = Y-Ybl
    Xec = X0-X0bl; Yec = Y0-Y0bl
    if GraphXYall and ibead==0:
        AFSG.MakeGraphXYall(tdms_track, imax, 100, 100, fname[ibead][:-5], path1+OutputFolder,
                            OutFormat, SaveGraph, True, X0bl, Y0bl, pprint=pprint);  # Xr0avg, Yr0avg)

    indnonan=(~np.isnan(Xc))*(~np.isnan(Yc))
    Xc=Xc[indnonan]; Yc=Yc[indnonan]; Zc=Zc[indnonan]
    P=P[indnonan]; T=T[indnonan]; MinLUT=MinLUT[indnonan]

    indlow=(P>=Plow*(1-tol))*(P<=Plow*(1+tol)); #indlow=(P<Phigh*0.9)+(P>Pmax)
    indhigh=(P>=Phigh*(1-tol))*(P<=Phigh*(1+tol))
    if forcecalib: indhigh=(P>=Phigh*(1-tol))
    ind0force=(P==0)
    indeforce0=(P0==0)
    mask_anchor=(T>=T[range_anchorpoint[0]-start[ibead]])&(T<T[range_anchorpoint[1]-start[ibead]])
    mask_spectrum=(T>=T[range_spectrum[0]]-start[ibead])&(T<T[range_spectrum[1]-start[ibead]])

    listresults=[]
    listresultsSpec=[]
    if export:    # exports raw data for third party analysis
        exporttrack(X,Y,Zu,T,P, path1+'Export_'+refname+'.csv')
        exporttrack(Xr,Yr,Zr,Tr,P, path1+'Export_'+refname+'_r'+str(bref)+'.csv')        
    
    p1Lo=p1Lo0; p1Lc=p1Lc0; p1Ll=p1Ll0; p1Zl=p1Zl0  # initialization of gaussian center for fit guess
    for NZlavg in RangeNZlavg:      # test loop for various choices of NZlavg
        if len(RangeNZlavg)>1: print('=============== NZlavg=', NZlavg, '===============')
        Zcavg = uniform_filter1d(Zc, size=NZlavg, mode=mode)
        Zecavg = uniform_filter1d(Zec, size=NZlavg, mode=mode)
         
        x = np.random.randn(len(P))/10
        Zch=Zc[indhigh]; Xch=Xc[indhigh]; Ych=Yc[indhigh] # ; dZh=dZ[indhigh]; Ph=P[indhigh]+x[indhigh]
        Zcl=Zcavg[indlow]; Xcl=Xc[indlow]; Ycl=Yc[indlow] # ; dZl=dZ[indlow]; Pl=P[indlow]+x[indlow]
        Zc0=Zcavg[ind0force]; Xc0=Xc[ind0force]; Yc0=Yc[ind0force] #; dZ0=dZ[ind0force]
        Zec0=Zecavg[indeforce0]; Xec0=Xec[indeforce0]; Yec0=Yec[indeforce0] # ; dZe0=dZe[indeforce0]
        Xec0_mean=np.median(Xec0); Yec0_mean=np.median(Yec0);
        Zec0custom=Zecavg[mask_anchor0]; Xec0custom=Xec[mask_anchor0]; Yec0custom=Yec[mask_anchor0] #; dZe0custom=dZe[mask_anchor]

        figname=refname+'Various intervals'
        plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
        ax.plot(T0, Zecavg, c='k', alpha=0.1, label='Zcavg')
        ax.scatter(T[indhigh], np.median(Zch)*np.ones(len(T))[indhigh], marker='.', s=3, alpha=0.5, label='high power')
        ax.scatter(T[indlow], np.median(Zcl)*np.ones(len(T))[indlow]*10, marker='.', s=5, alpha=0.5, label='low power')
        ax.scatter(T[ind0force], np.median(Zc0)*np.ones(len(T))[ind0force], marker='.', s=3, alpha=0.5, label='0_force')
        ax.scatter(T0[indeforce0], np.median(Zec0)*np.ones(len(T0))[indeforce0], marker='.', s=3, alpha=0.5, label='0_force full range')
        ax.scatter(T0[mask_anchor0], np.median(Zec0custom)*np.ones(len(T0))[mask_anchor0], marker='.', s=8, alpha=0.5, label='mask_anchor0')
        if Select_range_spectrum: ax.scatter(T0[mask_spectrum0], np.zeros(len(T0))[mask_spectrum0], marker='.', s=8, alpha=0.5, label='mask_spectrum0')
        ax.legend(fontsize=10, markerscale=3)
        if SaveGraph: plt.savefig(path1+OutputFolder+figname+OutFormat, transparent=True)

#        indstart=(T0<start[ibead]); Z0start=Z0c[indstart]
        Nl=len(Zc[indlow]); Nh=len(Zc[indhigh]); N0=len(Zc[ind0force]); Ntot=len(Zc)
        print('States Populations: low=', Nl, ' high=', Nh, ' 0_force=', N0, ' total=', Ntot)        
        print('States Fractions: low=',"%.3f" %  (Nl/Ntot), ' high=',"%.3f" %  (Nh/Ntot), ' 0_force=',"%.3f" %  (N0/Ntot), 'total=',"%.3f" %  ((Nl+Nh+N0)/Ntot))        
        for AnchorPointState in AnchorPointStateList:      # ['low_force','0_force','custom']     
            print('=============== Anchor Point Reference:', AnchorPointState, end=' ')
            if AnchorPointState=='None':
                print('No anchoring point'); Zs=Zc; Xs=Xc; Ys=Yc
            else:
                if AnchorPointState=='low_force':
                    Zc_=Zcl; Xc_=Xcl; Yc_=Ycl; print("") 
                elif AnchorPointState=='0_force':
                    NZlavg+=10; print("") 
                    Zc_=Zc0; Xc_=Xc0; Yc_=Yc0   
         #           Zc_=Zec0; Xc_=Xec0; Yc_=Yec0
                elif AnchorPointState=='custom':
                    Zc_=Zec0custom; Xc_=Xec0custom; Yc_=Yec0custom
                    print(" interval (min)=", range_anchorpoint_min[0], range_anchorpoint_min[1])
                elif AnchorPointState=='0_force_custom':
                    mask_anchor0_0_force = mask_anchor0 * indeforce0
                    Zc_=Zecavg[mask_anchor0_0_force]; Xc_=Xec[mask_anchor0_0_force]; Yc_=Yec[mask_anchor0_0_force]
                    print(" interval (min)=", range_anchorpoint_min[0], range_anchorpoint_min[1])
                Zc_sort=np.sort(Zc_)       #     cleanedZclsort = [x for x in Zclsort if str(x) != 'nan']
                if len(Zc_sort)>0: Za=np.median(Zc_sort[np.arange(sample)])  # Za=np.amin(Zcl) #-BeadRad
                else:
                    print('Zc_ null wave:  Bad Anchor Point Reference Choice')
                    break
                print('Z reference Anchoring point at ', AnchorPointState, '. Length calculation', 'min=', "%.3f" % np.min(Zc_), 'min_',sample, "  %.3f" %Za)
                Xc_ = Xc_[~np.isnan(Xc_)]; Yc_ = Yc_[~np.isnan(Yc_)]
                AnchorCov=np.cov(Xc_,Yc_); print('Anchor Point Covariance=', AnchorCov)
                w, v = np.linalg.eig(AnchorCov); SymmetryFactor = np.sqrt(np.amax(w)/np.amin(w))
                print('Eigenvalues=', w, 'SymmetryFactor=', SymmetryFactor)
                Xa=np.median(Xc_); Ya=np.median(Yc_)     
                CA=np.sqrt((Xc-Xa)**2+(Yc-Ya)**2+(Zc-Za)**2)
                Zs=Zc-Za; Xs=Xc-Xa; Ys=Yc-Ya
                Zs00=Zec-Za; Xs00=Xec-Xa; Ys00=Yec-Ya
                Xol=Xc[indlow]; Yol=Yc[indlow]      # Jim  XY fluctuation at low power without anchor point corrections for selected beads
            Xsl=Xs[indlow]; Ysl=Ys[indlow]; Zsl=Zs[indlow]; Tl=T[indlow]; MinLUTl=MinLUT[indlow]
            Xsh=Xs[indhigh]; Ysh=Ys[indhigh] ; Zsh=Zs[indhigh]; Th=T[indhigh]; MinLUTh=MinLUT[indhigh]
            Xs0=Xs[ind0force]; Ys0=Ys[ind0force]; Zs0=Zs[ind0force]     
            Xsl_md=np.median(Xsl); Ysl_md=np.median(Ysl); Zsl_md=np.median(Zsl)  # Jim When "0 force" is selected to set the anchoring point, this line gives the anchoring point difference between "0 force" and " low power".
            Length=np.sqrt(Xs**2+Ys**2+(Zs+BeadRad)**2)-BeadRad; Lengthl=Length[indlow]; Lengthh=Length[indhigh]; Length0=Length[ind0force]
            PullAngle=(180/np.pi)*np.arcsin( np.sqrt( (np.median(Xsh)-np.median(Xsl))**2 + (np.median(Ysh)-np.median(Ysl))**2  )/np.median(Lengthh))
            PullAngle2=(180/np.pi)*np.arcsin( np.sqrt( (np.median(Xsh)-np.median(Xs0))**2 + (np.median(Ysh)-np.median(Ys0))**2  )/np.median(Lengthh))            
            print('High length=', "%.3f" % np.median(Lengthh), 'Low length=', "%.3f" % np.median(Lengthl), ' Pulling Angle=', "%.3f" % PullAngle, ' Pulling Angle2=', "%.3f" % PullAngle2)
            PullPhi = np.arctan2(np.median(Ysh)-np.median(Ys0), np.median(Xsh)-np.median(Xs0)) * 180 / np.pi
            print('Pulling Angle (0° along X) Phi=', PullPhi)
            
            if Display3D: AFSG.Plot3D(Xs,Ys,Zs, 30, refname+'XsYsZs', path1+OutputFolder, OutFormat, SaveGraph)
            if HistoZ: MakeHistoZ(Lengthh, Lengthl, 'Lengthh', 'Lengthl', -500, MaxZHistoPlot, refname+'LengthhLengthl')
            if GraphXYZ:
                contourdataXYslXYsh = AFSG.MakeGraphXY(Xsl, Ysl, Xsh, Ysh, -dX, dX, -dY, dY, 'XYl', 'XYh', refname+'XYslXYsh', SaveGraph, path1+OutputFolder, OutFormat, levels=10)
                contourdataXY0XYsh = AFSG.MakeGraphXY(Xec0, Yec0, [], [], Xa-dX, Xa+dX, Ya-dY, Ya+dY, 'XY0', '__', refname+'XY0__', SaveGraph, path1+OutputFolder, OutFormat, levels=10)
                    # Jim Plot XY fluctuation of selected beads during 0 force without anchoring correction 
                contourdataXYoriginlXYsh = AFSG.MakeGraphXY(Xol, Yol, [], [], Xa-dX, Xa+dX, Ya-dY, Ya+dY, 'XYlow', '__', refname+'XYoriginl__', SaveGraph, path1+OutputFolder, OutFormat, levels=20)
                    # Jim Plot XY fluctuation of selected beads during low power without anchoring correction 
           #     AFSG.MakeGraphXY(Xs0, Ys0, Xsh, Ysh, -dX, dX, -dY, dY, 'XY0', 'XYh', refname+'XYs0XYsh', SaveGraph, path1+OutputFolder, OutFormat, levels=20)
                    # Jim Plot to test XY fluctuation of selected beads at 0 force            
                AFSG.MakeGraphXY(Tl, Zsl, Th, Zsh, 0, T[-1], -100, MaxZHistoPlot, 'Xs', 'Ys', refname+'TZslZsh', SaveGraph, path1+OutputFolder, OutFormat, contour=False)
                AFSG.SingleTrace(refname, tdms_track, -2, sampleGraph, Xs, Ys, Zs, T, MinLUT, SaveGraph, path1+OutputFolder, OutFormat )
                AFSG.MakeGraphXY(T, Length, T, Zs, 0, T[-1], 0, MaxZHistoPlot, 'Length', 'Zs', refname+'Length', SaveGraph, path1+OutputFolder, OutFormat, contour=False)
        
            # Jim the calculation of the length and pulling angle from the center of the contour    
     #       contourdataXY0XYsh = contourdataXYslXYsh
            (X0ca,Y0ca) = contourdataXY0XYsh[0][2]
            if manual_anchor_point:
                (Xca,Yca) = coordinates_anchor_point
                print("Manual Anchor point: ", (Xca,Yca))
            elif AnchorPointState=='custom' or AnchorPointState=='0_force_custom':
                contourdataXYcustom = AFSG.MakeGraphXY(Xec0custom, Yec0custom, [], [], -dX, dX, -dY, dY, 'XYcustom', '__', refname+'XYcustom__', SaveGraph, path1+OutputFolder, OutFormat, levels=20)
                (Xca,Yca) = contourdataXYcustom[0][2]
                print("Custom interval contour Anchor point: ", (Xca,Yca))
            else:
                (Xca,Yca) = contourdataXYoriginlXYsh[0][2] 
                print("Contour Anchor point: ", (Xca,Yca))
            Xsc=Xc-Xca; Ysc=Yc-Yca
            Xscl=Xsc[indlow]; Yscl=Ysc[indlow]
            Xsch=Xsc[indhigh]; Ysch=Ysc[indhigh] 
            Xsc0=Xsc[ind0force]; Ysc0=Ysc[ind0force]    
            Xscl_md=np.median(Xscl); Yscl_md=np.median(Yscl);   
            LengthC=np.sqrt(Xsc**2+Ysc**2+(Zs+BeadRad)**2)-BeadRad; LengthCl=LengthC[indlow]; LengthCh=LengthC[indhigh]; LengthC0=LengthC[ind0force]
            PullAngleC=(180/np.pi)*np.arcsin( np.sqrt( (np.median(Xsch)-np.median(Xscl))**2 + (np.median(Ysch)-np.median(Yscl))**2  )/np.median(LengthCh))
            PullAngleC2=(180/np.pi)*np.arcsin( np.sqrt( (np.median(Xsch)-np.median(Xsc0))**2 + (np.median(Ysch)-np.median(Ysc0))**2  )/np.median(LengthCh))
            PullPhiC = np.arctan2(np.median(Ysch)-np.median(Ysc0), np.median(Xsch)-np.median(Xsc0)) * 180 / np.pi
            
            if HistoZ: MakeHistoZ(LengthCh, LengthCl, 'LengthCh', 'LengthCl', -500, MaxZHistoPlot, refname+'LengthChLengthCl')
            if GraphXYZ:
               AFSG.SingleTrace(refname, tdms_track, -3, sampleGraph, Xsc, Ysc, Zs, T, MinLUT, SaveGraph, path1+OutputFolder, OutFormat )
               AFSG.MakeGraphXY(Xsc0, Ysc0, Xsh, Ysh, Xca-dX, Xca+dX, Yca-dY, Yca+dY, 'XY0', 'XYh', refname+'XYsc0XYsh', SaveGraph, path1+OutputFolder, OutFormat, levels=20)
                # Jim Plot to test XY fluctuation of selected beads at 0 force
               AFSG.MakeGraphXY(Xec0, Yec0, [], [], Xca-dX, Xca+dX, Yca-dY, Yca+dY, 'XY0', '__', refname+'XY0_center contour', SaveGraph, path1+OutputFolder, OutFormat, levels=10) 
               AFSG.MakeGraphXY(Xol, Yol, [], [], Xca-dX, Xca+dX, Yca-dY, Yca+dY, 'XYlow', '__', refname+'XYoriginl_center contour', SaveGraph, path1+OutputFolder, OutFormat, levels=10)
               AFSG.MakeGraphXY(Xol, Yol, Xec0, Yec0, Xca-dX, Xca+dX, Yca-dY, Yca+dY, 'XYlow', 'XY0', refname+'XYoriginlXY0', SaveGraph, path1+OutputFolder, OutFormat,contour=False)  
               AFSG.MakeGraphXY(Xol, Yol, Xec0, Yec0, Xca-dX, Xca+dX, Yca-dY, Yca+dY, 'XYlow', 'XY0', refname+'XYoriginlXY0_2', SaveGraph, path1+OutputFolder, OutFormat, disperse=False, contour=False)
                    # Jim Plot XY fluctuation of selected beads to compare between low power and 0 power without anchoring correction

      #      if AllenComputation:
     #           Allen(Zcl, T[indlow], nAllen, 'Zcl', [1,2,4,8])

            print('=============== Opening analysis ===============')       
     #       Zsavg=np.convolve(Zs, np.ones((NZavg,))/NZavg, mode='same')
            Zsavg = uniform_filter1d(Zs, size=NZavg, mode=mode)

            if not forcecalib:

                DeltaCOFit=0
                for i in [0,1]:         #  [0,1] first pass dzjump, 2nd pass DeltaCOFit
          #          dzjump_tmp = dzjump*(i==0) + DeltaCOFit*weightsecondpass_step*(i==1)
          #          dzjump_tmp = dzjump*(i==0) + DeltaCOFit*(i==1)
                    dzjump_tmp = dzjump*(i==0) + ( (1-weightsecondpass_step)*dzjump + DeltaCOFit*weightsecondpass_step )*(i==1)
                    print('Round ',i,' /0  dzjump=',dzjump_tmp)
                    Pdown=Plow
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
                    Xscclosed=Xsc[indclosed]; Yscclosed=Ysc[indclosed]; LengthCclosed=LengthC[indclosed]  #Jim the calculation with the center of the contour
                    Zsopen=Zs[indopen]; Topen=T[indopen]; Xsopen=Xs[indopen]; Ysopen=Ys[indopen]; Lengthopen=Length[indopen]
                    Xscopen=Xsc[indopen]; Yscopen=Ysc[indopen]; LengthCopen=LengthC[indopen]              #Jim the calculation with the center of the contour
                    Nc=len(Zsclosed); No=len(Zsopen)
                    print('States Population: close=', Nc, ' open=', No, ' high=', Nh)
                    Nh = max(1, Nh)
                    print('States Fractions: close=',"%.3f" %  (Nc/Nh), ' open=', "%.3f" %  (No/Nh), ' (Nc+No)/Nh=', "%.3f" %  ((Nc+No)/Nh))
        
                    print('=============== Gaussian fits ===============')
                    if i==0:
                        p0Lo, p1Lo, p2Lo, pLo = FitHistoGaussian(Lengthopen, 'Lengthopen nm ', DisplayGaussianFits , 100 )
                        p0Lc, p1Lc, p2Lc, pLc = FitHistoGaussian(Lengthclosed, 'Lengthclosed nm ', DisplayGaussianFits, 100 )
                        p0Ll, p1Ll, p2Ll, pLl = FitHistoGaussian(Lengthl, 'Lengthlow nm ', False, 300 )
                        
                        # Jim Gaussian analysis of the length and pulling angle from the center of the contour
                        Cp0Lo, Cp1Lo, Cp2Lo, CpLo = FitHistoGaussian(LengthCopen, 'LengthCopen nm ', DisplayGaussianFits , 100 )
                        Cp0Lc, Cp1Lc, Cp2Lc, CpLc = FitHistoGaussian(LengthCclosed, 'LengthCclosed nm ', DisplayGaussianFits, 100 )
                        Cp0Ll, Cp1Ll, Cp2Ll, CpLl = FitHistoGaussian(LengthCl, 'LengthClow nm ', False, 300 )
                        
                        
                        if len(Length0)>0: p0L0, p1L0, p2L0, pL0 = FitHistoGaussian(Length0, 'Length0 nm ', False, 300 )
                        p0Zl, p1Zl, p2Zl, pZl = FitHistoGaussian(Zsl, 'Zsl nm ', False, 300 )
                        p0Xo, p1Xo, p2Xo, pXo = FitHistoGaussian(Xsopen, 'Xsopen nm ', False , 300 )
                        p0Yo, p1Yo, p2Yo, pYo = FitHistoGaussian(Ysopen, 'Ysopen nm ', False , 300 )
                        
                        # Jim Gaussian analysis of the length and pulling angle from the center of the contour
                        Cp0Xo, Cp1Xo, Cp2Xo, CpXo = FitHistoGaussian(Xscopen, 'XsCopen nm ', False , 300 )
                        Cp0Yo, Cp1Yo, Cp2Yo, CpYo = FitHistoGaussian(Yscopen, 'YsCopen nm ', False , 300 )
                        
                    p0Zo, p1Zo, p2Zo, pZo = FitHistoGaussian(Zsopen, 'Zsopen nm ', False , 150 )
                    p0Zc, p1Zc, p2Zc, pZc = FitHistoGaussian(Zsclosed, 'Zsclosed nm ', False , 150 )
                    DeltaCOFit=p1Zo-p1Zc; MiddleOpenClose=(p1Zo+p1Zc)/2.
                    print('DeltaCOFit=', DeltaCOFit, 'MiddleOpenClose=', MiddleOpenClose)
                    
                if DisplayCloseOpen: AFSG.MakeDisplayCloseOpen(refname, Tclosed, Zsclosed, Topen, Zsopen, Tl, Zsl, T, Zsavg, P, Tup, SaveGraph, CloseAfterSave, path1+OutputFolder, OutFormat, top = True)
            
                MakeHistoZ(Zsclosed, Zsopen, 'Zsclosed', 'Zsopen', -500, MaxZHistoPlot, refname+'ZsclosedZsopen')
            elif forcecalib: print('Force calibration: no step detection')

            print('=============== Power Spectra fits ===============')
            if Select_range_spectrum or forcecalib:
                Zsopen_=Zs00[mask_spectrum0]; Xsopen_=Xs00[mask_spectrum0];  Ysopen_=Ys00[mask_spectrum0]
                Zsclosed_=Zs00[mask_spectrum0]; Xsclosed_=Xs00[mask_spectrum0];  Ysclosed_=Ys00[mask_spectrum0]
                print('Select_range_spectrum', range_spectrum[0], range_spectrum[1])
                if forcecalib:
                    print("Force calibration: measuring Zshigh and Lengthhigh")
                    indcalib = indhigh*mask_spectrum
                    p0Zo, p1Zo, p2Zo, pZo = FitHistoGaussian(Zs[indcalib], 'Zshigh nm ', False , 150 )
                    p0Lo, p1Lo, p2Lo, pLo = FitHistoGaussian(Length[indcalib], 'Lengthhigh nm ', DisplayGaussianFits , 100 )
                    Cp0Lo, Cp1Lo, Cp2Lo, CpLo = FitHistoGaussian(LengthC[indcalib], 'LengthCopen nm ', DisplayGaussianFits , 100 )
                    p0Xo, p1Xo, p2Xo, pXo = FitHistoGaussian(Xs[indcalib], 'Xshigh nm ', False , 300 )
                    p0Yo, p1Yo, p2Yo, pYo = FitHistoGaussian(Ys[indcalib], 'Yshigh nm ', False , 300 )
                    Cp0Xo, Cp1Xo, Cp2Xo, CpXo = FitHistoGaussian(Xsc[indcalib], 'XsChigh nm ', False , 300 )
                    Cp0Yo, Cp1Yo, Cp2Yo, CpYo = FitHistoGaussian(Ysc[indcalib], 'YsChigh nm ', False , 300 )
                    Cp0Lc, Cp1Lc, Cp2Lc, CpLc = Cp0Lo, Cp1Lo, Cp2Lo, CpLo
                    p0Lc, p1Lc, p2Lc, pLc = p0Lo, p1Lo, p2Lo, pLo
                    
            else:
                Zsopen_=Zs[indopen]; Xsopen_=Xsopen;  Ysopen_=Ysopen
                Zsclosed_=Zs[indclosed]; Xsclosed_=Xsclosed;  Ysclosed_=Ysclosed
                
            if nbead>1: axtable0=axAllSpec[ibead,0]; axtable1=axAllSpec[ibead,1]; axtable2=axAllSpec[ibead,2]; axtable3=axAllSpec[ibead,3]
            else: axtable0=None; axtable1=None; axtable2=None; axtable3=None
            exportspectrumname = path1+OutputFolder+refname+'ExportSpectrum_'
            kx, Dx, dkx, dDx, fc_x, frictionXY, ftx, Pxxt_specx, fx, Pxx_specx = Spectrum(Xsopen_, 'XY', p1Zo, fs,'XsopenRef'+AnchorPointState, DisplaySpectrum, fmin_Hz, fmax_Hz, exportname=exportspectrumname, axtable=axtable0) 
            ky, Dy, dky, dDy, fc_y, frictionXY, fty, Pxxt_specy, fy, Pxx_specy = Spectrum(Ysopen_, 'XY', p1Zo, fs,'YsopenRef'+AnchorPointState, DisplaySpectrum, fmin_Hz, fmax_Hz, exportname=exportspectrumname, axtable=axtable1)
            if forcecalib:
                kx2, Dx2, dkx2, dDx2, fc_x2, frictionXY2, ftx2, Pxxt_specx2, fx2, Pxx_specx2 = kx, Dx, dkx, dDx, fc_x, frictionXY, ftx, Pxxt_specx, fx, Pxx_specx
                ky2, Dy2, dky2, dDy2, fc_y2, frictionXY2, fty2, Pxxt_specy2, fy2, Pxx_specy2 = ky, Dy, dky, dDy, fc_y, frictionXY, fty, Pxxt_specy, fy, Pxx_specy
            else:
                kx2, Dx2, dkx2, dDx2, fc_x2, frictionXY2, ftx2, Pxxt_specx2, fx2, Pxx_specx2 = Spectrum(Xsclosed_, 'XY', p1Zo, fs,'XsclosedRef'+AnchorPointState, DisplaySpectrum, fmin_Hz, fmax_Hz, exportname=exportspectrumname, axtable=axtable0) 
                ky2, Dy2, dky2, dDy2, fc_y2, frictionXY2, fty2, Pxxt_specy2, fy2, Pxx_specy2 = Spectrum(Ysclosed_, 'XY', p1Zo, fs,'YsclosedRef'+AnchorPointState, DisplaySpectrum, fmin_Hz, fmax_Hz, exportname=exportspectrumname, axtable=axtable1)
            kz, Dz, dkz, dDz, fc_z, frictionZ, ftz, Pxxt_specz, fz, Pxx_specz = Spectrum(Zsopen_, 'Z', p1Zo, fs,'ZsopenRef'+AnchorPointState, DisplaySpectrum, fmin_Hz, fmax_Hz, exportname=exportspectrumname, axtable=axtable2)
            kr, Dr, dkr, dDr, fc_r, frictionR, ftr, Pxxt_specr, fz, Pxx_specr = Spectrum(np.sqrt(Xsopen_**2+Ysopen_**2), 'R', p1Zo, fs,'RsopenRef'+AnchorPointState, DisplaySpectrum, fmin_Hz, fmax_Hz, exportname=exportspectrumname, axtable=axtable3) 
            Fx=kx*(BeadRad+p1Lo); dFx=np.abs(dkx*(BeadRad+p1Lo))+np.abs(kx*p2Lo)
            Fy=ky*(BeadRad+p1Lo); dFy=np.abs(dky*(BeadRad+p1Lo))+np.abs(ky*p2Lo)

            Fx2=kx2*(BeadRad+p1Lc); dFx2=np.abs(dkx2*(BeadRad+p1Lc))+np.abs(kx2*p2Lc)
            Fy2=ky2*(BeadRad+p1Lc); dFy2=np.abs(dky2*(BeadRad+p1Lc))+np.abs(ky2*p2Lc)
            Fz=kz*(BeadRad+p1Lo); dFz=np.abs(dkz*(BeadRad+p1Lo))+np.abs(kz*p2Lo)
            Fr=kr*(BeadRad+p1Lo); dFr=np.abs(dkr*(BeadRad+p1Lo))+np.abs(kr*p2Lo)
            Fx3 = (p1Lc-p1Lo)*(kx*kx2)/(kx-kx2)
            Fy3 = (p1Lc-p1Lo)*(ky*ky2)/(ky-ky2)
            
            FxC=kx*(BeadRad+Cp1Lo); dFxC=np.abs(dkx*(BeadRad+Cp1Lo))+np.abs(kx*Cp2Lo)
            FyC=ky*(BeadRad+Cp1Lo); dFyC=np.abs(dky*(BeadRad+Cp1Lo))+np.abs(ky*Cp2Lo)
            Fx2C=kx2*(BeadRad+Cp1Lc); dFx2C=np.abs(dkx2*(BeadRad+Cp1Lc))+np.abs(kx2*Cp2Lc)
            Fy2C=ky2*(BeadRad+Cp1Lc); dFy2C=np.abs(dky2*(BeadRad+Cp1Lc))+np.abs(ky2*Cp2Lc)
            FzC=kz*(BeadRad+Cp1Lo); dFzC=np.abs(dkz*(BeadRad+Cp1Lo))+np.abs(kz*Cp2Lo)
            FrC=kr*(BeadRad+Cp1Lo); dFrC=np.abs(dkr*(BeadRad+Cp1Lo))+np.abs(kr*Cp2Lo)
            Fx3C = (Cp1Lc-Cp1Lo)*(kx*kx2)/(kx-kx2)
            Fy3C = (Cp1Lc-Cp1Lo)*(ky*ky2)/(ky-ky2)
            
            if nbead>1: 
                axAllSpec[ibead,0].set_title(str(b[ibead])+' Fx='+"%.3f" % Fx, fontsize=6)
                axAllSpec[ibead,1].set_title(str(b[ibead])+' Fy='+"%.3f" % Fy, fontsize=6)
                axAllSpec[ibead,2].set_title(str(b[ibead])+' Fz='+"%.3f" % Fz, fontsize=6)
                axAllSpec[ibead,3].set_title(str(b[ibead])+' Fr='+"%.3f" % Fr, fontsize=6)
            fcx = kx*Dx/(kBT_pN_nm*2*np.pi)  # pEq[0]/(2*np.pi*friction)
            fcy = ky*Dy/(kBT_pN_nm*2*np.pi)
            fcz = kz*Dz/(kBT_pN_nm*2*np.pi)
            fcr = kr*Dr/(kBT_pN_nm*2*np.pi)
            print("Corner frequencies (Hz) x,y, z, r: ", "%.3f" % fcx , "%.3f" %  fcy, "%.3f" %  fcz, "%.3f" %  fcr)
            print("Forces (pN) x,y, x2, y2, z, r, x3, y3: ", "%.3f" % Fx, "%.3f" % Fy, "%.3f" % Fx2, "%.3f" % Fy2, "%.3f" % Fz,  "%.3f" % Fr,  "%.3f" % Fx3,  "%.3f" % Fy3,
                  "Forces Errors (pN) x,y, x2, y2, z: ", "%.3f" % dFx, "%.3f" % dFy, "%.3f" % dFx2, "%.3f" % dFy2, "%.3f" % dFz, "%.3f" % dFr)
            Dxytheo=kBT_pN_nm/frictionXY; Dztheo=kBT_pN_nm/frictionZ; Drtheo=kBT_pN_nm/frictionR
            print("D_by_Dtheo x, y,z: ", "%.3f" %  (Dx/Dxytheo), "%.3f" %  (Dy/Dxytheo), "%.3f" %  (Dz/Dztheo), "%.3f" %  (Dr/Drtheo))
            FSDXFit=kBT_pN_nm*(p1Lo+BeadRad)/p2Xo**2; FSDXMod=kBT_pN_nm*(pLo+BeadRad)/p2Xo**2
            FSDYFit=kBT_pN_nm*(p1Lo+BeadRad)/p2Yo**2; FSDYMod=kBT_pN_nm*(pLo+BeadRad)/p2Yo**2
            print("Forces from SD (pN) xFit,xMod,yFit,yMod: ", "%.3f" % FSDXFit, "%.3f" % FSDXMod, "%.3f" % FSDYFit, "%.3f" % FSDYMod)

            columnsSpecRaw = 'fx,Pxx_specx,fy,Pxx_specy,fx'
            print('export spectrum', columnsSpecRaw, len(fx), len(Pxx_specx), len(fy), len(Pxx_specy))
            exportspectrum(fx, Pxx_specx, fy, Pxx_specy, np.append(ftx, np.zeros(len(fx)-len(ftx))), path1+OutputFolder+refname+'rawXYspectrum.csv', colname = columnsSpecRaw)
                  
            if forcecalib:
                results=(fmin_Hz, fmax_Hz, fcx , fcy, fcz, kx, ky, kx2, ky2, kz, kr, Dx, Dy, Dx2, Dy2, Dz, Dr, Dxytheo, Dztheo, Fx, Fy, Fx2, Fy2, Fx3, Fy3, Fz, dFx, dFy, dFx2, dFy2, dFz, FxC, dFxC, FyC, dFyC, Fx2C, Fy2C, Fx3C, Fy3C, FzC, dFx2C, dFy2C, dFzC, FSDXFit, FSDXMod, FSDYFit, FSDYMod,
                         Cp0Lo, Cp1Lo, Cp2Lo, CpLo, p0Lo, p1Lo, p2Lo, pLo, PullAngle, PullAngleC, PullAngleC2, PullPhi, PullPhiC)
                columns=['fmin_Hz', 'fmax_Hz', 'fcx_Hz' , 'fcy_Hz', 'fcz_Hz', 'kx_pN_per_nm', 'ky_pN_per_nm', 'kx2_pN_per_nm', 'ky2_pN_per_nm', 'kz_pN_per_nm', 'kr_pN_per_nm', 'Dx_nm2_per_s', 'Dy_nm2_per_s', 'Dx2_nm2_per_s', 'Dy2_nm2_per_s',
                         'Dz_nm2_per_s', 'Dr_nm2_per_s', 'Dxytheo_nm2_per_s',
                     'Dztheo', 'FSpecx_pN', 'FSpecy_pN', 'FSpec2x_pN', 'FSpec2y_pN', 'FSpec3x_pN', 'FSpec3y_pN', 'FSpecz_pN', 'dFSpecx_pN', 'dFSpecy_pN', 'dFSpec2x_pN', 'dFSpec2y_pN', 'dFSpecz_pN', 'FSpecx_Contour_pN', 'dFSpecx_Contour_pN', 'FSpecy_Contour_pN', 'dFSpecy_Contour_pN', 'FSpecz_Contour_pN',
                     'FSpecx2_Contour_pN', 'FSpecy2_Contour_pN', 'FSpecx3_Contour_pN', 'FSpecy3_Contour_pN',
                     'dFSpecx2_Contour_pN', 'dFSpecy2_Contour_pN', 'dFSpecz_Contour_pN','FSDXFit_pN', 'FSDXMod_pN', 'FSDYFit_pN', 'FSDYMod_pN',
                     'Cp0Lo', 'Cp1Lo', 'Cp2Lo', 'CpLo', 'p0Lo', 'p1Lo', 'p2Lo', 'pLo', 'PullAngle', 'PullAngleC', 'PullAngleC2', 'PullPhi', 'PullPhiC']



                listresults.append(results)
                listresultsTotal.append(results)
                df = pd.DataFrame(listresults, columns = columns)
                dfTotal = pd.DataFrame(listresultsTotal, columns =columns)
                dfTotal.to_csv(path1+OutputFolder+refname+'FitsResults.csv')             
                
                sys.exit('Force calibration: end of program')
    
            if GraphTZ_Power: AFSG.MakeGraphTZ_Power(refname, Tup, Zuplow, Zupmid, Zuphigh, Dup, T, P, SaveGraph, path1+OutputFolder, OutFormat)
            if GraphTupdown: AFSG.MakeGraphTupdown(refname, Tup, Tupdown, Tupjump, Tupdown_Tup, SaveGraph, path1+OutputFolder, OutFormat)
    
            mask_cycle = Tupjump_Tuphigh> -1000
            if bad_cycles!=[]:
                for c in bad_cycles: mask_cycle[c] = False
                
            Tupjump_Tuphigh_mask_cycle = Tupjump_Tuphigh[mask_cycle]
            Zupjump_Zupmid_mask_cycle = Zupjump_Zupmid[mask_cycle]
            print('Cycle Elimination:', len(Tupjump_Tuphigh), len(Tupjump_Tuphigh_mask_cycle))
                
            if CalculateSurvival and not forcecalib:
                nc = next(color)
           #     (pE0, pE1, x2, y2, y2_fit0, y2_fit)=survival(Tupjump_Tuphigh, TEhigh, countNoOpen, tshiftsurvival, refname, nc)
                exportsurvivalname = path1+OutputFolder+refname+'ExportSurvival_'
                (pEq0, dpEq0, pEq10, pEq11, x2, y2, y2_fit0, y2_fit)=survival(Tupjump_Tuphigh_mask_cycle, TEhigh, countNoOpen, tshiftsurvival, refname+'Mask_Cycle', nc, exportname=exportsurvivalname)
            else:
                pEq0=1.; dpEq0=0; pEq10=1; pEq11=0
                print('Skip Survival Curve')
     #       if CalculateSurvival: (pE0, pE1, x2, y2, y2_fit0, y2_fit)=survival(Tupjump_Tuphigh, TEhigh, 2., '', next(color))   # to plot all beads on th same figure
            if DisplayJumpGraphs:
                AFSG.JumpGraphs(refname, Zs, Zsavg, Zupjump_Zupmid, Tupjump_Tuphigh, Tup, Zuplow, Tupjump, Zupjump, Zupmid,
                           Tuphigh, Zuphigh, Zupzero, Dup, Phigh*(1+tol), T, P, Tclosed, Zsclosed, Topen, Zsopen,
                           Tl, Zsl, SaveGraph, CloseAfterSave, path1+OutputFolder, OutFormat)
 
            print(tdmsfilename, type(tdmsfilename))
            results=(tdmsfilename, b[ibead], bref, BeadRad, X0med, Y0med, NZavg, NZlavg, dzjump, NZrefavg, tol, sample, start[ibead], stop[ibead], fs, TEhigh, TElow, Phigh, Plow,
                     Xec0_mean, Yec0_mean, Xa, Ya, Za, Xsl_md, Ysl_md, Zsl_md, AnchorCov, SymmetryFactor, PullAngle, PullAngle2, PullPhi, p1Lo, p1Lc, p1Ll, p1Zl,p2Lo, p2Lc, p2Ll, p2Zl, pLo, pLc, pLl, pZl, 
                     X0ca, Y0ca, Xca, Yca, Xscl_md, Xscl_md, Cp0Lo, Cp1Lo, Cp2Lo, CpLo, Cp0Lc, Cp1Lc, Cp2Lc, CpLc, Cp0Ll, Cp1Ll, Cp2Ll, CpLl, PullAngleC, PullAngleC2, PullPhiC,
                     p0Xo, p1Xo, p2Xo, pXo, p0Yo, p1Yo, p2Yo, pYo, p0Zo, p1Zo, p2Zo, pZo, fmin_Hz, fmax_Hz,
                     kx, ky, kx2, ky2, kz, fcx, fcy, fcz, fcr, Dx, Dy, Dx2, Dy2, Dz, Dxytheo, Dztheo, 1/pEq0, dpEq0/pEq0**2, 1/pEq10, Fx, Fy, Fx2, Fy2, Fx3, Fy3, Fz, dFx, dFy, dFx2, dFy2, dFz, FxC, FyC,
                     Fx2C, Fy2C, Fx3C, Fy3C, FzC, dFxC, dFyC, dFx2C, dFy2C, dFzC, FSDXFit, FSDXMod, FSDYFit, FSDYMod,
                     step_detectON, thres_detect, windowsize_detect, countOpen, countNoClose, countNoOpen, Nc, No, Nh, Zupjump_Zupmid, Tupjump_Tuphigh,
                     Tupjump_Tuphigh_mask_cycle, Zupjump_Zupmid_mask_cycle, contourdataXYslXYsh, contourdataXY0XYsh, contourdataXYoriginlXYsh, Fx)  
            listresults.append(results)
            listresultsTotal.append(results)
            columns=['tdmsfilename', 'bead', 'refbeads', 'BeadRad_nm', 'X0med', 'Y0med', 'NZavg', 'NZlavg', 'dzjump_nm', 'NZrefavg', 'tol', 'sample', 'start_frame', 'stop_frame', 'fs_Hz', 'TEhigh_ms', 'TElow_ms', 'Phigh_%', 'Plow_%',
                     'X_mean@0_force', 'Y_mean@0_force','X_anchor', 'Y_anchor', 'Z_anchor', 'dX_anchor', 'dY_anchor', 'dZ_anchor', 'AnchorCovMatrix', 'SymmetryFactor', 'PullAngle_deg', 'PullAngle2_deg', 'PullPhi_deg', 'p1Lo_nm' ,'p1Lc_nm','p1Ll_nm','p1Zl_nm',
                     'p2Lo_nm' ,'p2Lc_nm','p2Ll_nm','p2Zl_nm','pLo_nm', 'pLc_nm', 'pLl_nm', 'pZl_nm',
                     'X0forceC_anchor', 'Y0forceC_anchor', 'XC_anchor', 'YC_anchor','dXC_low', 'dYC_low', 'Cp0Lo', 'Cp1Lo', 'Cp2Lo', 'CpLo', 'Cp0Lc', 'Cp1Lc', 'Cp2Lc', 'CpLc', 'Cp0Ll', 'Cp1Ll', 'Cp2Ll', 'CpLl', 'PullAngleC', 'PullAngleC2', 'PullPhiC',
                     'p0Xo_nm', 'p1Xo_nm', 'p2Xo_nm', 'pXo_nm', 'p0Yo_nm', 'p1Yo_nm', 'p2Yo_nm', 'pYo_nm', 'p0Zo_nm', 'p1Zo_nm', 'p2Zo_nm', 'pZo_nm', 'fmin_Hz', 'fmax_Hz',
                     'kx_pN_per_nm', 'ky_pN_per_nm', 'kx2_pN_per_nm', 'ky2_pN_per_nm', 'kz_pN_per_nm', 'fc_x_Hz', 'fc_y_Hz', 'fc_z_Hz', 'fc_r_Hz', 'Dx_nm2_per_s', 'Dy_nm2_per_s', 'Dx2_nm2_per_s', 'Dy2_nm2_per_s', 'Dz_nm2_per_s', 'Dxytheo_nm2_per_s',
                     'Dztheo', 'Offrate_s-1_1par', 'dOffrate_s-1_1par', 'Offrate_s-1_2par', 'FSpecx_pN', 'FSpecy_pN', 'FSpec2x_pN', 'FSpec2y_pN', 'FSpec2x_pN', 'FSpec2y_pN', 'FSpecz_pN', 'dFSpecx_pN', 'dFSpecy_pN', 'dFSpecx2_pN', 'dFSpecy2_pN', 'dFSpecz_pN', 'FSpecx_Contour_pN', 'FSpecy_Contour_pN',
                     'FSpecx2_Contour_pN', 'FSpecy2_Contour_pN', 'FSpecx3_Contour_pN', 'FSpecy3_Contour_pN', 'FSpecz_Contour_pN', 'dFSpecx_Contour_pN', 'dFSpecy_Contour_pN','dFSpecx2_Contour_pN', 'dFSpecy2_Contour_pN', 'dFSpecz_Contour_pN','FSDXFit_pN', 'FSDXMod_pN', 'FSDYFit_pN', 'FSDYMod_pN',
                     'step_detectON', 'thres_detect', 'windowsize_detect', 'countOpen', 'countNoClose', 'countNoOpen', 'Nc', 'No', 'Nh', 'Zupjump_Zupmid', 'Tupjump_Tuphigh',
                     'Tupjump_Tuphigh_mask_cycle', 'Zupjump_Zupmid_mask_cycle', 'contourdataXYslXYsh', 'contourdataXY0XYsh', 'contourdataXYoriginlXYsh', 'Input Force']    
    df = pd.DataFrame(listresults, columns = columns)
    AFSG.plotresults(refname, df, SaveGraph, path1+OutputFolder, OutFormat)
    dfTotal = pd.DataFrame(listresultsTotal, columns =columns)
    dfTotal.to_csv(path1+OutputFolder+fname[0][:-5]+'_FitsResults.csv')

    
if nbead>1 and SaveGraph: plt.figure('figAllSpec'); plt.savefig(path1+OutputFolder+fname[0]+'figAllSpec'+OutFormat, transparent=True)
# dfTotal = pd.DataFrame(listresultsTotal, columns =columns)
# dfTotal.to_csv(path1+OutputFolder+fname[0][:-5]+'_FitsResults.csv')

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
        if SaveGraph: plt.figure(listFig1[i]); plt.savefig(path1+OutputFolder+fname[0]+listFig1[i]+OutFormat, transparent=True)
    for i in range(nFig2):
        ax2[i].legend(fontsize=6)
        if SaveGraph: plt.figure(listFig2[i]); plt.savefig(path1+OutputFolder+fname[0]+listFig2[i]+OutFormat, transparent=True)
