#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 17:42:14 2020
updates: 08/01/2021
@author: laurent
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import pandas as pd

def MakeGraphXYall(d0, imax, n, refrms, name, outname,  OutFormat, SaveGraph, corr=False, Xr0avg=None, Yr0avg=None, pprint=False):   
# Displays all XY traj with number and SD on the same graph
    if corr: name=name+"_Corr"
    figname=name+'_XYall'
    plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
    xmin=-1000; xmax=350000; ymin=-1000; ymax=250000
    for i in range(1,imax+1):
   #     (X,Y,Z,MinLUT,SSIM255)=track(d0, i)
        X = d0['ROI{:04d}'.format(i)+' X (nm)'][:]
        Y = d0['ROI{:04d}'.format(i)+' Y (nm)'][:]
        if corr: X=X-Xr0avg; Y=Y-Yr0avg
        xx=X[::n]; yy=Y[::n]
        ds=np.sqrt((X-X.mean())**2+(Y-Y.mean())**2)
        dsmean=ds.mean()
        if np.isnan(dsmean): lab=str(i)+'/'+'Nan'
        else: lab=str(i)+'/'+str(int(dsmean))
        if dsmean<refrms: lab=lab+'ref'
        if pprint: print('plot xy traj:'+lab)
        ax.scatter(xx, yy, marker='.', alpha=0.5, label=str(i))
        plt.text(xx[-1], yy[-1], lab, fontsize=6)
    ax.axis([xmin, xmax, ymin, ymax])
    if corr: ax.axis([xmin-np.mean(Xr0avg), xmax-np.mean(Xr0avg), ymin-np.mean(Yr0avg), ymax-np.mean(Yr0avg)])
    print(outname,figname,OutFormat)
    if SaveGraph: plt.savefig(outname+figname+OutFormat, transparent=True)
    
def GraphSelfReference(T, Zu, Z1, Z1avg, Zref, label, name, outname,  OutFormat, SaveGraph):
    figname=name+label
    plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
    ax.scatter(T, Zu, marker='.', c='b', label=label, alpha=0.2)
    ax.plot(T, Z1, c='b', label=label+'1', alpha=0.5); ax.legend()
    if Zref!=np.array([]): ax.plot(T, Zref, c='g', label=label+'ref', alpha=0.5)
    ax.plot(T, Z1avg, c='r', label=label+'1avg', alpha=0.5); ax.legend()
    if SaveGraph: plt.savefig(outname+figname+OutFormat, transparent=True)

def Plot3D(x,y,z,n, refname, outname, OutFormat, SaveGraph):    # 3D plot
    fig = plt.figure(refname, figsize=(7,7), dpi=100)
    xx=x[::n]; yy=y[::n]; zz=z[::n]
    ax = fig.add_subplot(111, projection='3d')
    # For each set of style and range settings, plot n random points in the box
    # defined by x in [23, 32], y in [0, 100], z in [zlow, zhigh].
    for m, zlow, zhigh in [('.', -50, -25), ('.', -30, -5)]:
        ax.scatter(xx, yy, zz, marker=m)
    if SaveGraph: plt.savefig(outname+refname+OutFormat, transparent=True)

def MakeGraphXY(X1, Y1, X2, Y2, xmin, xmax, ymin, ymax, label1, label2, refname, SaveGraph, 
                outname, OutFormat, levels=10, disperse=True, contour=True, line1=False, line2=False, log=False):
    figname=refname+'_XY'
    n1=1; n2=30    # one point every n2 points for calculating contour plots
    if disperse: n1=30 # one point every n1 points in X, Y scatter plots (for which contour==True)
    plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
    if line1: 
        ax.plot(X1, Y1, alpha=0.5, label=label1)
    else:
        ax.scatter(X1[::n1], Y1[::n1], marker='.', alpha=0.1, label=label1)
    if line2:
        ax.plot(X2, Y2, alpha=0.5, label=label2)
    else:
        if X2!=[]: ax.scatter(X2[::n1], Y2[::n1], marker='.', alpha=0.1, label=label2)
    if log: ax.set_yscale('log')

    if contour:
  #      XX,YY = np.meshgrid(X1[::n2],Y1[::n2])
        sns.kdeplot(x=X1[::n2], y=Y1[::n2], color='k', levels=levels, alpha=0.3, ax=ax, linestyles="--") #    sns.heatmap(X1[::30], Y1[::30], ax=ax)
        if X2!=[]:
      #      XX,YY = np.meshgrid(X2[::n2],Y2[::n2])
            sns.kdeplot(x=X2[::n2], y=Y2[::n2], color='k', levels=levels, alpha=0.6, ax=ax)
        print("Contour", refname, "total plots=", len(ax.collections), "levels=", levels)
        def centeroidnp(arr):
            length = arr.shape[0]
            sum_x = np.sum(arr[:, 0])
            sum_y = np.sum(arr[:, 1])
            return sum_x/length, sum_y/length
        contourdata=[]
        for i,c in enumerate(ax.collections):     
            if isinstance(c, matplotlib.collections.LineCollection):
                if len(c.get_segments())>0:
                     v = c.get_segments()[0]; cv=centeroidnp(v)                  
                     print("contour:", i, " nb of points", len(v), "centroid=", cv)
                     ax.plot(list(zip(*v))[0], list(zip(*v))[1], color='k', alpha=0.3, label="Outer Contour", linewidth=2.0)
                     ax.scatter(cv[0], cv[1], color='k', alpha=0.3, label="Contour center")
                     contourdata.append((i, len(v), cv))
                     break
        ax.scatter(np.mean(X1), np.mean(Y1), color='r', alpha=0.3, label='center'+label1)
        if X2!=[]: ax.scatter(np.mean(X2), np.mean(Y2), color='r', alpha=0.6, label='center'+label2)
    ax.legend(fontsize=6); ax.axis([xmin, xmax, ymin, ymax])
    if SaveGraph: plt.savefig(outname+figname+OutFormat, transparent=True)
    if contour: return contourdata

def SingleTrace(refname, d0, n,m, X,Y,Z,T, MinLUT, SaveGraph, outname, OutFormat):    # SingleTrace(30,10, 0,0,0,0,0)
# Traces time traces and histograms for bead n (step m) or coordinates if n<=0
    if n>0:
        X = d0['ROI{:04d}'.format(n)+' X (nm)'][:]
        Y = d0['ROI{:04d}'.format(n)+' Y (nm)'][:]
        Z = d0['ROI{:04d}'.format(n)+' Z (nm)'][:]
        MinLUT = d0['ROI{:04d}'.format(n)+' MinLUT'][:]
        lab=str(n)
        T = d0['Time (ms)'][:]
    T_s=(1/1000.)*T
    if n==0: lab='Corr by ref'
    if n==-1: lab='ref'
    if n==-2: lab='After Standard  anchor point'
    if n==-3: lab='After Contour anchor point'
    figname=refname+'SingleTrace'+lab
    fig=plt.figure(figname, figsize=(9,6), dpi=100)
    spec = fig.add_gridspec(ncols=2, nrows=2, width_ratios=[3, 1], height_ratios=[1, 1])
    print('Single Trace Plot'+lab+': Averages: X=', "%.3f" % (X.mean()), 'Y=', "%.3f" % (Y.mean()), 'Z=', "%.3f" % (Z.mean()))
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
    if SaveGraph: plt.savefig(outname+figname+OutFormat, transparent=True)

def MakeDisplayCloseOpen(refname, Tclosed, Zsclosed, Topen, Zsopen, Tl, Zsl, T, Zsavg, P, Tup, SaveGraph, CloseAfterSave, outname, OutFormat, top=False):               
    figname=refname+'XY_'+'TZclosedopen'
  #  fig = plt.figure(figname, figsize=(8,6), dpi=100)
    fig = plt.figure(figname, figsize=(15,4), dpi=100)
    if top:
        spec = fig.add_gridspec(ncols=1, nrows=2, width_ratios=[1], height_ratios=[1, 0.2])
        ax = fig.add_subplot( spec[0, 0]); ax2 = fig.add_subplot( spec[1, 0])
    else:
        ax = plt.gca(); ax2=ax.twinx()
    timescale = 60000
    ax.scatter(Tclosed/timescale, Zsclosed, marker='.', alpha=0.2, label='Zsclosed')
    ax.scatter(Topen/timescale, Zsopen, marker='.', alpha=0.2, label='Zsopen')
    ax.scatter(Tl/timescale, Zsl, marker='.', alpha=0.2, label='Zsl')
    ax.plot(T/timescale, Zsavg, c='m', alpha=0.5, label='Zsavg')
    ax2.plot(T/timescale, P, c='k', alpha=0.3, label='P'), ax2.legend(fontsize=6); 
 #   ax.legend(fontsize=6); ax.axis([0, T[-1]/timescale, -100, 1100])
    ax.axis([0, 250, -100, 1000])   # modifs for final graph article
  #  ax.axis([174.2, 181, -100, 1000])   # modifs for final graph article
    if top:
        ax2.set_xlabel("Time (min)"); ax2.set_xlim(ax.get_xlim())
    else: 
        ax.set_xlabel("Time (min)")
    ax.set_ylabel("Height (nm)"); ax2.set_ylabel("Power (%)")
    plt.tight_layout()
 #   for ic, tic in enumerate(Tup): plt.text((tic+0.3e5)/timescale, np.amax(P)*0.8 , str(ic), fontsize=10)
    if SaveGraph: plt.savefig(outname+figname+OutFormat, transparent=False)
    if CloseAfterSave: plt.close()

def MakeGraphTZ_Power(refname, Tup, Zuplow, Zupmid, Zuphigh, Dup, T, P, SaveGraph, outname, OutFormat):
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
    if SaveGraph: plt.savefig(outname+figname+OutFormat, transparent=True)

def MakeGraphTupdown(refname, Tup, Tupdown, Tupjump, Tupdown_Tup, SaveGraph, outname, OutFormat):
    figname=refname+'Tupdown vs Tup_'
    plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
    ax.scatter(Tup, Tupdown, c='b', marker='o', alpha=0.5, label='Tupdown')
    ax.scatter(Tup, Tupjump, c='r', marker='o', alpha=0.5, label='Tupjump')
    ax.plot(Tup, Tup, c='k', alpha=0.5, label='Tup')
    ax.legend(fontsize=6)
    plt.figure('Tupdown_Tup histo_'+refname, figsize=(6,6), dpi=100); ax = plt.gca()
    ax=sns.distplot( Tupdown_Tup, kde=False, bins=np.linspace(-10.,200.,num=20) )
    plt.figure('Tup-Tupdown_'+refname, figsize=(6,6), dpi=100); ax = plt.gca()
    ax.scatter(np.arange(len(Tup)), Tup, c='r', marker='o', alpha=0.5, label='Tup')
    ax.scatter(np.arange(len(Tup)), Tupdown, c='b', marker='o', alpha=0.5, label='Tupdown')
    ax.scatter(np.arange(len(Tup)), Tupjump, c='k', marker='o', alpha=0.5, label='Tupjump')
    ax2 = ax.twinx(); plt.grid(); ax.legend(fontsize=6)
    ax2.scatter(np.arange(len(Tup)), Tupdown_Tup, c='g', marker='.', alpha=0.5, label='Tupdown_Tup')
    ax2.legend(fontsize=6)
    ax.set_xlabel("Event Number");  ax.set_ylabel("Time (ms)")
    if SaveGraph: plt.savefig(outname+figname+OutFormat, transparent=True)
    #ax2.set_ylim(0.,0.75)

def JumpGraphs(refname, Zs, Zsavg, Zupjump_Zupmid, Tupjump_Tuphigh, Tup, Zuplow,Tupjump, Zupjump, 
               Zupmid, Tuphigh, Zuphigh, Zupzero, Dup, Pmax, T, P, Tclosed, Zsclosed, Topen, Zsopen,
               Tl, Zsl, SaveGraph, CloseAfterSave, outname, OutFormat):
    figname=refname+'dZjumphisto_'
    plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
    ax=sns.distplot( Zupjump_Zupmid, kde=False, bins=np.linspace(-100.,500.,num=20) )
    ax.set_xlabel("Jump height (nm)");  ax.set_ylabel("Number"); plt.xlim(0, 500.);
    if SaveGraph: plt.savefig(outname+figname+OutFormat, transparent=True)
    
    figname=refname+'JumpHigh dZ vs dT'
    plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
    ax.scatter(Tupjump_Tuphigh, Zupjump_Zupmid, c='b', marker='o', alpha=0.5, label='dZ vs dT')
    ax.set_xlabel("dt jump (s)");  ax.set_ylabel("dZ jump (nm)")
    plt.figure('JumpHigh dZ vs Zzero_'+refname, figsize=(6,6), dpi=100); ax = plt.gca()
    ax.scatter(Zupzero, Zupjump_Zupmid, c='b', marker='o', alpha=0.5, label='dZ vs Zzero')
    ax.set_xlabel("Z zero (nm)");  ax.set_ylabel("dZ jump (nm)")    
    if SaveGraph: plt.savefig(outname+figname+OutFormat, transparent=True)

    figname=refname+'TZ'
    plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
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
    if SaveGraph: plt.savefig(outname+figname+OutFormat, transparent=True)
    if CloseAfterSave: plt.close()  

def MakeDisplayplotvariableNZlavg(DisplayplotvariableNZlavg, ax1, ax2, dfT, c, ib, DeltaCOFit, DeltaCOMod):    
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
    ax2[0].scatter(dfT['FSpecx_pN'], dfT['FSpecy_pN'], marker='.', c=c, alpha=0.5, label=str(ib)+' Angle='+str(dfT['PullAngle_deg'][0]))
    ax2[0].set_xlabel("ForceSpecX_pN");  ax2[0].set_ylabel("ForceSpecY_pN")
    ax2[1].scatter(dfT['FSpecx_pN'], dfT['FSDXFit_pN'], marker='.', c=c, alpha=0.5, label=str(ib))
    ax2[1].set_xlabel("ForceSpec_pN");  ax2[1].set_ylabel("ForceSD_pN")
    
def plotresults(refname, df, SaveGraph, outname, OutFormat):
    figname=refname+'Variation_NZlavg'
    plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
    for y in ['Lo_nm', 'Lc_nm', 'Ll_nm', 'Zl_nm']:
        ax.errorbar(df['NZlavg'], df['p1'+y], df['p2'+y], label=y, marker='o', ms=10)
    ax.legend(fontsize=6); ax.set_xlabel("NZlavg");  ax.set_ylabel("Peak length or Z_nm")
    ax.axis([0, 1000, -100, 1100])
    if SaveGraph: plt.savefig(outname+figname+OutFormat, transparent=True)

