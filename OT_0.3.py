# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 09:49:34 2020

@author: Laurent
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import dabest
import plotly.graph_objects as go
from plotly import plot
#from scipy.optimize import curve_fit
from scipy import stats
#from scipy.stats import wilcoxon

plt.close(fig='all')
path1="/home/laurent/DATA/Redaction/OTtubes/Fabio/Dec2020/"
SaveGraph=True; OutFormat=".jpg"; OutFold=path1+'Output20201215/'
if not os.path.exists(OutFold): os.makedirs(OutFold)
fused_CD3=True  # fusionb of cd3 and cd3_ucht1
countpopulations=True
plotswarm=False
plothist=False
plotcorr=False
test=True
testtable=True
estimation=False
SaveGraph=True
label=''

FF=path1+'fitted_data_all_curves.csv'
FF=path1+'fitted_data_all_curves_eta5.csv'
dF=pd.read_csv(FF, encoding = "ISO-8859-1")
print(dF.columns.tolist())
if fused_CD3: dF['Condition'] = dF['Condition'].str.replace(r'_ucht1', '')

def MultiSwarmPlots(wy, wymax, wymin, wpool, colhue, figname):
    fig = plt.figure(figname, figsize=(len(wpool)*2, len(wy)*2), dpi=100)
    for (iy,y) in enumerate(wy):
        for (ipool,pool) in enumerate(wpool):
            fig.add_subplot( len(wy), len(wpool),len(wpool)*iy+ipool+1)
            if colhue=="": ax = sns.swarmplot(x="Condition", y=y, data=pool,
                                              dodge=True, order=ordrepres)
            else: ax = sns.swarmplot(x="Condition", y=y, hue=colhue, data=pool,
                                     dodge=True, order=ordrepres)
            ax.set_ylim(wymin[iy],wymax[iy]); ax.set_ylabel(y)
            if iy!=len(wy)-1: ax.axes.get_xaxis().set_visible(False)
            if ipool!=0: ax.axes.get_yaxis().set_visible(False)
            if iy==0: ax.set_title(wpoolname[ipool]+' n='+str(len(y)), fontsize=6)
            #ax.set_yscale('log')
            ax.legend(loc="upper right", title=colhue, title_fontsize=5, prop={'size': 5})
            
def MultiHistPlots(axall, wy, wymax, wymin, wyRef, wpool, wpoolname, condCD, condLat, figname):
    for (iy,y) in enumerate(wy):
        bins=np.logspace(np.log10(wymin[iy]), np.log10(wymax[iy]), num=40)
   #     print('pool length=', len(wpool))
        for (ipool, upool) in enumerate(wpool):
            if condCD=='':
                pool=upool[upool['Latrunculine']==(condLat=='lat')]
            else:
                pool=upool[(upool['Condition']==condCD)&(upool['Latrunculine']==(condLat=='lat'))]
            if len(wpool)==1: ax = axall[iy]
            else: ax = axall[ipool, iy]
            sns.distplot(pool[y], kde=False, axlabel=y, bins=bins, 
                         label=condCD+condLat+':'+str(len(pool[y])), ax=ax)
            if iy==0: print(condCD, condLat, wpoolname[ipool], len(pool[y]) )
#            plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
            ax.set_xlim(wymin[iy],wymax[iy]); ax.set_xscale('log')
            ax.set_ylim(0,max(len(pool[y])/4.,30)); ax.set_ylabel(wpoolname[ipool])
            if ipool!=len(wpool)-1: ax.axes.get_xaxis().set_visible(False)       
            if iy!=0: ax.axes.get_yaxis().set_visible(False)
            #ax2=ax.twinx()
            if wyRef!=None: ax.plot([wyRef[iy],wyRef[iy]], [0,40], linestyle='--',
                                    c='k',  label='Schmitz '+ "%.4f" % wyRef[iy])
            ax.plot([np.median(pool[y]), np.median(pool[y])], [0,40], 
                         alpha=0.5, label='Median '+condLat+ " %.4f" % np.median(pool[y]))
            print(condLat, wpoolname[ipool], y, 'N=', len(pool[y]), 'Median=',
                  np.median(pool[y]))
            if iy==0: ax.legend(loc="upper right", title_fontsize=4, prop={'size': 5})
            
def MultiCorr(wy, wpool, wpoolname, condCD, condLat, figname0):
    ratio=0.5
    for (ipool, upool) in enumerate(wpool):
        #print('ipool=', ipool)
        figname=figname0+wpoolname[ipool]
        pool=upool[upool['Condition'].str.contains(condCD)&
                   (upool['Latrunculine']==(condLat=='lat'))]
#https://seaborn.pydata.org/examples/many_pairwise_correlations.html
        corr = pool[wy].corr(min_periods=5)  #      print(corr)
        mask = np.triu(np.ones_like(corr, dtype=np.bool))
        f, ax = plt.subplots(num=figname, figsize=(11*ratio,9*ratio))
        f.suptitle(figname, fontsize=14)
        cmap = sns.diverging_palette(220, 10, as_cmap=True)
        sns.heatmap(corr, mask=mask, cmap=cmap, vmin=-1., vmax=1., annot=True,
                    center=0, square=True, linewidths=.5,cbar_kws={"shrink": .5}) 
        ax.tick_params(axis='both', which='major', labelsize=5)
        if SaveGraph: plt.savefig(OutFold+figname+OutFormat)

wpool=[dF]; wpoolname=["All"]
            
wy1Schmitz=[1.5e-3, 0.2, 1.5e-3]; wy2Schmitz=[6e-3, 0., 0., 0., 0.]            
wy1=["E1 [pN/nm]", "E2 [pN/nm]", "E1N [pN/nm]"]; wy2=["eta [pN/nm.sec]", 
     "Slope [pN/nm]", "tBk1 [sec]", "tBkX [sec]", "tau [sec]"]

wy1=["E1 [pN/nm]", "E2 [pN/nm]", "E1N [pN/nm]", "eta [pN/nm.sec]"]; wy2=["tBk1 [sec]"]

wy1min=[0.0001, 0.001, 0.00001, 0.00001]; wy2min=[0.01]
wy1max=[2., 200., 0.1, 1.]; wy2max=[20.]
wdy1=[0.05, 0.2, 0.001, 0.01]; wdy2=[0.2]

ordrepres=['cd3', 'cd11a', 'cd45'] if fused_CD3 else ['cd3', 'cd3_ucht1', 'cd11a', 'cd45']

dF['ConditionFull'] = np.where(dF['Latrunculine']==True , dF['Condition']+'+Lat', dF['Condition']+'-Lat')

AD = dF[dF['category_string']=='adhesion']
TF = dF[dF['category_string']=='finTube']
TI = dF[dF['category_string']=='infinTube']
TU = dF[(dF['category_string']=='finTube')+(dF['category_string']=='infinTube')]

wpoolM3=[AD, TU, dF]; wpoolM3name=["Adhesions", "Tubes", "All"]
wpoolM4=[AD, TF, TI, dF]; wpoolM4name=["Adhesions", "Fin. Tubes", "Inf. Tubes", "All"]

if countpopulations:
    listgb=[['Condition','category_string','Latrunculine']]
    for (ipool,pool) in enumerate(wpool):
        for ngb in listgb:
            gb=pool.groupby(ngb); print(wpoolname[ipool]); print(gb.size())

if test:
    wpoolTU=[TU]; wpoolTUname=["Tubes"]
    listgb=[['Condition','Latrunculine']]
    ws1=[]; wn=[]; ws2=[]; wp=[]; wy0=[]
    diclat = {}
    for (ipool,pool) in enumerate(wpoolTU):
        for ngb in listgb:
            gb=pool.groupby(ngb); print(wpoolTUname[ipool]); print(gb.size())
            for y in wy1+wy2:
                print(y) # print( gb[y].size(), gb[y].median())
                values_per_group = [col for col_name, col in gb[y]]
                name_per_group = [col_name for col_name, col in gb[y]]
                for i1, v1 in enumerate(values_per_group):
                    n1=name_per_group[i1]
                    for i2, v2 in enumerate(values_per_group): 
                        n2=name_per_group[i2]
                        if i1<i2 and ( (n1[0]==n2[0]) or (n1[1]==n2[1]) ):
                            w, p = stats.ranksums(v1, v2);
                            if p<0.05:
                                if n1[0]==n2[0]: s1=n1[0]; s2='+/-Lat'
                                if n1[1]==n2[1]: s2='Lat'+str(n1[1]); s1=n1[0]+'/'+n2[0]; 
                                n=str(len(v1))+'/'+str(len(v2))                                
                                print(s1, s2, n, "%.4f" %p )
                          #      print(n1, len(v1), n2, len(v2), ' : p=' , "%.4f" %p)
                                wy0.append(y); ws1.append(s1); ws2.append(s2); wn.append(n); wp.append("%.4f" % p)
    if testtable:
        dftable = pd.DataFrame({'Param':wy0, 'Cond1':ws1,'Cond2':ws2,'N1/N2':wn,'p':wp})
   #     dftable.style.applymap(color, subset=['Date'])
   # https://stackoverflow.com/questions/56041337/how-to-draw-a-beautiful-colorful-table-with-pandas-or-other-package-in-python
        print(dftable)
    #    fig = go.Figure(data=[go.Table(header=dict(values=['Cond1', 'Cond2', 'N1/N2', 'p']),
      #           cells=dict(values=[ws1, ws2, wn, wp]))
       #              ])
    #    plot(fig)
#        if SaveGraph: plt.figure(fig); plt.savefig(OutFold+fig+OutFormat)

if plotswarm: MultiSwarmPlots(wy1, wy1max, wy1min, wpool, "Latrunculine", 'AllSwarmPlots1')

if plotcorr:
    listx = ["E2 [pN/nm]"]; listy = ["eta [pN/nm.sec]"]
    for ipool, upool in enumerate(wpoolM3):
        figname='Corr_'+listx[0][0:3]+listy[0][0:3]+wpoolM3name[ipool]
        fig = plt.figure(figname); ax = plt.gca(); fig.suptitle(figname, fontsize=14)
        sns.histplot(upool, x=listx[0], y=listy[0], log_scale=(True, True), ax=ax, label=figname)
        if SaveGraph: plt.figure(figname); plt.savefig(OutFold+figname+OutFormat)
    
    for listcondCD in [[''], ordrepres]:
       for condLat in ['nolat', 'lat']:
           for icondCD, condCD in enumerate(listcondCD):
                if plotcorr: MultiCorr(wy1+wy2, wpoolM3, wpoolM3name, condCD, condLat,
                                   'AllCorr'+condCD+'_'+condLat+'_'+label)

if plothist:   
    figname = 'CompareTU'
    fig , axall = plt.subplots(len(ordrepres), len(wy1+wy2), num=figname, 
                                       figsize=(len(wy1+wy2)*2, 6), dpi=100)
    fig.suptitle(figname, fontsize=14)
    for condLat in ['nolat', 'lat']:
        for (iy,y) in enumerate(wy1+wy2):
                bins=np.logspace(np.log10((wy1min+wy2min)[iy]), np.log10((wy1max+wy2max)[iy]), num=40)
                for (icondCD, condCD) in enumerate(ordrepres):
                    pool=TU[(TU['Condition']==condCD)&(TU['Latrunculine']==(condLat=='lat'))]
                    ax = axall[icondCD, iy]
                    sns.distplot(pool[y], kde=False, axlabel=y, bins=bins, 
                                 label=condCD+condLat+':'+str(len(pool[y])), ax=ax)
                    if iy==0: print(condCD, condLat, "Tubes", len(pool[y]) )
                    ax.set_xlim((wy1min+wy2min)[iy],(wy1max+wy2max)[iy]); ax.set_xscale('log')
                    ax.set_ylim(0,max(len(pool[y])/4.,40)); ax.set_ylabel(condCD)
                    if icondCD!=len(ordrepres)-1: ax.axes.get_xaxis().set_visible(False)
                    if iy!=0: ax.axes.get_yaxis().set_visible(False)
                    ax.plot([np.median(pool[y]), np.median(pool[y])], [0,40], 
                                 alpha=0.5, label='Median '+condLat+" "+ "%.4f" % np.median(pool[y]))
                    ax.legend(loc="upper right", title_fontsize=4, prop={'size': 5})
    if SaveGraph: plt.figure(figname); plt.savefig(OutFold+figname+OutFormat)
     
   
    for listcondCD in [[''], ordrepres]:          
        for icondCD, condCD in enumerate(listcondCD):
            fig , axall = plt.subplots(len(wpoolM4), len(wy1+wy2), num='AllHisto'+condCD, 
                                       figsize=(len(wy1+wy2)*2, 6), dpi=100)
            fig.suptitle(condCD, fontsize=14)
            MultiHistPlots(axall, wy1+wy2, wy1max+wy2max, wy1min+wy2min, 
                           wy1Schmitz+wy2Schmitz, wpoolM4, wpoolM4name, condCD, 'nolat',
                           'AllHisto'+condCD+'_'+'nolatNewCat'+'_'+label)
            MultiHistPlots(axall, wy1+wy2, wy1max+wy2max, wy1min+wy2min,
                           None, wpoolM4, wpoolM4name, condCD, 'lat', 
                           'AllHisto'+condCD+'_'+'latNewCat'+'_'+label)
            if SaveGraph: plt.figure('AllHisto'+condCD); plt.savefig(OutFold+'AllHisto'+condCD+OutFormat)

def MultiEstimate(wy, wymax, wymin, wdy, upool, condLat, xCond, idx, MeanMed, CI=0.95, figname='Default'):
    fig , axall = plt.subplots(1, len(wy), num=figname, figsize=(len(wy)*3.5, 6.),
                               gridspec_kw={'wspace': 0.25}, dpi=100)
    fig.suptitle(figname, fontsize=14)
    for (iy,y) in enumerate(wy):
  #      pool=upool[upool['Latrunculine']==(condLat=='lat')]
        pool=upool
        ax = axall[iy]
        multi_groups = dabest.load(pool, idx=idx, x=xCond, y=y)
  #      multi_groups = dabest.load(pool, idx=idx, x=xCond, y=y, ci=CI)
        if MeanMed=='Mean': multi_groups.mean_diff.plot(ax=ax)
        elif MeanMed=='Median': multi_groups.median_diff.plot(fig_size=(3.5, 6.), ax=ax)
        ax.set_title(y, fontsize=10)
        ax.set_ylabel('Parameter value') #
        if iy!=0: ax.set_ylabel(''); ax.contrast_axes.set_ylabel('')
        ax.xaxis.label.set_size(4)
        ax.tick_params(axis='x', which='major', labelsize=5)
        ax.contrast_axes.tick_params(axis='x', which='major', labelsize=7)
        ax.tick_params(axis="y", labelsize=8)
        ax.contrast_axes.tick_params(axis="y", labelsize=8)
        ax.set_ylim(wymin[iy],wymax[iy])
        ax.contrast_axes.set_ylim(-wdy[iy],wdy[iy])
        ax.set_yscale('log')
    if SaveGraph: plt.savefig(OutFold+figname+OutFormat)



wpoolEst=[TU]        
wpoolEstName=['Tubes']        
if estimation:
    for (ipool,pool) in enumerate(wpoolEst):
        MultiEstimate(wy1+wy2, wy1max+wy2max, wy1min+wy2min, wdy1+wdy2, pool,'', 
                      'Condition', ('cd3', 'cd11a', 'cd45'),
                      'Median', 0.95, 'Estimate'+wpoolEstName[ipool]+'_'+'+-Lat'+'_'+label)
        MultiEstimate(wy1+wy2, wy1max+wy2max, wy1min+wy2min, wdy1+wdy2, pool, '', 
                      'ConditionFull', ('cd3-Lat', 'cd11a-Lat', 'cd45-Lat'),
                      'Median', 0.95, 'Estimate'+wpoolEstName[ipool]+'_'+'-Lat'+'_'+label)
        MultiEstimate(wy1+wy2, wy1max+wy2max, wy1min+wy2min, wdy1+wdy2, pool, '',
                      'ConditionFull', ('cd3+Lat', 'cd11a+Lat', 'cd45+Lat'),
                      'Median', 0.95, 'Estimate'+wpoolEstName[ipool]+'_'+'+Lat'+'_'+label)
        MultiEstimate(wy1, wy1max, wy1min, wdy1+wdy2, pool, '',
                      'ConditionFull', (('cd3-Lat','cd3+Lat'), ('cd11a-Lat','cd11a+Lat'), ('cd45-Lat','cd45+Lat')),
                      'Median', 0.95, 'Estimate'+wpoolEstName[ipool]+'_'+'lat+-'+'_splitted'+label)