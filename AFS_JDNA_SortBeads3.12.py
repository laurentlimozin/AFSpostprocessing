# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 12:23:26 2020

@author: Laurent
"""
import numpy as np
import os, re, six
import matplotlib.pyplot as plt
import pandas as pd
import random
import math
import gc
import seaborn as sns
import statsmodels.api as sm
from matplotlib.pyplot import cm
from matplotlib.lines import Line2D
from itertools import cycle
from scipy import signal
from scipy.optimize import curve_fit
from scipy.stats import pearsonr, linregress
from lmfit import minimize, Parameters, fit_report, Model
import warnings; warnings.filterwarnings("ignore")
from matplotlib.axes._axes import _log as matplotlib_axes_logger
matplotlib_axes_logger.setLevel('ERROR')
plt.close('all')

dmax=1.e3; tagfontsize=7; nmaxcolor=10
DisplaySurvival = True
SortingBeads = False
SaveGraph = True
Domultiplots = False
SelectMolecule = True
SelectForce = False

# https://seaborn.pydata.org/generated/seaborn.lmplot.html

rapa = False

#outfolder = 'Out20220310_Rapa_AllMolec/'; rapa = True; nrow=11; ncol=17; nrow2=5; ncol2=7; deltaf = 2; minsize = 10; nmaxfit=50;  SelectForce = False
if rapa: 
    outfolder = 'Out20220519_Rapa/'
    nrow=9; ncol=16; nrow2=3; ncol2=7; deltaf = 2; MinForce = 1; MaxForce = 25; minsize = 10; nmaxfit=50; SelectForce = False
    HighPowerMax = 35; DeltaP=0.2
else:
    outfolder = 'Out20220519_Nef/'
    nrow=9; ncol=16; nrow2=3; ncol2=3; nf = 9; deltaf = 3.5; MinForce = 5; MaxForce = 30; minsize = 4
    HighPowerMax = 10.5; DeltaP=0.05
print('Outfolder: ', outfolder)
inputpath0="/home/laurent/DATA/Experiences Diverses/AFS/Data/SelectionJim/"

dfstage_nef = pd.read_excel(inputpath0+'Stage Summary_Nef_20220323.xlsx')
dfstage_nef = pd.read_excel(inputpath0+'Stage Summary_Nef_20220405_v5.xlsx')
#ListMolec_nef = np.array([1,2,3,5,6,7,8,9,10,11,12,13,14,15,16])
ListMolec_nef = np.array([1,2,5,6,12,13,14,15,16,18,19,20])
ListMolec_nef = np.array([1,2,5,13,14,16,18,19,20])

dfstage_rapa = pd.read_excel(inputpath0+'Summary_Rapamycin_sorted-by-molecule-ID_length_2022-03-04.xlsx')
ListMolec_rapa = dfstage_rapa['Molecule ID'].unique()
removeMolec = np.array([18,20,23,24,25,26,27,28,29,30,34,37,38])
ListMolec_rapa = np.setdiff1d(ListMolec_rapa,removeMolec)
removeMolec = np.array([22, 31, 32, 33, 39, 40, 41, 42, 43, 44, 45 ])
ListMolec_rapa = np.setdiff1d(ListMolec_rapa,removeMolec)
MolecIDa_fmin_Hz_set=14; MolecIDb_fmin_Hz_set=21; fmin_Hz_set=0.2
MolecIDa_fmin_Hz_set=22; MolecIDb_fmin_Hz_set=21; fmin_Hz_set=0.2

startstr = ''; endstr = '.csv'
if rapa:
    folder_cycle = "Output_csv file_cycle_2021-11-30 rapa update/"
    folder_calib = "Rapa_Calib_2022/" # previously  "Output_csv file_2022-Force calibration/"
    folder_spectra = "Rapa_Spectra_FOV1_8/"
    readfiles_calib = [f for f in os.listdir(inputpath0 + folder_calib) if ( f.endswith(endstr) & f.startswith(startstr) ) ]
    print(folder_calib , len(readfiles_calib),'loaded files')
    ziplist_calib = zip([readfiles_calib], [inputpath0 + folder_calib])
else:
    folder_cycle = "Output_csv files_2021-11-03 and 20 nef/"
 #   folder_spectra = "Output_csv files_raw PSD_force calibration Nef/"
    folder_spectra = "Output_csv files_4.57_2022-04-13 Raw Spectra Nef/"
    folder_spectra_gold = "20201230-144513_Bead1_goldenbeadrawspectra/"

sizefig = 1.5
mindur=0.5; mindz = 100; maxdz = 300
kres = 100; kBT = 4 # pN.nm
Temperature_C = 25; kBT_pN_nm= 1.38e-23*(Temperature_C+273)*1.e12*1.e9
deltaD = 0.4
#BeadRad=790. # Bead Radius in nm
kc = 0.5
p1Zo = 900 # typical open JDNA extension
def FaxenXY(R, Z):  return 1 / ( 1 - (9/16)*(R/(Z+R)) + (1/8)*(R/(Z+R))**3 )
def FaxenZ(R, Z):  return 1 / ( 1 - (9/8)*(R/(Z+R)) + (1/2)*(R/(Z+R))**3 )
#Temperature_C = 25; kBT_pN_nm= 1.38e-23*(Temperature_C+273)*1.e12*1.e9

if not rapa:
    koffSPR = 1.8e-4; tmaxsurvival = 1500; minfrac=0.1
    dfstage0 = dfstage_nef
    ListMolec = ListMolec_nef
    dfstage0['kx (pN/nm)']=0
else:
    koffSPR = 22e-3; tmaxsurvival = 150; minfrac=0.01; tmaxsurvival = 50; minfrac=0.1
    dfstage0 = dfstage_rapa
    ListMolec = ListMolec_rapa
Nmolec = len(ListMolec)

def f(u): return u.replace(' ',',').split(',')
dfstage0['Trace'] = dfstage0['Trace'].fillna(0)
dfstage0 = dfstage0.astype({'Tdms (cyclic measurement)': str, 'Tdms_force calibration': str, 'Trace':int, 'Trace calib':int})
dfstage0['TdmsList'] = [f(x) for x in dfstage0['Tdms (cyclic measurement)'].values]
dfstage0['TdmsList_calib'] = [f(x) for x in dfstage0['Tdms_force calibration'].values]
if rapa:
    dfstage0['Force_pN'] = [float(str(x).replace('pN','')) for x in dfstage0['Force_Python_NoR'].values]
    dfstage0['dForce_pN'] = [float(str(x).replace('pN','')) for x in dfstage0['Error bar_F_python'].values]
    dfstage0['ForceLumicks_pN'] = [float(str(x).replace('pN','')) for x in dfstage0['Force_Lumicks'].values]
    dfstage0['dForceLumicks_pN'] = [float(str(x).replace('pN','')) for x in dfstage0['Error bar_F_lumicks'].values]
else:
    dfstage0['Force_pN']=np.nan; dfstage0['dForce_pN']=np.nan; dfstage0['ForceLumicks_pN']=np.nan; dfstage0['dForceLumicks_pN']=np.nan
    dfstage0['Length_cyclic_open (nm)']=np.nan; dfstage0['error of length_cyclic_open']=np.nan   
    dfstage0['Length_calibration (nm)']=np.nan; dfstage0['error of length']=np.nan
    dfstage0['Length_cyclic_close (nm)']=np.nan; wdlcc = dfstage0['error of length_cyclic_close']=np.nan
    dfstage0['No Close']=np.nan; dfstage0['Total cycles']=np.nan
    dfstage0['Force_Python_NoR']=np.nan

if rapa: dfstage0['% High power'] = 100*dfstage0['% High power']

dfstage0['FOV'] = 0
ListFOV = np.array( [(1,5), (6,13), (14,23), (31,32), (33,33), (35,37), (40,42), (43,45)] )
if rapa:
    for iFOV, boundFOV in enumerate(ListFOV):
        dfstage0['FOV'] = np.where( (dfstage0['Molecule ID']>=boundFOV[0]) & (dfstage0['Molecule ID']<=boundFOV[1]) , iFOV+1, dfstage0['FOV'])
    
dfstage0['Intern'] = ''
if rapa: dfstage0['Intern'] = np.where( dfstage0['Molecule ID']<=33 , 'Maryne', 'Adrien')
if rapa: dfstage0['L0'] = np.where( dfstage0['Molecule ID'] <= 6 , 900, 1300)
dfstage0['ForceOk']=1

if SelectForce: dfstage0 = dfstage0[ (dfstage0['Intern']=='Adrien') | ( (dfstage0['Intern']=='Maryne') & (dfstage0['% High power']<=20)) ]

dfstage0 = dfstage0.sort_values(by=['Molecule ID', '% High power'], ignore_index=True)#, inplace=True)

if rapa: 
    dfstageForceOk = dfstage0[dfstage0['ForceOk']>0]
else:
    dfstageForceOk = dfstage0[dfstage0['Molecule ID'].isin(ListMolec)]

colormap = cm.hsv(np.linspace(0, 1, Nmolec+1))
colormapFOV = cm.rainbow(np.linspace(0, 1, len(ListFOV)))
SelecMarker = {'.': 'point','o': 'circle','s': 'square','p': 'pentagon','*': 'star','h': 'hexagon1','8': 'octagon','+': 'plus',
 'x': 'x','D': 'diamond', 'd': 'thin_diamond','|': 'vline','H': 'hexagon2','P': 'plus_filled','X': 'x_filled'}
poolmarker = cycle(SelecMarker) # poolmarker = cycle(Line2D.markers)
poolmarker2x = cycle(SelecMarker) 
poolmarker2y = cycle(SelecMarker) 
#if rapa: dFForceRapa = pd.read_excel(inputpath0+'Summary_Rapamycin_high sampling rate.xlsx')

#######

readfiles_cycle = [f for f in os.listdir(inputpath0 + folder_cycle) if ( f.endswith(endstr) & f.startswith(startstr) ) ]
ziplist = zip([readfiles_cycle], [inputpath0 + folder_cycle])    
print(folder_cycle , len(readfiles_cycle),'loaded files')
    
outpath = inputpath0+outfolder
if not os.path.exists(outpath): os.makedirs(outpath)     

df=[]; df_calib=[]; 
print(" =========== Cycle files:", folder_cycle)
for readfiles, inputpath in ziplist:
    for (b,file) in enumerate(readfiles):
        print(file)
        dfx = pd.read_csv(inputpath+file)   #   dfx['FileName'] = re.split('_', file)[0]
        dfx['FileName'] = file
        df.append((file,dfx))
if rapa:
    print(" =========== Calibration files:", folder_calib)
    for readfiles, inputpath in ziplist_calib:
        for (b,file) in enumerate(readfiles):
            print(file)
            dfx = pd.read_csv(inputpath+file)   #   dfx['FileName'] = re.split('_', file)[0]
            dfx['FileName'] = file
            df_calib.append((file,dfx))
        
dff=[]; dff_calib=[]; dff_spectra=[] 
for i, idf in enumerate(df): dff.append(idf[1])
df0=pd.concat(dff, ignore_index=True)  
if rapa:
    for i, idf in enumerate(df_calib): dff_calib.append(idf[1])
    df0_calib=pd.concat(dff_calib, ignore_index=True)  

def g(u): 
    u = u.replace('[','')
    u = u.replace(']','')
    u = u.replace('\n','')
    ' '.join(u.split())
    u = u.replace(' ',',')
    u = u.split(',')
    u = [i for j, i in enumerate(u) if i!='']
 #   print(u)
 #   if len(~np.isnan(np.array(u)))>0:
    u = np.array(list(map(float, u)))

    return u

def TextCorrel(wx, wy):
    nonan = ~np.isnan(wx) * ~np.isnan(wy)
    if sum(nonan)>=2:
        (correl, p ) = pearsonr(wx[nonan], wy[nonan])
    else:
        correl=0
    return 'Correlation '+"{:1.2f}".format(correl)
def TextLinear(wx, wy):
    nonan = ~np.isnan(wx) * ~np.isnan(wy)
    slope, intercept, r, p, se = linregress(wx[nonan], wy[nonan])
    return 'Slope '+"{:1.2f}".format(slope)

#    def FitExp1(x, *p): return p[1]+(1-p[1])*np.exp(-x/p[0])
def FitExp0(x, *p): return np.exp(-x/p[0])
    
def survival(dur, TEhigh, countNoOpen, shift, refname, tmaxsurvival = tmaxsurvival, mindur=0, ax=None, col='k', exportname=None, shortlegend=None):
# calculates, fits and plots survival curve until TEhigh/1000.-shift
    print('Survival:', len(dur),'countNoOpen=', countNoOpen )
    dur2s = np.sort(dur[dur>mindur])
    print('nb events:', len(dur2s), 'without Nan:', len(dur2s[~np.isnan(dur2s)]) )
    dur2s = dur2s[~np.isnan(dur2s)]
    dur2s_y = np.arange(len(dur2s))
    dur2s_y = 1-dur2s_y/len(dur2s)
    if DisplaySurvival:
        if shortlegend!=None: ax.scatter(dur2s, dur2s_y, marker=next(poolmarker), c= col, alpha=0.5, label=shortlegend)
        if shortlegend==None: ax.scatter(dur2s, dur2s_y, marker=next(poolmarker), c= col, alpha=0.5, label=refname)
        ax.set_xlabel("Time (s)");  ax.set_ylabel("Survival Fraction")
        ax.axis([0., TEhigh/1000., 0.01, 1.]); ax.set_yscale('log')
    maxtfit=0.95*TEhigh/1000-shift
    imax = min(np.sum(dur2s<maxtfit), len(dur2s)-2)
    if not rapa: imax = np.sum(dur2s<maxtfit)
    imax = np.sum(dur2s<maxtfit)-1
    x2 = dur2s[dur2s<maxtfit]
    y2 = dur2s_y[dur2s<maxtfit]   #    y2log = dur2s_logy[dur2s<maxtfit]
    if rapa: imax = min(imax, nmaxfit)
    x2 = dur2s[range(imax+1)]
    y2 = dur2s_y[range(imax+1)]
    if len(y2)>0:
        while True:
            try:
                initialguessEq=[20]
                pEq0, pcovEq0=curve_fit(FitExp0, x2, y2, p0=initialguessEq)
                break
            except RuntimeError:
                print("No convergence"); pEq0=[1.]
                break
            except TypeError:
                print("Improper input"); pEq0=[1.]
                break
    else: pEq0=[1.]; pcovEq0=[[0]]
    y2_left = dur2s_y[imax-1]; x2_left = dur2s[imax-1]; koffmin = (1 - y2_left) / x2_left
    y2_right = dur2s_y[imax]; x2_right = dur2s[imax]; koffmax = (1 - y2_right) / x2_right
    y2_fit0=FitExp0(x2, pEq0[0])
    pEq=[1.,0.]
    if DisplaySurvival:
        if shortlegend==None: ax.plot(x2, y2_fit0, c='k',alpha=0.5, label='t fit='+"%3.0f" % (pEq0[0]) )
        if shortlegend!=None: ax.plot(x2, y2_fit0, c='k',alpha=0.5)
    #    if shortlegend!=None: ax.legend(fontsize=6, loc="upper right")
        ax.set_xlim(0,tmaxsurvival); ax.set_ylim(minfrac,1)
    print('Survival 1param fit 0-', maxtfit, 's: offrate=', '{:.2E}'.format(1/pEq0[0]), '+/-', '{:.2E}'.format(np.sqrt(pcovEq0[0][0])/pEq0[0]**2) , 'tau (s)=', '{:.2E}'.format(pEq0[0]))
    return(pEq0[0], np.sqrt(pcovEq0[0][0]), min(koffmin, koffmax), max(koffmin, koffmax), x2, y2, y2_fit0)

#def FitSpectrum(x, *p): return 4*kBT_pN_nm**2/(p[1]*p[0]*p[0]) * 1/( 1+ (x*2*np.pi*kBT_pN_nm/(p[1]*p[0]))**2 )
        
def QuickSpectrumFit(f, Pxx_spec, axis, GuessDTheo = 1, pprint=True):
    if np.sum(np.isnan(Pxx_spec))>0:
        print("Array contains Nan"); pEq=[np.nan, np.nan]; eEq=[np.nan, np.nan]
    else:
        while True:
            try:
                pEq, pcovEq = curve_fit(FitSpectrumGen, f, Pxx_spec, p0=[1.e-3, GuessDTheo])
                eEq=np.sqrt(np.diag(pcovEq))
                break
            except RuntimeError:
                print("No convergence"); pEq=[np.nan, np.nan]; eEq=[np.nan, np.nan]; break 
    pEq[0]=np.abs(pEq[0])
    fc = pEq[0]*pEq[1]/(2*np.pi*kBT_pN_nm)   #  modification 26/01/2022    fc=pEq[0]/(2*np.pi*friction)
  #  if pprint: print('QuickSpectrumFit: friction / friction_bulk =', friction/friction0, ' Corner frequency fc=', fc)
    FitPxx_spec=FitSpectrumGen(f, pEq[0], pEq[1])
#    if pprint: print('Spectrum', ' k (pN/nm)=',"%.5f" %  (pEq[0]),' D (µm²/s)=',"%.3f" %  (pEq[1]*1.e-6), 'D/Dtheo=', pEq[1]*friction/kBT_pN_nm), 'fc (Hz)= ',"%.3f" %  fc)
    if pprint: print('Single Spectrum', ' k (pN/nm)=',"%.5f" %  (pEq[0]),' D (µm²/s)=',"%.3f" %  (pEq[1]*1.e-6), 'fc (Hz)= ',"%.3f" %  fc)
    return pEq[0], pEq[1], eEq[0], eEq[1], fc, FitPxx_spec

def FitSpectrumGen(f, *p): return 4*kBT_pN_nm**2/(p[1]*p[0]*p[0]) * 1/( 1+ (f*2*np.pi*kBT_pN_nm/(p[1]*p[0]))**2 )  # k = p[0]; D = p[1]
def FitSpectrumGenR(L, f, *p):
 #   friction0 = 6*np.pi*1.e-9*p[1]; frictionXY=friction0 * FaxenXY(BeadRad_i, p1Zo)
    D = kBT_pN_nm/( 6*np.pi*1.e-9*p[1] * FaxenXY(p[1], L) )
    return 4*kBT_pN_nm**2/(D*p[0]*p[0]) * 1/( 1+ (f*2*np.pi*kBT_pN_nm/(D*p[0]))**2 )

def GlobalSpecFit(dfstage, folder_spectra, MolecID, xy):   # D, listk = GlobalSpecFit(dfstage, 1)
# dfstage0['R'+xy] = kBT_pN_nm/FaxenXY(dfstage0['Bead radius (nm)'], p1Zo)/(6*np.pi*1e-9*dfstage0['D'+xy+'global'])   
    dfs = dfstage[dfstage['Molecule ID']==MolecID].reset_index()
    n = dfs['TdmsList']; n_calib = dfs['TdmsList_calib']
    npower = dfs.shape[0]; print('Global Fit Spec ID:', MolecID, ' #powers=', npower, xy)
    params = Parameters()
   # params.add('D',  value=0.15e6, min=0.01e6, max=1e6)
    params.add('R',  value=800, min=600, max=2000)
    listR = []; listdR = []
    listD = []; listdD = []
    datax = {}; listfile = []
    listk = []; listdk = []; listL = []
    for inb, nt in enumerate(n.values):
        prefixe = n_calib[inb][0]+'_Bead'+str(int(dfs['Trace calib'][inb]))
        if rapa and MolecID>=MolecIDa_fmin_Hz_set and MolecID<=MolecIDb_fmin_Hz_set:
            prefixe = n[inb][0]+'_Bead'+str(int(dfs['Trace'][inb]))
        if rapa:
            power = int(dfs['% High power'][inb]); suffixe = str(power)+'%.csv'
        else:
            power = dfs['% High power'][inb]; suffixe = f'{power:g}'+'%.csv'

        fname_fitresult = prefixe +'_FitsResults_'+ suffixe
        print(' FitResults file', fname_fitresult)
        fullnameRes = inputpath0 + folder_spectra + fname_fitresult
        if os.path.isfile(fullnameRes):
            dffit = pd.read_csv(fullnameRes, index_col=False)
            L = dffit['Cp1Lo'].values[0] if rapa else dffit['Cp1Lo'].values[0]
        else:
            print(' missing file'); L = 900 %np.nan
   #     L = 1000
        listL += [L]
        datax['L'+str(inb)] = L
        fname_spectrarawXY = prefixe +'_rawXYspectrum_'+ suffixe
        print(' Spectrum file', fname_spectrarawXY)#, end='')
        fullnameSpec = inputpath0 + folder_spectra + fname_spectrarawXY
        print(n[inb],n_calib[inb][0], end='')
        if os.path.isfile(fullnameSpec) and (n_calib[inb][0] not in n[inb] ):#n[inb][0]!=n_calib[inb][0]:
            print(' specific calibration file')
            listfile += [1]
            params.add('k'+str(inb), value=1e-2, min=1e-4, max=1e-1, vary=True)
            dfraw = pd.read_csv(fullnameSpec, index_col=False)                 
            fmin_Hz = np.min(dfraw['fx.1'][dfraw['fx.1']>0]); fmax_Hz = np.max(dfraw['fx.1'][dfraw['fx.1']>0])
            if power>=25: fmin_Hz=1           
            if rapa and MolecID>=MolecIDa_fmin_Hz_set and MolecID<=MolecIDb_fmin_Hz_set:
                fmin_Hz=fmin_Hz_set
                if power>=20: fmin_Hz=0.5
                if power>=25: fmin_Hz=1
                if power>=30: fmin_Hz=5
            if not rapa and fmin_Hz<1: fmin_Hz=0.6 #fmin_Hz=0.3
            if not rapa and fmax_Hz>100: fmax_Hz=150 #fmin_Hz=0.3
            fr = dfraw['fx']; intfr = (fr>=fmin_Hz)&(fr<=fmax_Hz); frt=fr[intfr]
            print(f' fmin_Hz= {fmin_Hz:3.2f} fmax_Hz={fmax_Hz:3.2f}')
            datax['frt'+str(inb)] = frt
            datax['Pxy'+str(inb)] = dfraw['Pxx_spec'+xy][intfr]    
        elif n_calib[inb][0] in n[inb]:
            print(' cyclic file for calibration'); listfile += [0]       
        else:
            print(' missing calibration file'); listfile += [0]         
    def TotalFit(param, x=None, aax=None):
        inb0 = -1
        for inb in range(npower):
            if listfile[inb] == 1:
                inb0+=1
            #    mi = FitSpectrumGen(datax['frt'+str(inb)], param['k'+str(inb)], param['D'])
                mi = FitSpectrumGenR(datax['L'+str(inb)], datax['frt'+str(inb)], param['k'+str(inb)], param['R'])
                ri_ = datax['Pxy'+str(inb)] - mi
                ri_ = ( datax['Pxy'+str(inb)] - mi ) / np.max(datax['Pxy'+str(inb)])
                if inb0==0: ri = ri_
                if inb0!=0: ri = np.concatenate((ri, ri_))
        return ri
    if sum(listfile)>0:
        out = minimize(TotalFit, params)#, nan_policy='omit')#, kws={"wt": wt, "df": df})
        print('-------------------------------')
        print('Parameter    Value       Stderr')
        for name, param in out.params.items():
            print(f'{name:7s} {param.value:11.5f} {param.stderr:11.5f}')
     #   D=out.params['D'].value; dD=out.params['D'].stderr
        R = out.params['R'].value; dR = out.params['R'].stderr
        for inb in range(npower):
            if listfile[inb] == 1:
                listk += [out.params['k'+str(inb)].value]
                listdk += [out.params['k'+str(inb)].stderr]
                D = kBT/(6*np.pi*1.e-9*R * FaxenXY(R, listL[inb]) ); dD = dR * D / R
       #         print(R, listL[inb], D)    
                listD += [D]; listdD += [dD]
                listR += [R]; listdR += [dR]
            else:
                listk += [np.nan]; listdk += [np.nan]
                listD += [np.nan]; listdD += [np.nan]   
                listR += [np.nan]; listdR += [np.nan]
    else:
        listD=[np.nan]*npower; listk=[np.nan]*npower; listdD=[np.nan]*npower; listdk=[np.nan]*npower; listL=[np.nan]*npower
        listR=[np.nan]*npower; listdR=[np.nan]*npower
    return listD, listk, listdD, listdk, listL, listR, listdR

def CalculatePSDBins(f, Pxx_spec):       
    nbins=101; fbins=np.logspace(-2,3,nbins); Pbins=np.ones(nbins); dPbins=np.zeros(nbins)
    for m in range(nbins-1):
        u=Pxx_spec[(f>=fbins[m])&(f<fbins[m+1])]
        Pbins[m]=np.mean(u)
        dPbins[m]=np.std(u)/np.sqrt(len((u[~np.isnan(u)])))
    return fbins, Pbins, dPbins

def ReadFromdf_Multisurvival(dfstage, nrow, ncol, nrow2, ncol2, figname, pprint=False, minsize=1, mindur=0.5, mindz=50, maxdz=500):
    fig, axall = plt.subplots(nrow, ncol, num=figname, figsize=(ncol*sizefig, nrow*sizefig), dpi=100)
    fig1, axall1 = plt.subplots(nrow, ncol, num='multispectrumdzjump', figsize=(ncol*sizefig, nrow*sizefig), dpi=100)
    fig2, axall2 = plt.subplots(nrow2, ncol2, num=figname+'Bis', figsize=(ncol*sizefig, nrow*sizefig), dpi=100)
    fig2A, axall2A = plt.subplots(nrow2, ncol2, num=figname+'BisA', figsize=(ncol*sizefig, nrow*sizefig), dpi=100)
    fig2B, axall2B = plt.subplots(nrow2, ncol2, num=figname+'BisB', figsize=(ncol*sizefig, nrow*sizefig), dpi=100)
    fig3x, axall3x = plt.subplots(nrow2, ncol2, num='multispectrumx', figsize=(ncol*sizefig, nrow*sizefig), dpi=100)
    fig3y, axall3y = plt.subplots(nrow2, ncol2, num='multispectrumy', figsize=(ncol*sizefig, nrow*sizefig), dpi=100)
    figSpectestname = 'Spectrum test'; fig = plt.figure(figSpectestname, figsize=(6,6), dpi=100); ax3t = plt.gca()

    n = dfstage['TdmsList']
    n_calib = dfstage['TdmsList_calib']
    woff = []; wdoff = []; wkmin=[]; wkmax=[]; wkminmax=[]
    woff2A = []; wdoff2A = []; woff2B = []; wdoff2B = []
    wDeltaZ=[]; wdDeltaz=[]; wangle=[]; wangleC=[]; wFCyclex=[]; wFCycley=[]; wDCyclex=[]; wDCycley=[]
    wkx_c=[]; wDx_c=[]; wky_c=[]; wDy_c=[]; wCp1Lo_c=[]; wp1Lo_c=[]
    wFSDx=[]; wFSDy=[]; wFSDx_c=[]; wFSDy_c=[]
    columnsfitspec = ['MolecID', 'power', 'fittype', 'GuessDTheo', 'fmin_Hz', 'fmax_Hz',
                      'kxraw', 'Dxraw', 'dkxraw', 'dDxraw',
                      'kyraw', 'Dyraw', 'dkyraw', 'dDyraw',
                      'kx', 'Dx', 'dkx', 'dDx',
                      'ky', 'Dy', 'dky', 'dDy',
                      'Cp1Lo_D', 'p1Lo_D', 'Cp2Lo_D', 'p2Lo_D', 'PullAngle_D', 'PullAngleC_D', 'PullPhi_D', 'PullPhiC_D']
    dffitspeclist = []# = pd.DataFrame(np.zeros(5), columns=columnsfitspec)
    columnsoffrate = ['newkoff', 'newdkoff', 'kmin', 'kmax', 'kminmax',
                      'newkoff2A', 'newdkoff2A', 'newkoff2B', 'newdkoff2B',
                      'DeltaZ', 'dDeltaZ', 'Angle', 'AngleC', 'kx_c', 'Dx_c', 'ky_c', 'Dy_c',
                      'Cp1Lo_c', 'p1Lo_c', 'Fcyclex', 'Fcycley', 'Dcyclex', 'Dcycley', 'FSDx',
                      'FSDy', 'FSDx_c', 'FSDy_c', 'p1Lo_nm', 'p2Lo_nm']
    dfoffratelist = []
    count=0
    files = [i for i in os.listdir(inputpath0 + folder_spectra)if os.path.isfile(os.path.join(inputpath0 + folder_spectra,i))]
    for inb, nt in enumerate(n.values):
        print('########################################################################')
        MolecID = dfstage['Molecule ID'][inb]
        BeadRad_i = dfstage['Bead radius (nm)'][inb]
        friction0 = 6*np.pi*1.e-9*BeadRad_i
        print('Molecule ID=', MolecID, 'Power(%)=', dfstage['% High power'][inb], ' Radius=', BeadRad_i, ' Cycle file',nt)
        print(' Calib file', n_calib[inb])

        if (SelectMolecule == False) or (MolecID in ListMolec):
            iID = np.where(ListMolec==MolecID)[0][0]
            col = colormap[iID]
            prefixe = n_calib[inb][0]+'_Bead'+str(int(dfstage['Trace calib'][inb]))
            if rapa and MolecID>=MolecIDa_fmin_Hz_set and MolecID<=MolecIDb_fmin_Hz_set:
                prefixe = n[inb][0]+'_Bead'+str(int(dfstage['Trace'][inb]))
            if rapa:
                power = int(dfstage['% High power'][inb])
                suffixe = str(power)+'%.csv'
            else:
                power = dfstage['% High power'][inb]
                suffixe = f'{power:g}'+'%.csv'
            fname_spectrarawXY = prefixe +'_rawXYspectrum_'+ suffixe
            start_fname_spectraX = prefixe +'_ExportSpectrum_XsopenRef'
            start_fname_spectraY = prefixe +'_ExportSpectrum_YsopenRef'
            fname_spectraresults = prefixe +'_FitsResults_'+ suffixe
            print(' Spectrum file', fname_spectrarawXY)
            ax3x = axall3x[iID//ncol2, iID%ncol2]
            ax3y = axall3y[iID//ncol2, iID%ncol2]
            fullfilename = inputpath0 + folder_spectra + fname_spectrarawXY
            IDshow=13
            if os.path.isfile(fullfilename):
                fnameXa = inputpath0 + folder_spectra + start_fname_spectraX + "custom_" + suffixe
                fnameYa = inputpath0 + folder_spectra + start_fname_spectraY + "custom_" + suffixe
                fnameXb = inputpath0 + folder_spectra + start_fname_spectraX + "0_force_" + suffixe
                fnameYb = inputpath0 + folder_spectra + start_fname_spectraY + "0_force_" + suffixe
                if os.path.isfile(fnameXa):
                    dfx = pd.read_csv(fnameXa, index_col=False); dfy = pd.read_csv(fnameYa, index_col=False)
                if os.path.isfile(fnameXb):
                    dfx = pd.read_csv(fnameXb, index_col=False); dfy = pd.read_csv(fnameYb, index_col=False) 
                dfraw = pd.read_csv(inputpath0 + folder_spectra + fname_spectrarawXY, index_col=False) #print(fname_spectrarawXY)
                fnameres = inputpath0 + folder_spectra + fname_spectraresults
                if os.path.isfile(fnameres):
                    dfres = pd.read_csv(fnameres, index_col=False)
                    if 'p1Lo' not in dfres.columns: dfres['p1Lo'] = dfres.Cp1Lo;# print('dfres.p1Lo = dfres.Cp1Lo')
                    if 'p1Lo_nm' in dfres.columns: dfres['p1Lo'] = dfres.p1Lo_nm;# print('dfres.p1Lo = dfres.p1Lo_nm')
                    if 'p2Lo' not in dfres.columns: dfres['p2Lo'] = dfres.Cp2Lo;# print('dfres.p1Lo = dfres.Cp1Lo')
                    if 'p2Lo_nm' in dfres.columns: dfres['p2Lo'] = dfres.p2Lo_nm;# print('dfres.p1Lo = dfres.p1Lo_nm')
                    if 'PullAngle_deg' in dfres.columns: dfres['PullAngle'] = dfres.PullAngle_deg;# print('dfres.PullAngle = dfres.PullAngle_deg')
                    if 'PullPhi' not in dfres.columns: dfres['PullPhi'] = dfres.PullPhiC;# print('dfres.PullPhi = dfres.PullPhiC')
                else:
         #           dfres = pd.DataFrame(data={'Cp1Lo_nm': [0], 'p1Lo_nm': [0], 'Cp2Lo_nm': [0], 'p2Lo': [0], 'PullAngle': [0], 'PullAngleC': [0], 'PullPhi': [0], 'PullPhiC': [0]} )
                    dfres = pd.DataFrame(data={'Cp1Lo': [0], 'p1Lo': [0], 'Cp2Lo': [0], 'p2Lo': [0], 'PullAngle': [0], 'PullAngleC': [0], 'PullPhi': [0], 'PullPhiC': [0]} )
                    print(fnameres, ' does not exist')
                fmin_Hz = np.min(dfraw['fx.1'][dfraw['fx.1']>0]); fmax_Hz = np.max(dfraw['fx.1'][dfraw['fx.1']>0])                   
         #       print('fmin_Hz, fmax_Hz:' ,fmin_Hz, fmax_Hz)  #        fmin_Hz = 0.1; fmax_Hz = 15
                if power>=25: fmin_Hz=1
                if rapa and MolecID>=MolecIDa_fmin_Hz_set and MolecID<=MolecIDb_fmin_Hz_set:
                    fmin_Hz=fmin_Hz_set
                    if power>=20: fmin_Hz=0.5
                    if power>=25: fmin_Hz=1
                    if power>=30: fmin_Hz=5
          #      if not rapa and fmin_Hz<1: fmin_Hz=0.6 #fmin_Hz=0.3
                if not rapa: fmin_Hz=0.5 #fmin_Hz=0.3
                if not rapa and fmax_Hz>100: fmax_Hz=150 #fmin_Hz=0.3
                print(f' fmin_Hz= {fmin_Hz:3.2f} fmax_Hz={fmax_Hz:3.2f}')
                fr = dfraw['fx']; intfr = (fr>=fmin_Hz)&(fr<=fmax_Hz); frt=fr[intfr]
                fbinsX, PbinsX, dPbinsX = CalculatePSDBins(dfraw['fx'], dfraw['Pxx_specx'])
                fbinsY, PbinsY, dPbinsY = CalculatePSDBins(dfraw['fx'], dfraw['Pxx_specy'])
                ax3x.loglog(fbinsX, PbinsX, 'k-', marker=next(poolmarker2x), ms=6 , c=col, label=str(power)+'%', alpha=0.3)
                ax3y.loglog(fbinsY, PbinsY, 'k-', marker=next(poolmarker2y), ms=6 , c=col, label=str(power)+'%', alpha=0.3)
                fx = dfx['Xfbins_Hz']; intfx = (fx>=fmin_Hz)&(fx<=fmax_Hz); fxt=fx[intfx]
                fy = dfy['Xfbins_Hz']; intfy = (fy>=fmin_Hz)&(fy<=fmax_Hz); fyt=fy[intfy]
                                 
                if MolecID==IDshow:
                    ax3t.loglog(fx, dfx['Pbins'], 'k.', c=colormap[iID+1], label='x '+str(power)+'%')
                if ( (inb==0) or ( (inb>0) and ( MolecID > dfstage['Molecule ID'][inb-1] ) ) ):
                    GuessDTheo = 0; Dx=0; frictionXY=0
                GuessDTheo = 0; Dx=0; frictionXY=0      # 20220415: caution introduced to avoid pb with empty powers
               
      #          for fittype, color in zip([0],['b']): # zip([0,1,2],['b','r','g']): #
                p1Zo = 900; GuessDTheo = 150000; Dx=0; frictionXY=friction0 * FaxenXY(BeadRad_i, p1Zo); # fmin_Hz=0.1; fmax_Hz=15
    #            if (fittype == 0) and ( (inb==0) or ( (inb>0) and ( MolecID > dfstage['Molecule ID'][inb-1] ) ) ):
    #                frictionXY = friction0 * FaxenXY(BeadRad_i, p1Zo)
    #                GuessDTheo = Dx; print('Ratio GuessDTheo/Dtheo=', GuessDTheo * frictionXY / kBT_pN_nm)
                print('frictionXY=', frictionXY, 'GuessDTheo =', GuessDTheo)
                kxr, Dxr, dkxr, dDxr, fc_xr, FitPxx_specxr = QuickSpectrumFit(frt, dfraw['Pxx_specx'][intfr], 'XY', GuessDTheo = GuessDTheo)                        
                kyr, Dyr, dkyr, dDyr, fc_yr, FitPxx_specyr = QuickSpectrumFit(frt, dfraw['Pxx_specy'][intfr], 'XY', GuessDTheo = GuessDTheo)                       
                kx, Dx, dkx, dDx, fc_x, FitPxx_specx = QuickSpectrumFit(fxt, dfx['Pbins'][intfx], 'XY', GuessDTheo = GuessDTheo)
                ky, Dy, dky, dDy, fc_y, FitPxx_specy = QuickSpectrumFit(fyt, dfy['Pbins'][intfy], 'XY', GuessDTheo = GuessDTheo)
                if MolecID==IDshow: ax3t.loglog(fxt, FitPxx_specx, 'b-', label='x fit', alpha=0.5)
                kxg = dfstage['kxglobal'][inb]; Dxg = dfstage['Dxglobal'][inb]
                kyg = dfstage['kyglobal'][inb]; Dyg = dfstage['Dyglobal'][inb]
                fcg = (kxg+kyg)*( Dxg+Dyg )/4/(2*np.pi*kBT_pN_nm)
                if fcg < fmax_Hz/2:
                    print('show individual fit: ',fcg,' < ', fmax_Hz/2)
                    ax3x.loglog(fxt, FitPxx_specx, 'r-')
                    ax3y.loglog(fyt, FitPxx_specy, 'r-')
                FitPxx_specxrg = FitSpectrumGen(frt, kxg, Dxg)
                FitPxx_specyrg = FitSpectrumGen(frt, kyg, Dyg)
                ax3x.loglog(frt, FitPxx_specxrg, 'k-') 
                ax3y.loglog(frt, FitPxx_specyrg, 'k-') 
                fitspec = [MolecID, power, 0, GuessDTheo, fmin_Hz, fmax_Hz,
                           kxr, Dxr, dkxr, dDxr, kyr, Dyr, dkyr, dDyr,
                           kx, Dx, dkx, dDx, ky, Dy, dky, dDy,
                           dfres['Cp1Lo'][0], dfres['p1Lo'][0], dfres['Cp2Lo'][0], dfres['p2Lo'][0],
                           dfres['PullAngle'][0], dfres['PullAngleC'][0], dfres['PullPhi'][0], dfres['PullPhiC'][0]]
                dffitspeclist.append(fitspec)
                    
        #     print(MolecID, iID)
            wdur = np.zeros(0); wdz = np.zeros(0)
            if len(n_calib[inb])>0 and rapa:
                if ~np.isnan(dfstage['Trace calib'][inb]):
                    fname_calib = n_calib[inb][0]+'_FitsResults_Bead'+str(int(dfstage['Trace calib'][inb]))+'_'+str(int(dfstage['% High power'][inb]))+'%.csv'
                    innt_calib = df0_calib['FileName'] == fname_calib
                    print('fname_calib ',fname_calib,':', innt_calib.sum(), ' file')
                    if innt_calib.sum()>0:
                        kx_c = df0_calib['kx_pN_per_nm'][innt_calib].values[0]
                        Dx_c = df0_calib['Dx_nm2_per_s'][innt_calib].values[0]
                        ky_c = df0_calib['ky_pN_per_nm'][innt_calib].values[0]
                        Dy_c = df0_calib['Dy_nm2_per_s'][innt_calib].values[0]
                        Cp1Lo_c = df0_calib['Cp1Lo'][innt_calib].values[0]
                        p1Lo_c = df0_calib['p1Lo'][innt_calib].values[0]
                        FSDx_c =  df0_calib['FSDXFit_pN'][innt_calib].values
                        FSDy_c =  df0_calib['FSDYFit_pN'][innt_calib].values
                        print('kx_c = ', kx_c, ' ky_c = ', ky_c)
                    else:
                        print('No file'); kx_c = np.nan; Dx_c = np.nan; ky_c = np.nan; Dy_c = np.nan; Cp1Lo_c = np.nan; p1Lo_c = np.nan
                else:
                    print('No file'); kx_c = np.nan; Dx_c = np.nan; ky_c = np.nan; Dy_c = np.nan; Cp1Lo_c = np.nan; p1Lo_c = np.nan
            else:
                kx_c = np.nan; Dx_c = np.nan; ky_c = np.nan; Dy_c = np.nan; Cp1Lo_c = np.nan; p1Lo_c = np.nan; FSDx_c=np.nan; FSDy_c=np.nan
            for nnt in nt:
                fname = nnt+'_FitsResults_Bead'+str(dfstage['Trace'][inb])+'.csv'       
                innt = df0['FileName'] == fname
           #     timehighpower = dfstage['Time at hight power (s)']
                if sum(innt)==1:
                    print('====== single trace csv file ======', fname, sum(innt))
                if sum(innt)==0:
                    innt = (df0['FileName'] == nnt+'_FitsResults.csv') & (df0['bead'] == dfstage['Trace'][inb])
                    print('====== multiple trace csv file ====== ', nnt+'_FitsResults.csv', sum(innt), 'Bead ', dfstage['Trace'][inb])
         #       u = df0['Tupjump_Tuphigh_mask_cycle'][innt].values
                u = df0['Tupjump_Tuphigh'][innt].values
                uz = df0['Zupjump_Zupmid'][innt].values
                v = df0['FileName'][innt].values
          #      print('nb dur= ', len(u));# print(uz); print(v); print(fname)
                angleC = df0['PullAngleC'][innt].values
                angleC = angleC[0] if len(angleC)>0 else np.nan
                angle = df0['PullAngle_deg'][innt].values
                angle = angle[0] if len(angle)>0 else np.nan
                p1Lo_nm = df0['p1Lo_nm'][innt].values
                p1Lo_nm = p1Lo_nm[0] if len(p1Lo_nm)>0 else np.nan
                p2Lo_nm = df0['p2Lo_nm'][innt].values
                p2Lo_nm = p2Lo_nm[0] if len(p2Lo_nm)>0 else np.nan
                Fcyclex = df0['FSpecx_pN'][innt].values
                Fcyclex = Fcyclex[0] if len(Fcyclex)>0 else np.nan
                Fcycley = df0['FSpecy_pN'][innt].values
                Fcycley = Fcycley[0] if len(Fcycley)>0 else np.nan  
                Dcyclex = df0['Dx_nm2_per_s'][innt].values
                Dcyclex = Dcyclex[0] if len(Dcyclex)>0 else np.nan
                Dcycley = df0['Dy_nm2_per_s'][innt].values
                Dcycley = Dcycley[0] if len(Dcycley)>0 else np.nan
                if MolecID==1 or MolecID==2 or MolecID==3:
                    if n_calib[inb][0]=='20200824-201201':
                        print('No Correction D cycle:'); Dcyclex = Dcyclex*2; Dcyclex = Dcyclex*2
                FSDx =  df0['FSDXFit_pN'][innt].values
                FSDx = FSDx[0] if len(FSDx)>0 else np.nan
                FSDy =  df0['FSDYFit_pN'][innt].values
                FSDy = FSDy[0] if len(FSDy)>0 else np.nan                
                print('#########',v)
                print('angle:', angle, 'angleC:', angleC)
                print('Fcyclex:', Fcyclex, 'Fcycley:', Fcycley)
                if len(u)>0:
                    dur = np.array(g(u[0]))
                    dz = np.array(g(uz[0]))
                    print('len(dur)=', len(dur), 'len(dz)=', len(dz))                 
                    ipos = dur>=mindur; ineg = dur<mindur
                    imed = (dz>=mindz) & (dz<=maxdz) ; ilow = dz<mindz; ihigh =dz>maxdz
                    durpos = dur[ipos]; durneg = dur[ineg]; dzpos = dz[ipos]
                    dzhigh = dz[ihigh]; dzlow = dz[ilow]
                    isel = ipos & imed; inosel = ineg | ilow
                    dursel = dur[isel]; dzsel = dz[isel]
                    wdur = np.append(durpos, wdur); wdz = np.append(dzpos, wdz)
                   
                    TEhigh = np.max(wdur)*1000 if len(wdur)>0 else 0   #      TEhigh = df0['TEhigh_ms'][innt].values[0]
                    countNoOpen = df0['countNoOpen'][innt].values[0]
                    countNoClose = df0['countNoClose'][innt].values[0]
                    print('len(wdur)=', len(wdur), 'len(wdz)=', len(wdz), 'TEhigh=', TEhigh, 'countNoOpen=', countNoOpen, 'countNoClose=', countNoClose)
            print(count, count//ncol, count%ncol)
           # wdur = wdur[wdur>mindur]
            ndur = len(wdur)
            wdurA = wdur[0:int(ndur/2)]
            wdurB = wdur[int(ndur/2):ndur]
            for ax3i in [ax3x, ax3y, ax3t]:
                ax3i.set_xlabel("Frequency (Hz)");  ax3i.set_ylabel("PSD"); ax3i.axis([0.05, 300, 0.1, 1e4])
                if iID%ncol2!=0: ax3i.axes.get_yaxis().set_visible(False)
                if iID//ncol2!=nrow2-1: ax3i.axes.get_xaxis().set_visible(False)
            ax3x.legend(title='ID '+str(MolecID)+' x', title_fontsize=5, prop={'size': 5}, fontsize=5, loc="upper right")
            ax3y.legend(title='ID '+str(MolecID)+' y', title_fontsize=5, prop={'size': 5}, fontsize=5, loc="upper right")
            if len(wdur) >= minsize:
                ax = axall[count//ncol, count%ncol]
                ax1 = axall1[count//ncol, count%ncol]
                ax2 = axall2[iID//ncol2, iID%ncol2]
                ax2A = axall2A[iID//ncol2, iID%ncol2]
                ax2B = axall2B[iID//ncol2, iID%ncol2]
                for axi in [ax, ax1]:
                    if count%ncol!=0: axi.axes.get_yaxis().set_visible(False)
                    if count//ncol!=nrow-1: axi.axes.get_xaxis().set_visible(False)
                for ax2i in [ax2, ax2A, ax2B]:
                    if iID%ncol2!=0: ax2i.axes.get_yaxis().set_visible(False)
                    if iID//ncol2!=nrow2-1: ax2i.axes.get_xaxis().set_visible(False)
                if rapa: tmax = min(tmaxsurvival*1000, TEhigh)
                else: tmax = TEhigh
                shortlegend = str(power)+'%'

                ax1.scatter(wdur, wdz, marker='.', c=col, label=shortlegend)
                ax1.axis([-20, 200, mindz, maxdz])
                ax1.set_xlabel("Time (s)");  ax1.set_ylabel("Delta Z jump (nm)")

              #  print(shortlegend)
                (pEq0, dpEq0, kmin, kmax, x2, y2, y2_fit0) = survival(wdur, tmax, countNoOpen, 0, '', mindur=mindur, ax=ax , col=col, shortlegend=shortlegend)
                (pEq0, dpEq0, kmin, kmax, x2, y2, y2_fit0) = survival(wdur, tmax, countNoOpen, 0, '', mindur=mindur, ax=ax2, col=col, shortlegend=shortlegend)
                koff = 1/pEq0; dkoff = dpEq0/pEq0**2; count+=1
                if rapa:
                    (pEq0, dpEq0, kmin, kmax, x2, y2, y2_fit0) = survival(wdurA, tmax, countNoOpen, 0, '', mindur=mindur, ax=ax2A , col=col, shortlegend=shortlegend)
                    koff2A = 1/pEq0; dkoff2A = dpEq0/pEq0**2
                else:
                    koff2A = np.nan; dkoff2A = np.nan
                if rapa:
                    (pEq0, dpEq0, kmin, kmax, x2, y2, y2_fit0) = survival(wdurB, tmax, countNoOpen, 0, '', mindur=mindur, ax=ax2B , col=col, shortlegend=shortlegend)
                    koff2B = 1/pEq0; dkoff2B = dpEq0/pEq0**2
                else:
                    koff2B = np.nan; dkoff2B = np.nan
                if kmin == 0: kmin = kmax/len(wdur)
                for ax2i in [ax, ax1, ax2, ax2A, ax2B]:
                    if not ( (ax2i==ax3t) & (MolecID!=IDshow) ):
                        ax2i.legend(title='ID '+str(MolecID), title_fontsize=5, prop={'size': 5}, fontsize=5, loc="upper right")
                DeltaZ = np.mean(np.array(wdz)); dDeltaZ = np.std(np.array(wdz))
            else: 
                print('insufficient or no data')
                koff2A = np.nan; dkoff2A = np.nan
                koff2B = np.nan; dkoff2B = np.nan
                koff = np.nan; dkoff = np.nan
                kmin = np.nan; kmax = np.nan
                DeltaZ = np.nan; dDeltaZ = np.nan
            if len(wdur) >= minsize and np.min(wdur) > 0.7*np.max(wdur):
                print('constant max duration')
                koff = 1/np.max(wdur)/len(wdur); dkoff = koff
                koff = np.nan; dkoff = np.nan
            #    ax.legend(title='toff='+str(int(1/koff)), title_fontsize=6, prop={'size': 6})
        else:
            print('Non selected Molecule ID')
            koff = np.nan; dkoff = np.nan
            koff2A = np.nan; dkoff2A = np.nan; 
            koff2B = np.nan; dkoff2B = np.nan; 
            kmin = np.nan; kmax = np.nan
            DeltaZ = np.nan; dDeltaZ=np.nan
            angle = np.nan; angleC = np.nan
            kx_c = np.nan; Dx_c = np.nan; ky_c = np.nan; Dy_c = np.nan; Cp1Lo_c = np.nan
        
        offrate = [koff, dkoff, kmin, kmax, (kmin*kmax)**0.5, koff2A, dkoff2A, koff2B, dkoff2B,
                   DeltaZ, dDeltaZ, angle, angleC, kx_c, Dx_c, ky_c, Dy_c, Cp1Lo_c, p1Lo_c,
                   Fcyclex, Fcycley, Dcyclex, Dcycley, FSDx, FSDy, FSDx_c, FSDy_c, p1Lo_nm, p2Lo_nm]
        dfoffratelist.append(offrate)
        
        print('toff=', 1/koff, 'koff=', koff, 'dkoff=', dkoff)
        
    dfoffrate = pd.DataFrame(dfoffratelist, columns = columnsoffrate)
    dffitspec = pd.DataFrame(dffitspeclist, columns = columnsfitspec)
   # palette=sns.color_palette("Set1", dffitspec.fittype.nunique())
    figfitspecname, axfitspecall = plt.subplots(2,2, num='Spectrum fit All', figsize=(2*2*sizefig, 2*2*sizefig), dpi=100)
    for ifit, icol in zip([0,1,2], ['r','b','g']):
        dfloc = dffitspec[ (dffitspec['fittype']==ifit) & (dffitspec['MolecID']==IDshow) ]
  #      axfitspecall[0,0].errorbar(dfloc["power"], dfloc['Dxraw/Dtheo'], dfloc['dDxraw/Dtheo'], fmt=icol+'o-', markerfacecolor='none', label=str(ifit)+' raw', alpha=0.5)
  #      axfitspecall[0,0].errorbar(dfloc["power"], dfloc['Dx/Dtheo'], dfloc['dDx/Dtheo'], fmt=icol+'s-', label=str(ifit)+' bin', alpha=0.5)
        axfitspecall[1,0].errorbar(dfloc["power"], dfloc['kxraw'], dfloc['dkxraw'], fmt=icol+'o-', markerfacecolor='none', label=str(ifit)+' raw', alpha=0.5)
        axfitspecall[1,0].errorbar(dfloc["power"], dfloc['kx'], dfloc['dkx'], fmt=icol+'s-', label=str(ifit)+' bin', alpha=0.5)
  #      axfitspecall[0,1].errorbar(dfloc["power"], dfloc['Dyraw/Dtheo'], dfloc['dDyraw/Dtheo'], fmt=icol+'o-',  markerfacecolor='none', label=str(ifit)+' raw', alpha=0.5)
  #      axfitspecall[0,1].errorbar(dfloc["power"], dfloc['Dy/Dtheo'], dfloc['dDy/Dtheo'], fmt=icol+'s-', label=str(ifit)+' bin', alpha=0.5)
        axfitspecall[1,1].errorbar(dfloc["power"], dfloc['kyraw'], dfloc['dkyraw'], fmt=icol+'o-',  markerfacecolor='none', label=str(ifit)+' raw', alpha=0.5)
        axfitspecall[1,1].errorbar(dfloc["power"], dfloc['ky'], dfloc['dky'], fmt=icol+'s-', label=str(ifit)+' bin', alpha=0.5)   
    for axii in [ axfitspecall[0,0], axfitspecall[1,0], axfitspecall[0,1], axfitspecall[1,1] ]:
        axii.legend(loc='upper left', fontsize=7); axii.set_xlabel('Power (%)'); axii.set_xlim(0,31)
    for axii in [axfitspecall[0,0], axfitspecall[0,1] ]: axii.set_ylim(0,2)
 #   axfitspecall[0,0].set_ylabel('Dx/Dtheo'); axfitspecall[0,1].set_ylabel('Dy/Dtheo')
    for axii in [axfitspecall[1,0], axfitspecall[1,1] ]: axii.set_ylim(0,0.01)
    axfitspecall[1,0].set_ylabel('kx (pN/nm)'); axfitspecall[1,1].set_ylabel('ky (pN/nm)')
    if SaveGraph: plt.figure('Spectrum fit All');plt.tight_layout(); plt.savefig(outpath+'Spectrum fit All', transparent=False)
    
    if SaveGraph:
        plt.figure(figname); plt.savefig(outpath+figname, transparent=False)
        plt.figure('multispectrumdzjump'); plt.savefig(outpath+'multispectrumdzjump', transparent=False)
        plt.figure(figname+'Bis'); plt.savefig(outpath+figname+'Bis', transparent=False)
        plt.figure('multispectrumx'); plt.savefig(outpath+'multispectrumx', transparent=False)
        plt.figure('multispectrumy'); plt.savefig(outpath+'multispectrumy', transparent=False)
        plt.figure('Spectrum test'); plt.savefig(outpath+'Spectrum test', transparent=False)
  #  return dffitspec, woff, wdoff, wkmin, wkmax, wkminmax, woff2A, wdoff2A, woff2B, wdoff2B, wDeltaZ, wdDeltaz, wangle, wangleC, wkx_c, wDx_c, wky_c, wDy_c, wCp1Lo_c, wp1Lo_c, wFCyclex, wFCycley, wDCyclex, wDCycley, wFSDx, wFSDy, wFSDx_c, wFSDy_c
    return dfoffrate, dffitspec

if not SortingBeads: df0['beadrenum'] = df0.index

def Bell(wf, k0, f0): return k0*np.exp(wf/f0)
def BellxB(wf, k0, xB): return k0*np.exp(wf*xB/kBT)
def logBell(wf, k0, f0): return np.log(k0) + wf/f0
def logBellxB(wf, k0, xB): return np.log(k0) + wf*xB/kBT
def Friddle(wf, k0, xB, kc): return k0*np.exp( wf*xB/kBT - 0.5*kc*xB**2/kBT )
def logFriddle(wf, k0, xB, kc): return np.log(k0) + wf*xB/kBT - 0.5*kc*xB**2/kBT

def Affine(x, *p): return (x - p[1])/p[0]
def AffineL(x, *p): return x/p[1] + p[0]
#def WLC0(x, *p): return p[0] * ( 0.25*(1-x/p[1])**-2 - 0.25 + x/p[1] ) # force vs length x; p[1] = L0; p[0] = kT/Lp

def WLC(x, *p):     # https://en.wikipedia.org/wiki/Worm-like_chain#cite_note-Petrosyan-12
    f = (kBT/p[0]) * ( 0.25*(1-x/p[1])**-2 - 0.25 + x/p[1] -0.8*(x/p[1])**2.15 ) # force vs length x; p[1] = L0; p[0] = Lp
   # f = (kBT/p[0]) * ( 0.25*(1-x/p[1])**-2 - 0.25 + x/p[1] ) # force vs length x; p[1] = L0; p[0] = Lp
    return f*(1-x/p[1]>0)     

def fitBellFriddle(func, wf, wk, weights=None):
    Hmodel = Model(func)  # Hmodel = Model(Bell)
    Hparams = Hmodel.make_params()
    xBinit = 1.; k0init = 1e-4
    if rapa: xBinit = 0.5; k0init = 1e-2
    Hparams.add('xB', value=xBinit, vary=True, min=1e-6)
    Hparams.add('k0', value=k0init, vary=True, min=1e-6)
    if func == logFriddle: Hparams.add('kc', value=kc, vary=False)
    while True:
        try:
            Hresult = Hmodel.fit(np.log(wk), Hparams, wf=wf, weights = weights)
            eEq=np.sqrt(np.diag(Hresult.covar))
            res=Hresult.redchi     #        res=Hresult.chisqr 
            pEq=[Hresult.best_values['k0'], Hresult.best_values['xB']]
            if func == logFriddle: pEq=[Hresult.best_values['k0'], Hresult.best_values['xB'], Hresult.best_values['kc']]
            break
        except ValueError:
            print("Non finite value");
            pEq=[0.,1., kc]; res=1
            eEq = [np.nan, np.nan] 
            break
        except TypeError:
            print("Improper input");
            pEq=[0.,1., kc ]; res=1
            eEq = [np.nan, np.nan] 
            break
    if len(wk)<=2:  #   print('2 data points only; set large error bars')
        eEq = [0.,1., kc]
    return (pEq,eEq,res)  

figname = 'offrate_vs_force'; fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
figname2 = 'koff0 vs BellForce'; fig = plt.figure(figname2, figsize=(6,6), dpi=100); ax2 = plt.gca()
figname2Bis = 'koff0_vs_xBeta'; fig = plt.figure(figname2Bis, figsize=(6,6), dpi=100); ax2Bis = plt.gca()
figname2Ter = 'kx vs xBeta'; fig = plt.figure(figname2Ter, figsize=(6,6), dpi=100); ax2Ter = plt.gca()
figname2Ter1 = 'Length vs xbeta'; fig = plt.figure(figname2Ter1, figsize=(6,6), dpi=100); ax2Ter1 = plt.gca()
figname2Ter2 = 'kz vs Bell Force'; fig = plt.figure(figname2Ter2, figsize=(6,6), dpi=100); ax2Ter2 = plt.gca()
figname2Ter3 = 'kz vs kx'; fig = plt.figure(figname2Ter3, figsize=(6,6), dpi=100); ax2Ter3 = plt.gca()
figname2Ter6 = 'Lp_vs_L0'; fig = plt.figure(figname2Ter6, figsize=(6,6), dpi=100); ax2Ter6 = plt.gca()
figname2Ter4 = 'kz vs Lp'; fig = plt.figure(figname2Ter4, figsize=(6,6), dpi=100); ax2Ter4 = plt.gca()
figname2Ter5 = 'Lp vs Bell Force'; fig = plt.figure(figname2Ter5, figsize=(6,6), dpi=100); ax2Ter5 = plt.gca()
figname2Ter7 = 'DZ vs Length'; fig = plt.figure(figname2Ter7, figsize=(6,6), dpi=100); ax2Ter7 = plt.gca()
figname3 = 'Lumicks vs Python force'; fig = plt.figure(figname3, figsize=(6,6), dpi=100); ax3 = plt.gca()

figsize = (ncol*sizefig, nrow*sizefig)
figBis, axallBis = plt.subplots(nrow2, ncol2, num=figname+'Bis', figsize=figsize, dpi=100)
figBisAB, axallBisAB = plt.subplots(nrow2, ncol2, num=figname+'BisAB', figsize=figsize, dpi=100)
figLength, axallLength = plt.subplots(nrow2, ncol2, num='Force Lengths', figsize=figsize, dpi=100)
figPLength, axallPLength = plt.subplots(nrow2, ncol2, num='Power_Lengths', figsize=figsize, dpi=100)
figLengthBis, axallLengthBis = plt.subplots(nrow2, ncol2, num='Force_Lengths_Bis', figsize=figsize, dpi=100)
figAngle, axallAngle = plt.subplots(nrow2, ncol2, num='Force Angles', figsize=figsize, dpi=100)
figPAngle, axallPAngle = plt.subplots(nrow2, ncol2, num='Power_Angles', figsize=figsize, dpi=100)
figPPlateau, axallPPlateau = plt.subplots(nrow2, ncol2, num='Power Plateau', figsize=figsize, dpi=100)
figkxky, axallkxky = plt.subplots(nrow2, ncol2, num='Power Stiffness', figsize=figsize, dpi=100)
figkxkyBis, axallkxkyBis = plt.subplots(nrow2, ncol2, num='Power Stiffness Bis', figsize=figsize, dpi=100)
figkxkyTer, axallkxkyTer = plt.subplots(nrow2, ncol2, num='Power_Stiffness_Ter', figsize=figsize, dpi=100)
figDxDy, axallDxDy = plt.subplots(nrow2, ncol2, num='Power Friction', figsize=figsize, dpi=100)
figDxDyBis, axallDxDyBis = plt.subplots(nrow2, ncol2, num='Power_Friction_Bis', figsize=figsize, dpi=100)
figNoClose, axallNoClose = plt.subplots(nrow2, ncol2, num='No_Close_Fraction', figsize=figsize, dpi=100)
figFLumicks, axallFLumicks = plt.subplots(nrow2, ncol2, num='Force Lumicks', figsize=figsize, dpi=100)
figPower, axallPower = plt.subplots(nrow2, ncol2, num='Power vs Force', figsize=figsize, dpi=100)
figPower, axallPowerBis = plt.subplots(nrow2, ncol2, num='Power_vs_Force_Bis', figsize=figsize, dpi=100)
figDeltaz, axallDeltaz = plt.subplots(nrow2, ncol2, num='Deltaz vs Force', figsize=figsize, dpi=100)
figLengthDeltaZ, axallLengthDeltaZ = plt.subplots(nrow2, ncol2, num='Deltaz vs Length', figsize=figsize, dpi=100)
figMat = 'Maturation'; fig = plt.figure(figMat, figsize=(6,6), dpi=100); axMat = plt.gca()
fig4 = 'Force Range'; fig = plt.figure(fig4, figsize=(6,6), dpi=100); ax4 = plt.gca()
figContour = 'LengthContour'; fig = plt.figure(figContour, figsize=(6,6), dpi=100); axContour = plt.gca()
figOtheFlyDx = 'OtheFlyDx'; fig = plt.figure(figOtheFlyDx, figsize=(6,6), dpi=100); axOtheFlyDx = plt.gca()
figOtheFlyDy = 'OtheFlyDy'; fig = plt.figure(figOtheFlyDy, figsize=(6,6), dpi=100); axOtheFlyDy = plt.gca()
figOtheFlykx = 'OtheFlykx'; fig = plt.figure(figOtheFlykx, figsize=(6,6), dpi=100); axOtheFlykx = plt.gca()
figOtheFlyky = 'OtheFlyky'; fig = plt.figure(figOtheFlyky, figsize=(6,6), dpi=100); axOtheFlyky = plt.gca()
figOtheFlySDx = 'OtheFlySDx'; fig = plt.figure(figOtheFlySDx, figsize=(6,6), dpi=100); axOtheFlySDx = plt.gca()
figOtheFlySDy = 'OtheFlySDy'; fig = plt.figure(figOtheFlySDy, figsize=(6,6), dpi=100); axOtheFlySDy = plt.gca()
figXYk = 'XvsYcalibk'; fig = plt.figure(figXYk, figsize=(6,6), dpi=100); axXYk = plt.gca()
figXYD = 'XvsYcalibD'; fig = plt.figure(figXYD, figsize=(6,6), dpi=100); axXYD = plt.gca()
figXYF = 'XvsYcalibF'; fig = plt.figure(figXYF, figsize=(6,6), dpi=100); axXYF = plt.gca()
figRR = 'Rretrieval'; fig = plt.figure(figRR, figsize=(6,6), dpi=100); axRR = plt.gca()
figRR1 = 'Rretrieval1'; fig = plt.figure(figRR1, figsize=(6,6), dpi=100); axRR1 = plt.gca()

#ax3.errorbar(wfold, wfL, wdfL, wdfold, c=c, label='ID='+str(x), fmt='o', alpha=0.5)
ax.axhline(y=koffSPR, color='r', linestyle='-.', label='SPR')

def findR(wL, wD):
    wRin = np.arange(500,2500)
    wRout = np.zeros(len(wD))
    for iD, Di in enumerate(wD):
        wy = np.abs(kBT_pN_nm/FaxenXY(wRin, wL[iD])/(6*np.pi*1e-9*Di) - wRin)      
    #    wy = np.abs(kBT_pN_nm/FaxenXY(wRin, 1000)/(6*np.pi*1e-9*Di) - wRin)      
        wRout[iD] = wRin[wy.argmin()]
        if wRout[iD]==500: wRout[iD]=np.nan
    return wRout
    
for xy in ['x', 'y']:
    listDglobal = []; listkglobal = []; listdDglobal = []; listdkglobal = []; listLglobal=[]; listRglobal = []; listdRglobal = []
    for ii, MolecID in enumerate(np.arange(np.max(ListMolec))+1):# enumerate(ListMolec_rapa): #:
        listD, listk, listdD, listdk, listL, listR, listdR = GlobalSpecFit(dfstage0, folder_spectra, MolecID, xy)
        listkglobal.append(listk); listdkglobal.append(listdk)
        listDglobal.append(listD); listdDglobal.append(listdD)
        listRglobal.append(listR); listdRglobal.append(listdR)
        listLglobal.append(listL)
    flk = [val for sublist in listkglobal for val in sublist]; fldk = [val for sublist in listdkglobal for val in sublist]
    flD = [val for sublist in listDglobal for val in sublist]; fldD = [val for sublist in listdDglobal for val in sublist]
    flR = [val for sublist in listRglobal for val in sublist]; fldR = [val for sublist in listdRglobal for val in sublist]
    flL = [val for sublist in listLglobal for val in sublist]
    ll = len(dfstage0['Molecule ID']) - len(flk)
    dfstage0['k'+xy+'global'] = np.array(flk+[np.nan]*ll); dfstage0['dk'+xy+'global'] = np.array(fldk+[np.nan]*ll)
    dfstage0['D'+xy+'global'] = np.array(flD+[np.nan]*ll); dfstage0['dD'+xy+'global'] = np.array(fldD+[np.nan]*ll)
    if xy=='x': dfstage0['Lglobal'] = np.array(flL+[np.nan]*ll)
    dfstage0['R'+xy] = np.array(flR+[np.nan]*ll); dfstage0['dR'+xy] = np.array(fldR+[np.nan]*ll)
#    dfstage0['R'+xy] = kBT_pN_nm/FaxenXY(dfstage0['Bead radius (nm)'], p1Zo)/(6*np.pi*1e-9*dfstage0['D'+xy+'global'])    

dfoffrate, dffitspec =  ReadFromdf_Multisurvival(dfstage0, nrow, ncol, nrow2, ncol2, 'multisurvival2', minsize=minsize, mindur=mindur, mindz=mindz, maxdz=maxdz, pprint=True)

dfstage1 = pd.concat([ dfstage0.reset_index(), dfoffrate.reset_index()], axis=1)
dffitspec['Molecule ID'] = dffitspec['MolecID']
dffitspec['% High power'] = dffitspec['power']
if rapa:
    dfstage1['% High power'] = dfstage1['% High power'].astype(int)
else:
    dfstage1['% High power'] = dfstage1['% High power'].astype(float)
dffitspec_ = (dffitspec[dffitspec['fittype']==0]).reset_index(drop=True)
dfstage = pd.merge(dfstage1, dffitspec_, on = ["Molecule ID", "% High power"], how="outer", indicator=True )

dfstage['newdkoff2AB'] = np.abs(dfstage['newkoff2A'] - dfstage['newkoff2B'])
dfstage['Ffit'] = np.zeros(len(dfstage['newkoff'])) 

for xy in ['x', 'y']: dfstage['Rprec'+xy] = findR(dfstage['p1Lo_D'], dfstage['D'+xy+'global'])

wlogk0 = np.ndarray(Nmolec); wdlogk0 = np.ndarray(Nmolec)   # log koff at 0 force Bell Fit
wlogk0p = np.ndarray(Nmolec); wdlogk0p = np.ndarray(Nmolec) # log koff at 0 force Bell Fit projected force
wlogk0F = np.ndarray(Nmolec); wdlogk0F = np.ndarray(Nmolec) # log koff at 0 force Friddle Fit
wfB = np.ndarray(Nmolec); wdfB = np.ndarray(Nmolec)      # Bell force Bell Fit
wfBF = np.ndarray(Nmolec); wdfBF = np.ndarray(Nmolec)    # Bell force Friddle Fit
wxB = np.ndarray(Nmolec); wdxB = np.ndarray(Nmolec)      # barrier distance Bell Fit
wxBp = np.ndarray(Nmolec); wdxBp = np.ndarray(Nmolec)    # barrier distance Bell Fit projected force
wxBF = np.ndarray(Nmolec); wdxBF = np.ndarray(Nmolec)    # barrier distance Friddle Fit
wfrange = np.ndarray(Nmolec); wdfrange = np.ndarray(Nmolec)
wkx = np.zeros(Nmolec); wdkx = np.zeros(Nmolec)
wkz = np.zeros(Nmolec); wdkz= np.zeros(Nmolec)
wL0range = np.zeros(Nmolec); wdL0range = np.zeros(Nmolec)
wLp = np.zeros(Nmolec); wdLp = np.zeros(Nmolec)
wL0 = np.zeros(Nmolec); wdL0 = np.zeros(Nmolec)
wdzrange = np.zeros(Nmolec); wddzrange = np.zeros(Nmolec)
wlrange = np.zeros(Nmolec); wdlrange = np.zeros(Nmolec)
warange = np.zeros(Nmolec); wdarange = np.zeros(Nmolec)
wlmax = np.zeros(Nmolec); wdlmax = np.zeros(Nmolec)

wCoefPow_F = np.zeros(Nmolec)
wCoefPow_Ffit = np.zeros(Nmolec)
wfov = np.zeros(Nmolec, dtype=int)

if SelectMolecule:
    for ix, x in enumerate(ListMolec):
        c = colormap[ix] ; #c = next(color)
        dfx = dfstage[dfstage['Molecule ID']==x].reset_index()
        BeadRad_i = dfx['Bead radius (nm)']
        symb='o' #if dfx['Intern'][0] == 'Maryne' else 's'
        wfov[ix] = dfx['FOV'].iloc[0]
        cFOV = colormapFOV[int(dfx['FOV'].iloc[0])-1]

   ###################  Plots as a function of Power  #############"    
        wP = dfx['% High power']; wfP = dfx['Force_pN']
        wkxP = dfx['kx_c']; wkyP = dfx['ky_c']
        wDxP = dfx['Dx_c']*1e-6; wDyP = dfx['Dy_c']*1e-6
        wlCcalP = dfx['Cp1Lo_c']; wlcalP = dfx['p1Lo_c']
   
        wfL = dfx['ForceLumicks_pN']; wdfL = dfx['dForceLumicks_pN']
        wFCyclex = dfx['Fcyclex']; wFCycley = dfx['Fcycley']
        wDCyclex = dfx['Dcyclex']*1e-6; wDCycley = dfx['Dcycley']*1e-6
        wlco = dfx['Length_cyclic_open (nm)']; wdlco = dfx['error of length_cyclic_open']
        wFSDx = dfx['FSDx']; wFSDy = dfx['FSDy']
        wFSDx_c = dfx['FSDx_c']; wFSDy_c = dfx['FSDy_c']
        wkxPB = dfx['kx']; wkyPB = dfx['ky']
        wDxPB = dfx['Dx']; wDyPB = dfx['Dy']
        wkxrPB = dfx['kxraw']; wkyrPB = dfx['kyraw']
        wDxrPB = dfx['Dxraw']; wDyrPB = dfx['Dyraw']
        wkxgP= dfx['kxglobal']; wDxgP = dfx['Dxglobal']; wdkxgP= dfx['dkxglobal']; wdDxgP = dfx['dDxglobal']
        wkygP = dfx['kyglobal']; wDygP = dfx['Dyglobal']; wdkygP = dfx['dkyglobal']; wdDygP = dfx['dDyglobal']   
   #     wlP = dfx['p1Lo_D']; wdlP = dfx['p2Lo_D']
        wlP = dfx['p1Lo_nm']; wdlP = dfx['p2Lo_nm']

        wfmin_Hz = dfx['fmin_Hz']; wfmax_Hz = dfx['fmax_Hz']
        wkCyclex = wFCyclex/(wlco+BeadRad_i)
        wkCycley = wFCycley/(wlco+BeadRad_i)
       
        wlgP = dfx['Lglobal']  # wlP #  wlgP = wlcalP
   #     wlcalP = wlgP
        wfxgP = wkxgP*(wlgP+BeadRad_i); wdfxgP = wdkxgP*(wlgP+BeadRad_i)
        wfygP = wkygP*(wlgP+BeadRad_i); wdfygP = wdkygP*(wlgP+BeadRad_i) 
  #      print(wkxgP); print(wlgP); print(BeadRad_i)
   #     print(list(zip(wkxgP,wlgP,BeadRad_i)))
        
        wfxcalP = wkxP*(wlcalP+BeadRad_i)
        wfycalP = wkyP*(wlcalP+BeadRad_i)
        waP = dfx['Angle']; waCP = dfx['AngleC']
        friction0 = 6*np.pi*1.e-9*BeadRad_i
        wDxytheoP = 1e-6*kBT/(friction0 * FaxenXY(BeadRad_i, wlcalP))  #   friction0 = 6*np.pi*1.e-9*BeadRad       # units pN.s/nm
        wfPproj = wfP/np.cos(waP*np.pi/180)
        # here the calculation take into account the problem in n the PSD spectrum fit
        wkxPcorr = wkxP*np.sqrt(wDxytheoP/wDxP)  # wkxP*np.sqrt(wDxP/wDxytheoP)
        wkyPcorr = wkyP*np.sqrt(wDxytheoP/wDyP)  # wkyP*np.sqrt(wDyP/wDxytheoP)
        wPlateaux = 4*kBT**2*wDxP/(wDxytheoP**2*wkxP**2)     # 2*kBT**2/(wDxP*wkxP**2)
        wPlateauy = 4*kBT**2*wDyP/(wDxytheoP**2*wkyP**2)     # 2*kBT**2/(wDyP*wkyP**2)
        wPlateauxS = 1e6*4*kBT**2/(wDxPB*wkxPB**2)   
        wPlateauyS = 1e6*4*kBT**2/(wDyPB*wkyPB**2) 
        wPlateauxrS = 1e6*4*kBT**2/(wDxrPB*wkxrPB**2) 
        wPlateauyrS = 1e6*4*kBT**2/(wDyrPB*wkyrPB**2)  
        wPlateauxg = 1e6*4*kBT**2/(wDxgP*wkxgP**2) 
        wPlateauyg = 1e6*4*kBT**2/(wDygP*wkygP**2)  
        
        axPf = axallPower[ix//ncol2, ix%ncol2]
        axPfB = axallPowerBis[ix//ncol2, ix%ncol2]
        axPk = axallkxky[ix//ncol2, ix%ncol2]
        axPkB = axallkxkyBis[ix//ncol2, ix%ncol2]
        axPkT = axallkxkyTer[ix//ncol2, ix%ncol2]
        axPD = axallDxDy[ix//ncol2, ix%ncol2]
        axPDB = axallDxDyBis[ix//ncol2, ix%ncol2]
        axPA = axallPAngle[ix//ncol2, ix%ncol2]
        axPP = axallPPlateau[ix//ncol2, ix%ncol2]
        axPL = axallPLength[ix//ncol2, ix%ncol2]
        
        axPf.errorbar(wP, wfP, c=c, fmt=symb+'-', label='F vertical', alpha=0.5) 
        axPf.errorbar(wP, wfPproj, c=c, fmt='s-', markerfacecolor='none', label='F projected', alpha=0.5) 
        axPf.errorbar(wP, wfxcalP, c='k', fmt='+', label='Fx direct cal') 
        axPf.errorbar(wP, wfycalP, c='k', fmt='x', label='Fy direct cal') 
        axPf.errorbar(wP, wfL, c=c, fmt='*--', label='F Lumicks', alpha=0.5) 
        axPf.errorbar(wP, wFCyclex, c='k', fmt='h-.', markerfacecolor='none', label='Fcyclex', alpha=0.5) 
        axPf.errorbar(wP, wFCycley, c='k', fmt='d-.', markerfacecolor='none', label='Fcycley', alpha=0.5) 
        axPf.errorbar(wP, wFSDx_c, c='k', fmt='o-.', label='FSDx_c', alpha=0.5) 
        axPf.errorbar(wP, wFSDy_c, c='k', fmt='s-.', label='FSDy_c', alpha=0.5) 
        axPfB.errorbar(wP-DeltaP, wfxgP, wdfxgP, c=c, fmt=symb+'-', label='Fx g') 
        axPfB.errorbar(wP+DeltaP, wfygP, wdfygP, c=c, fmt=symb+'-', markerfacecolor='none', label='Fy g')
        
        axPL.errorbar(wP, wlcalP, c=c, fmt=symb+'-', label='L median anchor point', alpha=0.5)
        if not rapa: axPL.errorbar(wP, wlgP, c=c, fmt=symb+'-', label='L cycle', alpha=0.5)
        axPL.errorbar(wP, wlCcalP, c=c, fmt=symb+'-', label='L contour anchor point', markerfacecolor='none', alpha=0.5)
                
        axPA.errorbar(wP, waP, c=c, fmt=symb+'-', label='From Median', alpha=0.5) 
        axPA.errorbar(wP, waCP, c=c, fmt='s-', markerfacecolor='none', label='From Contour', alpha=0.5) 
        axPA.axhline(y=waP.mean(), color='k', linestyle='--', label='Mean Angle')
        
        axPk.errorbar(wP, wkxP, c=c, fmt=symb+'-', label='kx', alpha=0.5) 
        axPk.errorbar(wP, wkyP, c=c, fmt='s-', markerfacecolor='none', label='ky', alpha=0.5)
        axPk.errorbar(wP, wkxPcorr, c='k', fmt=symb+'-', label='kxcorr', alpha=0.5) 
        axPk.errorbar(wP, wkyPcorr, c='k', fmt='s-', markerfacecolor='none', label='kycorr', alpha=0.5) 
        axPk.errorbar(wP, wkCyclex, c='k', fmt='+', label='kcyclex', alpha=0.5)
        axPk.errorbar(wP, wkCycley, c='k', fmt='x', label='kcycley', alpha=0.5)
        
        axPkB.errorbar(wP, wkxgP, wdkxgP, c=c, fmt=symb+'-', label='kxg', alpha=0.5)
        axPkB.errorbar(wP, wkygP, wdkygP, c=c, fmt=symb+'-', markerfacecolor='none', label='kyg', alpha=0.5)
        axPkB.errorbar(wP, wkxrPB, c=c, fmt='s--', label='kxrB', alpha=0.5)
        axPkB.errorbar(wP, wkyrPB, c=c, fmt='s--', markerfacecolor='none', label='kyrB', alpha=0.5)
        axPkB.errorbar(wP, wkxPB, c='k', fmt='x', label='kxB', alpha=0.5)
        axPkB.errorbar(wP, wkyPB, c='k', fmt='+', label='kyB', alpha=0.5)
        axPkT.errorbar(wP-DeltaP, wkxgP, wdkxgP, c=c, fmt=symb+'-', label='kxg', alpha=0.5)
        axPkT.errorbar(wP+DeltaP, wkygP, wdkygP, c=c, fmt=symb+'-', markerfacecolor='none', label='kyg', alpha=0.5)
        
        axPD.errorbar(wP, wDxP, c=c, fmt=symb+'-', label='Dx', alpha=0.5) 
        axPD.errorbar(wP, wDyP, c=c, fmt='s-', markerfacecolor='none', label='Dy', alpha=0.5) 
        axPD.errorbar(wP, wDxytheoP, c='k', fmt='-', label='Dtheo', alpha=0.5)
        axPD.errorbar(wP, wDCyclex, c='k', fmt='+', label='Dcyclex', alpha=0.5) 
        axPD.errorbar(wP, wDCycley, c='k', fmt='x', label='Dcycley', alpha=0.5)
    
        axPDB.axhline(y=wDxgP[0]*1e-6, color='m', linestyle='--', label='Dxg')
        axPDB.axhline(y=wDygP[0]*1e-6, color='k', linestyle='--', label='Dyg')
        Dxmax = ( wDxgP[0] + wdDxgP[0] )*1e-6; Dxmin = ( wDxgP[0] - wdDxgP[0] )*1e-6
        axPDB.fill_between(np.arange(35), Dxmin, Dxmax, alpha=0.2)
        Dymax = ( wDygP[0] + wdDygP[0] )*1e-6; Dymin = ( wDygP[0] - wdDygP[0] )*1e-6
        axPDB.fill_between(np.arange(35), Dymin, Dymax, alpha=0.2)
        axPDB.errorbar(wP, wDxgP*1e-6, c=c, fmt='+-', label='Dxg', alpha=0.5) 
        axPDB.errorbar(wP, wDygP*1e-6, c=c, fmt='d-', markerfacecolor='none', label='Dyg', alpha=0.5)    

        wfc = (wkxgP+wkygP)*( wDxgP[0] + wDygP[0] )/4/(2*np.pi*kBT_pN_nm)
  #       fcg = (kxg+kyg)*( Dxg+Dyg )/4/(2*np.pi*kBT_pN_nm)
  #      print(wfc, wfmax_Hz) 
        iPok = wfc < wfmax_Hz/2
  #      iPok = wDxPB>0
        # axPDB.errorbar(wP[iPok], wDxP[iPok], c=c, fmt=symb+'-', label='Dx', alpha=0.5) 
        # axPDB.errorbar(wP[iPok], wDyP[iPok], c=c, fmt='s-', markerfacecolor='none', label='Dy', alpha=0.5)
        axPDB.errorbar(wP[iPok], wDxrPB[iPok]*1e-6, c=c, fmt=symb+'-', label='Dx', alpha=0.5) 
        axPDB.errorbar(wP[iPok], wDyrPB[iPok]*1e-6, c=c, fmt='s-', markerfacecolor='none', label='Dy', alpha=0.5)    
        # axPDB.errorbar(wP[iPok], wDxPB[iPok]*1e-6, c=c, fmt='+-', label='Dxr', alpha=0.5) 
        # axPDB.errorbar(wP[iPok], wDyPB[iPok]*1e-6, c=c, fmt='d-', markerfacecolor='none', label='Dyr', alpha=0.5)    
        
        axPP.errorbar(wP, wPlateauxS, c=c, fmt=symb+'-', label='Plateau x S', alpha=0.5) 
        axPP.errorbar(wP, wPlateauyS, c=c, fmt=symb+'-', markerfacecolor='none', label='Plateau y S', alpha=0.5) 
        axPP.errorbar(wP, wPlateauxrS, c='k', fmt='+', label='Plateau xr S', alpha=0.5) 
        axPP.errorbar(wP, wPlateauyrS, c='k', fmt='x', label='Plateau yr S', alpha=0.5) 
        axPP.errorbar(wP, wPlateauxg, c='k', fmt='+-', label='Plateau xg', alpha=0.5) 
        axPP.errorbar(wP, wPlateauyg, c='k', fmt='x-', label='Plateau yg', alpha=0.5) 
        
        frictionXY = friction0 * FaxenXY(BeadRad_i, p1Zo)
        Dxytheo_mu2_per_s = 1e-6*kBT/frictionXY
        Dmax = Dxytheo_mu2_per_s*(1+deltaD); Dmin = Dxytheo_mu2_per_s*(1-deltaD)
    #    axPD.fill_between(np.arange(35), Dmin, Dmax, alpha=0.2)
        selDx = (Dmin < wDxP) * (wDxP < Dmax)
        selDy = (Dmin < wDyP) * (wDyP < Dmax)  
        def avgarray(selx, sely, wx, wy):            
            n = len(wx); avg = np.zeros(n)
            selavg = selx + sely
            for i in range(n):
                if selx[i] and sely[i]:
                    avg[i] = (wx[i] + wy[i]) / 2
                elif selx[i] and not sely[i]:
                    avg[i] = wx[i]
                elif sely[i] and not selx[i]:
                    avg[i] = wy[i]
            return selavg, avg                
        selDavg, wDavg = avgarray(selDx, selDy, wDxP, wDyP)
        selDavg, wkavg = avgarray(selDx, selDy, wkxP, wkyP)
        selDavg, wFavg = avgarray(selDx, selDy, wfxcalP, wfycalP)

        axPD.errorbar(wP[selDavg], wDavg[selDavg], c='k', fmt='s', label='avg D sel', alpha=0.5)       
        axPk.errorbar(wP[selDavg], wkavg[selDavg], c='k', fmt='s', label='from D sel', alpha=0.5)       
        axPf.errorbar(wP[selDavg], wFavg[selDavg], c='k', fmt='s', label='from D sel', alpha=0.5)
        
   ###################  Recalculation of force based on D interval and max force  #############"    

        if len(wfxgP)>=1:
            wf = (wfxgP+wfygP)/2; wdf = (wdfxgP+wdfygP)  # wf = wFavg
            selLowF0 = (0 < wf) & (wf < MaxForce)
            wP_ = np.insert(wP[selLowF0].values, 0, 0)
            wFavg_ = np.insert(wf[selLowF0].values, 0, 0)
      #      print(len(wP_), len(wFavg_))
            slopePf, interceptPf, r, p, se = linregress(wP_, wFavg_)
            wffit = slopePf*wP + interceptPf; wdffit = np.zeros(len(wffit))
            wffit = wf; 
            if rapa: wffit = wf; wdffit = wdf
            selLowF = (wffit < MaxForce)
            axPf.errorbar(wP[selLowF], wffit[selLowF], c='r', fmt='-', label='fit', alpha=1)
        else:
            wffit=[]; slopePf=np.nan; wdffit=[]
        for ii in range(len(wffit)):   # try dfstage.update() ?
            iix = dfstage[dfstage['Molecule ID']==x].index.values[ii]
            dfstage.loc[iix, 'Ffit'] = wffit[ii]
        wCoefPow_Ffit[ix] = slopePf
  
   ###################  Plots as a function of Force  #############"    
  #      wfv = dfx['Force_pN']; wdfv = dfx['dForce_pN']
        wf = wffit; wdf = wdffit # np.zeros(len(wf))
        wf = (wfxgP+wfygP)/2; wdf = (wdfxgP+wdfygP)
   #      wlgP = dfx['Lglobal']
        wk = dfx['newkoff']
        wl = dfx['Length_calibration (nm)']; wdl = dfx['error of length']
        
        wlD = dfx['p1Lo_D']; wdlD = dfx['p2Lo_D']
        wlCcal = dfx['Cp1Lo_c']; wlcal = dfx['p1Lo_c']
        'Cp1Lo_D'
        'Cp1Lo_c'
        wa = dfx['Angle']; waC = dfx['AngleC']
        wlcc = dfx['Length_cyclic_close (nm)']; wdlcc = dfx['error of length_cyclic_close']
   #     wlP = dfx['Length_calibration (nm)']; wdlP = dfx['error of length']
        wnc = dfx['No Close']/dfx['Total cycles']
        wdk = dfx['newdkoff2AB']; # wdk = dfx['newdkoff']; 
        wk2A = dfx['newkoff2A']; wdk2A = dfx['newdkoff2A']
        wk2B = dfx['newkoff2B']; wdk2B = dfx['newdkoff2B']        
        wDZ = dfx['DeltaZ']; wdDZ = dfx['dDeltaZ']
        wk2 = dfx['kminmax']; wkmin = dfx['kmin']; wkmax = dfx['kmax']
        wfproj = wf/np.cos(wa*np.pi/180)
        wDL = wlco - wlcc; wdDL = wdlcc + wdlco
         
 #       selkoff = ~np.isnan(wk)[np.append(selLowF, np.zeros(len(wk)-len(selLowF)),dtype='bool')]
        selkoff = selLowF * ~np.isnan(wk)
        selkoff = ~np.isnan(wf) * ~np.isnan(wk)
   
        axBis = axallBis[ix//ncol2, ix%ncol2]
        axL = axallLength[ix//ncol2, ix%ncol2]
        axLB = axallLengthBis[ix//ncol2, ix%ncol2]
        axA = axallAngle[ix//ncol2, ix%ncol2]
        axNC = axallNoClose[ix//ncol2, ix%ncol2]
        axFL = axallFLumicks[ix//ncol2, ix%ncol2]
        axDZ = axallDeltaz[ix//ncol2, ix%ncol2]
        axLZ = axallLengthDeltaZ[ix//ncol2, ix%ncol2]
        axBisAB = axallBisAB[ix//ncol2, ix%ncol2]
        
        wdzrange[ix] = wDZ.mean(); wddzrange[ix] = wDZ.std()
        wlrange[ix] = wl.mean(); wdlrange[ix] = wl.std()     
        warange[ix] = wa.mean(); wdarange[ix] = wa.std()     
        wfrange[ix] = wf.mean(); wdfrange[ix] = wf.std()
        wlmax[ix] = wl.max()
        
        ax.errorbar(wf[selkoff], wk[selkoff], wdk[selkoff], c=c, fmt=symb+'-', label='ID='+str(x), alpha=0.3) 
        axBis.errorbar(wf[selkoff], wk[selkoff], wdk[selkoff], c=c, fmt=symb+'-', label='F Dsel', alpha=0.3) 
        axBis.errorbar(wfproj[selkoff], wk[selkoff], wdk[selkoff], wdf[selkoff], c=c, fmt=symb+'-', markerfacecolor='none', label='F projected', alpha=0.3) 
        axBisAB.errorbar(wf[selkoff], wk2A[selkoff], wdk2A[selkoff], c='g', fmt='o-', label='ID='+str(x), alpha=0.5) 
        axBisAB.errorbar(wf[selkoff], wk2B[selkoff], wdk2B[selkoff], c='r', fmt='o-', label='ID='+str(x), alpha=0.5) 
        axMat.errorbar(wk2A[selkoff], wk2B[selkoff], c=c, fmt='o')
        axOtheFlyDx.scatter(wDxP[selkoff], wDCyclex[selkoff], c=np.log(wkxP[selkoff]/0.1), label='Dx', cmap='gray', alpha=0.5)
        axOtheFlyDy.scatter(wDyP[selkoff], wDCycley[selkoff], c=np.log(wkxP[selkoff]/0.1), label='Dy', cmap='gray', alpha=0.5)
        axOtheFlykx.scatter(wkxP[selkoff], wkCyclex[selkoff], c=np.log(wkxP[selkoff]/0.1), label='kx', cmap='gray', alpha=0.5)
        axOtheFlyky.scatter(wkyP[selkoff], wkCycley[selkoff], c=np.log(wkxP[selkoff]/0.1), label='ky', cmap='gray', alpha=0.5)
        axOtheFlySDx.scatter(wFSDx[selkoff], wFSDx_c[selkoff], c=c, label='FSDx', alpha=0.5)
        axOtheFlySDy.scatter(wFSDy[selkoff], wFSDy_c[selkoff], c=c, label='FSDy', alpha=0.5)
        axContour.scatter(wlcalP, wlCcalP, c=c, label='Open length', alpha=0.5)
        axXYk.scatter(wkxP[selkoff], wkyP[selkoff], c=c, label='XYkcal', alpha=0.5)
        axXYD.scatter(wDxP[selkoff], wDyP[selkoff], c=c, label='XYDcal', alpha=0.5)
        axXYF.scatter(wfxcalP[selkoff], wfycalP[selkoff], c=c, label='XYFcal', alpha=0.5)

        axA.errorbar(wf[~np.isnan(wf)], wa[~np.isnan(wf)], c=c, fmt=symb+'-', label='From Median', alpha=0.5) 
        axA.errorbar(wf[~np.isnan(wf)], waC[~np.isnan(wf)], c=c, fmt=symb+'-', markerfacecolor='none', label='From Contour', alpha=0.5) 
        axA.axhline(y=wf[~np.isnan(wf)].mean(), color='k', linestyle='--', label='Mean Angle')
        
        nonan_f_l = ~pd.isnull(wf) * ~pd.isnull(wl)
        nonan_f_lD = ~pd.isnull(wf) * ~pd.isnull(wlD)
        
        axL.errorbar(wf[nonan_f_l], wl[nonan_f_l], wdl[nonan_f_l], wdf[nonan_f_l], c=c, fmt=symb+'-', label='L calibration', alpha=0.5) 
   #     axL.errorbar(wfproj[nonan_f_l], wl[nonan_f_l], wdl[nonan_f_l], wdf[nonan_f_l], c=c, fmt=symb+'-', markerfacecolor='none', label='F projected', alpha=0.5) 
        axL.errorbar(wf[nonan_f_l], wlcal[nonan_f_l], c='k', fmt='+', label='L calib', alpha=0.5) 
        axL.errorbar(wf[nonan_f_l], wlCcal[nonan_f_l], c='k', fmt='x', label='LC calib', alpha=0.5) 
        axL.errorbar(wf[nonan_f_l], wlcc[nonan_f_l], wdlcc[nonan_f_l], wdf[nonan_f_l], c=c, fmt='s-', markerfacecolor='none', label='LC closed', alpha=0.5) 
        axL.errorbar(wf[nonan_f_l], wlco[nonan_f_l], wdlco[nonan_f_l], wdf[nonan_f_l], c=c, fmt='o-', markerfacecolor='none', label='LC open', alpha=0.5) 
        axLB.errorbar(wf, wlD, wdlD, wdf, c=c, fmt=symb+'-', label='L Direct cal', alpha=0.5) 
        axLB.errorbar(wf, dfx['Lglobal'], c=c, fmt=symb+'-', markerfacecolor='none', label='L spec cal', alpha=0.5) 
        
        axNC.errorbar(wf, wnc, c=c, fmt=symb+'-', label='ID='+str(x), alpha=0.5) 
        axFL.errorbar(wf, wfL, wdfL, wdf, c=c, fmt=symb+'-', label='ID='+str(x), alpha=0.5) 
        axFL.plot(wfL, wfL, c='k', linestyle='-', label='ID='+str(x), alpha=0.3) 

        axDZ.errorbar(wf, wDZ, wdDZ, wdf, c=c, fmt=symb+'-', label='ID='+str(x), alpha=0.5) 
        axDZ.errorbar(wf, wDL, wdDL, wdf, c=c, fmt=symb+'-', markerfacecolor='none', label='ID='+str(x), alpha=0.5) 
        axLZ.errorbar(wl, wDZ, wdDZ, wdlcc, c=c, fmt=symb+'-', label='ID='+str(x), alpha=0.5) 
        def WLCL0pen(x, *p):
            L0 = dfx['L0'][0]; penalization = abs(p[1]-L0)*10
            return WLC(x, p[0], p[1]) + penalization

   ###################  Fits for Length  #############"                   
        if sum(nonan_f_l)>3:
            pEqElast2, pcovEqElast2 = curve_fit(AffineL, wf[nonan_f_l], wl[nonan_f_l], p0=[1000, 1e-3])
            wffit = 0.15*(np.arange(100)+1)
            wlfit = AffineL(wffit, pEqElast2[0], pEqElast2[1])
            axL.plot(wffit, wlfit, linestyle='-.', c='k', label='Affine fit', alpha=0.5)
            while True:
                try:
           #         pEqWLC, pcovEqWLC = curve_fit(WLC, wl[nonan_f_l], wf[nonan_f_l], p0=[20, wl[nonan_f_l].max()+100])
                    pEqWLC, pcovEqWLC = curve_fit(WLC, wlD[nonan_f_lD], wf[nonan_f_lD], p0=[20, wlD[nonan_f_lD].max()+100])
            #        pEqWLC, pcovEqWLC = curve_fit(WLC, wlD, wf, p0=[20, wlD.max()])
                    break
                except RuntimeError:
                    print("No convergence"); pEqWLC=[np.nan,np.nan]; pcovEqWLC=[[np.nan, np.nan],[np.nan, np.nan]]
                    break
            wlfitWLC = 20*(np.arange(100)+1)
            wffitWLC = WLC(wlfitWLC, pEqWLC[0], pEqWLC[1])
       #     axL.plot(wffitWLC, wlfitWLC, linestyle='--', c='k', label='WLC fit', alpha=0.5)
            axLB.plot(wffitWLC, wlfitWLC, linestyle='--', c='k', label='WLC fit', alpha=0.5)
        else: 
            pEqElast2 =[np.nan,np.nan]; pEqWLC=[np.nan,np.nan]; pcovEqWLC=[[np.nan, np.nan],[np.nan, np.nan]]
            
        wElast = dfx['kx (pN/nm)'].dropna(); wkx[ix] = wElast.mean(); wdkx[ix] = wElast.std()
        wkz[ix] = pEqElast2[1]; wL0range[ix] = pEqElast2[0]; wdkz[ix]=0; wdL0range[ix]=0
        wLp[ix] = pEqWLC[0] ; wdLp[ix]=np.sqrt(np.diag(pcovEqWLC))[0]**0.5 
        wL0[ix]=pEqWLC[1] ; wdL0[ix]=np.sqrt(np.diag(pcovEqWLC))[1]**0.5

   ###################  Fits for Off-rates  #############"                           
        kc = 0.5; #wft = np.logspace(np.log(0.1), np.log(15), num=30)

        if len(wFavg[selkoff])>=3:
            (pEqB, eEqB, res) = fitBellFriddle(logBellxB, wf[selkoff], wk[selkoff])#, 1/sd)
            wft = np.linspace(0,wf[selkoff].max())
            wlEqB = BellxB(wft, pEqB[0], pEqB[1])
            axBis.errorbar(wft, wlEqB, c='k', fmt='-', label='fit Bell', alpha=0.5)
            (pEqBp, eEqBp, resp) = fitBellFriddle(logBellxB, wfproj[selkoff], wk[selkoff])#, 1/sd)
            wlEqBp = BellxB(wft, pEqBp[0], pEqBp[1])
            axBis.errorbar(wft, wlEqBp, c='k', fmt='--', label='fit Bell proj', alpha=0.5)
            (pEqF, eEqF, resF) = fitBellFriddle(logFriddle, wf[selkoff], wk[selkoff])#, 1/sd)
            wlEqF = Friddle(wft, pEqF[0], pEqF[1], pEqF[2])   
        else:
            pEqB = [np.nan, np.nan]; pEqF = [np.nan, np.nan, np.nan]
            eEqB =[np.nan, np.nan]; eEqF = [np.nan, np.nan, np.nan]
            pEqBp = [np.nan, np.nan]
            eEqBp = [np.nan, np.nan]
            
        xB = pEqB[1]; dxB = eEqB[1]; xBF = pEqF[1]; dxBF = eEqF[1]
        xBp = pEqBp[1]; dxBp = eEqBp[1]
        wxB[ix] = xB; wdxB[ix] = dxB; wxBF[ix] = xBF; wdxBF[ix] = dxBF
        wxBp[ix] = xBp; wdxBp[ix] = dxBp
        k0 = pEqB[0] ; dk0 = eEqB[0]; k0F = pEqF[0] ; dk0F = eEqF[0]
        k0p = pEqBp[0] ; dk0p = eEqBp[0]
        if k0>0:
            wlogk0[ix] = np.log(k0); wdlogk0[ix] = dk0/k0
        else:
            wlogk0[ix] = np.nan; wdlogk0[ix] = np.nan
        if k0p>0:    
            wlogk0p[ix] = np.log(k0p); wdlogk0p[ix] = dk0p/k0p
        else:
            wlogk0p[ix] = np.nan; wdlogk0p[ix] = np.nan
        if k0F>0:
            wlogk0F[ix] = np.log(k0F); wdlogk0F[ix] = dk0F/k0F
        else:
            wlogk0F[ix] = np.nan; wdlogk0F[ix] = np.nan
        if xB>0: 
            fB = kBT/xB; dfB = -dxB*kBT/xB**2
            fBF = kBT/xBF; dfBF = -dxBF*kBT/xBF**2
        else:
            fB = np.nan; dfB = np.nan
            fBF = np.nan; dfBF = np.nan
        wfB[ix] = fB; wdfB[ix] = dfB; wfBF[ix] = fBF; wdfBF[ix] = dfBF      
        fres = fB*np.log(kres/k0) if k0!=0 else 0
        fresp = kBT/xBp*np.log(kres/k0p) if k0p!=0 and xBp!=0 else 0

        print(x, 'Fitting Molec Bell' , '(', len(wk[selkoff]) , ' points):  k0=  %.6f' % k0,'+/-  %.6f' % dk0,  ' xB= %.2f'% xB,'+/- %.2f' % dxB, 'fres =%.2f' % fres)
        print(x, 'FitBis  Molec Bell' , '(', len(wk[selkoff]) , ' points):  k0=  %.6f' % k0p,'+/-  %.6f' % dk0p,  ' xB= %.2f'% xBp,'+/- %.2f' % dxBp, 'fres =%.2f' % fresp)
        print(x, 'Fitting Molec Frid' , '(', len(wk[selkoff]) , ' points):  k0=  %.6f' % k0F,'+/-  %.6f' % dk0F,  ' xB= %.2f'% xBF,'+/- %.2f' % dxBF, 'kc =%.2f' % pEqF[2])

     #   axRR.scatter(dfx['R'+xy][0], slopePf**(1/3), c=cFOV, alpha=0.5)
        axRR.scatter(dfx['Rprec'+xy][0], slopePf**(1/3), c=cFOV, marker='s', alpha=0.5)
        axRR1.scatter(dfx['Rprec'+xy][0], wl[nonan_f_l].median(), c=cFOV, alpha=0.5)

        if fB>0:
            ax2.errorbar(wfB[ix], np.exp(wlogk0[ix]), dk0, wdfB[ix], c=c, label='ID='+str(x), fmt=symb, alpha=0.5)
            ax2.errorbar(wfBF[ix], np.exp(wlogk0F[ix]), dk0F, wdfBF[ix], c=c, fmt=symb, markerfacecolor='none', alpha=0.5)
            ax2Bis.errorbar(wxB[ix], np.exp(wlogk0[ix]), dk0, wdxB[ix], c=c, label='ID='+str(x), fmt=symb, markerfacecolor='none', alpha=0.5)
            ax2Ter.errorbar(wxB[ix], wkx[ix], wdkx[ix], wdxB[ix], c=c, label='ID='+str(x), fmt=symb, alpha=0.5)
            ax2Ter1.errorbar(wxB[ix], wlmax[ix], wdlmax[ix], wdxBp[ix], c=c, label='ID='+str(x), fmt=symb, alpha=0.5)
            ax2Ter2.errorbar(wfB[ix], wkz[ix], wdkz[ix], wdfB[ix], c=c, label='ID='+str(x), fmt=symb, alpha=0.5)
            ax2Ter3.errorbar(wkx[ix], wkz[ix], wdkz[ix], wdkx[ix], c=c, label='ID='+str(x), fmt=symb, alpha=0.5)
            ax2Ter4.errorbar(wLp[ix], wkz[ix], wdkz[ix], wdLp[ix], c=c, label='ID='+str(x), fmt=symb, alpha=0.5)
            ax2Ter6.errorbar(wL0[ix], wLp[ix], wdLp[ix], wdL0[ix], c=c, label='ID='+str(x), fmt=symb, alpha=0.5)
            ax2Ter5.errorbar(wfB[ix], wLp[ix], wdLp[ix], wdfB[ix], c=c, label='ID='+str(x), fmt=symb, alpha=0.5)
            ax2Ter7.errorbar(wlrange[ix], wdzrange[ix], wddzrange[ix], wdlrange[ix], c=c, label='ID='+str(x), fmt=symb, alpha=0.5)
            ax4.errorbar(wfB[ix], wfrange[ix], wdfrange[ix], wdfB[ix], c=c, label='ID='+str(x), fmt=symb, alpha=0.5)
            ax3.errorbar(wf, wfL, wdfL, wdf, c=c, label='ID='+str(x), fmt=symb, alpha=0.5)
        axL.set_ylabel('Open JDNA Length (nm)'); axL.set_ylim(500,2000); axL.set_xlabel('Force')
        axLB.set_ylabel('Open JDNA Length (nm)'); axLB.set_ylim(500,2000); axLB.set_xlabel('Force')
        axA.set_ylabel('Pulling Angle (deg)'); axA.set_ylim(0,60); axA.set_xlabel('Force')
        axNC.set_ylabel('No close fraction'); axNC.set_ylim(0,1)
        axFL.set_ylabel('Force Lumicks (pN)'); axFL.set_ylim(0,15)
        axPf.set_ylabel('Force (pN)'); axPf.set_ylim(0,35)
        axPL.set_ylabel('Open JDNA length (nm)'); axPL.set_ylim(500,2000)
        axPfB.set_ylabel('Force (pN)'); axPfB.set_ylim(0,35)
        axPk.set_ylabel('kx, ky (pN/nm)'); axPk.set_ylim(0,1.2e-2)
        axPkB.set_ylabel('kx, ky (pN/nm)'); axPkB.set_ylim(0,1.2e-2)
        axPkT.set_ylabel('kx, ky (pN/nm)'); axPkT.set_ylim(0,1.5e-2)
        axPD.set_ylabel('Dx, Dy (µm²/s)'); axPD.set_ylim(0,0.4)
        axPDB.set_ylabel('Dx, Dy (µm²/s)'); axPDB.set_ylim(0,0.4)
        axPA.set_ylabel('Pulling Angle (deg)'); axPA.set_ylim(0,70)
        axPP.set_ylabel('Plateau PSD'); axPP.set_ylim(1e4, 1e10); axPP.set_yscale('log')
        if rapa: axPf.set_ylim(0,25)
        if rapa: axPfB.set_ylim(0,25)
        axDZ.set_ylabel('Delta Z jump/ Delta L (nm)'); axDZ.set_ylim(0,300)
        axLZ.set_ylabel('Delta Z jump (nm)'); axLZ.set_ylim(mindz, maxdz)
        for axBisi in [axBis, axBisAB]:
            axBisi.set_ylabel('koff (1/s)'); axBisi.set_yscale('log')
            axBisi.set_ylim(1e-5,0.2)
            if rapa: axBisi.set_ylim(0.01, 0.4)
        for axBisi in [axBis, axBisAB, axL, axLB, axA, axNC, axFL, axPf, axPL, axPfB, axPk, axPkB, axPkT, axPD, axPDB, axPA, axPP, axDZ, axLZ]:
            if axBisi == axLB:
                res = "L0={:1.0f}".format(wL0range[ix])+" kz={:1.3f}".format(wkz[ix])
                res = res+ "\n WLC L0={:1.0f}".format(wL0[ix])+" Lp={:1.0f}".format(wLp[ix])
                axBisi.text(0.5, 1400, 'ID '+str(x)+' FOV '+str(wfov[ix])+':'+res , size=8)
            elif axBisi == axNC or axBisi == axFL  or axBisi == axPf or axBisi == axPfB or axBisi == axDZ or axBisi == axA or axBisi == axPA:
                if axBisi == axPf and len(wfP[~np.isnan(wfP)])>2: 
                    nonan = ~np.isnan(wP) * ~np.isnan(wfP)
                    slope, intercept, r, p, se = linregress(wP[nonan], wfP[nonan])
                    text = 'Jim F/P '+"{:1.2f}".format(slope)+' F(0) '+"{:1.2f}".format(intercept)
                    text += '\nDsel F/P '+"{:1.2f}".format(wCoefPow_Ffit[ix])+' F(0) '+"{:1.2f}".format(interceptPf)
                    axBisi.text(1, 8, text, size=10)
                    wCoefPow_F[ix] = slope
              #      axBisi.text(80, 0.5, TextLinear(wfP, wP), size=10)
                axBisi.text(1, 5, 'ID '+str(x)+' FOV '+str(wfov[ix]) , size=8)
            elif axBisi == axLZ:
                axBisi.text(700, 110, 'ID '+str(x)+' FOV '+str(wfov[ix]) , size=8)
            elif axBisi == axBis or axBisi == axBisAB:
                res = " k0 ={:1.6f}".format(k0)+" xB ={:1.2f}".format(xB)
                res += "\n k0p={:1.6f}".format(k0p)+" xBp={:1.2f}".format(xBp)
                ypos = 0.01 if rapa else 0.02
                axBisi.text(0.5, ypos, 'ID '+str(x)+' FOV '+str(wfov[ix])+'\n'+res , size=10)
            elif axBisi == axPk or axBisi == axPkB or axBisi == axPkT:
                axBisi.text(0.5, 0.008, 'ID '+str(x)+' FOV '+str(wfov[ix]) , size=10)
            elif axBisi == axPD or axBisi == axPDB:
                axBisi.text(0.5, 0.05, 'ID '+str(x)+' FOV '+str(wfov[ix]) , size=10)
            elif axBisi == axPP:
                axBisi.text(0.5, 1e5, 'ID '+str(x)+' FOV '+str(wfov[ix]) , size=10)
            elif axBisi == axPL:
                axBisi.text(0.5, 1500, 'ID '+str(x)+' FOV '+str(wfov[ix]) , size=10)
            axBisi.set_xlabel('Force (pN)'); axBisi.set_xlim(0,35)
            if rapa: axBisi.set_xlim(0,25)
            if axBisi == axPf or axBisi == axPL or axBisi == axPfB or  axBisi == axPk or axBisi == axPkB or axBisi == axPkT or axBisi == axPD or axBisi == axPDB or axBisi == axPA or axBisi == axPP :
                axBisi.set_xlabel('% High power'); axBisi.set_xlim(0, HighPowerMax)
            if ix%ncol2!=0: axBisi.axes.get_yaxis().set_visible(False)
            if ix//ncol2!=nrow2-1: axBisi.axes.get_xaxis().set_visible(False)
        if ix==0:
            axBis.legend(loc="upper left", fontsize=10)
            for axii in [axPf, axPL, axPfB, axPA, axPP, axPD, axPDB, axPk, axPkB, axPkT, axL, axLB, axA]: axii.legend(fontsize=10)

        if rapa: axBis.set_xlim(0,25)#; axBis.set_xscale('log')
        axLZ.set_xlabel('Open JDNA Length (nm)'); axLZ.set_xlim(500,1800)
        axMat.plot(wk2A,wk2A,'k'); axMat.set_xscale('log'); axMat.set_yscale('log')
        axMat.set_xlabel('Early koff (1/s)'); axMat.set_ylabel('Late koff (1/s)')
        axOtheFlykx.plot(wkxP, wkxP,'k', label='x=y'); axOtheFlyky.plot(wkyP, wkyP,'k', label='x=y')
        axOtheFlyDx.plot(wDxP, wDxP,'k', label='x=y'); axOtheFlyDy.plot(wDyP, wDyP,'k', label='x=y')
        axOtheFlySDx.plot(wFSDx, wFSDx,'k', label='x=y'); axOtheFlySDy.plot(wFSDy, wFSDy,'k', label='x=y')
        axXYk.plot(wkxP, wkxP,'k', label='x=y')
        axXYD.plot(wDxP, wDxP,'k', label='x=y')
        axXYF.plot(wfxcalP, wfxcalP,'k', label='x=y')
        
        axOtheFlyDx.set_xlim(3e-2, 1); axOtheFlyDx.set_ylim(3e-2, 1)
        axOtheFlyDy.set_xlim(3e-2, 1); axOtheFlyDy.set_ylim(3e-2, 1)
        axXYD.set_xlim(3e-2, 1); axXYD.set_ylim(3e-2, 1)
        axOtheFlykx.set_xlim(3e-4, 2e-2); axOtheFlykx.set_ylim(3e-4, 2e-2)
        axOtheFlyky.set_xlim(3e-4, 2e-2); axOtheFlyky.set_ylim(3e-4, 2e-2)
        axXYk.set_xlim(3e-4, 2e-2); axXYk.set_ylim(3e-4, 2e-2)
        axOtheFlySDx.set_xlim(0, 15); axOtheFlySDx.set_ylim(0, 15)
        axOtheFlySDy.set_xlim(0, 15); axOtheFlySDy.set_ylim(0, 15)
        axXYF.set_xlim(0, 15); axXYF.set_ylim(0, 15)
        axContour.plot(wlcalP, wlcalP,'k', label='x=y')
        axContour.set_xlim(500,1500); axContour.set_ylim(500,1500)
        axContour.set_xlabel('Median anchor point'); axContour.set_ylabel('Contour anchor point')
        if x==1: axContour.legend()
        for axOtheFLy in [axOtheFlyDx, axOtheFlyDy, axOtheFlykx, axOtheFlyky, axOtheFlySDx, axOtheFlySDy]:
            axOtheFLy.set_xscale('log'); axOtheFLy.set_yscale('log')
            if x==1: axOtheFLy.legend()
            axOtheFLy.set_xlabel('Classical calibration'); axOtheFLy.set_ylabel('On the fly calibration')
        for axXY in [axXYD, axXYk, axXYF]:
            axXY.set_xscale('log'); axXY.set_yscale('log')
            if x==1: axXY.legend()
            axXY.set_xlabel('X parameter'); axXY.set_ylabel('Y parameter')

if rapa:
    ymin = 700; ymax = 1000
else:
    ymin = 600; ymax = 2000
for xy in ['precx', 'precy', 'x', 'y']:   
    figname4a = 'RDiffusion'+xy; fig = plt.figure(figname4a, figsize=(max(ListMolec)*0.25,6), dpi=100); ax4a = plt.gca()
 #   sns.boxplot(x="Molecule ID", y='R'+xy, data=dfstage, ax=ax4a)
    ax4a.set_xlabel('Molecule ID')#    ax4a.set_xticklabels(ax4a.get_xticklabels(), rotation=90)
    ax4a.set_ylabel('Bead Radius deduced from Diffusion '+xy+' (nm)')
    ax4a.set_ylim(ymin, ymax); #ax4a.set_xlim(0,max(ListMolec)+1)
    ListMolec_str = [str(ID)  for ID in ListMolec]
    ListRxy = [dfstage.loc[dfstage['Molecule ID'] == ID]['R'+xy].median()  for ID in ListMolec]
    ListBeadR = [dfstage.loc[dfstage['Molecule ID'] == ID]['Bead radius (nm)'].median()  for ID in ListMolec]
 #   ax4a.scatter(dfstage['Molecule ID'], dfstage['R'+xy], c='k', marker='x', label='Measured')
    ax4a.scatter(ListMolec_str, ListRxy, c='k', marker='x', label='Measured')
 #   ax4a.scatter(dfstage['Molecule ID'], dfstage['Bead radius (nm)'], c='r', marker='+', label="Nominal")
    ax4a.scatter(ListMolec_str, ListBeadR, c='r', marker='+', label="Nominal")
    ax4a.legend(loc="upper right", fontsize=5)
    if SaveGraph: plt.tight_layout(); plt.savefig(outpath+figname4a, transparent=False)

ax.set_xlabel('Force (pN)'); ax.set_ylabel('koff (1/s)'); ax.set_yscale('log')
ax.set_ylim(1e-5,0.2); ax.set_xlim(0,35)
if rapa: ax.set_ylim(1e-2, 0.4); ax.set_xlim(0,25)

nonan = ~np.isnan(wxB) * ~np.isnan(wlogk0)
if sum(nonan)>=2:
    pEq, pcovEq = curve_fit(Affine, wxB[nonan], wlogk0[nonan] - np.log(koffSPR), p0=[-1, 1])
    wxBfit = 0.05*(np.arange(100)+1)
    wlogk0_fit = Affine(wxBfit, pEq[0], pEq[1]) + np.log(koffSPR)
    (correl, p ) = pearsonr(wxB[nonan], wlogk0[nonan])
    text2Bis = 'ln(koff°/koff_SPR) = '+"({:1.2f}".format(abs(pEq[1]))+"-xB)/{:1.2f}".format(abs(pEq[0]))+'\nCorrelation '+"{:1.2f}".format(correl)
    nonanF = ~np.isnan(wxBF) * ~np.isnan(wlogk0F)
    pEqF, pcovEqF = curve_fit(Affine, wxBF[nonanF], wlogk0F[nonanF] - np.log(koffSPR), p0=[-1, 1])
    wlogk0F_fit = Affine(wxBfit, pEqF[0], pEqF[1]) + np.log(koffSPR)
    #ax2Bis.plot(wxBfit, np.exp(wlogk0F_fit), linestyle='-.', c='k', label='Fit', alpha=0.5)
    (correlF, pF ) = pearsonr(wxBF[nonanF], wlogk0F[nonanF])
    text2Bis = text2Bis+ '\nln(koff°/koff_SPR) = '+"({:1.2f}".format(abs(pEqF[1]))+"-xB)/{:1.2f}".format(abs(pEqF[0]))+'\nCorrelationF '+"{:1.2f}".format(correlF)
    #ax2Bis.text(1., 0.05, text2Bis, size=10)
    ax2.plot(kBT/wxBfit, np.exp(wlogk0_fit), linestyle='-.', c='k', label='Fit', alpha=0.5)
    ax4.text(1.5, 13, TextCorrel(wfrange, wfB), size=10)
    ax2Ter.text(1.5, 0.005, TextCorrel(wxB, wkx), size=10)
    #ax2Ter1.text(1.5, 700, TextCorrel(wxB, wlrange), size=10)
    ax2Ter1.text(1.5, 700, TextCorrel(wxB, wlmax), size=10)
    ax2Ter2.text(1.5, 0.04, TextCorrel(wfB, wkz), size=10)
    ax2Ter3.text(0.003, 0.01, TextCorrel(wkx, wkz), size=10)
    ax2Ter4.text(30, 0.01, TextCorrel(wLp, wkz), size=10)
    ax2Ter5.text(1.5, 30, TextCorrel(wfB, wLp), size=10)
    ax2Ter6.text(1400, 50, TextCorrel(wL0, wLp), size=10)
    ax2Ter7.text(800, 125, TextCorrel(wlrange, wdzrange), size=10)
else:
    print('No enough common data on xB and k0')


for ax2i in [ax2, ax2Bis]:
    if not rapa: ax2i.axhline(y=1.2e-4, color='r', linestyle='-.', label='BLI')
    ax2i.axhline(y=koffSPR, color='r', linestyle='--', label='SPR')
    ax2i.set_ylabel('koff at zero force (1/s)');
    ax2i.set_ylim(1e-7,1e-1); ax2i.set_yscale('log')
    ax2i.set_xlim(0,10)
    if rapa: ax2i.set_ylim(0.01,0.4); ax2.set_xlim(0,MaxForce+5)
ax2.set_xlabel('Bell Force (pN)')
ax2Bis.set_xlim(0,1.5); ax2Bis.set_xlabel('Distance Barrier (nm)')
ax2Ter.set_xlim(0,3); ax2Ter.set_xlabel('Distance Barrier (nm)'); ax2Ter.set_ylabel('kx DNA (pN/nm)')
#ax2Ter1.set_xlim(0,3); ax2Ter1.set_xlabel('Distance Barrier (nm)'); ax2Ter1.set_ylabel('JDNA max length (nm)'); ax2Ter1.set_ylim(600,2200)
ax2Ter1.set_xlim(0,1); ax2Ter1.set_xlabel('Distance Barrier (nm)'); ax2Ter1.set_ylabel('JDNA max closed length (nm)'); ax2Ter1.set_ylim(600,2200)
ax2Ter2.set_xlim(0,MaxForce+5); ax2Ter2.set_xlabel('Bell Force (pN)'); ax2Ter2.set_ylabel('kz DNA linear fit (pN/nm)')
ax2Ter3.set_xlim(0,0.006); ax2Ter3.set_xlabel('kx DNA (pN/nm)'); ax2Ter3.set_ylabel('kz DNA linear fit (pN/nm)')
ax2Ter4.set_xlim(0,60); ax2Ter4.set_xlabel('DNA persistence length (nm)'); ax2Ter4.set_ylabel('kz DNA linear fit (pN/nm)')
ax2Ter6.set_xlim(1000,2300); ax2Ter6.set_xlabel('DNA Full Length (nm)'); ax2Ter6.set_ylabel('DNA persistence length (nm)'); ax2Ter6.set_ylim(0,30)
ax2Ter5.set_xlim(0,MaxForce+5); ax2Ter5.set_xlabel('Bell Force (pN)'); ax2Ter5.set_ylabel('DNA persistence length (nm)'); ax2Ter5.set_ylim(0,60)
#ax2Ter7.set_xlim(650,1800); ax2Ter7.set_xlabel('Open JDNA length (nm)'); ax2Ter7.set_ylabel('Delta z jump (nm)'); ax2Ter7.set_ylim(100,300)
ax2Ter7.set_xlim(1000,2300); ax2Ter7.set_xlabel('Closed JDNA length (nm)'); ax2Ter7.set_ylabel('Delta z jump (nm)'); ax2Ter7.set_ylim(100,300)
ax4.set_ylim(0,MaxForce+5); ax4.set_xlabel('Bell Force (pN)'); ax4.set_ylabel('Force range (pN)'); ax4.set_xlim(0,MaxForce+5)
axRR.set_ylim(0,2); axRR.set_xlabel('Radius from diffusion'); axRR.set_ylabel('Slope F vs P ^(1/3)'); axRR.set_xlim(600,1100)
axRR1.set_ylim(500,1800); axRR1.set_xlabel('Radius from diffusion'); axRR1.set_ylabel('Average Length (nm)'); axRR1.set_xlim(600,1100)

if SaveGraph: plt.figure(figMat); plt.savefig(outpath+figMat, transparent=False)
if SaveGraph: plt.figure(figOtheFlyDx); plt.savefig(outpath+figOtheFlyDx, transparent=False)
if SaveGraph: plt.figure(figOtheFlyDy); plt.savefig(outpath+figOtheFlyDy, transparent=False)
if SaveGraph: plt.figure(figOtheFlykx); plt.savefig(outpath+figOtheFlykx, transparent=False)
if SaveGraph: plt.figure(figOtheFlyky); plt.savefig(outpath+figOtheFlyky, transparent=False)
if SaveGraph: plt.figure(figOtheFlySDx); plt.savefig(outpath+figOtheFlySDx, transparent=False)
if SaveGraph: plt.figure(figOtheFlySDy); plt.savefig(outpath+figOtheFlySDy, transparent=False)
if SaveGraph: plt.figure(figXYk); plt.savefig(outpath+figXYk, transparent=False)
if SaveGraph: plt.figure(figXYD); plt.savefig(outpath+figXYD, transparent=False)
if SaveGraph: plt.figure(figXYF); plt.savefig(outpath+figXYF, transparent=False)

if SaveGraph: plt.figure(figContour); plt.savefig(outpath+figContour, transparent=False)

if SaveGraph: plt.figure(fig4); plt.savefig(outpath+fig4, transparent=False)
if SaveGraph: plt.figure(figRR); plt.savefig(outpath+figRR, transparent=False)
if SaveGraph: plt.figure(figRR1); plt.savefig(outpath+figRR1, transparent=False)


##################    FORCE COMPARISON     ############################
Fchoice = 'Ffit'  # 'Force_pN'
if rapa:  
    zipped = list(zip(wfov, wCoefPow_F, wCoefPow_Ffit, wlmax, wLp, wL0, wxB, wkz, wxBp, warange, wlogk0, wlogk0p))
    listcol = ['Coef F_Pow','Coef Ffit_Pow', 'Lmax', 'Lp', 'L0', 'xB', 'kz', 'xBp', 'Avg angle', 'logk0', 'logk0p']
    listymin = [0, 0, 0, 0, 0, 0, 0, 0, 0, -5, -5]
    listymax = [2.5, 2.5, 2000, 50, 2000, 1, 0.05, 1, 90, -2, -2]
    dfres= pd.DataFrame(zipped, columns=['FOV']+listcol)
    dfres = dfres[(dfres['FOV']>0) & (dfres['Coef F_Pow']>0)]
    for icol, col in enumerate(listcol):
        fignamecol = col+' by FOV'; fig = plt.figure(fignamecol); axcol=plt.gca()
        sns.swarmplot(x="FOV", y=col, hue='Avg angle', data=dfres, ax=axcol)
        axcol.legend(loc="upper right", fontsize=3); axcol.get_legend().remove()
        axcol.set_ylim(listymin[icol], listymax[icol])
        if SaveGraph: plt.savefig(outpath+fignamecol, transparent=False)

gc.collect()
##################    POOL BY FORCE AND FITTING    ########################

avgf=[]; avgk=[]; sdk=[]
#listf = deltaf*(np.arange(int((MaxForce-5)/deltaf))+0.5)
listf = MinForce + deltaf*(np.arange(int((MaxForce-MinForce)/deltaf)))
for f in listf:
    dfstagenonull = dfstage[~dfstage.isnull()]
    dfx = dfstage[ (dfstage[Fchoice]<f+deltaf/2) & (dfstage[Fchoice]>f-deltaf/2) & (dfstage['newkoff']).notnull() ]
    wf = dfx[Fchoice].values; wk = dfx['newkoff'].values
    wf_ = wf.mean()
    wk_ = np.exp((np.log(wk)).mean())
    wdk_ = wk.std()
    avgf.append(wf_); avgk.append(wk_); sdk.append(wdk_)
    print('f=',f, 'wf_=', wf_, 'wk_=', wk_)
#(pEq, eEq, wlEq, res) = fitBellSingle(avgf, avgk)
(pEq, eEq, res) = fitBellFriddle(logBellxB, avgf, avgk)#, 1/sd)
wlEq = BellxB(np.array(avgf), pEq[0], pEq[1])
k0pool = pEq[0]; dk0pool = eEq[0]; xBpool = pEq[1]; dxBpool = eEq[1]
f0pool = kBT/xBpool; df0pool = kBT*dxBpool/xBpool**2
note = 'Fit Pool:  k0=  %.6f' % k0pool + '+/-  %.6f' % dk0pool + ' xB= %.3f'% xBpool + '+/- %.3f' % dxBpool
print(note)
if rapa: ax.annotate(note, xy=(0.5,0.25))
else: ax.annotate(note, xy=(0.5,0.1))
ax.errorbar(avgf, avgk, sdk, fmt='o', c='k')
ax.errorbar(avgf, wlEq, c='k', fmt='--', label='Fit Pool')
ax.legend(loc="lower right", fontsize=5)
ax2Bis.errorbar(xBpool, k0pool, dk0pool, dxBpool, fmt='o', c='k', label='Fit Pool')
ax2Bis.legend(loc="upper right", fontsize=5)
ax2.errorbar(f0pool, k0pool, dk0pool, df0pool, fmt='o', c='k', label='Fit Pool')
for ax2i in [ax2, ax4, ax2Ter, ax2Ter1, ax2Ter2, ax2Ter3, ax2Ter4, ax2Ter5, ax2Ter6, ax2Ter7, axRR, axRR1]: ax2i.legend(loc="lower right", fontsize=5)

for fignamei in [figname2, figname2Bis, figname2Ter, figname2Ter1, figname2Ter2, figname2Ter3, figname2Ter4, figname2Ter5, figname2Ter6,
                 figname2Ter7, figname, figname+'Bis', figname+'BisAB', 'Power_Lengths', 'Force Lengths', 'Force_Lengths_Bis', 'Power_Angles', 'Power Plateau', 'Force Angles', 'No_Close_Fraction',
                 'Force Lumicks', 'Power vs Force', 'Power_vs_Force_Bis', 'Deltaz vs Force', 'Deltaz vs Length', 'Power Stiffness', 'Power Stiffness Bis',
                 'Power_Stiffness_Ter', 'Power Friction', 'Power_Friction_Bis']:
    if SaveGraph: plt.figure(fignamei); plt.savefig(outpath+fignamei, transparent=False)


df0.to_csv(outpath+'Compiled.csv')
dfstage.to_csv(outpath+'CompletedStage.csv')