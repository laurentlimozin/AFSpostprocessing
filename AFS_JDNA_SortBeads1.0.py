# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 12:23:26 2020

@author: Laurent
"""
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.pyplot import cm
plt.close('all')

dmax=1.e3; tagfontsize=7; nmaxcolor=10
inputpathPos='C:/DATA/Experiences Diverses/AFS/Jim/FOV1/'
    
readfilesPos = [f for f in os.listdir(inputpathPos) if f.endswith('.csv')]
print(len(readfilesPos),'loaded files')
df=[]
for (b,file) in enumerate(readfilesPos):
    df.append((file,pd.read_csv(inputpathPos+file)))
    print(file)
#color=iter(cm.rainbow(np.linspace(0, 1, 20)))
color=cm.rainbow(np.linspace(0, 1, nmaxcolor))

## Sorting of beads
plt.figure('MapTraj', figsize=(6,6), dpi=100); ax = plt.gca()
for i, idf in enumerate(df):
    u=idf[1]
    x=u['X0med']; y=u['Y0med']; p=u['Phigh_%']; b=u['bead']
#    print(len(x))
    if i==0:
        x[np.isnan(x)]=-1e3; y[np.isnan(y)]=-1e3
        u['beadrenum']=np.arange(len(x)); ntotb=len(x)
        xp=np.copy(x); yp=np.copy(y); up=np.arange(len(x), dtype=int)
#        print('file', i, 'bead init', ntotb)
        for j in range(len(x)):
            print('file', int(p[j]),'%', 'bead', b[j], up[j], round(xp[j]), round(yp[j]))
 #           print('file', i, 'bead', j, b[j], up[j], round(xp[j]), round(yp[j]))
            print('   new bead', int(u['beadrenum'][j]))
            ax.scatter(x[j], y[j], s=10, c=color[j])
            lab=str(i)+' / '+str(int(u['beadrenum'][j]))+' / '+str(int(b[j]))+' / '+str(int(p[j]))+'%'
            lab=str(int(p[j]))+'%'+' / '+str(int(b[j]))+' / '+str(int(u['beadrenum'][j]))
            ax.text(x[j], y[j]+5000, lab, fontsize=tagfontsize)
    elif i>0:
#        print('file', i, len(x))
        u['beadrenum']=np.zeros(len(x))-1
        for j in range(len(x)):
            print('file', int(p[j]),'%', 'bead', b[j], up[j], round(xp[j]), round(yp[j]))
 #           print('file', i, 'bead', j, b[j], up[j], round(xp[j]), round(yp[j]))
            d=np.sqrt((x[j]-xp)**2+(y[j]-yp)**2)
            dnonan=d[~np.isnan(d)]
            if len(dnonan)>0:
                jpmin=np.argmin(dnonan); dmin=np.amin(dnonan)
     #           print('    ',j, dmin, dmax)
                if dmin<dmax:
                    u['beadrenum'][j]=up[jpmin]
                    print('   old bead', int(u['beadrenum'][j]))
                else: 
                    u['beadrenum'][j]=int(ntotb)
                    xp=np.append(xp,[x[j]])
                    yp=np.append(yp,[y[j]])
                    up=np.append(up,[ntotb])
                    ntotb+=1
                    print('   new bead', int(u['beadrenum'][j]))
                ax.scatter(x[j], y[j], s=10, c=color[int(u['beadrenum'][j])])
                lab=str(i)+' / '+str(int(u['beadrenum'][j]))+' / '+str(int(b[j]))+' / '+str(int(p[j]))+'%'
                lab=str(int(p[j]))+'%'+' / '+str(int(b[j]))+' / '+str(int(u['beadrenum'][j]))
                ax.text(x[j], y[j]-(i-1)*5000, lab, fontsize=tagfontsize)
            else:
                print('Empty array')

dff=[]
for i, idf in enumerate(df): dff.append(idf[1])
if not os.path.exists(inputpathPos+'/Out/'): os.makedirs(inputpathPos+'/Out/')     
df0=pd.concat(dff)      
df0.to_csv(inputpathPos+'/Out/Compiled.csv')

## Example of Plot

def multiplots(df0, bname, uname, vname, ax, umax, vmax):
    n=df0[bname]; u=df0[uname]; v=df0[vname]
#    print(int(np.amax(n)))
    for inb in range(int(np.amax(n))):
   #     ax.scatter(u[n==inb],v[n==inb], marker='+', label=str(inb)) 
        ax.errorbar(u[n==inb],v[n==inb], np.zeros(len(u[n==inb])), label=str(inb), fmt='+-')
        ax.set_xlim(0,umax); ax.set_ylim(0,vmax)
        ax.set_xlabel(uname); ax.set_ylabel(vname)
    ax.legend()
        
plt.figure('x_vs_power', figsize=(6,6), dpi=100); ax = plt.gca()
multiplots(df0, 'beadrenum', 'Phigh_%', 'FSpecx_pN', ax, 5, 5)
plt.figure('off_vs_force', figsize=(6,6), dpi=100); ax = plt.gca()
multiplots(df0, 'beadrenum', 'FSpecx_pN', 'Offrate_1_per_s', ax, 3, 0.1)