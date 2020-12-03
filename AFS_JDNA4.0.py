#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 15:40:32 2020

@author: laurent
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
20201203 v4.0   Full revision of code
                Adapt to last version of nptdms package
