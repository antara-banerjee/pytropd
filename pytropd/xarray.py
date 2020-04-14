#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 19:31:11 2020

@author: abanerjee
"""

import time as pytime

''' xarray, testing speed
   # try xarray
   t0 = pytime.time(); print(t0)
   
   #f_Z = xr.open_dataset('/Volumes/CESM-GLENS/GLENS/b.e15.B5505C5WCCML45BGCR.f09_g16.feedback.004/atm/proc/tseries/month_1/Combined/p.e15.B5505C5WCCML45BGCR.f09_g16.feedback.004.cam.h0.U.202001-209912.nc', cache=False)
   f_Z = xr.open_dataset('/Volumes/Data-Banerjee3TB/CESM-GLENS/GLENS/b.e15.B5505C5WCCML45BGCR.f09_g16.feedback.010/atm/proc/tseries/month_1/Combined/p.e15.B5505C5WCCML45BGCR.f09_g16.feedback.010.cam.h0.U.202001-209912.nc', cache=False)
   #Z = f_Z['U']
   #Z = Z.mean('lon')
   
   #print(Z)
   
   f_Z.close()
   t1 = pytime.time(); print(t1)
   
   dt = t1-t0; print(dt)
'''