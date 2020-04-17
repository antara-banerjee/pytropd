#from __future__ import division
import os
import numpy as np
import scipy.stats as ss
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.patches as patches
import sys
import TropD_GLENS

# Example codes for using the TropD package to calculate tropical width metrics
# The code assumes that the current directory is ... tropd/pytropd/
## Set display and meta parameters
y1 = 2020
y2 = 2099
seas = 'SON'
    
#*****************************************************************************
def boxstats(x):

   def confidence_interval(ensvals):

       alpha = 0.05
       ensvals = np.array(ensvals)
       mean = np.mean(ensvals)
       n = np.shape(ensvals)[0]
       se = np.std(ensvals)/np.sqrt(float(n))
       tcrit = ss.t.ppf(1-(alpha/2.), n-1)
       ci = tcrit*se

       return (mean, ci)

   def percentile(ensvals, perc=25):

       ensvals = np.array(ensvals)
       median = np.percentile(ensvals, 50)
       lperc = np.percentile(ensvals, perc)
       uperc = np.percentile(ensvals, 100-perc)

       IQR = np.percentile(ensvals,75) - np.percentile(ensvals,25)
       onehalfIQR = 1.5*IQR

       data_sorted = np.sort(ensvals)
       #print data_sorted
       #print len(data_sorted)
       for i in range(len(data_sorted)):
           if data_sorted[i] >= (np.percentile(ensvals,25)-onehalfIQR):
               lextreme = data_sorted[i]
               break
       for i in range(len(data_sorted)):
           if data_sorted[i] > (np.percentile(ensvals,75)+onehalfIQR):
               uextreme = data_sorted[i-1]
               break
           # check in case there are no top outliers
           if i==(len(data_sorted)-1):
               uextreme = data_sorted[i]
               break

       print(data_sorted)
       n = len(data_sorted)
       ilowerCI = int(round(n/2. - 1.96/2*np.sqrt(n)))-1
       iupperCI = int(round(n/2. + 1.96/2*np.sqrt(n)))-1
       lowerCI = data_sorted[ilowerCI] 
       upperCI = data_sorted[iupperCI] 
       print(' confidence intervals ', lowerCI, median, upperCI)

       print(' box statistics ', onehalfIQR, lextreme, uextreme, np.percentile(ensvals,25)-onehalfIQR, np.percentile(ensvals,75)+onehalfIQR)

       return (median, lperc, uperc, onehalfIQR, lextreme, uextreme, lowerCI, upperCI)

   def variance(ensvals, sigma=1.6):

       ensvals = np.array(ensvals)
       mean = np.mean(ensvals)
       sd = np.std(ensvals)
       lrange = mean - (sigma*sd) 
       urange = mean + (sigma*sd) 

       return (lrange, urange)

   def calc_outliers(ensvals, accept_range):

      outliers = []

      for ens in ensvals:
          if (ens < accept_range[0] or ens > accept_range[-1]):
             outliers.append(ens)

      return outliers

   def calc_nsign(ensvals):

      npos = 0 
      nneg = 0 

      for ens in ensvals:
          if (ens < 0):
             nneg+=1
          else:
             npos+=1

      return (nneg, npos)

   # confidence interval
   mean, CI = confidence_interval(x)

   median, lperc, uperc, onehalfIQR, lextreme, uextreme, lowerCI, upperCI = percentile(x)

   lrange, urange = variance(x)

   outliers = calc_outliers(x, [lextreme,uextreme]) 

   nneg, npos = calc_nsign(x) 

   boxstats = {'mean':mean,\
               'CI':CI,\
               'median':median,\
               'lperc':lperc,\
               'uperc':uperc,\
               'onehalfIQR':onehalfIQR,\
               'lextreme':lextreme,\
               'uextreme':uextreme,\
               'lowerCI':lowerCI,\
               'upperCI':upperCI,\
               'lrange':lrange,\
               'urange':urange,\
               'outliers':outliers,\
               'nneg':nneg,\
               'npos':npos}
               
   return boxstats
#*****************************************************************************

# FIGURE
fig, (ax1, ax2) = plt.subplots(nrows=2,ncols=1, sharex=True, figsize=(6,6))

def common_axis_attributes(*axis):
    
    for ax in axis:
        ax.set_ylim([-0.5,0.5]) 
        ax.set_xlim([0,10])
        ax.axhline(y=0, color='k')
        
    return

#def add_boxplot(data, xpos, axis):

#*****************************************************************************
# 5) EDJ
EDJ_slopes_NH = []
EDJ_slopes_SH = []
for i in range(1,21):
    print(i)
    fname = '/Volumes/Data-Banerjee3TB/CESM-GLENS/GLENS/b.e15.B5505C5WCCML45BGCR.f09_g16.feedback.0'+str(i).zfill(2)+'/atm/proc/tseries/month_1/Combined/p.e15.B5505C5WCCML45BGCR.f09_g16.feedback.0'+str(i).zfill(2)+'.cam.h0zm.U.202001-209912.nc'
    slope_NH, slope_SH = TropD_GLENS.EDJ(fname, y1, y2, seas)
    EDJ_slopes_NH.append(slope_NH*10)
    EDJ_slopes_SH.append(slope_SH*10)

EDJ_boxstats_NH = boxstats(EDJ_slopes_NH)
EDJ_boxstats_SH = boxstats(EDJ_slopes_SH)

#box
rect = patches.Rectangle((1,EDJ_boxstats_SH['lperc']), width=1, height=EDJ_boxstats_SH['uperc']-EDJ_boxstats_SH['lperc'], linewidth=1, facecolor='b', fill=True, alpha=0.8, edgecolor='k')
ax2.add_patch(rect)
ax2.plot([1,2], [EDJ_boxstats_SH['median'],EDJ_boxstats_SH['median']], color='k', linewidth=1.5) 
# whiskers
ax2.plot([1.5,1.5],[EDJ_boxstats_SH['uperc'],EDJ_boxstats_SH['uextreme']], linestyle='-', color='k', linewidth=0.5) 
ax2.plot([1.5,1.5],[EDJ_boxstats_SH['lperc'],EDJ_boxstats_SH['lextreme']], linestyle='-', color='k', linewidth=0.5) 
ax2.plot([1.4,1.6],[EDJ_boxstats_SH['uextreme'],EDJ_boxstats_SH['uextreme']], linestyle='-', color='k', linewidth=0.5) # whiskers bars
ax2.plot([1.4,1.6],[EDJ_boxstats_SH['lextreme'],EDJ_boxstats_SH['lextreme']], linestyle='-', color='k', linewidth=0.5) # whiskers bars
# outliers
outlier_points = EDJ_boxstats_SH['outliers']
ax2.scatter([1.5]*len(outlier_points), outlier_points, facecolors='none', edgecolors='k', marker='o', s=10)

#*****************************************************************************

common_axis_attributes(ax1, ax2)
plt.show()

sys.exit()
