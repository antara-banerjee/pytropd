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

EDJ_slopes_NH = []
EDJ_slopes_SH = []

for i in range(1,21):
    fname = '/Volumes/Data-Banerjee3TB/CESM-GLENS/GLENS/b.e15.B5505C5WCCML45BGCR.f09_g16.feedback.0'+str(i).zfill(2)+'/atm/proc/tseries/month_1/Combined/p.e15.B5505C5WCCML45BGCR.f09_g16.feedback.0'+str(i).zfill(2)+'.cam.h0zm.U.202001-209912.nc'
    slope_NH, slope_SH = TropD_GLENS.EDJ(fname, y1, y2, seas)
    EDJ_slopes_NH.append(slope_NH*10)
    EDJ_slopes_SH.append(slope_SH*10)

EDJ_boxstats_NH = boxstats(EDJ_slopes_NH)
EDJ_boxstats_SH = boxstats(EDJ_slopes_SH)

fig, ax = plt.subplots()
rect = patches.Rectangle((1,EDJ_boxstats_SH['lperc']), width=1, height=EDJ_boxstats_SH['uperc']-EDJ_boxstats_SH['lperc'], linewidth=1, facecolor='b', fill=True, alpha=0.8, edgecolor='k')
ax.add_patch(rect)

## one rectangle to rule them all
#rect = patches.Rectangle((ipos-width/2.,boxstats['lperc']), width=width, height=boxstats['uperc']-boxstats['lperc'], linewidth=1, facecolor=color,fill=fill, alpha=0.8, edgecolor=edgecolor)
#ax.add_patch(rect)

## median - dashed and solid based on significance
#if boxstats['lowerCI']<0 and boxstats['upperCI']>0: 
#   ax.plot([ipos-width/2,ipos+width/2], [boxstats['median'],boxstats['median']], color='k', linewidth=1.5, linestyle='--') 
#else:
#   ax.plot([ipos-width/2,ipos+width/2], [boxstats['median'],boxstats['median']], color='k', linewidth=1.5) 

ax.plot([1,2], [EDJ_boxstats_SH['median'],EDJ_boxstats_SH['median']], color='k', linewidth=1.5) 

#print 'CONFIDENCE INTERVALS ON MEDIAN: ', boxstats['lowerCI'], boxstats['upperCI']
#ax.plot([ipos-width/2,ipos+width/2], [boxstats['lowerCI'],boxstats['lowerCI']], color='k', linewidth=2) 
#ax.plot([ipos-width/2,ipos+width/2], [boxstats['upperCI'],boxstats['upperCI']], color='k', linewidth=2) 

# whiskers
ax.plot([1.5,1.5],[EDJ_boxstats_SH['uperc'],EDJ_boxstats_SH['uextreme']], linestyle='-', color='k', linewidth=0.5) 
ax.plot([1.5,1.5],[EDJ_boxstats_SH['lperc'],EDJ_boxstats_SH['lextreme']], linestyle='-', color='k', linewidth=0.5) 
#ax.plot([1.4,1.6],[boxstats['uextreme'],boxstats['uextreme']], linestyle='-', color='k', linewidth=0.5) # whiskers bars
#ax.plot([1.4,1.6],[boxstats['lextreme'],boxstats['lextreme']], linestyle='-', color='k', linewidth=0.5) # whiskers bars

# outliers
outlier_points = EDJ_boxstats_SH['outliers']
ax.scatter([1.5]*len(outlier_points), outlier_points, facecolors='none', edgecolors='k', marker='o', s=10)
ax.set_ylim([-0.5,0.5]) 
ax.set_xlim([0,10])
ax.axhline(y=0, color='k')
plt.show()

# number in each rectangle
#ax.text(ipos+width/2., -0.4, str(boxstats['nsign'][0]), fontsize=14, color='b')
#ax.text(ipos+width/2., 0.1, str(boxstats['nsign'][1]), fontsize=14, color='r') 

sys.exit()
