#from __future__ import division
import numpy as np
import scipy.stats as ss
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as patches
import sys
import TropD_GLENS

mpl.rcParams['pdf.fonttype']=42
mpl.rcParams['ps.fonttype']=42
mpl.rcParams['font.sans-serif'] = "Arial"
mpl.rcParams['font.family'] = "Arial"
mpl.rcParams['mathtext.fontset'] = 'stixsans'

# Example codes for using the TropD package to calculate tropical width metrics
# The code assumes that the current directory is ... tropd/pytropd/
## Set display and meta parameters
y1 = 2020
y2 = 2099
seas = 'MAM'
    
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

       #print(' box statistics ', onehalfIQR, lextreme, uextreme, np.percentile(ensvals,25)-onehalfIQR, np.percentile(ensvals,75)+onehalfIQR)

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
fig, (ax1, ax2) = plt.subplots(nrows=2,ncols=1, sharex=True, figsize=(7,10))

#***************************************
def common_axis_attributes(*axis):
    
    for ax in axis:
        ax.set_ylim([-0.5,0.5]) 
        ax.set_xlim([0,9])
        ax.axhline(y=0, color='k')
        ax.set_xticks(range(1,9))
        ax.set_xticklabels(['PSI','TPB','OLR','STJ','EDJ','P-E','UAS','PSL'])
        ax.tick_params(axis='both', which='major', labelsize=16)
        
    return

#***************************************
def add_boxplot(data, xpos, ax, metric):
    
    #box
    rect = patches.Rectangle((xpos-0.4,data['lperc']), width=0.8, height=data['uperc']-data['lperc'], linewidth=1, facecolor='b', fill=True, alpha=0.8, edgecolor='k')
    ax.add_patch(rect)
    ax.plot([xpos-0.4,xpos+0.4], [data['median'],data['median']], color='k', linewidth=1.5) 
    # whiskers
    ax.plot([xpos,xpos],[data['uperc'],data['uextreme']], linestyle='-', color='k', linewidth=0.5) 
    ax.plot([xpos,xpos],[data['lperc'],data['lextreme']], linestyle='-', color='k', linewidth=0.5) 
    ax.plot([xpos-0.2,xpos+0.2],[data['uextreme'],data['uextreme']], linestyle='-', color='k', linewidth=0.5) # whiskers bars
    ax.plot([xpos-0.2,xpos+0.2],[data['lextreme'],data['lextreme']], linestyle='-', color='k', linewidth=0.5) # whiskers bars
    # outliers
    outlier_points = data['outliers']
    ax.scatter([xpos]*len(outlier_points), outlier_points, facecolors='none', edgecolors='k', marker='o', s=10)
    #ax.set_xticks([xpos])
    #ax.set_xticklabels(xpos, metric)
    
#*****************************************************************************
# 1) EDJ
PSI_slopes_NH = []
PSI_slopes_SH = []
for i in range(1,21):
    fname = '/Volumes/Data-Banerjee3TB/CESM-GLENS/GLENS/b.e15.B5505C5WCCML45BGCR.f09_g16.feedback.0'+str(i).zfill(2)+'/atm/proc/tseries/month_1/Combined/p.e15.B5505C5WCCML45BGCR.f09_g16.feedback.0'+str(i).zfill(2)+'.cam.h0zm.V.202001-209912.netCDF3.nc'
    slope_NH, slope_SH = TropD_GLENS.PSI(fname, y1, y2, seas)
    PSI_slopes_NH.append(slope_NH*10)
    PSI_slopes_SH.append(slope_SH*10)

PSI_boxstats_NH = boxstats(PSI_slopes_NH)
PSI_boxstats_SH = boxstats(PSI_slopes_SH)

add_boxplot(PSI_boxstats_NH, 1, ax1, 'PSI')
add_boxplot(PSI_boxstats_SH, 1, ax2, 'PSI')

#*****************************************************************************
# 3) UAS
OLR_slopes_NH = []
OLR_slopes_SH = []
for i in range(1,21):
    fnameOLR = '/Volumes/Data-Banerjee3TB/CESM-GLENS/GLENS/b.e15.B5505C5WCCML45BGCR.f09_g16.feedback.0'+str(i).zfill(2)+'/atm/proc/tseries/month_1/Combined/b.e15.B5505C5WCCML45BGCR.f09_g16.feedback.0'+str(i).zfill(2)+'.cam.h0zm.FLNT.202001-209912.netCDF3.nc'
    fnameOLRCS = '/Volumes/Data-Banerjee3TB/CESM-GLENS/GLENS/b.e15.B5505C5WCCML45BGCR.f09_g16.feedback.0'+str(i).zfill(2)+'/atm/proc/tseries/month_1/Combined/b.e15.B5505C5WCCML45BGCR.f09_g16.feedback.0'+str(i).zfill(2)+'.cam.h0zm.FLNTC.202001-209912.netCDF3.nc'
    slope_NH, slope_SH = TropD_GLENS.OLR(fnameOLR, fnameOLRCS, y1, y2, seas)
    OLR_slopes_NH.append(slope_NH*10)
    OLR_slopes_SH.append(slope_SH*10)

OLR_boxstats_NH = boxstats(OLR_slopes_NH)
OLR_boxstats_SH = boxstats(OLR_slopes_SH)

add_boxplot(OLR_boxstats_NH, 3, ax1, 'OLR')
add_boxplot(OLR_boxstats_SH, 3, ax2, 'OLR')

#*****************************************************************************
# 4) STJ
STJ_slopes_NH = []
STJ_slopes_SH = []
for i in range(1,21):
    fname = '/Volumes/Data-Banerjee3TB/CESM-GLENS/GLENS/b.e15.B5505C5WCCML45BGCR.f09_g16.feedback.0'+str(i).zfill(2)+'/atm/proc/tseries/month_1/Combined/p.e15.B5505C5WCCML45BGCR.f09_g16.feedback.0'+str(i).zfill(2)+'.cam.h0zm.U.202001-209912.nc'
    slope_NH, slope_SH = TropD_GLENS.STJ(fname, y1, y2, seas)
    STJ_slopes_NH.append(slope_NH*10)
    STJ_slopes_SH.append(slope_SH*10)

STJ_boxstats_NH = boxstats(STJ_slopes_NH)
STJ_boxstats_SH = boxstats(STJ_slopes_SH)

add_boxplot(STJ_boxstats_NH, 4, ax1, 'STJ')
add_boxplot(STJ_boxstats_SH, 4, ax2, 'STJ')

#*****************************************************************************
# 5) EDJ
EDJ_slopes_NH = []
EDJ_slopes_SH = []
for i in range(1,21):
    fname = '/Volumes/Data-Banerjee3TB/CESM-GLENS/GLENS/b.e15.B5505C5WCCML45BGCR.f09_g16.feedback.0'+str(i).zfill(2)+'/atm/proc/tseries/month_1/Combined/p.e15.B5505C5WCCML45BGCR.f09_g16.feedback.0'+str(i).zfill(2)+'.cam.h0zm.U.202001-209912.nc'
    slope_NH, slope_SH = TropD_GLENS.EDJ(fname, y1, y2, seas)
    EDJ_slopes_NH.append(slope_NH*10)
    EDJ_slopes_SH.append(slope_SH*10)

EDJ_boxstats_NH = boxstats(EDJ_slopes_NH)
EDJ_boxstats_SH = boxstats(EDJ_slopes_SH)

add_boxplot(EDJ_boxstats_NH, 5, ax1, 'EDJ')
add_boxplot(EDJ_boxstats_SH, 5, ax2, 'EDJ')

#*****************************************************************************
# 7) UAS
UAS_slopes_NH = []
UAS_slopes_SH = []
for i in range(1,21):
    fnameU = '/Volumes/Data-Banerjee3TB/CESM-GLENS/GLENS/b.e15.B5505C5WCCML45BGCR.f09_g16.feedback.0'+str(i).zfill(2)+'/atm/proc/tseries/month_1/Combined/p.e15.B5505C5WCCML45BGCR.f09_g16.feedback.0'+str(i).zfill(2)+'.cam.h0zm.U.202001-209912.nc'
    fnameUAS = '/Volumes/Data-Banerjee3TB/CESM-GLENS/GLENS/b.e15.B5505C5WCCML45BGCR.f09_g16.feedback.0'+str(i).zfill(2)+'/atm/proc/tseries/month_1/Combined/b.e15.B5505C5WCCML45BGCR.f09_g16.feedback.0'+str(i).zfill(2)+'.cam.h0zm.U10.202001-209912.netCDF3.nc'
    slope_NH, slope_SH = TropD_GLENS.UAS(fnameU, fnameUAS, y1, y2, seas)
    UAS_slopes_NH.append(slope_NH*10)
    UAS_slopes_SH.append(slope_SH*10)

UAS_boxstats_NH = boxstats(UAS_slopes_NH)
UAS_boxstats_SH = boxstats(UAS_slopes_SH)

add_boxplot(UAS_boxstats_NH, 7, ax1, 'UAS')
add_boxplot(UAS_boxstats_SH, 7, ax2, 'UAS')

#*****************************************************************************
# 8) PSL
PSL_slopes_NH = []
PSL_slopes_SH = []
for i in range(1,21):
    fname = '/Volumes/Data-Banerjee3TB/CESM-GLENS/GLENS/b.e15.B5505C5WCCML45BGCR.f09_g16.feedback.0'+str(i).zfill(2)+'/atm/proc/tseries/month_1/Combined/b.e15.B5505C5WCCML45BGCR.f09_g16.feedback.0'+str(i).zfill(2)+'.cam.h0zm.PSL.202001-209912.netCDF3.nc'
    slope_NH, slope_SH = TropD_GLENS.PSL(fname, y1, y2, seas)
    PSL_slopes_NH.append(slope_NH*10)
    PSL_slopes_SH.append(slope_SH*10)

PSL_boxstats_NH = boxstats(PSL_slopes_NH)
PSL_boxstats_SH = boxstats(PSL_slopes_SH)

add_boxplot(PSL_boxstats_NH, 8, ax1, 'PSL')
add_boxplot(PSL_boxstats_SH, 8, ax2, 'PSL')

#*****************************************************************************

common_axis_attributes(ax1, ax2)
plt.suptitle(seas, weight='bold', fontsize=18, y=0.95, color='r')
ax1.set_title('Northern Hemisphere', fontsize=18)
ax2.set_title('Southern Hemisphere', fontsize=18)
plt.savefig('GLENS_metrics_'+seas+'.png', dpi=300)
#plt.tight_layout()
plt.show()

sys.exit()
