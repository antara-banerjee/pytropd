from __future__ import division
import numpy as np
from functions import *

def TropD_Metric_EDJ(U, lat, lev=np.array([1]), method='max', n=0):
  
  '''TropD EDJ metric
  
  Latitude of the eddy driven jet (EDJ) 
  Latitude of maximum of the zonal wind at the level closest to the 850 hPa level
  Written by Ori Adam Mar.20.2017
  Edited by Alison Ming Jul.4.2017
  
  Positional arguments:
  U(lat,lev) or U (lat,)-- Zonal mean zonal wind. Also takes surface wind 
  lat -- latitude vector
  lev -- vertical level vector in hPa units

  Keyword arguments:
  method (optional) -- 'max' (default) |  'peak'
  n (optional, default = 6) -- rank of moment used to calculate the position of max value. n = 1,2,4,6,8,...  

  Outputs:
  PhiSH -- latitude of EDJ in the SH
  PhiNH -- latitude of EDJ in the NH
  '''

  try:
    assert (not hasattr(n, "__len__") and n >= 0)  
  except AssertionError:
   print 'TropD_Metric_EDJ: ERROR : the smoothing parameter n must be >= 0'
   
  try:
    assert(method in ['max','peak'])
  except AssertionError:
    print 'TropD_Metric_EDJ: ERROR : unrecognized method ',method

  eq_boundary = 15
  polar_boundary = 60
  
  if len(lev) > 1:
    u = U[:,find_nearest(lev, 850)]
  else:
    u = np.copy(U)
    
  if method=='max':
    if n:
      PhiNH = TropD_Calculate_MaxLat(u[(lat > eq_boundary) & (lat < polar_boundary)],\
              lat[(lat > eq_boundary) & (lat < polar_boundary)],n)
      PhiSH = TropD_Calculate_MaxLat(u[(lat > -polar_boundary) & (lat < -eq_boundary)],\
              lat[(lat > -polar_boundary) & (lat < -eq_boundary)],n)

    else:
      PhiNH = TropD_Calculate_MaxLat(u[(lat > eq_boundary) & (lat < polar_boundary)],\
              lat[(lat > eq_boundary) & (lat < polar_boundary)])
      PhiSH = TropD_Calculate_MaxLat(u[(lat > -polar_boundary) & (lat < -eq_boundary)],\
              lat[(lat > -polar_boundary) & (lat < -eq_boundary)])
  elif method=='peak':
    if n:
      PhiNH = TropD_Calculate_MaxLat(u[(lat > eq_boundary) & (lat < polar_boundary)],\
              lat[(lat > eq_boundary) & (lat < polar_boundary)],n)
      PhiSH = TropD_Calculate_MaxLat(u[(lat > -polar_boundary) & (lat < -eq_boundary)],\
              lat[(lat > -polar_boundary) & (lat < -eq_boundary)],n)
    else:
      PhiNH = TropD_Calculate_MaxLat(u[(lat > eq_boundary) & (lat < polar_boundary)],\
              lat[(lat > eq_boundary) & (lat < polar_boundary)],30)
      PhiSH = TropD_Calculate_MaxLat(u[(lat > -polar_boundary) & (lat < -eq_boundary)],\
              lat[(lat > -polar_boundary) & (lat < -eq_boundary)],30)
  
  else:
    print 'TropD_Metric_EDJ: ERROR: unrecognized method ',method

  return PhiSH, PhiNH



# Written by Ori Adam Mar.21.2017
# Edited by Alison Ming Jul.4.2017
"""Returns TropD Outgoing Longwave Radiation (OLR) metric
   
Parameters
----------
olr(lat,) :  zonal mean TOA olr (positive)
lat -- equally spaced latitude column vector
   
method (optional) : '250W' (default) | 'cutoff' | '10Perc' | '20W' | 'max' | 'peak'

'250W'(Default): the first latitude poleward of the tropical OLR maximum in each hemisphere where OLR crosses 250W/m^2

'20W': the first latitude poleward of the tropical OLR maximum in each hemisphere where OLR crosses the tropical OLR max minus 20W/m^2

'cutoff': the first latitude poleward of the tropical OLR maximum in each hemisphere where OLR crosses a specified cutoff value

'10Perc': the first latitude poleward of the tropical OLR maximum in each hemisphere where OLR is 10# smaller than the tropical OLR maximum

'max': the latitude of maximum of tropical olr in each hemisphere with the smoothing paramerer n=6 in TropD_Calculate_MaxLat

'peak': the latitude of maximum of tropical olr in each hemisphere with the smoothing parameter n=30 in TropD_Calculate_MaxLat


Cutoff (optional) : Scalar. For the method 'cutoff', Cutoff specifies the OLR cutoff value. 

n (optional) : For the 'max' method, n is the smoothing parameter in TropD_Calculate_MaxLat

Returns
-------

PhiSH : latitude of near equator OLR threshold crossing in the SH
PhiNH : latitude of near equator OLR threshold crossing in the NH
"""

def TropD_Metric_OLR(olr, lat, method='250W', Cutoff=50, n=int(6)):

  try:
    assert(isinstance(n, int)) 
  except AssertionError:
    print 'TropD_Metric_OLR: ERROR: the smoothing parameter n must be an integer'
  
  try:
    assert(n>=1) 
  except AssertionError:
    print 'TropD_Metric_OLR: ERROR: the smoothing parameter n must be >= 1'

  
  # make latitude vector monotonically increasing
  if lat[-1] < lat[0]:
    olr = np.flip(olr,0)
    lat = np.flip(lat,0)
    
  eq_boundary = 5
  subpolar_boundary = 40
  polar_boundary = 60
  # NH
  olr_max_lat_NH = TropD_Calculate_MaxLat(olr[(lat > eq_boundary) & (lat < subpolar_boundary)],\
                    lat[(lat > eq_boundary) & (lat < subpolar_boundary)])
  olr_max_NH = max(olr[(lat > eq_boundary) & (lat < subpolar_boundary)])

  # SH
  olr_max_lat_SH = TropD_Calculate_MaxLat(olr[(lat > -subpolar_boundary) & (lat < -eq_boundary)],\
                    lat[(lat > -subpolar_boundary) & (lat < -eq_boundary)])
  olr_max_SH = max(olr[(lat > -subpolar_boundary) & (lat < -eq_boundary)])

  if method=='20W':
    PhiNH = TropD_Calculate_ZeroCrossing(olr[(lat > olr_max_lat_NH) & (lat < polar_boundary)] - olr_max_NH + 20,\
                    lat[(lat > olr_max_lat_NH) & (lat < polar_boundary)])
    PhiSH = TropD_Calculate_ZeroCrossing(np.flip(olr[(lat < olr_max_lat_SH) & \
                    (lat > -polar_boundary)],0) - olr_max_SH + 20,\
                    np.flip(lat[(lat < olr_max_lat_SH) & (lat > -polar_boundary)],0))

  elif method=='250W':
    PhiNH = TropD_Calculate_ZeroCrossing(olr[(lat > olr_max_lat_NH) & (lat < polar_boundary)] - 250,\
                    lat[(lat > olr_max_lat_NH) & (lat < polar_boundary)])
    PhiSH = TropD_Calculate_ZeroCrossing(np.flip(olr[(lat < olr_max_lat_SH) & (lat > -polar_boundary)],0) - 250,\
                    np.flip(lat[(lat < olr_max_lat_SH) & (lat > -polar_boundary)],0))

  elif method=='cutoff':
    PhiNH = TropD_Calculate_ZeroCrossing(olr[(lat > olr_max_lat_NH) & (lat < polar_boundary)] - Cutoff,\
                    lat[(lat > olr_max_lat_NH) & (lat < polar_boundary)])
    PhiSH = TropD_Calculate_ZeroCrossing(np.flip(olr[(lat < olr_max_lat_SH) & (lat > -polar_boundary)],0) - Cutoff,\
                    np.flip(lat[(lat < olr_max_lat_SH) & (lat > -polar_boundary)],0))
  
  elif method=='10Perc':
    PhiNH = TropD_Calculate_ZeroCrossing(olr[(lat > olr_max_lat_NH) & (lat < polar_boundary)] / olr_max_NH - 0.9,\
                    lat[(lat > olr_max_lat_NH) & (lat < polar_boundary)])
    PhiSH = TropD_Calculate_ZeroCrossing(np.flip(olr[(lat < olr_max_lat_SH) & (lat > -polar_boundary)],0) \
                    / olr_max_SH - 0.9, np.flip(lat[(lat < olr_max_lat_SH) & (lat > -polar_boundary)],0))

  elif method=='max':
    if Cutoff_is_set:
      PhiNH = TropD_Calculate_MaxLat(olr[(lat > eq_boundary) & (lat < subpolar_boundary)],\
                    lat[(lat > eq_boundary) & (lat < subpolar_boundary)], n=n)
      PhiSH = TropD_Calculate_MaxLat(olr[(lat > -subpolar_boundary) & (lat < -eq_boundary)],\
                    lat[(lat > -subpolar_boundary) & (lat < -eq_boundary)], n=n)
    else:
      PhiNH = np.copy(olr_max_lat_NH)
      PhiSH = np.copy(olr_max_lat_SH)
 
  elif method=='peak':
    PhiNH = TropD_Calculate_MaxLat(olr[(lat > eq_boundary) & (lat < subpolar_boundary)],\
                    lat[(lat > eq_boundary) & (lat < subpolar_boundary)],30)
    PhiSH = TropD_Calculate_MaxLat(olr[(lat > -subpolar_boundary) & (lat < -eq_boundary)],\
                    lat[(lat > -subpolar_boundary) & (lat < -eq_boundary)],30)

  else:
    print 'TropD_Metric_OLR: unrecognized method ', method

    PhiNH = np.empty(0)
    PhiSH = np.empty(0)
  
  return PhiSH, PhiNH
    
    
def TropD_Metric_PE(pe,lat,method='zero_crossing',Lat_Uncertainty=0.0):

  '''TropD Precipitation minus Evaporation (PE) metric
  Written by Ori Adam Mar.21.2017
  Edited by Alison Ming Jul.4.2017
     
  Positional arguments:
  pe(lat,) -- zonal-mean precipitation minus evaporation
  lat -- equally spaced latitude column vector

  Keyword arguments:
  method -- 'zero_crossing': the first latitude poleward of the subtropical minimum where P-E changes from negative to positive values. Only one method so far.
  Lat_Uncertainty (optional) -- The minimal distance allowed between the first and second zero crossings along lat

  Output:
  PhiSH -- latitude of first subtropical P-E zero crossing in the SH
  PhiNH -- latitude of first subtropical P-E zero crossing in the NH
  '''    
  try:
    assert(method in ['zero_crossing'])
  except AssertionError:
    print 'TropD_Metric_PE: ERROR : unrecognized method ',method
    
  # make latitude vector monotonically increasing
  if lat[-1] < lat[0]:
      pe = np.flip(pe)
      lat = np.flip(lat)
    
  # The gradient of PE is used to determine whether PE becomes positive at the zero crossing
  ped = np.interp(lat, (lat[:-1] + lat[1:])/2.0, np.diff(pe))
    
  # define latitudes of boundaries certain regions 
  eq_boundary=5
  subpolar_boundary=50
  polar_boundary=70

    
  # NH
  M1 = TropD_Calculate_MaxLat(-pe[(lat > eq_boundary) & (lat < subpolar_boundary)],\
                 lat[(lat > eq_boundary) & (lat < subpolar_boundary)], 30)
  ZC1 = TropD_Calculate_ZeroCrossing(pe[(lat > M1) & (lat < polar_boundary)], \
                 lat[(lat > M1) & (lat < polar_boundary)], Lat_Uncertainty)
  if np.interp(ZC1, lat, ped) > 0:
    PhiNH = ZC1
  else:
    PhiNH = TropD_Calculate_ZeroCrossing(pe[(lat > ZC1) & (lat < polar_boundary)], \
                  lat[(lat > ZC1) & (lat < polar_boundary)], Lat_Uncertainty)
  
  # SH
  # flip arrays to find the most equatorward zero crossing
  M1 = TropD_Calculate_MaxLat(np.flip(-pe[(lat < -eq_boundary) & (lat > -subpolar_boundary)],0),\
                 np.flip(lat[(lat < -eq_boundary) & (lat > -subpolar_boundary)],0), 30)               
  ZC1 = TropD_Calculate_ZeroCrossing(np.flip(pe[(lat < M1) & (lat > -polar_boundary)],0), \
                 np.flip(lat[(lat < M1) & (lat > -polar_boundary)],0), Lat_Uncertainty)

  if np.interp(ZC1, lat, ped) < 0:
    PhiSH = ZC1
  else:
    PhiSH = TropD_Calculate_ZeroCrossing(np.flip(pe[(lat < ZC1) & (lat > -polar_boundary)],0), \
                  np.flip(lat[(lat < ZC1) & (lat > -polar_boundary)],0), Lat_Uncertainty)

  return PhiSH, PhiNH

def TropD_Metric_PSI(V, lat, lev, method='Psi_500', Lat_Uncertainty=0):
  '''TropD PSI metric 
  Latitude of the meridional mass streamfunction subtropical zero crossing
  Written by Ori Adam Mar.20.2017
  Edited by Alison Ming Jul.4.2017
     
  Positional arguments:
  V(lat,lev) -- zonal-mean meridional wind
  lat -- latitude vector
  lev -- vertical level vector in hPa units
  
  Keyword arguments:  
  method (optional) -- 'Psi_500' (default) |  'Psi_500_10Perc'  |  'Psi_300_700' |  'Psi_500_Int'  |  'Psi_Int'
  
  'Psi_500'{default}: Zero crossing of the stream function (Psi) at the 500hPa level
  'Psi_500_10Perc': Crossing of 10# of the extremum value of Psi in each hemisphre at the 500hPa level
  'Psi_300_700': Zero crossing of Psi vertically averaged between the 300hPa and 700 hPa levels
  'Psi_500_Int': Zero crossing of the vertically-integrated Psi at the 500 hPa level
  'Psi_Int'    : Zero crossing of the column-averaged Psi
    
  Lat_Uncertainty (optional) -- The minimal distance allowed between the first and second zero crossings. For example, for Lat_Uncertainty = 10, the function will return a NaN value if a second zero crossings is found within 10 degrees of the most equatorward zero crossing.   
  
  Outputs:
  PhiSH -- latitude of Psi zero crossing in the SH
  PhiNH -- latitude of Psi zero crossing in the NH
  '''


  try:
    assert (Lat_Uncertainty >= 0)  
  except AssertionError:
    print 'TropD_Metric_PSI: ERROR : Lat_Uncertainty must be >= 0'
  
  try:
    assert(method in ['Psi_500','Psi_500_10Perc','Psi_300_700','Psi_500_Int','Psi_Int'])
  except AssertionError:
    print 'TropD_Metric_PSI: ERROR : unrecognized method ',method
    
  subpolar_boundary=30
  polar_boundary=60
    
  Psi = TropD_Calculate_StreamFunction(V, lat, lev)
  Psi[np.isnan(Psi)]=0
  # make latitude vector monotonically increasing
  if lat[-1] < lat[0]:
      Psi = np.flip(Psi, 0)
      lat = np.flip(lat, 0)
    
  COS = np.repeat(np.cos(lat*np.pi/180), len(lev), axis=0).reshape(len(lat),len(lev))
    
  if ( method=='Psi_500' or method=='Psi_500_10Perc'):
    # Use Psi at the level nearest to 500 hPa
    P = Psi[:,find_nearest(lev, 500)]

  elif method == 'Psi_300_700':
    # Use Psi averaged between the 300 and 700 hPa level
    P = np.trapz(Psi[:,(lev <= 700) & (lev >= 300)] * COS[:,(lev <= 700) & (lev >= 300)],\
                  lev[(lev <= 700) & (lev >= 300)]*100, axis=1)

  elif method == 'Psi_500_Int':
    # Use integrated Psi from p=0 to level mearest to 500 hPa
    PPsi_temp = sp.integrate.cumtrapz(Psi*COS, lev, axis=1)
    PPsi = np.zeros(np.shape(Psi))
    PPsi[:,1:] = PPsi_temp
    P = PPsi[:,find_nearest(lev, 500)]
     
  elif method == 'Psi_Int':
    # Use vertical mean of Psi 
    P = np.trapz(Psi*COS, lev, axis=1)
  
  else:
    print 'TropD_Metric_PSI: ERROR : Unrecognized method ', method
  
    
  # 1. Find latitude of maximal (minimal) tropical Psi in the NH (SH)
  # 2. Find latitude of minimal (maximal) subtropical Psi in the NH (SH)
  # 3. Find the zero crossing between the above latitudes

  # NH
  Lmax = TropD_Calculate_MaxLat(P[(lat > 0) & (lat < subpolar_boundary)],\
                                lat[(lat > 0) & (lat < subpolar_boundary)])

  Lmin = TropD_Calculate_MaxLat(-P[(lat > Lmax) & (lat < polar_boundary)],\
                                lat[(lat > Lmax) & (lat < polar_boundary)])
  if method=='Psi_500_10Perc':
    Pmax = max(P[(lat > 0) & (lat < subpolar_boundary)])
    PhiNH = TropD_Calculate_ZeroCrossing(P[(lat > Lmax) & (lat < Lmin)] - 0.1*Pmax,\
            lat[(lat > Lmax) & (lat < Lmin)])

  else:
    PhiNH = TropD_Calculate_ZeroCrossing(P[(lat > Lmax) & (lat < Lmin)],\
            lat[(lat > Lmax) & (lat < Lmin)], Lat_Uncertainty=Lat_Uncertainty)
  
  # SH
  Lmax = TropD_Calculate_MaxLat(-P[(lat < 0) & (lat > -subpolar_boundary)],\
         lat[(lat < 0) & (lat > -subpolar_boundary)])

  Lmin = TropD_Calculate_MaxLat(P[(lat < Lmax) & (lat > -polar_boundary)],\
         lat[(lat < Lmax) & (lat > -polar_boundary)])

  if method=='Psi_500_10Perc':
    Pmin = min(P[(lat < 0) & (lat > -subpolar_boundary)])
    PhiSH = TropD_Calculate_ZeroCrossing(np.flip(P[(lat < Lmax) & (lat > Lmin)], 0) + 0.1*Pmin,\
            np.flip(lat[(lat < Lmax) & (lat > Lmin)], 0))
  else:
    PhiSH = TropD_Calculate_ZeroCrossing(np.flip(P[(lat < Lmax) & (lat > Lmin)], 0),\
            np.flip(lat[(lat < Lmax) & (lat > Lmin)], 0), Lat_Uncertainty=Lat_Uncertainty)
  return PhiSH, PhiNH

    
def TropD_Metric_PSL(ps, lat, method='max'):

  '''TropD PSL metric
  Latitude of maximum of the subtropical sea-level pressure
  Written by Ori Adam Mar.17.2017
  Edited by Alison Ming Jul.4.2017
  
  Positional arguments:
  ps(lat,) -- sea-level pressure
  lat -- equally spaced latitude column vector

  Keyword arguments:
  method (optional) -- 'max' (default) | 'peak'
  
  Outputs:
  PhiSH -- latitude of subtropical sea-level pressure maximum in the SH
  PhiNH -- latitude of subtropical sea-level pressure maximum in the NH
  '''
  try:
    assert(method in ['max','peak'])
  except AssertionError:
    print 'TropD_Metric_PSL: ERROR : unrecognized method ',method

  eq_boundary = 15
  polar_boundary = 60
    
  if method == 'max':
    PhiNH = TropD_Calculate_MaxLat(ps[(lat > eq_boundary) & (lat < polar_boundary)],\
            lat[(lat > eq_boundary) & (lat < polar_boundary)])
    PhiSH = TropD_Calculate_MaxLat(ps[(lat > -polar_boundary) & (lat < -eq_boundary)],\
            lat[(lat > -polar_boundary) & (lat < -eq_boundary)])

  elif method == 'peak':
    PhiNH = TropD_Calculate_MaxLat(ps[(lat > eq_boundary) & (lat < polar_boundary)],\
            lat[(lat > eq_boundary) & (lat < polar_boundary)], 30)
    PhiSH = TropD_Calculate_MaxLat(ps[(lat > -polar_boundary) & (lat < -eq_boundary)],\
            lat[(lat > -polar_boundary) & (lat < -eq_boundary)], 30)
  else:
    print 'TropD_Metric_PSL: ERROR: unrecognized method ',method
  
  return PhiSH, PhiNH
    

def TropD_Metric_STJ(U, lat, lev, method='adjusted', n=0):

  '''TropD Subtropical Jet (STJ) metric
  Written by Ori Adam Mar.20.2017
  Edited by Alison Ming Jul.4.2017
  
  Positional arguments:
  U(lat,lev) -- zonal mean zonal wind
  lat -- latitude vector
  lev -- vertical level vector in hPa units
  
  Keyword arguments:
  method (optional) -- 'adjusted' (default) | 'core' | 'adjusted_peak' | 'core_peak'

  'adjusted' : Latitude of maximum (smoothing parameter n=6) of the zonal wind averaged between the 100 and 400 hPa levels minus the zonal mean zonal wind at the level closes to the 850 hPa level, poleward of 10 degrees and equatorward of the Eddy Driven Jet latitude

  'adjusted_peak': Latitude of maximum (smoothing parameter n=30) of the zonal wind averaged between the 100 and 400 hPa levels minus the zonal mean zonal wind at the level closes to the 850 hPa level, poleward of 10 degrees and equatorward of the Eddy Driven Jet latitude

  'core': Latitude of maximum of the zonal wind (smoothing parameter n=6) averaged between the 100 and 400 hPa levels, poleward of 10 degrees and equatorward of 70 degrees
  
  'core_peak': Latitude of maximum of the zonal wind (smoothing parameter n=30) averaged between the 100 and 400 hPa levels, poleward of 10 degrees and equatorward of 70 degrees
  
  Outputs:
  PhiSH -- latitude of STJ in the SH
  PhiNH -- latitude of STJ in the NH
  '''

  try:
    assert (not hasattr(n, "__len__") and n >= 0)  
  except AssertionError:
    print 'TropD_Metric_STJ: ERROR : the smoothing parameter n must be >= 0'
  
  try:
    assert(method in ['adjusted','core','adjusted_peak','core_peak'])
  except AssertionError:
    print 'TropD_Metric_STJ: ERROR : unrecognized method ',method

  eq_boundary=10
  polar_boundary=70

  lev_int = lev[(lev >= 100) & (lev <= 400)]

  if (method=='adjusted' or method=='adjusted_peak'): 
    idx_850 = find_nearest(lev, 850)

    # Pressure weighted vertical mean of U minus near surface U
    if len(lev_int) > 1:
      u = np.trapz(U[:, (lev >= 100) & (lev <= 400)], lev_int, axis=1) \
          / (lev_int[-1] - lev_int[0]) - U[:,idx_850]

    else:
      u = np.mean(U[:,(lev >= 100) & (lev <= 400)], axis=1) - U[:,idx_850]

  elif (method=='core' or method=='core_peak'):
    # Pressure weighted vertical mean of U
    if len(lev_int) > 1:
      u = np.trapz(U[:, (lev >= 100) & (lev <= 400)], lev_int, axis=1) \
          / (lev_int[-1] - lev_int[0])

    else:
      u = np.mean(U[:, (lev >= 100) & (lev <= 400)], axis=1)

  else:
    print 'TropD_Metric_STJ: unrecognized method ',method
    print 'TropD_Metric_STJ: optional methods are: adjusted (default), adjusted_peak, core, core_peak'

  if method=='core':
    if n:
      PhiNH = TropD_Calculate_MaxLat(u[(lat > eq_boundary) & (lat < polar_boundary)],\
          lat[(lat > eq_boundary) & (lat < polar_boundary)], n)
      PhiSH = TropD_Calculate_MaxLat(u[(lat > -polar_boundary) & (lat < -eq_boundary)],\
          lat[(lat > -polar_boundary) & (lat < -eq_boundary)], n)
    else:
      PhiNH = TropD_Calculate_MaxLat(u[(lat > eq_boundary) & (lat < polar_boundary)],\
          lat[(lat > eq_boundary) & (lat < polar_boundary)])
      PhiSH = TropD_Calculate_MaxLat(u[(lat > -polar_boundary) & (lat < -eq_boundary)],\
          lat[(lat > -polar_boundary) & (lat < -eq_boundary)])

  elif method=='core_peak':
    if n:
      PhiNH = TropD_Calculate_MaxLat(u[(lat > eq_boundary) & (lat < polar_boundary)],\
          lat[(lat > eq_boundary) & (lat < polar_boundary)], n)
      PhiSH = TropD_Calculate_MaxLat(u[(lat > -polar_boundary) & (lat < - eq_boundary)],\
          lat[(lat > -polar_boundary) & (lat < -eq_boundary)], n)
    else:
      PhiNH = TropD_Calculate_MaxLat(u[(lat > eq_boundary) & (lat < polar_boundary)],\
          lat[(lat > eq_boundary) & (lat < polar_boundary)], 30)
      PhiSH = TropD_Calculate_MaxLat(u[(lat > -polar_boundary) & (lat < -eq_boundary)],\
          lat[(lat > -polar_boundary) & (lat < -eq_boundary)], 30)

  elif method=='adjusted':
    PhiSH_EDJ, PhiNH_EDJ = TropD_Metric_EDJ(U,lat,lev)
    if n:
      PhiNH = TropD_Calculate_MaxLat(u[(lat > eq_boundary) & (lat < PhiNH_EDJ)],\
          lat[(lat > eq_boundary) & (lat < PhiNH_EDJ)], n)
      PhiSH = TropD_Calculate_MaxLat(u[(lat > PhiSH_EDJ) & (lat < -eq_boundary)],\
          lat[(lat > PhiSH_EDJ) & (lat < -eq_boundary)], n)

    else:
      PhiNH = TropD_Calculate_MaxLat(u[(lat > eq_boundary) & (lat < PhiNH_EDJ)],\
          lat[(lat > eq_boundary) & (lat < PhiNH_EDJ)])
      PhiSH = TropD_Calculate_MaxLat(u[(lat > PhiSH_EDJ) & (lat < -eq_boundary)],\
          lat[(lat > PhiSH_EDJ) & (lat < -eq_boundary)])

  elif method=='adjusted_peak':
    PhiSH_EDJ,PhiNH_EDJ = TropD_Metric_EDJ(U,lat,lev)
    if n:
      PhiNH = TropD_Calculate_MaxLat(u[(lat > eq_boundary) & (lat < PhiNH_EDJ)],\
          lat[(lat > eq_boundary) & (lat < PhiNH_EDJ)], n)
      PhiSH = TropD_Calculate_MaxLat(u[(lat > PhiSH_EDJ) & (lat < -eq_boundary)],\
          lat[(lat > PhiSH_EDJ) & (lat < -eq_boundary)], n)
    else:
      PhiNH = TropD_Calculate_MaxLat(u[(lat > eq_boundary) & (lat < PhiNH_EDJ)],\
          lat[(lat > eq_boundary) & (lat < PhiNH_EDJ)], 30)
      PhiSH = TropD_Calculate_MaxLat(u[(lat > PhiSH_EDJ) & (lat < -eq_boundary)],\
          lat[(lat > PhiSH_EDJ) & (lat < -eq_boundary)], 30)

  return PhiSH, PhiNH    

def TropD_Metric_TPB(T, lat, lev, method='max_gradient', n=0, Z=None, Cutoff=15*1000):

  '''TropD Tropopause break (TPB) metric
  Written by Ori Adam Mar.17.2017
  Edited by Alison Ming Jul.4.2017
  
  Positional arguments:
  T(lat,lev) -- temperature (K)
  lat -- latitude vector
  lev -- pressure levels column vector in hPa

  Keyword arguments:
  method (optional) -- 'max_gradient' (default) |  'max_potemp'  |  'cutoff' 
  
  'max_gradient': The latitude of maximal poleward gradient of the tropopause height
  
  'cutoff': The most equatorward latitude where the tropopause crosses a prescribed cutoff value
  
  'max_potemp': The latitude of maximal difference between the potential temperature at the tropopause and at the surface
  
  Z(lat,lev) (optional) -- geopotential height (m)
  Cutoff (optional, scalar) -- geopotential height (m) cutoff that marks the location of the tropopause break

  Outputs:
  PhiSH -- latitude of tropopause break in the SH
  PhiNH -- latitude of tropopause break in the NH
  '''


  Rd = 287.04
  Cpd = 1005.7
  k = Rd / Cpd
  try:
    assert (not hasattr(n, "__len__") and n >= 0)  
  except AssertionError:
    print 'TropD_Metric_TPB: ERROR : the smoothing parameter n must be >= 0'

  try:
    assert(method in ['max_gradient','max_potemp','cutoff'])
  except AssertionError:
    print 'TropD_Metric_TPB: ERROR : unrecognized method ',method
  
  polar_boundary=70

  if method=='max_gradient':
    Pt = TropD_Calculate_TropopauseHeight(T,lev)
    Ptd = np.diff(Pt) / (lat[1] - lat[0])
    lat2 = (lat[1:] + lat[:-1]) / 2
    
    if (n >= 1):
      PhiNH = TropD_Calculate_MaxLat(Ptd[:,(lat2 > 0) & (lat2 < polar_boundary)],\
              lat2[(lat2 > 0) & (lat2 < polar_boundary)], n=n)
      PhiSH = TropD_Calculate_MaxLat(-Ptd[:,(lat2 > -polar_boundary) & (lat2 < 0)],\
              lat2[(lat2 > -polar_boundary) & (lat2 < 0)], n=n)
    
    else:
      PhiNH = TropD_Calculate_MaxLat(Ptd[:,(lat2 > 0) & (lat2 < polar_boundary)],\
              lat2[(lat2 > 0) & (lat2 < polar_boundary)])
      PhiSH = TropD_Calculate_MaxLat(-Ptd[:,(lat2 > -polar_boundary) & (lat2 < 0)],\
              lat2[(lat2 > -polar_boundary) & (lat2 < 0)])
     
  elif method=='max_potemp':
    XF = np.tile((lev / 1000) ** k, (len(lat), 1))
    PT = T / XF
    Pt, PTt = TropD_Calculate_TropopauseHeight(T, lev, Z=PT)
    PTdif = PTt - np.nanmin(PT, axis = 1)
    
    if (n >= 1):
      PhiNH = TropD_Calculate_MaxLat(PTdif[:,(lat > 0) & (lat < polar_boundary)],\
              lat[(lat > 0) & (lat < polar_boundary)], n=n)
      PhiSH = TropD_Calculate_MaxLat(PTdif[:,(lat > - polar_boundary) & (lat < 0)],\
              lat[(lat > -polar_boundary) & (lat < 0)], n=n)
    
    else:
      PhiNH = TropD_Calculate_MaxLat(PTdif[:,(lat > 0) & (lat < polar_boundary)],\
              lat[(lat > 0) & (lat < polar_boundary)])
      PhiSH = TropD_Calculate_MaxLat(PTdif[:,(lat > - polar_boundary) & (lat < 0)],\
              lat[(lat > -polar_boundary) & (lat < 0)])
   
  elif method=='cutoff':
    Pt, Ht = TropD_Calculate_TropopauseHeight(T, lev, Z)
    
    # make latitude vector monotonically increasing
    if lat[-1] < lat[0]:
      Ht = np.flip(np.squeeze(Ht),0)
      lat = np.flip(lat,0)
    
    polar_boundary = 60
      
    PhiNH = TropD_Calculate_ZeroCrossing(Ht[(lat > 0) & (lat < polar_boundary)] - Cutoff,
              lat[(lat > 0) & (lat < polar_boundary)])
    PhiSH = TropD_Calculate_ZeroCrossing(np.flip(Ht[(lat < 0) & (lat > -polar_boundary)], 0) - Cutoff,
              np.flip(lat[(lat < 0) & (lat > -polar_boundary)], 0))
  
  else:
    print 'TropD_Metric_TPB: ERROR : Unrecognized method ', method

  return PhiSH, PhiNH
  

def TropD_Metric_UAS(U, lat, lev=np.array([1]), method='zero_crossing', Lat_Uncertainty = 0):
  
  '''TropD near-surface zonal wind metric
  Written by Ori Adam Mar.21.2017
  Edited by Alison Ming Jul.4.2017
  
  Positional arguments:
  U(lat,lev) or U (lat,)-- Zonal mean zonal wind. Also takes surface wind 
  lat -- latitude vector
  lev -- vertical level vector in hPa units. lev=np.array([1]) for single-level input zonal wind U(lat,)

  method (optional) -- 'zero_crossing' (default)
  'zero_crossing': the first subtropical latitude where near-surface zonal wind changes from negative to positive

  Keyword arguments:
  Lat_Uncertainty (optional) -- the minimal distance allowed between the first and second zero crossings
  
  Outputs:
  PhiSH -- latitude of first subtropical zero crossing in the SH
  PhiNH -- latitude of first subtropical zero crossing in the NH
  '''

  try:
    assert (Lat_Uncertainty >= 0)  
  except AssertionError:
    print 'TropD_Metric_PSI: ERROR : Lat_Uncertainty must be >= 0'
    
  try:
    assert(method in ['zero_crossing'])
  except AssertionError:
    print 'TropD_Metric_PSI: ERROR : unrecognized method ',method
    
  if len(lev) > 1:
    uas = U[:,find_nearest(lev, 850)]
  else:
    uas = np.copy(U)
    
  # make latitude vector monotonically increasing
  if lat[-1] < lat[0]:
      uas = np.flip(uas)
      lat = np.flip(lat)

  # define latitudes of boundaries certain regions 
  eq_boundary=5
  subpolar_boundary=30
  polar_boundary=60

  # NH
  uas_min_lat_NH = TropD_Calculate_MaxLat(-uas[(lat > eq_boundary) & (lat < subpolar_boundary)],\
                   lat[(lat > eq_boundary) & (lat < subpolar_boundary)])
  # SH
  uas_min_lat_SH = TropD_Calculate_MaxLat(-uas[(lat > -subpolar_boundary) & (lat < -eq_boundary)],\
      lat[(lat > -subpolar_boundary) & (lat < -eq_boundary)])
  try:
    assert(method=='zero_crossing')
    PhiNH = TropD_Calculate_ZeroCrossing(uas[(lat > uas_min_lat_NH) & (lat < polar_boundary)],\
            lat[(lat > uas_min_lat_NH) & (lat < polar_boundary)], Lat_Uncertainty)
    # flip arrays to find the most equatorward zero crossing
    PhiSH = TropD_Calculate_ZeroCrossing(np.flip(uas[(lat < uas_min_lat_SH) & (lat > -polar_boundary)],0),\
            np.flip(lat[(lat < uas_min_lat_SH) & (lat > -polar_boundary)],0), Lat_Uncertainty)

    return PhiSH, PhiNH
  except AssertionError:
    print 'TropD_Metric_UAS: ERROR : unrecognized method ', method

  
