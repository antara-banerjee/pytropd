## Set display and meta parameters
    #y1=1979
    #y2=2016
    #time=linspace(y1,y2 + 1,dot(12,(y2 - y1 + 1)) + 1)
    #time=(time(arange(1,end() - 1)) + time(arange(2,end()))) / 2
    #red_color=concat([1,0.3,0.4])
    #orange_color=concat([255,140,0]) / 256
    #blue_color=concat([0,0.447,0.741])
    #purple_color=concat([0.494,0.184,0.556])
    #green_color=concat([0.466,0.674,0.188])
    #lightblue_color=concat([0.301,0.745,0.933])
    #maroon_color=concat([0.635,0.078,0.184])
    ## Ticks
    #Ytick1=arange(- 60,60,1)
    #Ytick2=arange(- 60,60,2)
    #Ytick5=arange(- 60,60,5)
    #YtickLabels1=cellarray([])
    #YtickLabels2=cellarray([])
    #YtickLabels5=cellarray([])
    #S0N='S0N'
    #for y in arange(1,length(Ytick1),2).reshape(-1):
    #    if Ytick1(y) == 0:
    #        YtickLabels1[y]=cellarray(['0'])
    #    else:
    #        YtickLabels1[y]=cellarray([concat([int2str(abs(Ytick1(y))),S0N(sign(Ytick1(y)) + 2)])])
    #    if y < length(Ytick1):
    #        YtickLabels1[y + 1]=cellarray([''])
    #
    #for y in arange(1,length(Ytick2),2).reshape(-1):
    #    if Ytick2(y) == 0:
    #        YtickLabels2[y]=cellarray(['0'])
    #    else:
    #        YtickLabels2[y]=cellarray([concat([int2str(abs(Ytick2(y))),S0N(sign(Ytick2(y)) + 2)])])
    #    if y < length(Ytick2):
    #        YtickLabels2[y + 1]=cellarray([''])
    #
    #for y in arange(1,length(Ytick5),2).reshape(-1):
    #    if Ytick5(y) == 0:
    #        YtickLabels5[y]=cellarray(['0'])
    #    else:
    #        YtickLabels5[y]=cellarray([concat([int2str(abs(Ytick5(y))),S0N(sign(Ytick5(y)) + 2)])])
    #    if y < length(Ytick5):
    #        YtickLabels5[y + 1]=cellarray([''])

from __future__ import division
import numpy as np
import scipy as sp
from scipy.io import netcdf
from TropD_Metric_PSI import TropD_Metric_PSI 
from TropD_Metric_TPB import TropD_Metric_TPB 
from TropD_Calculate_Mon2Season import TropD_Calculate_Mon2Season

## 1) PSI -- Streamfunction zero crossing
f_V = netcdf.netcdf_file('../ValidationData/va.nc','r')
V = f_V.variables['va'][:]
#Change axes of V to be [time, lat, lev]
V = np.transpose(V, (2,1,0))
lat = f_V.variables['lat'][:]
lev = f_V.variables['lev'][:]

Phi_psi_nh = np.zeros((np.shape(V)[0],1))
Phi_psi_sh = np.zeros((np.shape(V)[0],1))

for j in range(np.shape(V)[0]):
  Phi_psi_sh[j], Phi_psi_nh[j] = TropD_Metric_PSI(V[j,:,:], lat, lev)


# Calculate metric from annual mean
V_ANN=TropD_Calculate_Mon2Season(V,np.arange(12))

Phi_psi_nh_ANN = np.zeros((np.shape(V_ANN)[0],1))
Phi_psi_sh_ANN = np.zeros((np.shape(V_ANN)[0],1))

for j in range(np.shape(V_ANN)[0]):
  Phi_psi_sh_ANN[j], Phi_psi_nh_ANN[j] = TropD_Metric_PSI(V_ANN[j,:,:], lat, lev)


#figure
#subplot('211')
#plot(time,Phi_psi_nh,'-','linewidth',1,'color',green_color)
#hold('on')
#plot(concat([arange(y1,y2)]) + 0.5,Phi_psi_nh_ANN,'-','linewidth',2,'color',blue_color)
#plot(concat([arange(y1,y2)]) + 0.5,TropD_Calculate_Mon2Season(Phi_psi_nh,concat([arange(1,12)])),'-k','linewidth',2)
#set(gca,'fontsize',12,'linewidth',2,'tickdir','out','box','off','xtick',concat([arange(1980,2020,5)]),'xticklabels',cellarray(['']),'ytick',Ytick5,'yticklabels',YtickLabels5)
#ylabel(cellarray([['NH \\Psi_{500}'],['latitude']]))
#xlim(concat([y1,y2 + 1]))
#l=legend('Latitude of \\Psi_{500} zero crossing from monthly mean V','Latitude of \\Psi_{500} zero crossing from annual mean V','Latitude of \\Psi_{500} zero crossing from annual means of monthly metric values')
#set(l,'box','off','location','north')
#subplot('212')
#plot(time,Phi_psi_sh,'-','linewidth',1,'color',green_color)
#hold('on')
#plot(concat([arange(y1,y2)]) + 0.5,Phi_psi_sh_ANN,'-','linewidth',2,'color',blue_color)
#plot(concat([arange(y1,y2)]) + 0.5,TropD_Calculate_Mon2Season(Phi_psi_sh,concat([arange(1,12)])),'-k','linewidth',2)
#set(gca,'fontsize',12,'linewidth',2,'tickdir','out','box','off','xtick',concat([arange(1980,2020,5)]),'ytick',Ytick5,'yticklabels',YtickLabels5)
#xlim(concat([y1,y2 + 1]))
#xlabel('Year','fontsize',14)
#ylabel(cellarray([['SH \\Psi zero crossing'],['latitude']]))
# Introduce latitude unertainty condition: no additional zero crossing is allowed within 10 degrees

Phi_psi_nh_L = np.zeros((np.shape(V)[0],1))
Phi_psi_sh_L = np.zeros((np.shape(V)[0],1))

for j in range(np.shape(V)[0]):
  Phi_psi_sh_L[j], Phi_psi_nh_L[j] = TropD_Metric_PSI(V[j,:,:], lat, lev, 'Psi_500', 10)

#figure
#plot(time,Phi_psi_nh,'-','linewidth',2,'color',green_color)
#hold('on')
#plot(time(isnan(Phi_psi_nh_L)),Phi_psi_nh(isnan(Phi_psi_nh_L)),'*','markersize',10,'color',red_color)
#hold('on')
#set(gca,'fontsize',12,'linewidth',2,'tickdir','out','box','off','xtick',concat([arange(1980,2020,5)]),'ytick',Ytick5,'yticklabels',YtickLabels5)
#ylabel(cellarray([['NH \\Psi_{500}'],['latitude']]))
#xlim(concat([y1,y2 + 1]))
#l=legend('\\Psi_{500} zero crossing from monthly mean V that qualify uncertainty criterion','\\Psi_{500} zero crossing from monthly mean V that fail uncertainty criterion')
#set(l,'box','off','location','north')
#xlabel('Year','fontsize',14)

## 2) TPB -- Tropopause break latitude
f_T = netcdf.netcdf_file('../ValidationData/ta.nc','r')
f_Z = netcdf.netcdf_file('../ValidationData/zg.nc','r')
T = f_T.variables['ta'][:]
Z = f_Z.variables['zg'][:]
#Change axes of T and Z to be [time, lat, lev]
T = np.transpose(T, (2,1,0))
Z = np.transpose(Z, (2,1,0))

lat = f_T.variables['lat'][:]
lev = f_T.variables['lev'][:]

Phi_tpb_nh = np.zeros((np.shape(T)[0],1))
Phi_tpb_sh = np.zeros((np.shape(T)[0],1))

for j in range(np.shape(V)[0]):
  Phi_tpb_sh[j],Phi_tpb_nh[j] = TropD_Metric_TPB(T[j,:,:], lat, lev)

# Calculate tropopause break from annual mean
T_ANN = TropD_Calculate_Mon2Season(T, np.arange(12))

Z_ANN = TropD_Calculate_Mon2Season(Z, np.arange(12))

Phi_tpb_nh_ANN = np.zeros((np.shape(T_ANN)[0],1))

Phi_tpb_sh_ANN = np.zeros((np.shape(T_ANN)[0],1))

Phi_tpbZ_nh_ANN = np.zeros((np.shape(T_ANN)[0],1))
                                                   
Phi_tpbZ_sh_ANN = np.zeros((np.shape(T_ANN)[0],1))

Phi_tpbT_nh_ANN = np.zeros((np.shape(T_ANN)[0],1))
                                                   
Phi_tpbT_sh_ANN = np.zeros((np.shape(T_ANN)[0],1))

for j in range(np.shape(T_ANN)[0]):
  Phi_tpb_sh_ANN[j], Phi_tpb_nh_ANN[j] = TropD_Metric_TPB(T_ANN[j,:,:], lat, lev, method='max_gradient')
  Phi_tpbT_sh_ANN[j],Phi_tpbT_nh_ANN[j] = TropD_Metric_TPB(T_ANN[j,:,:], lat, lev, method='max_potemp')
  Phi_tpbZ_sh_ANN[j],Phi_tpbZ_nh_ANN[j] = TropD_Metric_TPB(T_ANN[j,:,:], lat, lev, method='cutoff',\
                                                           Z=Z_ANN[j,:,:,], Cutoff=15*1000)

#Working up to here

figure
subplot('211')
plot(time,Phi_tpb_nh,'-','linewidth',1,'color',green_color)
hold('on')
plot(concat([arange(y1,y2)]) + 0.5,Phi_tpb_nh_ANN,'-','linewidth',2,'color',blue_color)
plot(concat([arange(y1,y2)]) + 0.5,Phi_tpbZ_nh_ANN,'--','linewidth',1,'color',blue_color)
plot(concat([arange(y1,y2)]) + 0.5,Phi_tpbT_nh_ANN,'--','linewidth',1,'color',red_color)
plot(concat([arange(y1,y2)]) + 0.5,TropD_Calculate_Mon2Season(Phi_tpb_nh,concat([arange(1,12)])),'-k','linewidth',2)
set(gca,'fontsize',12,'linewidth',2,'tickdir','out','box','off','xtick',concat([arange(1980,2020,5)]),'xticklabels',cellarray(['']),'ytick',Ytick5,'yticklabels',YtickLabels5)
ylabel(cellarray([['NH tropopause break'],['latitude']]))
xlim(concat([y1,y2 + 1]))
l=legend('Latitude of tropopause break from monthly mean T -- potential temperature difference','Latitude of tropopause break from annual mean T -- maximal gradient','Latitude of tropopause break from annual mean T -- 15km cutoff height','Latitude of tropopause break from annual mean T -- potential temperature difference','Latitude of tropopause break from annual mean of monthly metric values -- potential temperature difference')
set(l,'box','off','location','north')
subplot('212')
plot(time,Phi_tpb_sh,'-','linewidth',1,'color',green_color)
hold('on')
plot(concat([arange(y1,y2)]) + 0.5,Phi_tpb_sh_ANN,'-','linewidth',2,'color',blue_color)
plot(concat([arange(y1,y2)]) + 0.5,Phi_tpbZ_sh_ANN,'--','linewidth',1,'color',blue_color)
plot(concat([arange(y1,y2)]) + 0.5,Phi_tpbT_sh_ANN,'--','linewidth',1,'color',red_color)
plot(concat([arange(y1,y2)]) + 0.5,TropD_Calculate_Mon2Season(Phi_tpb_sh,concat([arange(1,12)])),'-k','linewidth',2)
set(gca,'fontsize',12,'linewidth',2,'tickdir','out','box','off','xtick',concat([arange(1980,2020,5)]),'ytick',Ytick5,'yticklabels',YtickLabels5)
xlim(concat([y1,y2 + 1]))
xlabel('Year','fontsize',14)
ylabel(cellarray([['SH tropopause break'],['latitude']]))
## 3) OLR -- OLR cutoff
#te: OLR is assumed to be positive upwards and in units of W/m^2
olr=- ncread('../ValidationData/rlnt.nc','rlnt')

olrcs=- ncread('../ValidationData/rlntcs.nc','rlntcs')

lat=ncread('../ValidationData/rlnt.nc','lat')

olr_ANN=TropD_Calculate_Mon2Season(olr,concat([arange(1,12)]))

olrcs_ANN=TropD_Calculate_Mon2Season(olrcs,concat([arange(1,12)]))

Phi_olr_nh=zeros(size(olr,1),1)

Phi_olr_sh=zeros(size(olr,1),1)

Phi_olr_nh_ANN=zeros(size(olr_ANN,1),1)

Phi_olr_sh_ANN=zeros(size(olr_ANN,1),1)

Phi_olrcs_nh=zeros(size(olr,1),1)

Phi_olrcs_sh=zeros(size(olr,1),1)

Phi_olrcs_nh_ANN=zeros(size(olr_ANN,1),1)

Phi_olrcs_sh_ANN=zeros(size(olr_ANN,1),1)

Phi_olr20_nh_ANN=zeros(size(olr_ANN,1),1)

Phi_olr20_sh_ANN=zeros(size(olr_ANN,1),1)

Phi_olr240_nh_ANN=zeros(size(olr_ANN,1),1)

Phi_olr240_sh_ANN=zeros(size(olr_ANN,1),1)

for j in arange(1,size(olr,1)).reshape(-1):
  Phi_olr_sh[j],Phi_olr_nh[j]=TropD_Metric_OLR(squeeze(olr(j,arange(),arange())),lat,nargout=2)
  Phi_olrcs_sh[j],Phi_olrcs_nh[j]=TropD_Metric_OLR(squeeze(olrcs(j,arange(),arange())),lat,nargout=2)

for j in arange(1,size(olr_ANN,1)).reshape(-1):
  Phi_olr_sh_ANN[j],Phi_olr_nh_ANN[j]=TropD_Metric_OLR(squeeze(olr_ANN(j,arange(),arange())),lat,nargout=2)
  Phi_olrcs_sh_ANN[j],Phi_olrcs_nh_ANN[j]=TropD_Metric_OLR(squeeze(olrcs_ANN(j,arange(),arange())),lat,nargout=2)
  Phi_olr20_sh_ANN[j],Phi_olr20_nh_ANN[j]=TropD_Metric_OLR(squeeze(olr_ANN(j,arange(),arange())),lat,'20W',nargout=2)
  Phi_olr240_sh_ANN[j],Phi_olr240_nh_ANN[j]=TropD_Metric_OLR(squeeze(olr_ANN(j,arange(),arange())),lat,'cutoff',240,nargout=2)

figure
subplot('211')
plot(time,Phi_olr_nh,'-','linewidth',3,'color',green_color)
hold('on')
plot(time,Phi_olrcs_nh,'-','linewidth',1,'color',multiply(green_color,0.5))
hold('on')
plot(concat([arange(y1,y2)]) + 0.5,Phi_olr_nh_ANN,'-','linewidth',3,'color',blue_color)
plot(concat([arange(y1,y2)]) + 0.5,Phi_olrcs_nh_ANN,'-','linewidth',1,'color',multiply(blue_color,0.5))
set(gca,'fontsize',12,'linewidth',2,'tickdir','out','box','off','xtick',concat([arange(1980,2020,5)]),'xticklabels',cellarray(['']),'ytick',Ytick2,'yticklabels',YtickLabels2)
ylabel('NH OLR cutoff latitude')
xlim(concat([y1,y2 + 1]))
l=legend('Latitude of OLR 250W/m^2 cutoff latitude from monthly OLR','Latitude of OLR 250W/m^2 cutoff latitude from monthly clear-sky OLR','Latitude of OLR 250W/m^2 cutoff latitude from annual mean OLR','Latitude of OLR 250W/m^2 cutoff latitude from annual mean clear-sky OLR')
set(l,'box','off','location','north')
subplot('212')
plot(time,Phi_olr_sh,'-','linewidth',3,'color',green_color)
hold('on')
plot(time,Phi_olrcs_sh,'-','linewidth',1,'color',dot(green_color,0.5))
hold('on')
plot(concat([arange(y1,y2)]) + 0.5,Phi_olr_sh_ANN,'-','linewidth',3,'color',blue_color)
plot(concat([arange(y1,y2)]) + 0.5,Phi_olrcs_sh_ANN,'-','linewidth',1,'color',dot(blue_color,0.5))
set(gca,'fontsize',12,'linewidth',2,'tickdir','out','box','off','xtick',concat([arange(1980,2020,5)]),'ytick',Ytick2,'yticklabels',YtickLabels2)
xlim(concat([y1,y2 + 1]))
xlabel('Year','fontsize',14)
ylabel('SH OLR cutoff latitude')
figure
subplot('211')
hold('on')
plot(concat([arange(y1,y2)]) + 0.5,Phi_olr_nh_ANN,'-','linewidth',3,'color',multiply(blue_color,0.5))
plot(concat([arange(y1,y2)]) + 0.5,Phi_olr240_nh_ANN,'-','linewidth',3,'color',blue_color)
plot(concat([arange(y1,y2)]) + 0.5,Phi_olr20_nh_ANN,'-','linewidth',3,'color',green_color)
set(gca,'fontsize',12,'linewidth',2,'tickdir','out','box','off','xtick',concat([arange(1980,2020,5)]),'xticklabels',cellarray(['']),'ytick',Ytick2,'yticklabels',YtickLabels2)
ylabel('NH OLR cutoff latitude')
xlim(concat([y1,y2 + 1]))
l=legend('Latitude of OLR 250W/m^2 {default} cutoff latitude from annual-mean OLR','Latitude of OLR 240W/m^2 cutoff latitude from annual-mean OLR','Latitude of OLR -20W/m^2 cutoff latitude from annual-mean OLR')
set(l,'box','off','location','north')
subplot('212')
hold('on')
plot(concat([arange(y1,y2)]) + 0.5,Phi_olr_sh_ANN,'-','linewidth',3,'color',multiply(blue_color,0.5))
plot(concat([arange(y1,y2)]) + 0.5,Phi_olr240_sh_ANN,'-','linewidth',3,'color',blue_color)
plot(concat([arange(y1,y2)]) + 0.5,Phi_olr20_sh_ANN,'-','linewidth',3,'color',green_color)
set(gca,'fontsize',12,'linewidth',2,'tickdir','out','box','off','xtick',concat([arange(1980,2020,5)]),'ytick',Ytick2,'yticklabels',YtickLabels2)
xlim(concat([y1,y2 + 1]))
xlabel('Year','fontsize',14)
ylabel('SH OLR cutoff latitude')
## 4) STJ -- Subtropical Jet (STJ) latitude
U=ncread('../ValidationData/ua.nc','ua')

lat=ncread('../ValidationData/ua.nc','lat')

lev=ncread('../ValidationData/ua.nc','lev')

# Calculate STJ latitude from annual mean
U_ANN=TropD_Calculate_Mon2Season(U,concat([arange(1,12)]))

Phi_stj_nh_ANN_adj=zeros(size(U_ANN,1),1)

Phi_stj_sh_ANN_adj=zeros(size(U_ANN,1),1)

Phi_stj_nh_ANN_core=zeros(size(U_ANN,1),1)

Phi_stj_sh_ANN_core=zeros(size(U_ANN,1),1)

for j in arange(1,size(U_ANN,1)).reshape(-1):
  Phi_stj_sh_ANN_adj[j],Phi_stj_nh_ANN_adj[j]=TropD_Metric_STJ(squeeze(U_ANN(j,arange(),arange())),lat,lev,nargout=2)

  Phi_stj_sh_ANN_core[j],Phi_stj_nh_ANN_core[j]=TropD_Metric_STJ(squeeze(U_ANN(j,arange(),arange())),lat,lev,'core',nargout=2)

figure
subplot('211')
plot(concat([arange(y1,y2)]) + 0.5,Phi_stj_nh_ANN_adj,'-','linewidth',2,'color',green_color)
hold('on')
plot(concat([arange(y1,y2)]) + 0.5,Phi_stj_nh_ANN_core,'-','linewidth',2,'color',blue_color)
set(gca,'fontsize',12,'linewidth',2,'tickdir','out','box','off','xtick',concat([arange(1980,2020,5)]),'xticklabels',cellarray(['']),'ytick',Ytick2,'yticklabels',YtickLabels2)
ylabel('NH STJ latitude')
xlim(concat([y1,y2 + 1]))
l=legend('Latitude of STJ from anual mean U, using \'adjusted\' method','Latitude of STJ from anual mean U, using \'core\' method')
set(l,'box','off','location','north')
subplot('212')
plot(concat([arange(y1,y2)]) + 0.5,Phi_stj_sh_ANN_adj,'-','linewidth',2,'color',green_color)
hold('on')
plot(concat([arange(y1,y2)]) + 0.5,Phi_stj_sh_ANN_core,'-','linewidth',2,'color',blue_color)
set(gca,'fontsize',12,'linewidth',2,'tickdir','out','box','off','xtick',concat([arange(1980,2020,5)]),'ytick',Ytick2,'yticklabels',YtickLabels2)
xlim(concat([y1,y2 + 1]))
xlabel('Year','fontsize',14)
ylabel('SH STJ latitude')
## 5) EDJ -- Eddy Driven Jet (EDJ) latitude
U=ncread('../ValidationData/ua.nc','ua')

lat=ncread('../ValidationData/ua.nc','lat')

lev=ncread('../ValidationData/ua.nc','lev')

Phi_edj_nh=zeros(size(U,1),1)

Phi_edj_sh=zeros(size(U,1),1)

for j in arange(1,size(U,1)).reshape(-1):
  Phi_edj_sh[j],Phi_edj_nh[j]=TropD_Metric_EDJ(squeeze(U(j,arange(),arange())),lat,lev,'max',nargout=2)

# Calculate EDJ latitude from annual mean
U_ANN=TropD_Calculate_Mon2Season(U,concat([arange(1,12)]))

Phi_edj_nh_ANN=zeros(size(U_ANN,1),1)

Phi_edj_sh_ANN=zeros(size(U_ANN,1),1)

for j in arange(1,size(U_ANN,1)).reshape(-1):
  Phi_edj_sh_ANN[j],Phi_edj_nh_ANN[j]=TropD_Metric_EDJ(squeeze(U_ANN(j,arange(),arange())),lat,lev,nargout=2)

figure
subplot('211')
plot(time,Phi_edj_nh,'-','linewidth',1,'color',green_color)
hold('on')
plot(concat([arange(y1,y2)]) + 0.5,Phi_edj_nh_ANN,'-','linewidth',2,'color',blue_color)
plot(concat([arange(y1,y2)]) + 0.5,TropD_Calculate_Mon2Season(Phi_edj_nh,concat([arange(1,12)])),'-k','linewidth',2)
set(gca,'fontsize',12,'linewidth',2,'tickdir','out','box','off','xtick',concat([arange(1980,2020,5)]),'xticklabels',cellarray(['']),'ytick',Ytick2,'yticklabels',YtickLabels2)
ylabel('NH EDJ latitude')
xlim(concat([y1,y2 + 1]))
l=legend('Latitude of EDJ from monthly mean U','Latitude of EDJ from annual mean U','Latitude of EDJ from annual mean of monthly metric values')
set(l,'box','off','location','north')
subplot('212')
plot(time,Phi_edj_sh,'-','linewidth',1,'color',green_color)
hold('on')
plot(concat([arange(y1,y2)]) + 0.5,Phi_edj_sh_ANN,'-','linewidth',2,'color',blue_color)
plot(concat([arange(y1,y2)]) + 0.5,TropD_Calculate_Mon2Season(Phi_edj_sh,concat([arange(1,12)])),'-k','linewidth',2)
set(gca,'fontsize',12,'linewidth',2,'tickdir','out','box','off','xtick',concat([arange(1980,2020,5)]),'ytick',Ytick2,'yticklabels',YtickLabels2)
xlim(concat([y1,y2 + 1]))
xlabel('Year','fontsize',14)
ylabel('SH EDJ latitude')
## 6) PE -- Precipitation minus evaporation subtropical zero crossing latitude
pr=ncread('../ValidationData/pr.nc','pr')

L=2510400.0

er=- ncread('../ValidationData/hfls.nc','hfls') / L

PE=pr - er
lat=ncread('../ValidationData/pr.nc','lat')

PE_ANN=TropD_Calculate_Mon2Season(PE,concat([arange(1,12)]))

Phi_pe_nh=zeros(size(pr,1),1)

Phi_pe_sh=zeros(size(pr,1),1)

Phi_pe_nh_ANN=zeros(size(PE_ANN,1),1)

Phi_pe_sh_ANN=zeros(size(PE_ANN,1),1)

for j in arange(1,size(PE,1)).reshape(-1):
  Phi_pe_sh[j],Phi_pe_nh[j]=TropD_Metric_PE(squeeze(PE(j,arange(),arange())),lat,nargout=2)

for j in arange(1,size(PE_ANN,1)).reshape(-1):
  Phi_pe_sh_ANN[j],Phi_pe_nh_ANN[j]=TropD_Metric_PE(squeeze(PE_ANN(j,arange(),arange())),lat,nargout=2)

figure
subplot('211')
plot(time,Phi_pe_nh,'-','linewidth',2,'color',green_color)
hold('on')
plot(concat([arange(y1,y2)]) + 0.5,Phi_pe_nh_ANN,'-','linewidth',2,'color',blue_color)
plot(concat([arange(y1,y2)]) + 0.5,TropD_Calculate_Mon2Season(Phi_pe_nh,concat([arange(1,12)])),'-k','linewidth',2)
set(gca,'fontsize',12,'linewidth',2,'tickdir','out','box','off','xtick',concat([arange(1980,2020,5)]),'xticklabels',cellarray(['']),'ytick',Ytick2,'yticklabels',YtickLabels2)
ylabel('NH P - E zero-crossing')
xlim(concat([y1,y2 + 1]))
l=legend('Latitude of P minus E zero-crossing','Latitude of P minus E zero-crossing from annual mean field','Latitude of P minus E zero-crossing from annual mean of monthly metric')
set(l,'box','off','location','north')
subplot('212')
plot(time,Phi_pe_sh,'-','linewidth',2,'color',green_color)
hold('on')
plot(concat([arange(y1,y2)]) + 0.5,Phi_pe_sh_ANN,'-','linewidth',2,'color',blue_color)
plot(concat([arange(y1,y2)]) + 0.5,TropD_Calculate_Mon2Season(Phi_pe_sh,concat([arange(1,12)])),'-k','linewidth',2)
set(gca,'fontsize',12,'linewidth',2,'tickdir','out','box','off','xtick',concat([arange(1980,2020,5)]),'ytick',Ytick2,'yticklabels',YtickLabels2)
xlim(concat([y1,y2 + 1]))
xlabel('Year','fontsize',14)
ylabel('SH P - E zero-crossing')
## 7) UAS -- Zonal surface wind subtropical zero crossing latitude
U=ncread('../ValidationData/ua.nc','ua')

lat=ncread('../ValidationData/ua.nc','lat')

lev=ncread('../ValidationData/ua.nc','lev')

uas=ncread('../ValidationData/uas.nc','uas')

uas_ANN=TropD_Calculate_Mon2Season(uas,concat([arange(1,12)]))

U_ANN=TropD_Calculate_Mon2Season(U,concat([arange(1,12)]))

Phi_uas_nh=zeros(size(uas,1),1)

Phi_uas_sh=zeros(size(uas,1),1)

Phi_uas_nh_ANN=zeros(size(uas_ANN,1),1)

Phi_uas_sh_ANN=zeros(size(uas_ANN,1),1)

Phi_Uas_nh_ANN=zeros(size(uas_ANN,1),1)

Phi_Uas_sh_ANN=zeros(size(uas_ANN,1),1)

# The TropD_Metric_UAS metric accepts both 2D and 1D zonal wind. 
#r surface 1D wind, both of the following syntaxes are allowed and equivalent
[Phi_sh , Phi_nh] = TropD_Metric_UAS(uas,lat) ;   
[Phi_sh , Phi_nh] = TropD_Metric_UAS(uas,lat,1) ;

for j in arange(1,size(uas,1)).reshape(-1):
  Phi_uas_sh[j],Phi_uas_nh[j]=TropD_Metric_UAS(squeeze(uas(j,arange(),arange())),lat,nargout=2)

for j in arange(1,size(uas_ANN,1)).reshape(-1):
  Phi_uas_sh_ANN[j],Phi_uas_nh_ANN[j]=TropD_Metric_UAS(squeeze(uas_ANN(j,arange(),arange())),lat,nargout=2)
  Phi_Uas_sh_ANN[j],Phi_Uas_nh_ANN[j]=TropD_Metric_UAS(squeeze(U_ANN(j,arange(),arange())),lat,lev,nargout=2)

figure
subplot('211')
plot(time,Phi_uas_nh,'-','linewidth',2,'color',green_color)
hold('on')
plot(concat([arange(y1,y2)]) + 0.5,Phi_uas_nh_ANN,'-','linewidth',2,'color',blue_color)
plot(concat([arange(y1,y2)]) + 0.5,Phi_Uas_nh_ANN,'-','linewidth',2,'color',red_color)
plot(concat([arange(y1,y2)]) + 0.5,TropD_Calculate_Mon2Season(Phi_uas_nh,concat([arange(1,12)])),'-k','linewidth',2)
set(gca,'fontsize',12,'linewidth',2,'tickdir','out','box','off','xtick',concat([arange(1980,2020,5)]),'xticklabels',cellarray(['']),'ytick',Ytick2,'yticklabels',YtickLabels2)
ylabel('NH uas zero-crossing')
xlim(concat([y1,y2 + 1]))
l=legend('Latitude of surface zonal wind zero crossing','Latitude of surface zonal wind zero crossing from annual mean field','Latitude of 850 hPa zonal wind zero crossing from annual mean field','Latitude of surface zonal wind zero crossing from annual mean of monthly metric')
set(l,'box','off','location','north')
subplot('212')
plot(time,Phi_uas_sh,'-','linewidth',2,'color',green_color)
hold('on')
plot(concat([arange(y1,y2)]) + 0.5,Phi_uas_sh_ANN,'-','linewidth',2,'color',blue_color)
plot(concat([arange(y1,y2)]) + 0.5,Phi_Uas_sh_ANN,'-','linewidth',2,'color',red_color)
plot(concat([arange(y1,y2)]) + 0.5,TropD_Calculate_Mon2Season(Phi_uas_sh,concat([arange(1,12)])),'-k','linewidth',2)
set(gca,'fontsize',12,'linewidth',2,'tickdir','out','box','off','xtick',concat([arange(1980,2020,5)]),'ytick',Ytick2,'yticklabels',YtickLabels2)
xlim(concat([y1,y2 + 1]))
xlabel('Year','fontsize',14)
ylabel('SH uas zero-crossing')
## 8) PSL -- Sea-level Pressure Maximum
ps=ncread('../ValidationData/psl.nc','psl')

lat=ncread('../ValidationData/psl.nc','lat')

ps_DJF=TropD_Calculate_Mon2Season(ps,concat([1,2,12]))

ps_JJA=TropD_Calculate_Mon2Season(ps,concat([6,7,8]))

Phi_ps_DJF_nh=zeros(size(ps_DJF,1),1)

Phi_ps_JJA_nh=zeros(size(ps_JJA,1),1)

Phi_ps_DJF_sh=zeros(size(ps_DJF,1),1)

Phi_ps_JJA_sh=zeros(size(ps_JJA,1),1)

for j in arange(1,size(ps_DJF,1)).reshape(-1):
  Phi_ps_DJF_sh[j],Phi_ps_DJF_nh[j]=TropD_Metric_PSL(squeeze(ps_DJF(j,arange(),arange())).T,lat,nargout=2)

for j in arange(1,size(ps_JJA,1)).reshape(-1):
  Phi_ps_JJA_sh[j],Phi_ps_JJA_nh[j]=TropD_Metric_PSL(squeeze(ps_JJA(j,arange(),arange())).T,lat,nargout=2)

figure
subplot('211')
plot(concat([arange(y1,y2)]) + 0.5,Phi_ps_DJF_nh,'-','linewidth',2,'color',green_color)
hold('on')
plot(concat([arange(y1,y2)]) + 0.5,Phi_ps_JJA_nh,'-','linewidth',2,'color',blue_color)
set(gca,'fontsize',12,'linewidth',2,'tickdir','out','box','off','xtick',concat([arange(1980,2020,5)]),'ytick',Ytick1,'yticklabels',YtickLabels1)
ylabel('NH max psl latitude')
xlim(concat([y1,y2 + 1]))
l=legend('Latitude of max sea-level pressure during DJF','Latitude of max sea-level pressure during JJA')
set(l,'box','off','location','south')
subplot('212')
plot(concat([arange(y1,y2)]) + 0.5,Phi_ps_DJF_sh,'-','linewidth',2,'color',green_color)
hold('on')
plot(concat([arange(y1,y2)]) + 0.5,Phi_ps_JJA_sh,'-','linewidth',2,'color',blue_color)
set(gca,'fontsize',12,'linewidth',2,'tickdir','out','box','off','xtick',concat([arange(1980,2020,5)]),'ytick',Ytick1,'yticklabels',YtickLabels1)
ylabel('SH max psl latitude')
xlim(concat([y1,y2 + 1]))
## 9) Compare annual mean metrics
i500
V=ncread('../ValidationData/va.nc','va')

lat=ncread('../ValidationData/va.nc','lat')

lev=ncread('../ValidationData/va.nc','lev')

V_ANN=TropD_Calculate_Mon2Season(V,concat([arange(1,12)]))

Phi_psi_nh_ANN=zeros(size(V_ANN,1),1)

Phi_psi_sh_ANN=zeros(size(V_ANN,1),1)

for j in arange(1,size(V_ANN,1)).reshape(-1):
  Phi_psi_sh_ANN[j],Phi_psi_nh_ANN[j]=TropD_Metric_PSI(squeeze(V_ANN(j,arange(),arange())),lat,lev,nargout=2)

# Tropopause break
T=ncread('../ValidationData/ta.nc','ta')

T_ANN=TropD_Calculate_Mon2Season(T,concat([arange(1,12)]))

Phi_tpb_nh_ANN=zeros(size(T_ANN,1),1)

Phi_tpb_sh_ANN=zeros(size(T_ANN,1),1)

for j in arange(1,size(T_ANN,1)).reshape(-1):
  Phi_tpb_sh_ANN[j],Phi_tpb_nh_ANN[j]=TropD_Metric_TPB(squeeze(T_ANN(j,arange(),arange())),lat,lev,nargout=2)

# Surface pressure max
ps=ncread('../ValidationData/psl.nc','psl')

ps_ANN=TropD_Calculate_Mon2Season(ps,concat([arange(1,12)]))

Phi_ps_nh_ANN=zeros(size(ps_ANN,1),1)

Phi_ps_sh_ANN=zeros(size(ps_ANN,1),1)

for j in arange(1,size(ps_ANN,1)).reshape(-1):
  Phi_ps_sh_ANN[j],Phi_ps_nh_ANN[j]=TropD_Metric_PSL(squeeze(ps_ANN(j,arange(),arange())).T,lat,nargout=2)

# Eddy driven jet
U=ncread('../ValidationData/ua.nc','ua')

U_ANN=TropD_Calculate_Mon2Season(U,concat([arange(1,12)]))

Phi_edj_nh_ANN=zeros(size(U_ANN,1),1)

Phi_edj_sh_ANN=zeros(size(U_ANN,1),1)

for j in arange(1,size(U_ANN,1)).reshape(-1):
  Phi_edj_sh_ANN[j],Phi_edj_nh_ANN[j]=TropD_Metric_EDJ(squeeze(U_ANN(j,arange(),arange())),lat,lev,nargout=2)

# Subtropical jet
Phi_stj_nh_ANN=zeros(size(U_ANN,1),1)

Phi_stj_sh_ANN=zeros(size(U_ANN,1),1)

for j in arange(1,size(U_ANN,1)).reshape(-1):
  Phi_stj_sh_ANN[j],Phi_stj_nh_ANN[j]=TropD_Metric_STJ(squeeze(U_ANN(j,arange(),arange())),lat,lev,nargout=2)

# OLR
olr=- ncread('../ValidationData/rlnt.nc','rlnt')

olr_ANN=TropD_Calculate_Mon2Season(olr,concat([arange(1,12)]))

Phi_olr_nh_ANN=zeros(size(olr_ANN,1),1)

Phi_olr_sh_ANN=zeros(size(olr_ANN,1),1)

for j in arange(1,size(olr_ANN,1)).reshape(-1):
  Phi_olr_sh_ANN[j],Phi_olr_nh_ANN[j]=TropD_Metric_OLR(squeeze(olr_ANN(j,arange(),arange())),lat,nargout=2)

# P minus E
pr=ncread('../ValidationData/pr.nc','pr')

L=2510400.0

er=- ncread('../ValidationData/hfls.nc','hfls') / L

PE=pr - er
PE_ANN=TropD_Calculate_Mon2Season(PE,concat([arange(1,12)]))

Phi_pe_nh_ANN=zeros(size(PE_ANN,1),1)

Phi_pe_sh_ANN=zeros(size(PE_ANN,1),1)

for j in arange(1,size(PE_ANN,1)).reshape(-1):
  Phi_pe_sh_ANN[j],Phi_pe_nh_ANN[j]=TropD_Metric_PE(squeeze(PE_ANN(j,arange(),arange())),lat,nargout=2)

# Surface winds
uas=ncread('../ValidationData/uas.nc','uas')

uas_ANN=TropD_Calculate_Mon2Season(uas,concat([arange(1,12)]))

Phi_uas_nh_ANN=zeros(size(uas_ANN,1),1)

Phi_uas_sh_ANN=zeros(size(uas_ANN,1),1)

for j in arange(1,size(uas_ANN,1)).reshape(-1):
  Phi_uas_sh_ANN[j],Phi_uas_nh_ANN[j]=TropD_Metric_UAS(squeeze(uas_ANN(j,arange(),arange())),lat,nargout=2)


figure
subplot('211')
plot(concat([arange(y1,y2)]),Phi_psi_nh_ANN,'-','linewidth',2,'color',concat([0,0,0]))
hold('on')
plot(concat([arange(y1,y2)]),Phi_tpb_nh_ANN,'-','linewidth',2,'color',green_color)
plot(concat([arange(y1,y2)]),Phi_edj_nh_ANN,'-','linewidth',2,'color',blue_color)
plot(concat([arange(y1,y2)]),Phi_stj_nh_ANN,'-','linewidth',2,'color',red_color)
plot(concat([arange(y1,y2)]),Phi_olr_nh_ANN,'-','linewidth',2,'color',lightblue_color)
plot(concat([arange(y1,y2)]),Phi_pe_nh_ANN,'-','linewidth',2,'color',orange_color)
plot(concat([arange(y1,y2)]),Phi_uas_nh_ANN,'-','linewidth',2,'color',purple_color)
plot(concat([arange(y1,y2)]),Phi_ps_nh_ANN,'--','linewidth',2,'color',maroon_color)
set(gca,'fontsize',12,'linewidth',2,'tickdir','out','box','off','xtick',concat([arange(1980,2020,5)]),'xticklabels',cellarray(['']),'ytick',Ytick5,'yticklabels',YtickLabels5)
l=legend('PSI','TPB','EDJ','STJ','OLR','P-E','UAS','PSL')
set(l,'box','off','location','eastoutside')
ylabel('NH HC edge')
xlim(concat([y1,y2 + 1]))
ylim(concat([25,50]))
subplot('212')
plot(concat([arange(y1,y2)]),Phi_psi_sh_ANN,'-','linewidth',2,'color',concat([0,0,0]))
hold('on')
plot(concat([arange(y1,y2)]),Phi_tpb_sh_ANN,'-','linewidth',2,'color',green_color)
plot(concat([arange(y1,y2)]),Phi_edj_sh_ANN,'-','linewidth',2,'color',blue_color)
plot(concat([arange(y1,y2)]),Phi_stj_sh_ANN,'-','linewidth',2,'color',red_color)
plot(concat([arange(y1,y2)]),Phi_olr_sh_ANN,'-','linewidth',2,'color',lightblue_color)
plot(concat([arange(y1,y2)]),Phi_pe_sh_ANN,'-','linewidth',2,'color',orange_color)
plot(concat([arange(y1,y2)]),Phi_uas_sh_ANN,'-','linewidth',2,'color',purple_color)
plot(concat([arange(y1,y2)]),Phi_ps_sh_ANN,'--','linewidth',2,'color',maroon_color)
set(gca,'fontsize',12,'linewidth',2,'tickdir','out','box','off','xtick',concat([arange(1980,2020,5)]),'ytick',Ytick5,'yticklabels',YtickLabels5)
l=legend('PSI','TPB','EDJ','STJ','OLR','P-E','UAS','PSL')
set(l,'box','off','location','eastoutside')
xlim(concat([y1,y2 + 1]))
xlabel('Year','fontsize',14)
ylabel('SH HC edge')
## 10) Validate metrics
i500
V=ncread('../ValidationData/va.nc','va')

lat=ncread('../ValidationData/va.nc','lat')

lev=ncread('../ValidationData/va.nc','lev')

V_ANN=TropD_Calculate_Mon2Season(V,concat([arange(1,12)]))

Phi_psi_nh=zeros(size(V,1),1)

Phi_psi_sh=zeros(size(V,1),1)

Phi_psi_nh_ANN=zeros(size(V_ANN,1),1)

Phi_psi_sh_ANN=zeros(size(V_ANN,1),1)

for j in arange(1,size(V,1)).reshape(-1):
  Phi_psi_sh[j],Phi_psi_nh[j]=TropD_Metric_PSI(squeeze(V(j,arange(),arange())),lat,lev,nargout=2)

for j in arange(1,size(V_ANN,1)).reshape(-1):
  Phi_psi_sh_ANN[j],Phi_psi_nh_ANN[j]=TropD_Metric_PSI(squeeze(V_ANN(j,arange(),arange())),lat,lev,nargout=2)

Phi_nh=ncread('../ValidationMetrics/PSI_ANN.nc','PSI_NH')
Phi_sh=ncread('../ValidationMetrics/PSI_ANN.nc','PSI_SH')
if logical_not((std(Phi_nh - Phi_psi_nh_ANN) < 1e-10 and std(Phi_sh - Phi_psi_sh_ANN) < 1e-10)):
  disp('Warning: annual-mean Validation and calculated PSI metrics are NOT equal!')
else:
  disp('OK. Annual-mean Validation and calculated PSI metrics are the same!')

Phi_nh=ncread('../ValidationMetrics/PSI.nc','PSI_NH')
Phi_sh=ncread('../ValidationMetrics/PSI.nc','PSI_SH')
if logical_not((std(Phi_nh - Phi_psi_nh) < 1e-10 and std(Phi_sh - Phi_psi_sh) < 1e-10)):
  disp('Warning: monthly Validation and calculated PSI metrics are NOT equal!')
else:
  disp('OK. Monthly Validation and calculated PSI metrics are the same!')

# Tropopause break
T=ncread('../ValidationData/ta.nc','ta')

T_ANN=TropD_Calculate_Mon2Season(T,concat([arange(1,12)]))

Phi_tpb_nh=zeros(size(T,1),1)

Phi_tpb_sh=zeros(size(T,1),1)

Phi_tpb_nh_ANN=zeros(size(T_ANN,1),1)

Phi_tpb_sh_ANN=zeros(size(T_ANN,1),1)

for j in arange(1,size(T,1)).reshape(-1):
  Phi_tpb_sh[j],Phi_tpb_nh[j]=TropD_Metric_TPB(squeeze(T(j,arange(),arange())),lat,lev,nargout=2)

for j in arange(1,size(T_ANN,1)).reshape(-1):
  Phi_tpb_sh_ANN[j],Phi_tpb_nh_ANN[j]=TropD_Metric_TPB(squeeze(T_ANN(j,arange(),arange())),lat,lev,nargout=2)

Phi_nh=ncread('../ValidationMetrics/TPB_ANN.nc','TPB_NH')
Phi_sh=ncread('../ValidationMetrics/TPB_ANN.nc','TPB_SH')
if logical_not((std(Phi_nh - Phi_tpb_nh_ANN) < 1e-10 and std(Phi_sh - Phi_tpb_sh_ANN) < 1e-10)):
  disp('Warning: annual-mean Validation and calculated TPB metrics are NOT equal!')
else:
  disp('OK. Annual-mean Validation and calculated TPB metrics are the same!')

Phi_nh=ncread('../ValidationMetrics/TPB.nc','TPB_NH')
Phi_sh=ncread('../ValidationMetrics/TPB.nc','TPB_SH')
if logical_not((std(Phi_nh - Phi_tpb_nh) < 1e-10 and std(Phi_sh - Phi_tpb_sh) < 1e-10)):
  disp('Warning: monthly Validation and calculated TPB metrics are NOT equal!')
else:
  disp('OK. Monthly Validation and calculated TPB metrics are the same!')

# Surface pressure max (Invalid in NH)
ps=ncread('../ValidationData/psl.nc','psl')

ps_DJF=TropD_Calculate_Mon2Season(ps,concat([1,2,12]))

ps_JJA=TropD_Calculate_Mon2Season(ps,concat([arange(6,8)]))

ps_MAM=TropD_Calculate_Mon2Season(ps,concat([arange(3,5)]))

ps_SON=TropD_Calculate_Mon2Season(ps,concat([arange(9,11)]))

Phi_ps_sh_DJF=zeros(size(ps_DJF,1),1)
Phi_ps_sh_JJA=zeros(size(ps_JJA,1),1)
Phi_ps_sh_MAM=zeros(size(ps_MAM,1),1)
Phi_ps_sh_SON=zeros(size(ps_SON,1),1)
Phi_ps_nh_DJF=zeros(size(ps_DJF,1),1)
Phi_ps_nh_JJA=zeros(size(ps_JJA,1),1)
Phi_ps_nh_MAM=zeros(size(ps_MAM,1),1)
Phi_ps_nh_SON=zeros(size(ps_SON,1),1)
for j in arange(1,size(ps_DJF,1)).reshape(-1):
  Phi_ps_sh_DJF[j],Phi_ps_nh_DJF[j]=TropD_Metric_PSL(squeeze(ps_DJF(j,arange(),arange())).T,lat,nargout=2)

for j in arange(1,size(ps_JJA,1)).reshape(-1):
  Phi_ps_sh_JJA[j],Phi_ps_nh_JJA[j]=TropD_Metric_PSL(squeeze(ps_JJA(j,arange(),arange())).T,lat,nargout=2)

for j in arange(1,size(ps_MAM,1)).reshape(-1):
  Phi_ps_sh_MAM[j],Phi_ps_nh_MAM[j]=TropD_Metric_PSL(squeeze(ps_MAM(j,arange(),arange())).T,lat,nargout=2)

for j in arange(1,size(ps_SON,1)).reshape(-1):
  Phi_ps_sh_SON[j],Phi_ps_nh_SON[j]=TropD_Metric_PSL(squeeze(ps_SON(j,arange(),arange())).T,lat,nargout=2)

Phi_sh=ncread('../ValidationMetrics/PSL_DJF.nc','PSL_SH')
Phi_nh=ncread('../ValidationMetrics/PSL_DJF.nc','PSL_NH')
if logical_not((std(Phi_sh - Phi_ps_sh_DJF) < 1e-10)) or logical_not((std(Phi_nh - Phi_ps_nh_DJF) < 1e-10)):
  disp('Warning: DJF Validation and calculated PSL metrics are NOT equal!')
else:
  disp('OK. DJF Validation and calculated PSL metrics are the same!')

Phi_sh=ncread('../ValidationMetrics/PSL_JJA.nc','PSL_SH')
Phi_nh=ncread('../ValidationMetrics/PSL_JJA.nc','PSL_NH')
if logical_not((std(Phi_sh - Phi_ps_sh_JJA) < 1e-10)) or logical_not((std(Phi_nh - Phi_ps_nh_JJA) < 1e-10)):
  disp('Warning: JJA Validation and calculated PSL metrics are NOT equal!')
else:
  disp('OK. JJA Validation and calculated PSL metrics are the same!')

Phi_sh=ncread('../ValidationMetrics/PSL_MAM.nc','PSL_SH')
Phi_nh=ncread('../ValidationMetrics/PSL_MAM.nc','PSL_NH')
if logical_not((std(Phi_sh - Phi_ps_sh_MAM) < 1e-10)) or logical_not((std(Phi_nh - Phi_ps_nh_MAM) < 1e-10)):
  disp('Warning: MAM Validation and calculated PSL metrics are NOT equal!')
else:
  disp('OK. MAM Validation and calculated PSL metrics are the same!')

Phi_sh=ncread('../ValidationMetrics/PSL_SON.nc','PSL_SH')
Phi_nh=ncread('../ValidationMetrics/PSL_SON.nc','PSL_NH')
if logical_not((std(Phi_sh - Phi_ps_sh_SON) < 1e-10)) or logical_not((std(Phi_nh - Phi_ps_nh_SON) < 1e-10)):
  disp('Warning: SON Validation and calculated PSL metrics are NOT equal!')
else:
  disp('OK. SON Validation and calculated PSL metrics are the same!')

# Eddy driven jet
U=ncread('../ValidationData/ua.nc','ua')

U_ANN=TropD_Calculate_Mon2Season(U,concat([arange(1,12)]))

Phi_edj_nh=zeros(size(U,1),1)

Phi_edj_sh=zeros(size(U,1),1)

Phi_edj_nh_ANN=zeros(size(U_ANN,1),1)

Phi_edj_sh_ANN=zeros(size(U_ANN,1),1)

for j in arange(1,size(U,1)).reshape(-1):
  Phi_edj_sh[j],Phi_edj_nh[j]=TropD_Metric_EDJ(squeeze(U(j,arange(),arange())),lat,lev,nargout=2)

for j in arange(1,size(U_ANN,1)).reshape(-1):
  Phi_edj_sh_ANN[j],Phi_edj_nh_ANN[j]=TropD_Metric_EDJ(squeeze(U_ANN(j,arange(),arange())),lat,lev,nargout=2)

Phi_nh=ncread('../ValidationMetrics/EDJ_ANN.nc','EDJ_NH')
Phi_sh=ncread('../ValidationMetrics/EDJ_ANN.nc','EDJ_SH')
if logical_not((std(Phi_nh - Phi_edj_nh_ANN) < 1e-10 and std(Phi_sh - Phi_edj_sh_ANN) < 1e-10)):
  disp('Warning: annual-mean Validation and calculated EDJ metrics are NOT equal!')
else:
  disp('OK. Annual-mean Validation and calculated EDJ metrics are the same!')

Phi_nh=ncread('../ValidationMetrics/EDJ.nc','EDJ_NH')
Phi_sh=ncread('../ValidationMetrics/EDJ.nc','EDJ_SH')
if logical_not((std(Phi_nh - Phi_edj_nh) < 1e-10 and std(Phi_sh - Phi_edj_sh) < 1e-10)):
  disp('Warning: monthly Validation and calculated EDJ metrics are NOT equal!')
else:
  disp('OK. Monthly Validation and calculated EDJ metrics are the same!')

# Subtropical jet
Phi_stj_nh=zeros(size(U,1),1)

Phi_stj_sh=zeros(size(U,1),1)

Phi_stj_nh_ANN=zeros(size(U_ANN,1),1)

Phi_stj_sh_ANN=zeros(size(U_ANN,1),1)

for j in arange(1,size(U,1)).reshape(-1):
  Phi_stj_sh[j],Phi_stj_nh[j]=TropD_Metric_STJ(squeeze(U(j,arange(),arange())),lat,lev,nargout=2)

for j in arange(1,size(U_ANN,1)).reshape(-1):
  Phi_stj_sh_ANN[j],Phi_stj_nh_ANN[j]=TropD_Metric_STJ(squeeze(U_ANN(j,arange(),arange())),lat,lev,nargout=2)

Phi_nh=ncread('../ValidationMetrics/STJ_ANN.nc','STJ_NH')
Phi_sh=ncread('../ValidationMetrics/STJ_ANN.nc','STJ_SH')
if logical_not((std(Phi_nh - Phi_stj_nh_ANN) < 1e-10 and std(Phi_sh - Phi_stj_sh_ANN) < 1e-10)):
  disp('Warning: annual-mean Validation and calculated STJ metrics are NOT equal!')
else:
  disp('OK. Annual-mean Validation and calculated STJ metrics are the same!')

Phi_nh=ncread('../ValidationMetrics/STJ.nc','STJ_NH')
Phi_sh=ncread('../ValidationMetrics/STJ.nc','STJ_SH')
if logical_not((std(Phi_nh - Phi_stj_nh) < 1e-10 and std(Phi_sh - Phi_stj_sh) < 1e-10)):
  disp('Warning: monthly Validation and calculated STJ metrics are NOT equal!')
else:
  disp('OK. Monthly Validation and calculated STJ metrics are the same!')

# OLR
olr=- ncread('../ValidationData/rlnt.nc','rlnt')

olr_ANN=TropD_Calculate_Mon2Season(olr,concat([arange(1,12)]))

Phi_olr_nh=zeros(size(olr,1),1)

Phi_olr_sh=zeros(size(olr,1),1)

Phi_olr_nh_ANN=zeros(size(olr_ANN,1),1)

Phi_olr_sh_ANN=zeros(size(olr_ANN,1),1)

for j in arange(1,size(olr,1)).reshape(-1):
  Phi_olr_sh[j],Phi_olr_nh[j]=TropD_Metric_OLR(squeeze(olr(j,arange(),arange())),lat,nargout=2)

for j in arange(1,size(olr_ANN,1)).reshape(-1):
  Phi_olr_sh_ANN[j],Phi_olr_nh_ANN[j]=TropD_Metric_OLR(squeeze(olr_ANN(j,arange(),arange())),lat,nargout=2)

Phi_nh=ncread('../ValidationMetrics/OLR_ANN.nc','OLR_NH')
Phi_sh=ncread('../ValidationMetrics/OLR_ANN.nc','OLR_SH')
if logical_not((std(Phi_nh - Phi_olr_nh_ANN) < 1e-10 and std(Phi_sh - Phi_olr_sh_ANN) < 1e-10)):
  disp('Warning: annual-mean Validation and calculated OLR metrics are NOT equal!')
else:
  disp('OK. Annual-mean Validation and calculated OLR metrics are the same!')

Phi_nh=ncread('../ValidationMetrics/OLR.nc','OLR_NH')
Phi_sh=ncread('../ValidationMetrics/OLR.nc','OLR_SH')
if logical_not((std(Phi_nh - Phi_olr_nh) < 1e-10 and std(Phi_sh - Phi_olr_sh) < 1e-10)):
  disp('Warning: monthly Validation and calculated OLR metrics are NOT equal!')
else:
  disp('OK. Monthly Validation and calculated OLR metrics are the same!')

# P minus E
pr=ncread('../ValidationData/pr.nc','pr')

L=2510400.0

er=- ncread('../ValidationData/hfls.nc','hfls') / L

PE=pr - er
PE_ANN=TropD_Calculate_Mon2Season(PE,concat([arange(1,12)]))

Phi_pe_nh=zeros(size(PE,1),1)

Phi_pe_sh=zeros(size(PE,1),1)

Phi_pe_nh_ANN=zeros(size(PE_ANN,1),1)

Phi_pe_sh_ANN=zeros(size(PE_ANN,1),1)

for j in arange(1,size(PE,1)).reshape(-1):
  Phi_pe_sh[j],Phi_pe_nh[j]=TropD_Metric_PE(squeeze(PE(j,arange(),arange())),lat,nargout=2)

for j in arange(1,size(PE_ANN,1)).reshape(-1):
  Phi_pe_sh_ANN[j],Phi_pe_nh_ANN[j]=TropD_Metric_PE(squeeze(PE_ANN(j,arange(),arange())),lat,nargout=2)

Phi_nh=ncread('../ValidationMetrics/PE_ANN.nc','PE_NH')
Phi_sh=ncread('../ValidationMetrics/PE_ANN.nc','PE_SH')
if logical_not((std(Phi_nh - Phi_pe_nh_ANN) < 1e-10 and std(Phi_sh - Phi_pe_sh_ANN) < 1e-10)):
  disp('Warning: annual-mean Validation and calculated P-E metrics are NOT equal!')
else:
  disp('OK. Annual-mean Validation and calculated P-E metrics are the same!')

Phi_nh=ncread('../ValidationMetrics/PE.nc','PE_NH')
Phi_sh=ncread('../ValidationMetrics/PE.nc','PE_SH')
if logical_not((std(Phi_nh - Phi_pe_nh) < 1e-10 and std(Phi_sh - Phi_pe_sh) < 1e-10)):
  disp('Warning: monthly Validation and calculated P-E metrics are NOT equal!')
else:
  disp('OK. Monthly Validation and calculated P-E metrics are the same!')

# Surface winds
uas=ncread('../ValidationData/uas.nc','uas')

uas_ANN=TropD_Calculate_Mon2Season(uas,concat([arange(1,12)]))

Phi_uas_nh=zeros(size(uas,1),1)

Phi_uas_sh=zeros(size(uas,1),1)

Phi_uas_nh_ANN=zeros(size(uas_ANN,1),1)

Phi_uas_sh_ANN=zeros(size(uas_ANN,1),1)

for j in arange(1,size(uas,1)).reshape(-1):
  Phi_uas_sh[j],Phi_uas_nh[j]=TropD_Metric_UAS(squeeze(uas(j,arange(),arange())),lat,nargout=2)

for j in arange(1,size(uas_ANN,1)).reshape(-1):
  Phi_uas_sh_ANN[j],Phi_uas_nh_ANN[j]=TropD_Metric_UAS(squeeze(uas_ANN(j,arange(),arange())),lat,nargout=2)

Phi_nh=ncread('../ValidationMetrics/UAS_ANN.nc','UAS_NH')
Phi_sh=ncread('../ValidationMetrics/UAS_ANN.nc','UAS_SH')
if logical_not((std(Phi_nh - Phi_uas_nh_ANN) < 1e-10 and std(Phi_sh - Phi_uas_sh_ANN) < 1e-10)):
  disp('Warning: annual-mean Validation and calculated UAS metrics are NOT equal!')
else:
  disp('OK. Annual-mean Validation and calculated UAS metrics are the same!')

Phi_nh=ncread('../ValidationMetrics/UAS.nc','UAS_NH')
Phi_sh=ncread('../ValidationMetrics/UAS.nc','UAS_SH')
if logical_not((std(Phi_nh - Phi_uas_nh) < 1e-10 and std(Phi_sh - Phi_uas_sh) < 1e-10)):
  disp('Warning: monthly Validation and calculated UAS metrics are NOT equal!')
else:
  disp('OK. Monthly Validation and calculated UAS metrics are the same!')
