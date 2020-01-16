########################################
####    CMIP5 - Climate Analysis    ####
########################################
####     Behzad Asadieh, Ph.D.      ####
####  University of Pennsylvania    ####
####    basadieh@sas.upenn.edu      ####
####     github.com/behzadasd       ####
########################################

from Behzadlib import func_latlon_regrid, func_regrid, func_oceanlandmask, func_oceanindex, func_gridcell_area, func_plotmap_contourf, func_EOF
########################################
import numpy as np
import xarray as xr
import numpy.ma as ma
from netCDF4 import MFDataset, Dataset, num2date, date2num, date2index
import os
import matplotlib
import matplotlib.mlab as ml
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm, maskoceans
from scipy.interpolate import griddata
import math
import copy
###############################################################################
dir_pwd = os.getcwd() # Gets the current directory (and in which the code is placed)

### Regrdridding calculations ###
# creating new coordinate grid, same which was used in interpolation in data processing code
lat_n_regrid, lon_n_regrid = 180, 360 # Number of Lat and Lon elements in the regridded data
lon_min_regrid, lon_max_regrid = 0, 360 # Min and Max value of Lon in the regridded data
lat_min_regrid, lat_max_regrid = -90, 90 # Min and Max value of Lat in the regridded data

# This function for creating new Lat-Lon fields is saved in Behzadlib code in this directory - imported at the begenning
Lat_regrid_1D, Lon_regrid_1D, Lat_bound_regrid, Lon_bound_regrid = func_latlon_regrid(lat_n_regrid, lon_n_regrid, lat_min_regrid, lat_max_regrid, lon_min_regrid, lon_max_regrid)
Lon_regrid_2D, Lat_regrid_2D = np.meshgrid(Lon_regrid_1D, Lat_regrid_1D)

# Land/Ocean mask - The function is saved in Behzadlib code in this directory - imported at the begenning
Ocean_Land_mask = func_oceanlandmask(Lat_regrid_2D, Lon_regrid_2D) # 1= ocean, 0= land

###############################################################################
GCM_Names = ['GFDL-ESM2M', 'GFDL-ESM2G', 'IPSL-CM5A-MR', 'IPSL-CM5A-LR', 'MIROC-ESM', 'MIROC-ESM-CHEM', 'CESM1-BGC', 'CMCC-CESM', 'CanESM2', 'GISS-E2-H-CC', 'GISS-E2-R-CC', 'MPI-ESM-MR', 'MPI-ESM-LR', 'NorESM1-ME']
GCM = 'GFDL-ESM2G'

lat_t='lat'
lon_t='lon'
time_t='time'

###############################################################################
###  Empirical Orthogonal Functions (EOFs) of Sea Surface Temperature  ########
###############################################################################
dir_data_in1 = ('/data2/scratch/cabre/CMIP5/CMIP5_models/ocean_physics/') # Directory to raed raw data from
dir_data_in2=(dir_data_in1+ GCM + '/historical/mo/')

year_start=1901
year_end=2000

Var_name='thetao' # The variable name to be read from .nc files
dset = xr.open_mfdataset(dir_data_in2+Var_name+'*12.nc')
Data_all = dset[Var_name].sel(lev=0,method='nearest') # data at lev=0   

Lat_orig = dset[lat_t]
Lat_orig = Lat_orig.values
Lon_orig = dset[lon_t]
Lon_orig = Lon_orig.values     

Data_monthly = Data_all[ (Data_all.time.dt.year >= year_start ) & (Data_all.time.dt.year <= year_end)]
Data_monthly = Data_monthly.values # Converts to a Numpy array
Data_monthly [ Data_monthly > 1e19 ] = np.nan
Data_rannual = Data_monthly.reshape( np.int(Data_monthly.shape[0]/12) ,12,Data_monthly.shape[1],Data_monthly.shape[2])
Data_rannual = np.nanmean(Data_rannual,axis=1)

# Regriding data into 1degree by 1degree fields -  The regridding function is saved in Behzadlib code in this directory - imported at the begenning
# The SST data is regridded to coarser 4x4 degree fields to reduce the calculation time and memory needed
Lat_regrid_1D_2, Lon_regrid_1D_2, Lat_bound_regrid_2, Lon_bound_regrid_2 = func_latlon_regrid(90,180, lat_min_regrid, lat_max_regrid, lon_min_regrid, lon_max_regrid)
Lon_regrid_2D_2, Lat_regrid_2D_2 = np.meshgrid(Lon_regrid_1D_2, Lat_regrid_1D_2)
SST_annual = func_regrid(Data_rannual, Lat_orig, Lon_orig, Lat_regrid_2D_2, Lon_regrid_2D_2)

# Keeping only Pacific Ocean grid-cells
# Ocean Index mask - The function is saved in Behzadlib code in this directory - imported at the begenning
Ocean_Index = func_oceanindex (Lat_regrid_2D_2, Lon_regrid_2D_2) # [0=land] [2=Pacific] [3=Indian Ocean] [6=Atlantic] [10=Arctic] [8=Baffin Bay (west of Greenland)] [9=Norwegian Sea (east of Greenland)] [11=Hudson Bay (Canada)] 
for ii in range(SST_annual.shape[0]):
    SST_annual[ii,:,:][ Ocean_Index != 2 ] = np.nan

# Empirical Orthogonal Functions (EOFs) - The function is saved in Behzadlib code in this directory - imported at the begenning
EOF_spatial_pattern, EOF_time_series, EOF_variance_prcnt = func_EOF (SST_annual, Lat_regrid_2D_2)

## ButterWorth Filtering of high frequency (less than 3 yr time-period) signals in the EOF time series
from scipy import stats, signal
CutOff_T = 3 # Cut-off period 
n_order = 4 # Order of filtering
fs = 1  # Sampling frequency, equal to 1 year in our case
fc = 1/CutOff_T  # Cut-off frequency of the filter
ww = fc / (fs / 2) # Normalize the frequency
bb, aa = signal.butter(n_order, ww, 'low')

EOF_time_series_BWfilt = copy.deepcopy(EOF_time_series)
for ii in range(EOF_time_series_BWfilt.shape[0]):
    EOF_time_series_BWfilt[ii,:] = signal.filtfilt(bb, aa, EOF_time_series_BWfilt[ii,:])

### Ploting the Spatial Patterns and Indices
Plot_Var = EOF_spatial_pattern
Plot_prcnt=EOF_variance_prcnt
P_cmap=plt.cm.seismic; P_proj='cyl'; P_lon0=210.; P_latN=90.; P_latS=-90.; P_range=np.linspace(-0.08,0.08,41); P_Lon=Lon_regrid_2D_2; P_Lat=Lat_regrid_2D_2;
n_r=2 ; n_c=2 ; n_t=4
fig=plt.figure()
for M_i in range(n_t):
    ax = fig.add_subplot(n_r,n_c,M_i+1)     
    Var_plot_ii=Plot_Var[M_i,:,:]    
    
    m = Basemap( projection=P_proj, lon_0=P_lon0, llcrnrlon=P_lon0-180, llcrnrlat=P_latS, urcrnrlon=P_lon0+180, urcrnrlat=P_latN)    
    if M_i == (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[True,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,True], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes    
    elif M_i==0 or M_i==n_c or M_i==n_c*2 or M_i==n_c*3 or M_i==n_c*4 or M_i==n_c*5 or M_i==n_c*6 or M_i==n_c*7 or M_i==n_c*8:
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[True,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
    elif M_i >= n_t-n_c and M_i != (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,True], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes
    else:
        m.drawparallels(np.arange(P_latS, P_latN+0.001, 30.),labels=[False,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Latitutes
        m.drawmeridians(np.arange(0,360,90.),labels=[False,False,False,False], linewidth=0.01, color='k', fontsize=16) # labels = [left,right,top,bottom] # Longitudes

    m.fillcontinents(color='0')
    m.drawcoastlines(linewidth=1.0, linestyle='solid', antialiased=1, ax=None, zorder=None)
    m.fillcontinents(color='0.95')
    im=m.contourf(P_Lon, P_Lat, Var_plot_ii, P_range, latlon=True, cmap=P_cmap, extend='both')
    plt.title('EOF #'+str(M_i+1)+'  ,  '+str(round(Plot_prcnt[M_i], 2))+' % of the variance', fontsize=14)
        
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, hspace=0.05, wspace=0.1) # the amount of height/width reserved for space between subplots
cbar_ax = fig.add_axes([0.93, 0.1, 0.015, 0.8]) # [right,bottom,width,height] 
fig.colorbar(im, cax=cbar_ax)
plt.suptitle('Empirical Orthogonal Functions (EOFs) of Pacific Ocean sea surface temperature'+'\n'+'Spatial Pattern maps - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=20)    
plt.show()
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full   
fig.savefig(dir_pwd+'/'+'Fig_EOF_SST_SpatialPattern_'+str(GCM)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


Plot_Var = EOF_time_series_BWfilt
Plot_prcnt=EOF_variance_prcnt
fig, ax = plt.subplots(nrows=2, ncols=2)
for ii in range(4):
    EOF_time_series_plot_norm=(Plot_Var[ii,:]-np.nanmean(Plot_Var[ii,:]))/np.std(Plot_Var[ii,:])    
    #EOF_time_series_plot_norm_rm=runningMeanFast(EOF_time_series_plot_norm, 10)
    EOF_time_series_plot_norm_rm=EOF_time_series_plot_norm
    plt.subplot(2, 2, ii+1)
    n_l=EOF_time_series_plot_norm.shape[0]
    years=np.linspace(0, n_l, n_l)
    plt.plot(years,EOF_time_series_plot_norm_rm, 'k') 
    y2=np.zeros(len(EOF_time_series_plot_norm_rm))
    plt.fill_between(years, EOF_time_series_plot_norm_rm, y2, where=EOF_time_series_plot_norm_rm >= y2, color = 'r', interpolate=True)
    plt.fill_between(years, EOF_time_series_plot_norm_rm, y2, where=EOF_time_series_plot_norm_rm <= y2, color = 'b', interpolate=True)
    plt.axhline(linewidth=1, color='k')
    plt.title('EOF # '+str(ii+1)+'  ,  '+str(round(Plot_prcnt[ii], 2))+' % of the variance', fontsize=18)
    plt.xticks(fontsize = 18); plt.yticks(fontsize = 18)
plt.suptitle('Empirical Orthogonal Functions (EOFs) of Pacific Ocean sea surface temperature'+'\n'+'Time Series (Butterworth filtered, Normalized) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=20)    
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_pwd+'/'+'Fig_EOF_SST_Indices_'+str(GCM)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

#################################################################################################
#%%### North Atlantic Oscillation (NAO) calculated as the 1st EOF of Sea-Level Air Presure ######
#################################################################################################
dir_data_in1 = ('/data2/scratch/cabre/CMIP5/CMIP5_models/atmosphere_physics/') # Directory to raed raw data from
dir_data_in2=(dir_data_in1+ GCM + '/historical/mo/')

year_start=1901
year_end=2000

Var_name='psl' # The variable name to be read from .nc files
dset = xr.open_mfdataset(dir_data_in2+Var_name+'*12.nc')
Data_all = dset[Var_name] 

Lat_orig = dset[lat_t]
Lon_orig = dset[lon_t]  

#Data_NAO = Data_all[ (Data_all.time.dt.year >= year_start ) & (Data_all.time.dt.year <= year_end)]
#Data_NAO = Data_NAO.isel(lat=slice(49, 90), lon=slice(112, 144))  # Latitude range [10,85]  # Longitude range [280,360]
#Lat_NAO = Lat_orig.isel(lat=slice(49, 90))
#Lat_NAO=Lat_NAO.values
#Lon_NAO = Lon_orig.isel(lon=slice(112, 144))
#Lon_NAO=Lon_NAO.values
Data_NAO = Data_all[ (Data_all.time.dt.year >= year_start ) & (Data_all.time.dt.year <= year_end)]
Data_NAO = Data_NAO.sel(lat=np.arange(10,85,2), lon=np.arange(280,360,2), method='nearest')
Lat_NAO = Lat_orig.sel(lat=np.arange(10,85,2), method='nearest') # Latitude range [10,85]
Lat_NAO=Lat_NAO.values
Lon_NAO = Lon_orig.sel(lon=np.arange(280,360,2), method='nearest') # Longitude range [280,360]
Lon_NAO=Lon_NAO.values

Data_NAO = Data_NAO.values # Converts to a Numpy array
Data_NAO [ Data_NAO > 1e19 ] = np.nan
Data_NAO_rannual = Data_NAO.reshape( np.int(Data_NAO.shape[0]/12) ,12,Data_NAO.shape[1],Data_NAO.shape[2]) # Calculating annual average values, from monthly data
Data_NAO_rannual = np.nanmean(Data_NAO_rannual,axis=1)

Lon_NAO_2D, Lat_NAO_2D = np.meshgrid(Lon_NAO, Lat_NAO)
# Empirical Orthogonal Functions (EOFs) - The function is saved in Behzadlib code in this directory - imported at the begenning
EOF_spatial_pattern, EOF_time_series, EOF_variance_prcnt = func_EOF (Data_NAO_rannual, Lat_NAO_2D)
# NAO calculated as the 1st EOF of Sea-Level Air Presure
NAO_spatial_pattern = EOF_spatial_pattern[0,:,:]
NAO_index = EOF_time_series[0,:]


Plot_Var = NAO_spatial_pattern
cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99)))
Plot_range=np.linspace(-cmap_limit,cmap_limit,27)

fig=plt.figure()
P_lon0=360.
m = Basemap( projection='cyl', lon_0=P_lon0, llcrnrlon=P_lon0-120, llcrnrlat=-10., urcrnrlon=P_lon0+60, urcrnrlat=80.)
m.drawparallels(np.arange(-20., 80.+0.001, 20.),labels=[True,False,False,False], linewidth=0.01, color='k', fontsize=20) # labels = [left,right,top,bottom] # Latitutes
m.drawmeridians(np.arange(-180,+180,30.),labels=[False,False,False,True], linewidth=0.01, color='k', fontsize=20) # labels = [left,right,top,bottom] # Longitudes        
m.drawcoastlines(linewidth=1.0, linestyle='solid', color='k', antialiased=1, ax=None, zorder=None)
im=m.contourf(Lon_NAO_2D, Lat_NAO_2D, Plot_Var,Plot_range,latlon=True, cmap=plt.cm.seismic, extend='both')
cbar = m.colorbar(im,"right", size="3%", pad="2%")
cbar.ax.tick_params(labelsize=20) 
plt.show()
plt.title('North Atlantic Oscillation (NAO) calculated as the 1st EOF of Sea-Level Air Presure'+'\n'+'Spatial Pattern map - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_pwd+'/'+'Fig_NAO_SpatialPattern_'+str(GCM)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


Plot_Var = (NAO_index - np.nanmean(NAO_index))/np.std(NAO_index) # NAO index normalized   

fig=plt.figure()
n_l=Plot_Var.shape[0]
years=np.linspace(year_start, year_start+n_l-1, n_l)
plt.plot(years,Plot_Var, 'k') 
y2=np.zeros(len(Plot_Var))
plt.fill_between(years, Plot_Var, y2, where=Plot_Var >= y2, color = 'r', interpolate=True)
plt.fill_between(years, Plot_Var, y2, where=Plot_Var <= y2, color = 'b', interpolate=True)
plt.axhline(linewidth=1, color='k')
plt.xticks(fontsize = 18); plt.yticks(fontsize = 18)
plt.title('North Atlantic Oscillation (NAO) calculated as the 1st EOF of Sea-Level Air Presure'+'\n'+'Time Series - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM), fontsize=20)    
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_pwd+'/'+'Fig_NAO_Indices_'+str(GCM)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


##################################################################################################
###  Calculating Wind Curls using wind stress in latitudinal and longitudinal directions  ########
##################################################################################################
dir_data_in1 = ('/data2/scratch/cabre/CMIP5/CMIP5_models/atmosphere_physics/') # Directory to raed raw data from
dir_data_in2=(dir_data_in1+ GCM + '/historical/mo/')

year_start=1901
year_end=2000

Var_name='tauu' # The variable name to be read from .nc files
dset = xr.open_mfdataset(dir_data_in2+Var_name+'*12.nc')
Data_all = dset[Var_name]

Lat_orig = dset[lat_t]
Lat_orig = Lat_orig.values
Lon_orig = dset[lon_t]
Lon_orig = Lon_orig.values     

Data_monthly = Data_all[ (Data_all.time.dt.year >= year_start ) & (Data_all.time.dt.year <= year_end)]
Data_monthly = Data_monthly.values # Converts to a Numpy array
Data_rannual = Data_monthly.reshape( np.int(Data_monthly.shape[0]/12) ,12,Data_monthly.shape[1],Data_monthly.shape[2])
Data_rannual = np.nanmean(Data_rannual,axis=1)
# Regriding data into 1degree by 1degree fields -  The regridding function is saved in Behzadlib code in this directory - imported at the begenning
Tau_X = func_regrid(Data_rannual, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D)

Var_name='tauv' # The variable name to be read from .nc files
dset = xr.open_mfdataset(dir_data_in2+Var_name+'*12.nc')
Data_all = dset[Var_name] 
Data_monthly = Data_all[ (Data_all.time.dt.year >= year_start ) & (Data_all.time.dt.year <= year_end)]
Data_monthly = Data_monthly.values # Converts to a Numpy array
Data_rannual = Data_monthly.reshape( np.int(Data_monthly.shape[0]/12) ,12,Data_monthly.shape[1],Data_monthly.shape[2])
Data_rannual = np.nanmean(Data_rannual,axis=1)
Tau_Y = func_regrid(Data_rannual, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D)

# Wind_Curl = ( D_Tau_Y / D_X ) - ( D_Tau_X / D_Y ) # D_X = (Lon_1 - Lon_2) * COS(Lat)
Wind_Curl = np.zeros(( Tau_X.shape[0], Tau_X.shape[1], Tau_X.shape[2]))  
for tt in range (0,Tau_X.shape[0]):  
    for ii in range (1,Tau_X.shape[1]-1):
        for jj in range (1,Tau_X.shape[2]-1): # Wind_Curl = ( D_Tau_Y / D_X ) - ( D_Tau_X / D_Y ) # D_X = (Lon_1 - Lon_2) * COS(Lat)
            
            Wind_Curl[tt,ii,jj] = (  ( Tau_Y[tt, ii,jj+1] - Tau_Y[tt, ii,jj-1] ) /  np.absolute(  ( ( Lon_regrid_2D[ii,jj+1] -  Lon_regrid_2D[ii,jj-1] ) * 111321  *  math.cos(math.radians( ( Lat_regrid_2D[ii,jj])))   )  )     )   -   (  ( Tau_X[tt, ii+1,jj] - Tau_X[tt, ii-1,jj] ) / np.absolute( ( ( Lat_regrid_2D[ii+1,jj] -  Lat_regrid_2D[ii-1,jj] ) * 111321 ) )  )

        Wind_Curl[tt,ii,0] = (  ( Tau_Y[tt, ii,1] - Tau_Y[tt, ii,-1] ) /  np.absolute(  ( ( Lon_regrid_2D[ii,1] -  Lon_regrid_2D[ii,-1] ) * 111321  *  math.cos(math.radians( ( Lat_regrid_2D[ii,0])))   )  )     )   -   (  ( Tau_X[tt, ii+1,0] - Tau_X[tt, ii-1,0] ) / np.absolute( ( ( Lat_regrid_2D[ii+1,0] -  Lat_regrid_2D[ii-1,0] ) * 111321 ) )  )
        Wind_Curl[tt,ii,-1] = (  ( Tau_Y[tt, ii,0] - Tau_Y[tt, ii,-2] ) /  np.absolute(  ( ( Lon_regrid_2D[ii,0] -  Lon_regrid_2D[ii,-2] ) * 111321  *  math.cos(math.radians( ( Lat_regrid_2D[ii,-1])))   )  )     )   -   (  ( Tau_X[tt, ii+1,-1] - Tau_X[tt, ii-1,-1] ) / np.absolute( ( ( Lat_regrid_2D[ii+1,-1] -  Lat_regrid_2D[ii-1,-1] ) * 111321 ) )  )

# Wind_Crul / f # f = coriolis parameter = 2Wsin(LAT) , W = 7.292E-5 rad/s
Wind_Curl_f = np.zeros(( Tau_X.shape[0], Tau_X.shape[1], Tau_X.shape[2])) 
for tt in range (0,Tau_X.shape[0]):   
    for ii in range (1,Tau_X.shape[1]-1):
        if np.absolute( Lat_regrid_2D[ii,0] ) >= 5: # Only calulate for Lats > 5N and Lats < 5S, to avoid infinit numbers in equator where f is zero
            for jj in range (1,Tau_X.shape[2]-1): # Wind_Curl = ( D_Tau_Y / D_X ) - ( D_Tau_X / D_Y ) # D_X = (Lon_1 - Lon_2) * COS(Lat)
            
                Wind_Curl_f[tt,ii,jj] = (  ( ( Tau_Y[tt, ii,jj+1] - Tau_Y[tt, ii,jj-1] ) / ( 2*7.292E-5 *  math.sin(math.radians( ( Lat_regrid_2D[ii,jj]))) ) ) /  np.absolute(  ( ( Lon_regrid_2D[ii,jj+1] -  Lon_regrid_2D[ii,jj-1] ) * 111321  *  math.cos(math.radians( ( Lat_regrid_2D[ii,jj])))   )  )     )   -   (  ( ( Tau_X[tt, ii+1,jj] - Tau_X[tt, ii-1,jj] ) / ( 2*7.292E-5 *  math.sin(math.radians( ( Lat_regrid_2D[ii,jj]))) ) ) / np.absolute( ( ( Lat_regrid_2D[ii+1,jj] -  Lat_regrid_2D[ii-1,jj] ) * 111321 ) )  )

            Wind_Curl_f[tt,ii,0] = (  ( ( Tau_Y[tt, ii,1] - Tau_Y[tt, ii,-1] ) / ( 2*7.292E-5 *  math.sin(math.radians( ( Lat_regrid_2D[ii,0]))) ) ) /  np.absolute(  ( ( Lon_regrid_2D[ii,1] -  Lon_regrid_2D[ii,-1] ) * 111321  *  math.cos(math.radians( ( Lat_regrid_2D[ii,0])))   )  )     )   -   (  ( ( Tau_X[tt, ii+1,jj] - Tau_X[tt, ii-1,0] ) / ( 2*7.292E-5 *  math.sin(math.radians( ( Lat_regrid_2D[ii,0]))) ) ) / np.absolute( ( ( Lat_regrid_2D[ii+1,0] -  Lat_regrid_2D[ii-1,0] ) * 111321 ) )  )
            Wind_Curl_f[tt,ii,-1] = (  ( ( Tau_Y[tt, ii,0] - Tau_Y[tt, ii,-2] ) / ( 2*7.292E-5 *  math.sin(math.radians( ( Lat_regrid_2D[ii,-1]))) ) ) /  np.absolute(  ( ( Lon_regrid_2D[ii,0] -  Lon_regrid_2D[ii,-2] ) * 111321  *  math.cos(math.radians( ( Lat_regrid_2D[ii,-1])))   )  )     )   -   (  ( ( Tau_X[tt, ii+1,-1] - Tau_X[tt, ii-1,-1] ) / ( 2*7.292E-5 *  math.sin(math.radians( ( Lat_regrid_2D[ii,-1]))) ) ) / np.absolute( ( ( Lat_regrid_2D[ii+1,-1] -  Lat_regrid_2D[ii-1,-1] ) * 111321 ) )  )

###############################################################################
Plot_Var = np.nanmean(Wind_Curl,axis=0) * 1E7
Plot_Var[ Ocean_Land_mask==0 ]=np.nan # masking over land, so grid cells that fall on land area (value=0) will be deleted
Plot_Var2 = np.nanmean(Tau_X,axis=0)
Plot_Var2 [ Ocean_Land_mask==0 ]=np.nan # masking over land, so grid cells that fall on land area (value=0) will be deleted

cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99)))
Plot_range=np.linspace(-cmap_limit,cmap_limit,27)
Plot_unit='(1E-7 N/m3)'; Plot_title= 'Wind Curl (1E-7 N/m3) - (contour lines = Tau_x) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)

fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.seismic, 'cyl', 210., 80., -80., '-')
im2=m.contour(Lon_regrid_2D, Lat_regrid_2D,Plot_Var2, 20, latlon=True, colors='k')
plt.clabel(im2, fontsize=8, inline=1)
plt.show()
fig.savefig(dir_pwd+'/'+'Fig_Wind_Curl_'+str(GCM)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')

###############################################################################
Lat_regrid_1D_4, Lon_regrid_1D_4, Lat_bound_regrid_4, Lon_bound_regrid_4 = func_latlon_regrid(45, 90, lat_min_regrid, lat_max_regrid, lon_min_regrid, lon_max_regrid)
Lon_regrid_2D_4, Lat_regrid_2D_4 = np.meshgrid(Lon_regrid_1D_4, Lat_regrid_1D_4)
Tau_X_4 = func_regrid(np.nanmean(Tau_X,axis=0), Lat_regrid_2D, Lon_regrid_2D, Lat_regrid_2D_4, Lon_regrid_2D_4)
Tau_Y_4 = func_regrid(np.nanmean(Tau_Y,axis=0), Lat_regrid_2D, Lon_regrid_2D, Lat_regrid_2D_4, Lon_regrid_2D_4)

Plot_Var_f = np.nanmean(Wind_Curl_f,axis=0) * 1E3

cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var_f, 99)))
Plot_range=np.linspace(-cmap_limit,cmap_limit,27)
Plot_unit='(1E-3 N.S/m3.rad)'; Plot_title= 'Curl of (Wind/f) (Ekman upwelling) (1E-3 N.S/m3.rad) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)+'\n(Arrows: wind direction) (contour line: Curl(wind/f)=0)'

fig, m = func_plotmap_contourf(Plot_Var_f, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.seismic, 'cyl', 210., 80., -80., '-')
im2=m.quiver(Lon_regrid_2D_4, Lat_regrid_2D_4, Tau_X_4, Tau_Y_4, latlon=True, pivot='middle')
plt.show()
im3=m.contour(Lon_regrid_2D[25:50,:], Lat_regrid_2D[25:50,:],Plot_Var_f[25:50,:], levels = [0], latlon=True, colors='darkgreen')
plt.clabel(im3, fontsize=8, inline=1)
fig.savefig(dir_pwd+'/'+'Fig_Wind_Curl_f_WQuiver_'+str(GCM)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')


###############################################################
#### Average Arctic Sea Ice Concentration for each month  #####
###############################################################
dir_data_in1 = ('/data2/scratch/cabre/CMIP5/CMIP5_models/ocean_physics/') # Directory to raed raw data from
dir_data_in2=(dir_data_in1+ GCM + '/historical/mo/')

year_start=1991
year_end=2000

Var_name='sic' # The variable name to be read from .nc files
dset = xr.open_mfdataset(dir_data_in2+Var_name+'*12.nc')
Data_all = dset[Var_name]

Lat_orig = dset[lat_t]
Lat_orig = Lat_orig.values
Lon_orig = dset[lon_t]
Lon_orig = Lon_orig.values  

Data_monthly = Data_all[ (Data_all.time.dt.year >= year_start ) & (Data_all.time.dt.year <= year_end)]
Data_monthly = Data_monthly.values # Converts to a Numpy array
Data_monthly_ave = Data_monthly.reshape( np.int(Data_monthly.shape[0]/12) ,12,Data_monthly.shape[1],Data_monthly.shape[2])
Data_monthly_ave = np.nanmean(Data_monthly_ave,axis=0)
# Regriding data into 1degree by 1degree fields -  The regridding function is saved in Behzadlib code in this directory - imported at the begenning
SIC_monthly = func_regrid(Data_monthly_ave, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D)
SIC_monthly [SIC_monthly==0] = np.nan # To mask the ice-free ocean in the map

Time_months = ('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec');
n_r=3; n_c=4 ; n_t=12
fig=plt.figure()
for ii in range(0,n_t):
    ax = fig.add_subplot(n_r,n_c,ii+1)
    Plot_Var=SIC_monthly[ii,:,:]
    m = Basemap( projection='npstere',lon_0=0, boundinglat=30)
    if ii==0 or ii==n_c or ii==n_c*2 or ii==n_c*3 or ii==n_c*4 or ii==n_c*5 or ii==n_c*6 or ii==n_c*7 or ii==n_c*8:
        m.drawmeridians(np.arange(0,360,30), labels=[1,0,0,0], linewidth=1, color='k', fontsize=12) # labels = [left,right,top,bottom] # Longitudes
    elif ii == (n_c*(n_r-1)): # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawmeridians(np.arange(0,360,30), labels=[1,0,0,1], linewidth=1, color='k', fontsize=12) # labels = [left,right,top,bottom] # Longitudes      
    elif ii >= n_t-n_c: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawmeridians(np.arange(0,360,30), labels=[0,0,0,1], linewidth=1, color='k', fontsize=12) # labels = [left,right,top,bottom] # Longitudes      
    else:
        m.drawmeridians(np.arange(0,360,30), labels=[0,0,0,0], linewidth=1, color='k', fontsize=12) # labels = [left,right,top,bottom] # Longitudes
    m.drawparallels(np.arange(-90,90,20), labels=[1,0,0,0], linewidth=1, color='k', fontsize=12) # labels = [left,right,top,bottom] # Latitutes  
    m.drawcoastlines(linewidth=1.0, linestyle='solid', color='k', antialiased=1, ax=None, zorder=None)
    m.fillcontinents(color='0.8')
    im=m.contourf(Lon_regrid_2D, Lat_regrid_2D, Plot_Var,np.linspace(0,100,51) ,latlon=True, cmap=plt.cm.jet, extend='max')      
    plt.title(Time_months[ii])
plt.suptitle( ( 'Arctic monthly Sea Ice concentration - average of '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)), fontsize=18)
plt.subplots_adjust(left=0.15, bottom=0.1, right=0.85, top=0.92, hspace=0.2, wspace=0.1) # the amount of height/width reserved for space between subplots
cbar_ax = fig.add_axes([0.87, 0.1, 0.015, 0.82]) # [right,bottom,width,height] 
fig.colorbar(im, cax=cbar_ax)
plt.show()
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_pwd+'/'+'Fig_SeaIce_Arctic_monthly_'+str(GCM)+'.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
#plt.close()















