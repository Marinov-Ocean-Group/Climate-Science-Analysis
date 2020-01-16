#############################################
####   CMIP5 - Sea Surface Temperature   ####
#############################################
####     Behzad Asadieh, Ph.D.      ####
####  University of Pennsylvania    ####
####    basadieh@sas.upenn.edu      ####
########################################

from Behzadlib import func_latlon_regrid, func_regrid, func_oceanlandmask, func_gridcell_area
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
import copy
########################################

GCM_Names = ['GFDL-ESM2M', 'GFDL-ESM2G', 'IPSL-CM5A-MR', 'IPSL-CM5A-LR', 'MIROC-ESM', 'MIROC-ESM-CHEM', 'CESM1-BGC', 'CMCC-CESM', 'CanESM2', 'GISS-E2-H-CC', 'GISS-E2-R-CC', 'MPI-ESM-MR', 'MPI-ESM-LR', 'NorESM1-ME']

dir_pwd = os.getcwd() # Gets the current directory (and in which the code is placed)
dir_data_in1 = ('/data2/scratch/cabre/CMIP5/CMIP5_models/ocean_physics/') # Directory to raed raw data from
dir_figs = (dir_pwd + '/Figures_SST/') # Directory to save processed data

### Variables to be edited by user ###
start_date_hist=1980 # Start Date for Calculations - Historical
end_date_hist=1999 # End Date for Calculations - Historical
start_date_rcp8p5=2080 # Start Date for Calculations - RCP8.5
end_date_rcp8p5=2099 # End Date for Calculations - RCP8.5
Var_name='thetao' # The variable name to be read from .nc files

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

Multimodel_Variable_Surface_Ave_hist=np.zeros((len(GCM_Names), lat_n_regrid, lon_n_regrid))# Multimodel surface average of specified variable, regridded - Historical
Multimodel_Variable_Surface_Ave_rcp8p5=np.zeros((len(GCM_Names), lat_n_regrid, lon_n_regrid))# Multimodel surface average of specified variable, regridded - RCP8.5

######################################
### Historical Period Calculations ###
######################################
for M_i in range(len(GCM_Names)): # M_i=1
    
    GCM=GCM_Names[M_i]
    dir_data_in2=(dir_data_in1+ GCM + '/historical/mo/')

    if GCM=='HadGEM2-ES' or GCM=='MIROC-ESM' : # These two models have only concated data (the file name ends with concat.nc) - and the coordinates are upper case
        lat_t='LAT'
        lon_t='LON'
        time_t='TIME'
        Var_t='THETAO'
        dset = xr.open_mfdataset(dir_data_in2+Var_name+'*concat.nc')
        Data_all = dset[Var_t].sel(LEV=0,method='nearest') # data at lev=0        
        
    else:
        lat_t='lat'
        lon_t='lon'
        time_t='time'
        Var_t='thetao'  
        dset = xr.open_mfdataset(dir_data_in2+Var_name+'*12.nc')
        Data_all = dset[Var_t].sel(lev=0,method='nearest') # data at lev=0        

    Lat_orig = dset[lat_t]
    Lat_orig = Lat_orig.values
    Lon_orig = dset[lon_t]
    Lon_orig = Lon_orig.values    
    
    Data_regrid = Data_all[ (Data_all[time_t].dt.year >= start_date_hist ) & (Data_all[time_t].dt.year <= end_date_hist)]
    Data_regrid = Data_regrid.values # Converts to a Numpy array
    Data_regrid = np.nanmean(Data_regrid,axis=0)
    Data_regrid [ Data_regrid > 1e19 ] = np.nan
    dset.close()
    
    # Regriding data into 1degree by 1degree fields
    # The regridding function is saved in Behzadlib code in this directory - imported at the begenning
    Data_regrid = func_regrid(Data_regrid, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D)  
    
    Multimodel_Variable_Surface_Ave_hist[M_i,:,:]=Data_regrid

    print('Model '+str(GCM)+' - hist - processed')
    
##########################################
### 21st century (RCP8.5) Calculations ###
##########################################
for M_i in range(len(GCM_Names)): # M_i=1
    
    GCM=GCM_Names[M_i]
    dir_data_in2=(dir_data_in1+ GCM + '/rcp85/mo/')

    if GCM=='MIROC-ESM-CHEM': # These two models have only concated data (the file name ends with concat.nc) - and the coordinates are upper case
        lat_t='LAT'
        lon_t='LON'
        time_t='TIME'
        Var_t='THETAO'
        dset = xr.open_mfdataset(dir_data_in2+Var_name+'*concat.nc')
        Data_all = dset[Var_t].sel(LEV=0,method='nearest') # data at lev=0        
        
    else:
        lat_t='lat'
        lon_t='lon'
        time_t='time'
        Var_t='thetao'  
        dset = xr.open_mfdataset(dir_data_in2+Var_name+'*12.nc')
        Data_all = dset[Var_t].sel(lev=0,method='nearest') # data at lev=0        

    Lat_orig = dset[lat_t]
    Lat_orig = Lat_orig.values
    Lon_orig = dset[lon_t]
    Lon_orig = Lon_orig.values    
    
    Data_regrid = Data_all[ (Data_all[time_t].dt.year >= start_date_rcp8p5 ) & (Data_all[time_t].dt.year <= end_date_rcp8p5)]
    Data_regrid = Data_regrid.values # Converts to a Numpy array
    Data_regrid = np.nanmean(Data_regrid,axis=0)
    Data_regrid [ Data_regrid > 1e19 ] = np.nan
    dset.close()
    
    # Regriding data into 1degree by 1degree fields
    # The regridding function is saved in Behzadlib code in this directory - imported at the begenning
    Data_regrid = func_regrid(Data_regrid, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D)  
    
    Multimodel_Variable_Surface_Ave_rcp8p5[M_i,:,:]=Data_regrid

    print('Model '+str(GCM)+' - rcp8.5 - processed')

###################################################
#### Ploting for all models - Historical Period ###
n_r=4 # Number of rows for subplot
n_c=4 # Number of columns for subplot
n_range=list(range(len(GCM_Names)))
bounds_max=33
bounds = np.arange(0, bounds_max, bounds_max/33)
norm = matplotlib.colors.BoundaryNorm(boundaries=bounds, ncolors=256)
Var_plot_unit='Unit = °C'

fig=plt.figure()
for ii in n_range:
    ax = fig.add_subplot(n_r,n_c,ii+1)
    Plot_Var=Multimodel_Variable_Surface_Ave_hist[ii,:,:] - 273.15
    m = Basemap(projection='cyl', lat_0=0, lon_0=0)
    m.drawcoastlines(linewidth=1.25)
    m.fillcontinents(color='0.95')
    #m.drawmapboundary(fill_color='0.9')
    m.drawparallels(np.arange(-90.,90.001,30.),labels=[True,False,False,False], linewidth=0.01) # labels = [left,right,top,bottom]
    if ii+1 >= n_range[-n_c+1]: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawmeridians(np.arange(0.,360.,60.),labels=[False,False,False,True], linewidth=0.01) # labels = [left,right,top,bottom]
    im1 = m.pcolormesh(Lon_regrid_2D, Lat_regrid_2D, Plot_Var, norm=norm, shading='flat', cmap=plt.cm.jet, latlon=True) # Choose colormap: https://matplotlib.org/users/colormaps.html
    plt.title(GCM_Names[ii])
plt.suptitle( ( 'CMIP5 Sea Surface Temperature - hist - average of '+str(start_date_hist)+'-'+str(end_date_hist) ), fontsize=18)
cbar = plt.colorbar(cax=plt.axes([0.93, 0.1, 0.015, 0.8]), extend='max') # cax = [left position, bottom postion, width, height] 
cbar.set_label(Var_plot_unit)
plt.show()
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_pwd+'/'+'Fig_CMIP5_SST_hist_' +str(start_date_hist)+'_'+str(end_date_hist) + '.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
#plt.close()

#######################################################
#### Ploting for all models - 21st century (RCP8.5) ###
n_r=4 # Number of rows for subplot
n_c=4 # Number of columns for subplot
n_range=list(range(len(GCM_Names)))
bounds_max=33
bounds = np.arange(0, bounds_max, bounds_max/33)
norm = matplotlib.colors.BoundaryNorm(boundaries=bounds, ncolors=256)
Var_plot_unit='Unit = °C'

fig=plt.figure()
for ii in n_range:
    ax = fig.add_subplot(n_r,n_c,ii+1)
    Plot_Var=Multimodel_Variable_Surface_Ave_rcp8p5[ii,:,:] - 273.15
    m = Basemap(projection='cyl', lat_0=0, lon_0=0)
    m.drawcoastlines(linewidth=1.25)
    m.fillcontinents(color='0.95')
    #m.drawmapboundary(fill_color='0.9')
    m.drawparallels(np.arange(-90.,90.001,30.),labels=[True,False,False,False], linewidth=0.01) # labels = [left,right,top,bottom]
    if ii+1 >= n_range[-n_c+1]: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawmeridians(np.arange(0.,360.,60.),labels=[False,False,False,True], linewidth=0.01) # labels = [left,right,top,bottom]
    im1 = m.pcolormesh(Lon_regrid_2D, Lat_regrid_2D, Plot_Var, norm=norm, shading='flat', cmap=plt.cm.jet, latlon=True) # Choose colormap: https://matplotlib.org/users/colormaps.html
    plt.title(GCM_Names[ii])
plt.suptitle( ( 'CMIP5 Sea Surface Temperature - RCP8.5 - average of '+str(start_date_rcp8p5)+'-'+str(end_date_rcp8p5) ), fontsize=18)
cbar = plt.colorbar(cax=plt.axes([0.93, 0.1, 0.015, 0.8]), extend='max') # cax = [left position, bottom postion, width, height] 
cbar.set_label(Var_plot_unit)
plt.show()
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_pwd+'/'+'Fig_CMIP5_SST_rcp8p5_' +str(start_date_rcp8p5)+'_'+str(end_date_rcp8p5) + '.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
#plt.close()

#######################################################
### Ploting change in SST - rcp85 minus hist ###
n_r=4 # Number of rows for subplot
n_c=4 # Number of columns for subplot
n_range=list(range(len(GCM_Names)))
bounds_max=6
bound_ranges=bounds_max/20
bounds = np.arange(-1*bounds_max, bounds_max+bound_ranges, bound_ranges)
norm = matplotlib.colors.BoundaryNorm(boundaries=bounds, ncolors=256)
Var_plot_unit='Unit = °C'

fig=plt.figure()
for ii in n_range:
    ax = fig.add_subplot(n_r,n_c,ii+1)
    Plot_Var=Multimodel_Variable_Surface_Ave_rcp8p5[ii,:,:] - Multimodel_Variable_Surface_Ave_hist[ii,:,:]
    m = Basemap(projection='cyl', lat_0=0, lon_0=0)
    m.drawcoastlines(linewidth=1.25)
    m.fillcontinents(color='0.95')
    #m.drawmapboundary(fill_color='0.9')
    m.drawparallels(np.arange(-90.,90.001,30.),labels=[True,False,False,False], linewidth=0.01) # labels = [left,right,top,bottom]
    if ii+1 >= n_range[-n_c+1]: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawmeridians(np.arange(0.,360.,60.),labels=[False,False,False,True], linewidth=0.01) # labels = [left,right,top,bottom]
    im1 = m.pcolormesh(Lon_regrid_2D, Lat_regrid_2D, Plot_Var, norm=norm, shading='flat', cmap=plt.cm.RdBu_r, latlon=True) # Choose colormap: https://matplotlib.org/users/colormaps.html
    
    plt.title(GCM_Names[ii])
plt.suptitle( ('Climate Change impact on sea surface temperature under RCP8.5 Scenario - '+str(start_date_rcp8p5)+'-'+str(end_date_rcp8p5)+' minus '+ str(start_date_hist)+'-'+str(end_date_hist)), fontsize=18)
cbar = plt.colorbar(cax=plt.axes([0.93, 0.1, 0.015, 0.8]), extend='both') # cax = [left position, bottom postion, width, height] 
cbar.set_label(Var_plot_unit)
plt.show()
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_pwd+'/'+'Fig_CMIP5_SST_climate_change_Impact_RCP8p5.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
#plt.close()

####################################################################
### Change in global SST in 21st century vs 20st centiru average ###

# Calculate Grid Cell areas in million km2 - This function is saved in Behzadlib code in this directory - imported at the begenning
GridCell_Area = func_gridcell_area(Lat_bound_regrid, Lon_bound_regrid)
GridCell_Area [ Ocean_Land_mask == 0] = np.nan # masking over land, so grid cells that fall on land area will be deleted

# Calculating grid-area-weighted global averages in SST in historical period and 21st century period, and their difference
SST_hist=np.zeros((len(GCM_Names))) # Historical average of SST
Delta_SST=np.zeros((len(GCM_Names))) # Historical average of SST minus 21st century average of SST
for ii in range(len(GCM_Names)):
    SST_hist[ii]= np.nansum( np.multiply(Multimodel_Variable_Surface_Ave_hist[ii,:,:]- 273.15, GridCell_Area) ) / np.nansum(GridCell_Area)
    Delta_SST[ii]= (np.nansum( np.multiply(Multimodel_Variable_Surface_Ave_rcp8p5[ii,:,:], GridCell_Area) ) / np.nansum(GridCell_Area)) - (np.nansum( np.multiply(Multimodel_Variable_Surface_Ave_hist[ii,:,:], GridCell_Area) ) / np.nansum(GridCell_Area))

fig, ax = plt.subplots()
ax.scatter(SST_hist, Delta_SST, s=200, marker='d', c='r')
ax.scatter(SST_hist, Delta_SST, s=20, marker='d', c='b')
for ii, txt in enumerate(GCM_Names):
    ax.annotate(txt, (SST_hist[ii],Delta_SST[ii]), fontsize=14)
plt.xlabel('Global SST (hist ave) [°C]', fontsize=18)
plt.xlim(17, 19.5)
plt.xticks( fontsize = 18)
plt.ylabel('Δ SST (rcp8.5 minus hist) [°C]', fontsize=18)
plt.ylim(1.5, 4)
plt.yticks( fontsize = 18)
plt.title( ('Global Change in SST under climate change (rcp8.5 '+str(start_date_rcp8p5)+'-'+str(end_date_rcp8p5)+' minus '+ str(start_date_hist)+'-'+str(end_date_hist))+') VS. '+str(start_date_hist)+'-'+str(end_date_hist)+ ' average', fontsize=18)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
fig.savefig(dir_pwd+'/'+'Fig_CMIP5_SST_climate_change_Impact_RCP8p5_scatter.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
#plt.close()




