#############################################
####   CMIP5 - Sea Surface Temperature   ####
#############################################
####     Behzad Asadieh, Ph.D.      ####
####  University of Pennsylvania    ####
####    basadieh@sas.upenn.edu      ####
########################################

from Behzadlib import func_latlon_regrid, func_regrid
########################################
import numpy as np
from numpy import zeros, ones, empty, nan, shape
from numpy import isnan, nanmean, nanmax, nanmin
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
########################################

GCM_Names = ['GFDL-ESM2M', 'GFDL-ESM2G', 'HadGEM2-ES','IPSL-CM5A-MR', 'IPSL-CM5A-LR', 'MIROC-ESM', 'MIROC-ESM-CHEM', 'CESM1-BGC', 'CMCC-CESM', 'CanESM2', 'GISS-E2-H-CC', 'GISS-E2-R-CC', 'MPI-ESM-MR', 'MPI-ESM-LR', 'NorESM1-ME']

dir_pwd = os.getcwd() # Gets the current directory (and in which the code is placed)
dir_data_in1 = ('/data2/scratch/cabre/CMIP5/CMIP5_models/ocean_physics/') # Directory to raed raw data from
dir_figs = (dir_pwd + '/Figures_SST/') # Directory to save processed data

### Variables to be edited by user ###
start_date_cal_hist='1980' # Start Date for Calculations # Define as string, not number
end_date_cal_hist='1999' # End Date for Calculations
Var_name='thetao' # The variable name to be read from .nc files

### Regrdridding calculations ###
# creating new coordinate grid, same which was used in interpolation in data processing code
lat_n_regrid, lon_n_regrid = 180, 360 # Number of Lat and Lon elements in the regridded data
lon_min_regrid, lon_max_regrid = 0, 360 # Min and Max value of Lon in the regridded data
lat_min_regrid, lat_max_regrid = -90, 90 # Min and Max value of Lat in the regridded data

Lat_regrid_1D, Lon_regrid_1D, Lat_bound_regrid, Lon_bound_regrid = func_latlon_regrid(lat_n_regrid, lon_n_regrid, lat_min_regrid, lat_max_regrid, lon_min_regrid, lon_max_regrid)
Lon_regrid_2D, Lat_regrid_2D = np.meshgrid(Lon_regrid_1D, Lat_regrid_1D)

Multimodel_Variable_Surface_Ave_Regrid_hist_plt=zeros((len(GCM_Names), lat_n_regrid, lon_n_regrid))# Multimodel surface average of specified variable, regridded
Multimodel_Time_AllYears_hist=[]
Multimodel_InputFileNames_hist=[]

for M_i in range(len(GCM_Names)):
    
    GCM=GCM_Names[M_i]
    
    dir_data_in2=(dir_data_in1+ GCM + '/historical/mo/')
    Input_File_Names = [xx for xx in sorted(os.listdir(dir_data_in2)) if xx.startswith(Var_name) and xx.endswith(".nc")] # List all the files in the directory that are .nc and end with year number (to avoid concated files)
    if GCM=='HadGEM2-ES' or GCM=='MIROC-ESM': # Some models have only concated data which the file name ends with concat.nc, so the year characters in their file name is different
        Input_File_Names = [xx for xx in Input_File_Names if ( int(xx[-23:-19])>=int(start_date_cal_hist) and int(xx[-23:-19])<=int(end_date_cal_hist) ) or ( int(xx[-16:-12])>=int(start_date_cal_hist) and int(xx[-16:-12])<=int(end_date_cal_hist) ) or ( int(xx[-23:-19])<=int(start_date_cal_hist) and int(xx[-16:-12])>=int(end_date_cal_hist) )] # Keep only the files that the time range is in the specified time interval
    else:
        Input_File_Names = [xx for xx in Input_File_Names if not xx.endswith("concat.nc") ] # Some models have both decadal files and a concated file, which may result in duplication
        Input_File_Names = [xx for xx in Input_File_Names if ( int(xx[-16:-12])>=int(start_date_cal_hist) and int(xx[-16:-12])<=int(end_date_cal_hist) ) or ( int(xx[-9:-5])>=int(start_date_cal_hist) and int(xx[-9:-5])<=int(end_date_cal_hist) ) or ( int(xx[-16:-12])<=int(start_date_cal_hist) and int(xx[-9:-5])>=int(end_date_cal_hist) )] # Keep only the files that the time range is in the specified time interval

    dir_data_in_file=(dir_data_in2 +Input_File_Names[0])
    dset = Dataset(dir_data_in_file)
    ## dset.variables  # Shows the variables in the .nc file
    
    if GCM=='HadGEM2-ES' or GCM=='MIROC-ESM' :
        lat_char='LAT'
        lon_char='LON'
        time_char='TIME'
        Var_char='THETAO'
    else:
        lat_char='lat'
        lon_char='lon'
        time_char='time'
        Var_char='thetao'
    
    # Reading lat and lon values from the first file
    Lat=np.asarray(dset.variables[lat_char][:])
    Lon=np.asarray(dset.variables[lon_char][:])
    
    Variable_Surface_AllYears=[]
    Time_AllYears=[]   
    for F_i in  range(len(Input_File_Names[:])): # F_i=0
        
        try:#open netcdf file #MFDataset function can read multiple netcdf files, note * - wildcard which is used to do it                
            dir_data_in_file=(dir_data_in2 +Input_File_Names[F_i])
            dset = Dataset(dir_data_in_file)
            ## dset.variables  # Shows the variables in the .nc file
        except:
            print ('There is no such'+ GCM)
            continue
        
        Time_dset=dset.variables[time_char] # append time to variable
        Var_dset = dset.variables[Var_char] # append thetao data to variable
        # next lines of code are looking for specified dates
        start_date_i=0 # Conter of start_date (to be calculated below)
        end_date_i=0 # Conter of end_date (to be calculated below)
        text = Time_dset.units.split(' ') # Splits the Time_dset.UNITS wherever there's a ' ' separator
        Data_start_date=text[2][:4] # From TEXT, the (2+1)th row and the first 4 chars show the begenning year of data 
        # the loop goes through all the time indeces and checks if the time equals desired value
        ## This loop finds which row in the Time_dset variable is in the desired start_date_cal and end_date_cal range
        for ii in range(len(Time_dset[:])): 
            # converting time to date 
            date=str(num2date(Time_dset[ii],units=Time_dset.units,calendar=Time_dset.calendar))
            date_sp = date.split('-') # Splits the TIMES.UNITS wherever there's a '-' separator
            # if date_sp[0] which is year equals specified start year, append index of the year to start_date variable
            if (int(date_sp[0])>=int(start_date_cal_hist) and int(date_sp[0])<=int(end_date_cal_hist)): # Splits the DATE wherever there's a '-' separator and saves it as date_yr . Now the first row of WORDS would show the year
                Variable_Surface_AllYears.append(np.asarray(Var_dset[ii,0,:,:])) # Variable(time, lev, rlat, rlon) - Variable of all the desired years, only on the surface
                Time_AllYears.append(date)
    dset.close() # closing netcdf file    
    Variable_Surface_AllYears = np.asarray(Variable_Surface_AllYears) # Converts the appended list to an array
    
    Variable_Surface_Ave=nanmean(Variable_Surface_AllYears,0) # Averages over the index 0 (over time)
    Variable_Surface_Ave_mask=ma.masked_where(Variable_Surface_Ave >= 1e19, Variable_Surface_Ave) # Masking grids with nodata(presented by value=1E20 by GCM)
    
    # Regriding data into 1degree by 1degree fields
    Data_regrid = func_regrid(Variable_Surface_Ave, Lat, Lon, Lat_regrid_2D, Lon_regrid_2D)
    
    Variable_Surface_Ave_Regrid=np.asarray(Data_regrid) # Converts to numpy array - fills empty (masked) gridcells with nan
    Multimodel_Variable_Surface_Ave_Regrid_hist_plt[M_i,:,:]=Variable_Surface_Ave_Regrid

    #Multimodel_Time_AllYears[M_i,:]=Time_AllYears
    Multimodel_Time_AllYears_hist.append(Time_AllYears)
    Multimodel_InputFileNames_hist.append(Input_File_Names)

    print('Model '+str(GCM)+' - hist - processed successfully')

#################################################
#### Ploting for all models (with pcolormesh) ###
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
    Plot_Var=Multimodel_Variable_Surface_Ave_Regrid_hist_plt[ii,:,:] - 273.15
    m = Basemap(projection='cyl', lat_0=0, lon_0=0)
    m.drawcoastlines(linewidth=1.25)
    m.fillcontinents(color='0.95')
    #m.drawmapboundary(fill_color='0.9')
    m.drawparallels(np.arange(-90.,90.001,30.),labels=[True,False,False,False], linewidth=0.01) # labels = [left,right,top,bottom]
    if ii+1 >= n_range[-n_c+1]: # Adds longitude ranges only to the last subplots that appear at the bottom of plot
        m.drawmeridians(np.arange(0.,360.,60.),labels=[False,False,False,True], linewidth=0.01) # labels = [left,right,top,bottom]
    im1 = m.pcolormesh(Lon_regrid_2D, Lat_regrid_2D, Plot_Var, norm=norm, shading='flat', cmap=plt.cm.jet, latlon=True) # Choose colormap: https://matplotlib.org/users/colormaps.html
    plt.title(GCM_Names[ii])
plt.suptitle( ( 'CMIP5 Sea Surface Temperature - hist - average of '+start_date_cal_hist+'-'+end_date_cal_hist ), fontsize=18)
cbar = plt.colorbar(cax=plt.axes([0.93, 0.1, 0.015, 0.8]), extend='max') # cax = [left position, bottom postion, width, height] 
cbar.set_label(Var_plot_unit)
plt.show()
mng = plt.get_current_fig_manager()
mng.window.showMaximized() # Maximizes the plot window to save figures in full
#fig.savefig(dir_pwd+Var_name+'_AllGCMs_surface_hist-' +start_date_cal_hist+'-'+end_date_cal_hist + '.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
#plt.close()




