############################################
####   Data Analysis/Plotting Library   ####
############################################
####     Behzad Asadieh, Ph.D.      ####
####  University of Pennsylvania    ####
####    basadieh@sas.upenn.edu      ####
####     github.com/behzadasd       ####
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


def func_latlon_regrid(lat_n_regrid, lon_n_regrid, lat_min_regrid, lat_max_regrid, lon_min_regrid, lon_max_regrid): 
    # This function centers lat at [... -1.5, -0.5, +0.5, +1.5 ...] - No specific equator lat cell
    ###lat_n_regrid, lon_n_regrid = 180, 360 # Number of Lat and Lon elements in the regridded data
    ###lon_min_regrid, lon_max_regrid = 0, 360 # Min and Max value of Lon in the regridded data
    ###lat_min_regrid, lat_max_regrid = -90, 90 # Min and Max value of Lat in the regridded data
    ####creating arrays of regridded lats and lons ###
    #### Latitude Bounds ####
    Lat_regrid_1D=zeros ((lat_n_regrid));
    Lat_bound_regrid = zeros ((lat_n_regrid,2)); Lat_bound_regrid[0,0]=-90;  Lat_bound_regrid[0,1]=Lat_bound_regrid[0,0] + (180/lat_n_regrid); Lat_regrid_1D[0]=(Lat_bound_regrid[0,0]+Lat_bound_regrid[0,1])/2
    for ii in range(1,lat_n_regrid):
        Lat_bound_regrid[ii,0]=Lat_bound_regrid[ii-1,1]
        Lat_bound_regrid[ii,1]=Lat_bound_regrid[ii,0] +  (180/lat_n_regrid)
        Lat_regrid_1D[ii]=(Lat_bound_regrid[ii,0]+Lat_bound_regrid[ii,1])/2
    #### Longitude Bounds ####
    Lon_regrid_1D=zeros ((lon_n_regrid));
    Lon_bound_regrid = zeros ((lon_n_regrid,2)); Lon_bound_regrid[0,0]=0;  Lon_bound_regrid[0,1]=Lon_bound_regrid[0,0] + (360/lon_n_regrid); Lon_regrid_1D[0]=(Lon_bound_regrid[0,0]+Lon_bound_regrid[0,1])/2
    for ii in range(1,lon_n_regrid):
        Lon_bound_regrid[ii,0]=Lon_bound_regrid[ii-1,1]
        Lon_bound_regrid[ii,1]=Lon_bound_regrid[ii,0] +  (360/lon_n_regrid)
        Lon_regrid_1D[ii]=(Lon_bound_regrid[ii,0]+Lon_bound_regrid[ii,1])/2
    
    return Lat_regrid_1D, Lon_regrid_1D, Lat_bound_regrid, Lon_bound_regrid


def func_gridcell_area(Lat_bound_regrid, Lon_bound_regrid): 
    #### Calculate Grid Cell areas in million km2 ####
    earth_R = 6378 # Earth Radius - Unit is kilometer (km)
    GridCell_Area = empty((Lat_bound_regrid.shape[0], Lon_bound_regrid.shape[0] )) *nan
    for ii in range(Lat_bound_regrid.shape[0]):
        for jj in range(Lon_bound_regrid.shape[0]):
            GridCell_Area [ii,jj] = math.fabs( (earth_R**2) * (math.pi/180) * (Lon_bound_regrid[jj,1] - Lon_bound_regrid[jj,0])  * ( math.sin(math.radians(Lat_bound_regrid[ii,1])) - math.sin(math.radians(Lat_bound_regrid[ii,0]))) )
    GridCell_Area = GridCell_Area / 1e6 # to convert the area to million km2
    
    return GridCell_Area


def func_oceanlandmask(Lat_regrid_2D, Lon_regrid_2D):
    lat_n_regrid, lon_n_regrid =Lat_regrid_2D.shape[0], Lat_regrid_2D.shape[1]
    Ocean_Land_mask = empty ((lat_n_regrid, lon_n_regrid)) * nan
    ocean_mask= maskoceans(Lon_regrid_2D-180, Lat_regrid_2D, Ocean_Land_mask)
    for ii in range(lat_n_regrid):
        for jj in range(lon_n_regrid):
            if ma.is_masked(ocean_mask[ii,jj]):
                Ocean_Land_mask[ii,jj]=1 # Land_Ocean_mask=1 means grid cell is ocean (not on land)
            else:
                Ocean_Land_mask[ii,jj]=0 # Land_Ocean_mask=0 means grid cell is land
    land_mask2 = copy.deepcopy ( Ocean_Land_mask ) # The created land_mask's longitude is from -180-180 - following lines transfer it to 0-360
    Ocean_Land_mask=empty((Lat_regrid_2D.shape[0], Lat_regrid_2D.shape[1])) *nan
    Ocean_Land_mask[:,0:int(Ocean_Land_mask.shape[1]/2)]=land_mask2[:,int(Ocean_Land_mask.shape[1]/2):]
    Ocean_Land_mask[:,int(Ocean_Land_mask.shape[1]/2):]=land_mask2[:,0:int(Ocean_Land_mask.shape[1]/2)]
    
    return Ocean_Land_mask # 1= ocean, 0= land

def func_oceanindex (Lat_regrid_2D, Lon_regrid_2D):
    
    Ocean_Land_mask = func_oceanlandmask(Lat_regrid_2D, Lon_regrid_2D) # 1= ocean, 0= land
    
    #directory= '/data1/home/basadieh/behzadcodes/behzadlibrary/'
    directory = os.path.dirname(os.path.realpath(__file__)) # Gets the directory where the code is located - The gx3v5_OceanIndex.nc should be placed in the same directory
    file_name='gx3v5_OceanIndex.nc'
    dset_n = Dataset(directory+'/'+file_name)
    
    REGION_MASK=np.asarray(dset_n.variables['REGION_MASK'][:])
    TLAT=np.asarray(dset_n.variables['TLAT'][:])
    TLONG=np.asarray(dset_n.variables['TLONG'][:])
    
    REGION_MASK_regrid = func_regrid(REGION_MASK, TLAT, TLONG, Lat_regrid_2D, Lon_regrid_2D)
    Ocean_Index = copy.deepcopy(REGION_MASK_regrid)    
    for tt in range(0,6): # Smoothing the coastal gridcells - If a cell in the regrid has fallen on land but in Ocean_Land_mask it's in ocean, a neighboring Ocean_Index value will be assigned to it
        for ii in range(Ocean_Index.shape[0]):
            for jj in range(Ocean_Index.shape[1]):
                
                if Ocean_Index[ii,jj] == 0 and Ocean_Land_mask[ii,jj] == 1:
                    if ii>2 and jj>2:
                        Ocean_Index[ii,jj] = np.max(Ocean_Index[ii-1:ii+2,jj-1:jj+2])

    Ocean_Index[ np.where( np.logical_and(   np.logical_and( 0 <= Lon_regrid_2D , Lon_regrid_2D < 20 ) , np.logical_and( Lat_regrid_2D <= -30 , Ocean_Index != 0 )   ) ) ] = 6 ## Assigning Atlantic South of 30S to Atlantic Ocean Index
    Ocean_Index[ np.where( np.logical_and(   np.logical_and( 290 <= Lon_regrid_2D , Lon_regrid_2D <= 360 ) , np.logical_and( Lat_regrid_2D <= -30 , Ocean_Index != 0 )   ) ) ] = 6   
    Ocean_Index[ np.where( np.logical_and(   np.logical_and( 20 <= Lon_regrid_2D , Lon_regrid_2D < 150 ) , np.logical_and( Lat_regrid_2D <= -30 , Ocean_Index != 0 )   ) ) ] = 3 ## Assigning Pacifi South of 30S to Atlantic Ocean Index   
    Ocean_Index[ np.where( np.logical_and(   np.logical_and( 150 <= Lon_regrid_2D , Lon_regrid_2D < 290 ) , np.logical_and( Lat_regrid_2D <= -30 , Ocean_Index != 0 )   ) ) ] = 2 ## Assigning Pacifi South of 30S to Atlantic Ocean Index
    
    return Ocean_Index # [0=land] [2=Pacific] [3=Indian Ocean] [6=Atlantic] [10=Arctic] [8=Baffin Bay (west of Greenland)] [9=Norwegian Sea (east of Greenland)] [11=Hudson Bay (Canada)] 
                       # [-7=Mediterranean] [-12=Baltic Sea] [-13=Black Sea] [-5=Red Sea] [-4=Persian Gulf] [-14=Caspian Sea]

def func_regrid(Data_orig, Lat_orig, Lon_orig, Lat_regrid_2D, Lon_regrid_2D):    
    
    Lon_orig[Lon_orig < 0] +=360
    
    if np.ndim(Lon_orig)==1: # If the GCM grid is not curvlinear
        Lon_orig,Lat_orig=np.meshgrid(Lon_orig, Lat_orig)
        
    lon_vec = np.asarray(Lon_orig)
    lat_vec = np.asarray(Lat_orig)
    lon_vec = lon_vec.flatten()
    lat_vec = lat_vec.flatten()
    coords=np.squeeze(np.dstack((lon_vec,lat_vec)))

    Data_orig=np.squeeze(Data_orig)
    if Data_orig.ndim==2:#this is for 2d regridding
        data_vec = np.asarray(Data_orig)
        if np.ndim(data_vec)>1:
            data_vec = data_vec.flatten()
        Data_regrid = griddata(coords, data_vec, (Lon_regrid_2D, Lat_regrid_2D), method='nearest')
        return np.asarray(Data_regrid)
    if Data_orig.ndim==3:#this is for 3d regridding
        Data_regrid=[]
        for d in range(len(Data_orig)):
            z = np.asarray(Data_orig[d,:,:])
            if np.ndim(z)>1:
                z = z.flatten()
            zi = griddata(coords, z, (Lon_regrid_2D, Lat_regrid_2D), method='nearest')
            Data_regrid.append(zi)
        return np.asarray(Data_regrid)


def func_EOF (Calc_Var, Calc_Lat): # Empirical Orthogonal Functions maps and indices
#%% Example :
#Calc_Var = data_set_regrid [:,10:61,300:]
#Calc_Lat = Lat_regrid_2D [10:61,300:]
#Calc_Lon = Lon_regrid_2D [10:61,300:]
#
#EOF_spatial_pattern, EOF_time_series, EOF_variance_prcnt = func_EOF (Calc_Var, Calc_Lat)
#%%
    EOF_all=[]
    for i in range(Calc_Var.shape[0]):
        
        print ('EOF calc - Year: ', i+1)
        data_i=Calc_Var[i,:,:]
        data_i=np.squeeze(data_i)        

        data_EOF=[]
        if i==0:
            [lat_ii,lon_jj] = np.where(~np.isnan(data_i))

        for kk in range(len(lat_ii)):
            EOF_i=data_i[lat_ii[kk],lon_jj[kk]]*np.sqrt(np.cos(np.deg2rad(Calc_Lat[lat_ii[kk],lon_jj[kk]])))
            data_EOF.append(EOF_i)
    
        EOF_all.append(data_EOF)    
    
    EOF_all=np.asarray(EOF_all)
    
    C=np.cov(np.transpose(EOF_all))
    #C= np.array(C, dtype=np.float32)
    eigval,eigvec=np.linalg.eig(C)
    eigval=np.real(eigval)
    eigvec=np.real(eigvec)
    
    EOF_spatial_pattern = empty((10,Calc_Var.shape[1],Calc_Var.shape[2]))*nan # Stores first 10 EOFs for spatial pattern map
    for ss in range(EOF_spatial_pattern.shape[0]):
        for kk in range(len(lat_ii)):
            EOF_spatial_pattern[ss,lat_ii[kk],lon_jj[kk]] = eigvec[kk,ss]

    EOF_time_series = empty((10,Calc_Var.shape[0]))*nan # Stores first 10 EOFs times series
    for ss in range(EOF_time_series.shape[0]):
        EOF_time_series[ss,:] = np.dot(np.transpose(eigvec[:,ss]),np.transpose(EOF_all))
        
    EOF_variance_prcnt = empty((10))*nan # Stores first 10 EOFs variance percentage
    for ss in range(EOF_variance_prcnt.shape[0]):
        EOF_variance_prcnt[ss]=( eigval[ss]/np.nansum(eigval,axis=0) ) * 100        

    return EOF_spatial_pattern, EOF_time_series, EOF_variance_prcnt


def func_plotmap_contourf(P_Var, P_Lon, P_Lat, P_range, P_title, P_unit, P_cmap, P_proj, P_lon0, P_latN, P_latS, P_c_fill):
### P_Var= Plotting variable, 2D(lat,lon) || P_Lon=Longitude, 2D || P_range=range of plotted values, can be vector or number || P_title=Plot title || P_unit=Plot colorbar unit
### P_cmap= plt.cm.seismic , plt.cm.jet || P_proj= 'cyl', 'npstere', 'spstere' || P_lon0=middle longitude of plot || P_latN=upper lat bound of plot || P_latS=lower lat bound of plot || P_c_fill= 'fill' fills the continets with grey color
    
#%% Example :
    
#Plot_Var = np.nanmean(Wind_Curl,axis=0)
#cmap_limit=np.nanmax(np.abs( np.nanpercentile(Plot_Var, 99)))
#Plot_range=np.linspace(-cmap_limit,cmap_limit,101) ### Or:  Plot_range=100
#Plot_unit='(N/m3)'; Plot_title= 'Wind Curl (1E-7 N/m3) - (contour lines = Tau_x) - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
#fig, m = func_plotmap_contourf(Plot_Var, Lon_regrid_2D, Lat_regrid_2D, Plot_range, Plot_title, Plot_unit, plt.cm.seismic, 'cyl', 210., 80., -80., 'fill')
#fig.savefig(dir_figs+'name.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
#%%  
    fig=plt.figure()
    
    if P_proj=='npstere':
        m = Basemap( projection='npstere',lon_0=0, boundinglat=30)
        m.drawparallels(np.arange(-90,90,20), labels=[1,1,0,1], linewidth=1, color='k', fontsize=20)
        m.drawmeridians(np.arange(0,360,30), labels=[1,1,0,1], linewidth=1, color='k', fontsize=20)
    elif P_proj=='spstere':
        m = Basemap( projection='spstere',lon_0=180, boundinglat=-30)
        m.drawparallels(np.arange(-90,90,20), labels=[1,1,0,1], linewidth=1, color='k', fontsize=20)
        m.drawmeridians(np.arange(0,360,30), labels=[1,1,0,1], linewidth=1, color='k', fontsize=20)        
    else:
         m = Basemap( projection=P_proj, lon_0=P_lon0, llcrnrlon=P_lon0-180, llcrnrlat=P_latS, urcrnrlon=P_lon0+180, urcrnrlat=P_latN)
         m.drawparallels(np.arange(P_latS, P_latN+0.001, 40.),labels=[True,False,False,False], linewidth=0.01, color='k', fontsize=20) # labels = [left,right,top,bottom] # Latitutes
         m.drawmeridians(np.arange(P_lon0-180,P_lon0+180,60.),labels=[False,False,False,True], linewidth=0.01, color='k', fontsize=20) # labels = [left,right,top,bottom] # Longitudes        
    if P_c_fill=='fill':
        m.fillcontinents(color='0.8')
    m.drawcoastlines(linewidth=1.0, linestyle='solid', color='k', antialiased=1, ax=None, zorder=None)
    im=m.contourf(P_Lon, P_Lat, P_Var,P_range,latlon=True, cmap=P_cmap, extend='both')
    if P_proj=='npstere' or P_proj=='spstere':
        cbar = m.colorbar(im,"right", size="4%", pad="14%")
    else:
        cbar = m.colorbar(im,"right", size="3%", pad="2%")
    cbar.ax.tick_params(labelsize=20) 
    cbar.set_label(P_unit)
    plt.show()
    plt.title(P_title, fontsize=18)
    mng = plt.get_current_fig_manager()
    mng.window.showMaximized() # Maximizes the plot window to save figures in full
    
    #m = Basemap( projection='cyl',lon_0=210., llcrnrlon=30.,llcrnrlat=-80.,urcrnrlon=390.,urcrnrlat=80.)    
    #m.drawparallels(np.arange(-90.,90.001,30.),labels=[True,False,False,False], linewidth=0.01, color='k', fontsize=20) # labels = [left,right,top,bottom] # Latitutes
    #m.drawmeridians(np.arange(0.,360.,60.),labels=[False,False,False,True], linewidth=0.01, color='k', fontsize=18) # labels = [left,right,top,bottom] # Longitudes
    #plt.close()
    
    return fig, m


#### Compute and plot lagged correlation plot####
def func_plot_lagcor_sig(P_Var_x, P_Var_y, lag, P_title, P_color, P_legend, P_legend_loc, P_v_i):
#%% Example :
    
#fig=plt.figure()
#func_plot_lagcor_sig(WS_index_surface_rm[:-9], LAB_index_rm[:-9], 40, P_title, 'y', 'LAB index (Temp at 500m depth)', 'best', '-')
#func_plot_lagcor_sig(WS_index_surface_rm[:-9], Winds_40S60S_0W30W_m[:-9], 40, P_title, 'g', 'Zonal Winds [40S60S, 0W30W]', 'best', '-')
#func_plot_lagcor_sig(WS_index_surface_rm[:-9], AMOC_max_50S_m[:-9], 40, P_title, 'b', 'AMOC max 50S', 'best', '-')
#func_plot_lagcor_sig(WS_index_surface_rm[:-9], AMOC_max_30S_m[:-9], 40, P_title, 'r', 'AMOC max 30S', 'best', '-')
#plt.ylim(-0.7,0.7)
#fig.savefig(dir_figs+'name.png', format='png', dpi=300, transparent=True, bbox_inches='tight')
###  P_v_i='yes' shows the P-value plot of the corrolation
### P_legend_loc = 'best' or 'lower left' or 'right' or 'center' or ...
#%%    
    R_val=[]
    P_val=[]
    from scipy import stats
    for i in range(2*lag):
        slope, intercept, r_value, p_value, std_err = stats.linregress(P_Var_x[lag:len(P_Var_x)-lag], P_Var_y[i:len(P_Var_y)-2*lag+i])
        R_val.append(r_value)
        P_val.append(p_value)
    xx=np.linspace(-lag,lag+1, 2*lag)
    plt.grid(True,which="both",ls="-", color='0.65')
    plt.plot(xx, R_val, P_color, label=P_legend, linewidth=3.0)
    if P_v_i=='yes':
        plt.plot(xx, P_val, P_color, ls='--')#, label='Significance (P-value)')
    plt.xlabel('Years lag', fontsize=18)
    plt.ylabel('Correlation coefficient (R)', fontsize=18)
    plt.title(P_title, fontsize=18)
    plt.xticks(fontsize = 20); plt.yticks(fontsize = 24)
    plt.legend()
    #plt.show()
    plt.legend(prop={'size': 15}, loc=P_legend_loc, fancybox=True, framealpha=0.8)
    l = plt.axhline(y=0, color='k')
    mng = plt.get_current_fig_manager()
    mng.window.showMaximized() # Maximizes the plot window to save figures in full

    #return fig


#### Compute and plot Power Spectral Density plots using welch method####
def plot_PSD_welch_conf(P_Var, P_fq, P_c_probability, P_rho, P_smooth, P_title, P_color, P_legend, P_legend_loc):

### P_c_probability= 0.95 or 0.05 or '-' = confidence interval levels, if the input is not a number then teh confidence intervals won't be plptted
### P_rho = 'yes' or 0.7 or '-' = Red noise significance line, if 'yes' then the Rho and line will be calculated, if P_rho is a number then it will be given as the Rho value, else NO Red noise significance line will be plotted
### P_smooth = 9 (should be even number) or '-' = Number of years for smoothing the PSD line and confidence intervals, if no number is given then NO smoothing will be applied
### P_Var= Plotting variable at X-axis, 1D(X) || P_title=Plot title || P_color=Plot line color || P_legend=Plot variable name to be shown in Plot Legend 
### P_fq=PSD Welch method's sampling frequency of the x time series in units of Hz (Defaults to 1.0) || P_legend_loc=location of legend, 'best' or 'lower left' or 'right' or 'center' or ...
#%% Example :

#P_Var=WS_index ; P_legend='WS Convection Index' ; P_color='r'
#P_title='Power Spectral Density - WS Convection Index with 95% confidence intervals - '+str(year_start)+'-'+str(year_end)+' - '+str(GCM)
#fig=plt.figure()
#ff,PSD=plot_PSD_welch_conf(P_Var, 1, 0.95, 'yes', 9, P_title, P_color, P_legend, 'lower right')
#plt.ylim(1e-5,1e3) #;plt.xlim(1,500)
#fig.savefig(dir_figs+'name.png', format='png', dpi=300, transparent=True, bbox_inches='tight') 
### 

    from scipy.stats import chi2    
    from scipy import signal

    ff,PSD = signal.welch(P_Var,P_fq) # ff= Array of sample frequencies  ,  PSD = Pxx, Power spectral density or power spectrum of x (which is P_Var)
    X_var=np.linspace(1,len(P_Var)/2+1,len(P_Var)/2+1)
    X_var=ff**(-1); X_var=X_var
    
    if P_rho=='yes':
        Rho = np.nansum( (P_Var[:-1] - np.nanmean(P_Var,axis=0)) * (P_Var[1:] - np.nanmean(P_Var,axis=0)) ,axis=0) / np.nansum( (P_Var[:-1] - np.nanmean(P_Var,axis=0)) * (P_Var[:-1] - np.nanmean(P_Var,axis=0)) ,axis=0)
    elif type(P_rho)==float:
        Rho=P_rho # Rho is the memory parameter
    
    if type(P_smooth)==float or type(P_smooth)==int:
        P_smooth=np.int(P_smooth)
        v = 2*P_smooth # P is the number of estimates in welch function and also the degree of freedom.
        
        sm= np.int( (P_smooth-1)/2)
        P_smooth = np.int( (sm*2)+1 ) # P_smoothhs to be even number for centered smoothing; in case odd number was given, it's changed to an even number by subtracting 1
        PSD_m=copy.deepcopy(PSD) ## Smoothing the variable
        for ii in range(sm,PSD.shape[0]-sm+1):
            PSD_m[ii]=np.nanmean(PSD[ii-sm:ii+sm])
        
        PSD_m=PSD_m[sm:-sm]
        X_var_m=X_var[sm:-sm]
        P_legend=P_legend+' ('+str(np.int(P_smooth))+'yr smoothed)'
        
    else:
        v=2

        PSD_m=copy.deepcopy(PSD)
        X_var_m=copy.deepcopy(X_var)
        
    if type(P_c_probability)==float:
        if P_c_probability < 0.5:
            P_c_probability=1-P_c_probability # In case the P_c_probability is input 0.05 instead of 0.95 for example
        alfa = 1 - P_c_probability

    if P_rho=='yes' or type(P_rho)==float or type(P_rho)==int:
        if type(P_c_probability)!=float: # In case the P_c_probability is not given since confidence interval calculation is not necessary, but red noise significance line is needed
            alfa=0.05
            
        F_x_v = (1-Rho**2) / (1 + Rho**2 - 2*Rho*np.cos(2*np.pi*ff ) )  #  F_x_v is the power spectraum   
        F_x_v_star=np.float( np.real( np.nanmean(PSD,axis=0) / np.nanmean(F_x_v,axis=0) ) ) * F_x_v 
        Pr_alpha = (1/v) * F_x_v_star * np.float( chi2.ppf([1 - alfa], v) )
    
    plt.grid(True,which="both",ls="-", color='0.65')
    plt.loglog(X_var_m,PSD_m, color=P_color, label=P_legend)
    plt.legend(loc='best')
    plt.xlabel('Period (years)', fontsize=18)
    plt.ylabel('Spectral Density', fontsize=18) 
    plt.xticks(fontsize = 20); plt.yticks(fontsize = 20)
    #plt.gca().invert_xaxis()
    if type(P_c_probability)==float:
        Chi = chi2.ppf([1 - alfa / 2, alfa / 2], v)
        
        PSDc_lower = PSD_m * ( v / Chi[0] )
        PSDc_upper = PSD_m * ( v / Chi[1] ) 
        plt.loglog(X_var_m,PSDc_lower, color='g', ls='--', label=str(np.int( (1 - alfa) *100))+'% confidence intervals')
        plt.loglog(X_var_m,PSDc_upper, color='g', ls='--')
    if P_rho=='yes' or type(P_rho)==float or type(P_rho)==int:
        plt.loglog(X_var,Pr_alpha , color='b', ls='--', label=str(np.int( (1 - alfa) *100))+'% Red Noise Significance Level')
    plt.legend(prop={'size': 20}, loc=P_legend_loc, fancybox=True, framealpha=0.8)
    plt.title(P_title, fontsize=18) 
    plt.show()
    mng = plt.get_current_fig_manager()
    mng.window.showMaximized() # Maximizes the plot window to save figures in full
    return ff,PSD











