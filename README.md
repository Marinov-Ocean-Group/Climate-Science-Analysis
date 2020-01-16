## Climate-Science-Analysis

---------------------------------------------------------------------------------------------------------
# Climate Change impact on global Sea Surface Temperatures

Main code #1: CMIP5_SST.py

* The code reads 4-Dimensional sea surface temperature (SST) data from 14 CMIP5 climate models for two periods of historical 1980-1999 and future 2080-2099 (under RCP8.5 scenario). It regrids them into 1 degree by 1 degree fields, and plots all models' SST averages for each period. It also plots the difference between 2080-2099 and 1980-1999 which shows the impact of climate change on global SSTs under RCP8.5 scenario. It also plots the global average SST of each model versus the changes in SSTs between two periods as a scatter plot. The averages are grid-cell-area-weighted meaning the lager tropical gridcells have higher impact on global averages than the smaller polar cells

* The climate model data are stored at UPenn's local server

Functions code: Behzadlib.py

* This code contains various analysis/plotting functions that are imported in the main code as needed

Code #2: CMIP5_SST_2.py
* Does same job as the CMIP5_SST.py code, but reads .nc files using a different method by listing all file names and reading the files one by one (instead of loading them all using MFDATSET) - This is useful when the data are huge and the computer memory is low

---------------------------------------------------------------------------------------------------------

Final plotting products:

* Fig_CMIP5_SST_hist_1980_1999.png   = Global SST maps of 14 CMIP5 models for the historical 1980-1999 period
* Fig_CMIP5_SST_rcp8p5_2080_2099.png = Global SST maps of 14 CMIP5 models for the 2080-2099 under RCP8.5 scenario
* Fig_CMIP5_SST_climate_change_Impact_RCP8p5.png = Climate Change impact - 2080-2099 average minus 1980-1999 average
* Fig_CMIP5_SST_climate_change_Impact_RCP8p5_scatter.png = Scatter plot of global average SST versus the changes in SSTs between two periods, for each model


![Alt text](https://raw.githubusercontent.com/Marinov-Ocean-Group/Climate-Science-Analysis/master/Fig_CMIP5_SST_climate_change_Impact_RCP8p5.png)


---------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------
# Other Climate Science analysis and plots

Main code #2: CMIP5_Climate.py

Various climate science analysis and plots:
* Calculating and plotting Empirical Orthogonal Functions (EOFs) of sea surface temperature over Pacific Ocean, where the EOF indices are Butter-worth filtered to smooth out high frequency noises
* Calculating North Atlantic Oscillation (NAO) defined as the 1st EOF of Sea-Level Air Presure over North Atlantic
* Calculating Curl of the Wind using wind stress in X and Y directions
* Plotting Arctic Sea Ice Concentration average over each month of the year

* The climate model data are stored at UPenn's local server

Functions code: Behzadlib.py

* This code contains various analysis/plotting functions that are imported in the main code as needed

---------------------------------------------------------------------------------------------------------

Final plotting products:

* Fig_EOF_SST_SpatialPattern_GFDL-ESM2G.png = Spatial Pattern of EOFs of sea surface temperature over Pacific Ocean

* Fig_EOF_SST_Indices_GFDL-ESM2G.png = Indices of EOFs of sea surface temperature over Pacific Ocean
![Alt text](https://raw.githubusercontent.com/Marinov-Ocean-Group/Climate-Science-Analysis/master/Fig_EOF_SST_Indices_GFDL-ESM2G.png)

* Fig_NAO_SpatialPattern_GFDL-ESM2G.png = Spatial Pattern of North Atlantic Oscillation (NAO)

* Fig_NAO_Indices_GFDL-ESM2G.png = Indices of North Atlantic Oscillation (NAO)
![Alt text](https://raw.githubusercontent.com/Marinov-Ocean-Group/Climate-Science-Analysis/master/Fig_NAO_SpatialPattern_GFDL-ESM2G.png)

* Fig_Wind_Curl_GFDL-ESM2G.png = Curl of the wind, calculated as

Wind_Curl = ( D_Tau_Y / D_X ) - ( D_Tau_X / D_Y ) # D_X = (Lon_1 - Lon_2) * COS(Lat)

* Fig_Wind_Curl_f_WQuiver_GFDL-ESM2G.png = Ekman transport, equal to wind curl divided by coriolis parameter. The quivers are the wind direction - Wind_Crul / f , f = coriolis parameter = 2Wsin(LAT) , W = 7.292E-5 rad/s
![Alt text](https://raw.githubusercontent.com/Marinov-Ocean-Group/Climate-Science-Analysis/master/Fig_Wind_Curl_f_WQuiver_GFDL-ESM2G.png)

* Fig_SeaIce_Arctic_monthly_GFDL-ESM2G.png = Arctic Sea Ice concentration average for each month - average of 1991-2000
![Alt text](https://raw.githubusercontent.com/Marinov-Ocean-Group/Climate-Science-Analysis/master/Fig_SeaIce_Arctic_monthly_GFDL-ESM2G.png)



