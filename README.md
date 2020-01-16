# Climate-Science-Analysis

Various climate science analysis and plots:
* Calculating and plotting Empirical Orthogonal Functions (EOFs) of sea surface temperature over Pacific Ocean, where the EOF indices are Butter-worth filtered to smooth out high frequency noises
* Calculating North Atlantic Oscillation (NAO) defined as the 1st EOF of Sea-Level Air Presure over North Atlantic
* Calculating Curl of the Wind using wind stress in X and Y directions
* Plotting Arctic Sea Ice Concentration average over each month of the year

---------------------------------------------------------------------------------------------------------

Main code: CMIP5_Climate.py

* The climate model data are stored at UPenn's local server

Functions code: Behzadlib.py

* This code contains various analysis/plotting functions that are imported in the main code as needed

---------------------------------------------------------------------------------------------------------

Final plotting products:

* Fig_EOF_SST_SpatialPattern_GFDL-ESM2G.png = Spatial Pattern of EOFs of sea surface temperature over Pacific Ocean

* Fig_EOF_SST_Indices_GFDL-ESM2G.png = Indices of EOFs of sea surface temperature over Pacific Ocean
![Alt text](https://raw.githubusercontent.com/behzadasd/Climate-Science-Analysis/master/Fig_EOF_SST_Indices_GFDL-ESM2G.png)

* Fig_NAO_SpatialPattern_GFDL-ESM2G.png = Spatial Pattern of North Atlantic Oscillation (NAO)

* Fig_NAO_Indices_GFDL-ESM2G.png = Indices of North Atlantic Oscillation (NAO)
![Alt text](https://raw.githubusercontent.com/behzadasd/Climate-Science-Analysis/master/Fig_NAO_SpatialPattern_GFDL-ESM2G.png)

* Fig_Wind_Curl_GFDL-ESM2G.png = Curl of the wind, calculated as

Wind_Curl = ( D_Tau_Y / D_X ) - ( D_Tau_X / D_Y ) # D_X = (Lon_1 - Lon_2) * COS(Lat)

* Fig_Wind_Curl_f_WQuiver_GFDL-ESM2G.png = Ekman transport, equal to wind curl divided by coriolis parameter. The quivers are the wind direction - Wind_Crul / f , f = coriolis parameter = 2Wsin(LAT) , W = 7.292E-5 rad/s
![Alt text](https://raw.githubusercontent.com/behzadasd/Climate-Science-Analysis/master/Fig_Wind_Curl_f_WQuiver_GFDL-ESM2G.png)

* Fig_SeaIce_Arctic_monthly_GFDL-ESM2G.png = Arctic Sea Ice concentration average for each month - average of 1991-2000
![Alt text](https://raw.githubusercontent.com/behzadasd/Climate-Science-Analysis/master/Fig_SeaIce_Arctic_monthly_GFDL-ESM2G.png)



