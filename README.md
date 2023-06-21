# KILN
Ksn from Inverse LiNear method or KILN

This code calculates channel steepness index (ksn) using the inverse scheme set out by Fox (2019) and Smith et al. (2022). 

The code requires topotoolbox (Schwanghart and Scherler 2014) to run, and runs in MATLAB.

Plot_Ksn reads in a DEM, extracts a stream network, and allows the user to set a number of parameters, before calling the KsnInversion function. 
KsnInversion calculates the channel steepness index values of a discretised grid over a DEM. The grid_size function determines the size of the
discretised grid. 

Once KsnInversion has completed the channel steepness index calculations, a plot is made using the viridis colour map and the fuction cmapscale (courtesy of Prof. Mark Brandon, Yale), then MATLAB matrices containing the x coordinates, y coordinates and channel steepness values are saved. 


