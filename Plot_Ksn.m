%%% Plotting routine to calculate ksn from the linear inverse method
%%% described by Fox 2019 and Smith et al. 2022
clear all
close all


addpath('PATH TO TOPOTOOLBOX/topotoolbox-master')
addpath('PATH TO TOPOTOOLBOX/topotoolbox-master/utilities')

%read in DEM 
DEM = GRIDobj('LOCATION OF YOUR DEM');
DEM = reproject2utm(DEM, 30); %this function may not be needed
FD = FLOWobj(DEM,'preprocess','carve');
DEM = imposemin(FD,DEM,0.0001);
A = flowacc(FD);

%extract stream network
S = STREAMobj(FD, 'minarea', 1e6, 'unit','map'); %extract stream network that you desire

grid_size = 500; %controls size of discretised grid
alpha = 1000; %controls the weighting of smoothness
mn = 0.45; %m/n ratio or concavity ratio

%call the function to calculate ksn
inversion = KsnInversion(S, DEM, A, grid_size, alpha, mn);


%plotting routine to plot ksn values with topotoolbox
%x, y and ksn values are also saved as matlab matrixes
figure
pcolor(inversion.xgrid, inversion.ygrid, inversion.ksn);
shading interp
cMap = viridis(800);
cMap = cmapscale((inversion.ksn(inversion.ksn>0)), cMap, 0.5);
cMap = [[0.7 0.7 0.7]; cMap];
colormap( cMap);
hold on
plot(S, 'w')
box on
h = gca;
xticklabels({})
h.TickDir = 'in';
h.Layer = 'top';
h.FontSize = 16;
h.LineWidth = 1;
%... Create colorbar
hCB = colorbar;
hCB.Label.String = 'k_{sn} [m]';
hCB.Label.FontSize = 16;
hCB.Label.LineWidth = 1;
%... Write labels
%hT = title('U*', 'FontSize', 16);
%hT.Units = 'normalized';
%hT.Position(2) = hT.Position(2) + 0.02;
%xlabel('Longitude (deg)', 'FontSize', 16)
ylabel('Northing [m]', 'FontSize', 16)
xlabel('Easting [m]', 'FontSize', 16)


