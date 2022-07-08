%%10Be Production Rate Calculator Code

%Author: Chia-Yu Chen
%Modified by: Richard Ott and Erica Erlanger
% Remodified by Yanyan Wang to make it a function on March 5, 2018
% Input:
%       DEM: GRIDobj
%       ZONE: the UTM zone number, for southern hemisphere points, use negative values for zone -ZONE.
% Output:
%        GRIDobj of production rates
% Last modified by Yanyan Wang on Aug. 8, 2019, comment shielding calculation, see DiBiase,
% 2018 paper about shiedling and attenuation length change

function production = CNP(DEM,ZONE)

%Convert GRIDobj to matrix and coordinate vectors
[~,X,Y] = GRIDobj2mat(DEM);

% Find max and min X and Y for DEM, and convert them to lat lon
[lat_min,lon_min]=utm2ll(min(X),min(Y),ZONE,'wgs84');
[lat_max,lon_max]=utm2ll(max(X),max(Y),ZONE,'wgs84');

lat_in=(lat_min:0.1:lat_max+0.1)';
lon_in=(lon_min:0.1:lon_max+0.1)';

% DEM pixel size
pixelsize = DEM.cellsize; % pixcel size of DEM

%% Calculate topographic shielding 
% Tshd= toposhielding(DEM,10,5);

%% Import Scaling Equation Constants and interpolate for each pixel in DEM

% Scaling Equation Constants by Latitude (Stone,2000)
Lat = [0,10,20,30,40,50,60,90];
a = [31.8518,34.3699,40.3153,42.0983,56.7733,69.0720,71.8733,71.8733];
b = [250.3193,258.4759,308.9894,512.6857,649.1343,832.4566,863.1927,863.1927];
c = [-0.083393,-0.089807,-0.106248,-0.120551,-0.160859,-0.199252,-0.207069,-0.207069];
d = [7.4260e-5,7.9457e-5,9.4508e-5, 1.1752e-4,1.5463e-4,1.9391e-4,2.0127e-4,2.0127e-4];
e = [-2.2397e-8,-2.3697e-8,-2.8234e-8,-3.8809e-8,-5.0330e-8,-6.3653e-8,-6.6043e-8,-6.6043e-8];
m = [0.587,0.600,0.678,0.833,0.933,1.000,1.000,1.000];

% Convert latitude and longitude to UTM
% [x,~]=ll2utm(ones(length(lon_in),1)*lat_in(1),lon_in); 
[~,y]=ll2utm(lat_in,ones(length(lat_in),1)*lon_in(1)); 

% Create x and y array for interpolation of production rate within drainage
%basin
x1 = min(X):pixelsize:max(X); 
y1 =(min(Y):pixelsize:max(Y))'; 
% to deal with south hemisphere negative points in order to interpolate
if lat_in(1)<0
   lat_in = abs(lat_in); 
end
% Interpolation of scaling fator - to find values at [x1(:),y1(:)]
A = GRIDobj(x1,y1,repmat(interp1(y,interp1(Lat,a,lat_in),y1),1,length(x1)));
B = GRIDobj(x1,y1,repmat(interp1(y,interp1(Lat,b,lat_in),y1),1,length(x1)));
C = GRIDobj(x1,y1,repmat(interp1(y,interp1(Lat,c,lat_in),y1),1,length(x1)));
D = GRIDobj(x1,y1,repmat(interp1(y,interp1(Lat,d,lat_in),y1),1,length(x1)));
E = GRIDobj(x1,y1,repmat(interp1(y,interp1(Lat,e,lat_in),y1),1,length(x1)));
M = GRIDobj(x1,y1,repmat(interp1(y,interp1(Lat,m,lat_in),y1),1,length(x1)));

%% Calculate Production Rates for Spallation and Muons

% Spallation-induced Production Rate of 10Be (at/g/yr) assuming 'Lm" scaling framework (Borchers et al., 2015)
Pn_SLHL = 4.00; 

% Slow Muon-induced Production Rate of 10Be (at/g/yr) scaled to sea level (Braucher et al. 2011)
Pms_SLHL = 0.012; 

% Fast Muon-induced Production Rate of 10Be (at/g/yr) scaled to sea level (Braucher et al. 2011)
Pmf_SLHL = 0.039; 

% Perform air pressure correction for elevation (Eq. 1 from Stone, 2000)to account for diminished cosmic ray attenuation with decreased atmospheric pressue 
pres=1013.25*exp(((-0.03417)/6.5e-3)*(log(288.15)-log(288.15-(6.5e-3*DEM.Z)))); 

%% Calculation 10Be production from spallation and muons

% 10Be production from Spallation, assuming neutron attenuation length in air of 150 g/cm2
% Pn= Tshd.Z*Pn_SLHL.*(A.Z+B.Z.*exp(-pres/150)+C.Z.*pres+D.Z.*pres.^2+E.Z.*pres.^3); % Stone, 2000
Pn= Pn_SLHL.*(A.Z+B.Z.*exp(-pres/150)+C.Z.*pres+D.Z.*pres.^2+E.Z.*pres.^3); % Stone, 2000

% 10Be production from Slow Muons, assuming (1) sea level pressure of 1013.25 mbar and 
%(2) muon attentuation length in air of 260 g/cm2 (Braucher et al., 2011)
% Pms=Tshd.Z*Pms_SLHL.*exp((1013.25-pres)/260); 
Pms=Pms_SLHL.*exp((1013.25-pres)/260); 

% 10Be production from Fast Muons, assuming (1) sea level pressure of 1013.25 mbar and 
%(2) muon attentuation length in air of 510 g/cm2 (Braucher et al., 2011)
% Pmf=Tshd.Z*Pmf_SLHL.*exp((1013.25-pres)/510); 
Pmf=Pmf_SLHL.*exp((1013.25-pres)/510); 

%Convert production to GRIDobjs
Pn_obj=GRIDobj(X,Y,Pn);
Pms_obj=GRIDobj(X,Y,Pms);
Pmf_obj=GRIDobj(X,Y,Pmf);

production.neutron = Pn_obj;
production.smuon = Pms_obj;
production.fmuon = Pmf_obj;



% % Convert output to text file
% GRIDobj2ascii(Pn_obj,['Pn_' appd_name '.txt'])
% GRIDobj2ascii(Pms_obj,['Pms_' appd_name '.txt'])
% GRIDobj2ascii(Pmf_obj,['Pmf_' appd_name '.txt'])

% PnAvg = mean(Pn(:),'omitnan');
% PmsAvg = mean(Pms(:),'omitnan');
% PmfAvg = mean(Pmf(:),'omitnan');
% thres_ele = 120;
% tag = DEM.Z>thres_ele;
% Pn_effec = Pn.*tag;
% Pms_effec = Pms.*tag;
% Pmf_effec = Pmf.*tag;
% numeffec = sum(tag(:));
% PnAvg = sum(Pn_effec(:),'omitnan')/numeffec;
% PmsAvg = sum(Pms_effec(:),'omitnan')/numeffec;
% PmfAvg = sum(Pmf_effec(:),'omitnan')/numeffec;

% 
% imageschs(Pmf_obj)
% figure
% imageschs(DEM)
% 
% figure
% imageschs(Tshd)



