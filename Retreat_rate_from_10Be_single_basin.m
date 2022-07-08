% This script use cosmogenic nuclide concentration to calculate 3 kinds of
%     rate with unit of m/Ma: 
%       (1) the conventional "erosion rate" (Rate_convention)
%       (2) the retreat rate from the dot product method (Rate_dot)
%       (3) the retreat rate from the shadow method (Rate_shadow)
% 
% Geographic direction is defined here as : North-0 degree, East-90 degree,
%                                          South-180 degree, West-270 degree
%
% 
% Inputs:
%      1) measured 10Be concentration
%      2) error of 10Be concentration 
%      3) .tif file of the basin
%      4) regional retreat direction
%      5) UTM zone of the basin
%      6) geographic heimisphere of the basin
%
% Outputs: A .mat file that stores a matrix of 10 columns. 
%          1-retreat direction; 2 to 4-retreat rate from Local scalar product 
%          method and the 34% uncertainty; 5 to 7-retreat rate from Basin 
%          projection method and the 34% uncertainty; 8 to 10-conventional erosion 
%          rate and the 34% uncertainty.
%          File is named with the keyword in DEM tif file.
%          File is saved to a folder named <Retreat_rate_from_10Be>
%

% 
% Author: yanyan.wang@erdw.ethz.ch
%

%% clean Workspace
clear; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 0_0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% addpath 
% put path of where you store the TopoToolBox in the bracket
addpath('/Users/yanywang/Dropbox/MATLAB/topotoolbox-master')

% put path of where you store the mass_flux_calculator package
addpath('/Users/yanywang/Downloads/mass_flux_from_cosmo-master')

% put path of where you store your DEM tif files
addpath('/Users/yanywang/Document/DEMs')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Give the inputs, follow the units when giving numbers to input parameters. 

% Input 1) measured 10Be concentration, [atoms/g] 
conc = 123345; 

% Input 2) error of 10Be concentration, [atoms/g]
conc_std = 6789;

% Input 3) name of the .tif file of your basin DEM
demtxt = 'ghats_0_4.TIF' ; %  DEM file name

% Input 4) this is the estimated regional divide/escarpment retreat
% direction azimuth, [degree]. Minimum is 0, maximum is 360. 
% Multiple numbers are supported, 
% e.g. regional_dir = [57 100 130 200];
% e.g. regional_dir = linspace(0 50, 10); 
regional_dir = [57 60 100]; 

% Input 5) UTM zone numbe, For southern hemisphere, use negative values.
% e.g. Madagascar (-38), western India (43). Please refer to the World UTM
% grid zone map for your study area. 
ZONE = 43; 

% Input 6) UTM heimisphere, northern hemisphere, use 'N'; southern
% hemisphere, use 'S'.
utmhemi = 'N';



%%%%%%%%%%%%%%% process DEM before calculations
%%

% deal with DEM 
DEM = GRIDobj(demtxt);
DEM = fillsinks(DEM);
% grid size of DEM
cs = DEM.cellsize;
% index of non-nan data point
IDnan =find(~isnan(DEM.Z));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate surface production rate and basin-wide nuclides flux
production = CNP(DEM, ZONE); % this step will take some time, please be patient

% Calculate total Be10 nuclides atoms, don't change these values unless you
% have other reliable reference sources. 
rho_rock = 2.65; % density of rock, [g/cm^3]
attN = 160; %Spallation attenuation length (Braucher et al., 2011) [g/cm2]
attMs = 1500; % Slow muons attenuation length (Braucher et al., 2011) [g/cm2]
attMf = 4320;% Fast muons attenuation length (Braucher et al., 2011) [g/cm2]

muN=attN/rho_rock;
muMs=attMs/rho_rock;
muMf=attMf/rho_rock;

% basin-wide nuclide flux
M_Be = muN*sum(production.neutron.Z(:),'omitnan')*cs^2+muMs*sum(production.smuon.Z(:),'omitnan')*cs^2+ ...
       muMf*sum(production.fmuon.Z(:),'omitnan')*cs^2; 

   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate conventional erosion rate (Rate_convention)  
Rate_convention = (1/conc*M_Be./(length(IDnan)*cs^2)*1e+4); % conventional "erosion rate", [m/Ma]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  use Local-Scalar-Product method to calculate retreat rate (Rate_dot)
% to define the inclination flux vector with the horizontal surface,
% 90-vertical, 0-purely horizontal. No need to change this if you only want
% to calculate horizontal retreat rate and conventional erosion rate.
phi = 0; % the declination of the flux vector with the horizontal surface, [degree]

Vasp = regional_dir+180;
% discretize basin surface into elemental surfaces and find their norms and areas
tri3D = DEM_3D_triangle_surfaces(DEM);

fluxNormx = cosd(phi)*sind(Vasp);
fluxNormy = cosd(phi)*cosd(Vasp);

% flux aspect unit vector
fluxNorm = [fluxNormx.' , fluxNormy.', ones(size(fluxNormx.'))*sind(phi)];

% volume of rocks 
Vrock = zeros(length(fluxNormx),1);
for kk = 1:length(fluxNormx)
    
    VV = repmat(fluxNorm(kk,:), length(tri3D.unitnorm),1);
    Vrock(kk) = sum(dot(VV,tri3D.unitnorm,2).*tri3D.area);  % sum of all positive and neagative volume
    
end

%Retreat rate from Local scaler product method, [m/Myr]
Rate_dot = 1/conc*M_Be./Vrock*1e+4; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% use Bains-projection method to calculate retreat rate (Rate_shadow)
sita = Vasp+180+90; % angle of coordinate clockwise rotation 
% to avoid singling angles
check = mod(sita,45)==0;
if sum(check)
   sita =  sita + 0.05;    
end

% xy is the geographic x-y coordinates data points
[demx,demy]=getcoordinates(DEM);
xydem1 = repmat(demx,length(demy),1);
xydem2 = repmat(demy,1,length(demx));
Parea = zeros(size(sita));
discs = cs; % discretize cellsize
  for k = 1:length(sita)
        xnew=xydem1; % x coordinates for the whole rectagular domain that enclose the basin
        ynew=xydem2; % y coordinates for the whole rectagular domain that enclose the basin
    % to rotate x-o-y surface with angle sita? in counterclockwise direction and use the centroid as rotate original point
        xco = mean(xnew(IDnan)); % x-coordinates of centroid
        yco = mean(ynew(IDnan)); % y-coordinates of centroid
        xabs = xnew-xco; % New coordinates in respect of the centroid
        yabs= ynew-yco;   
        xnewp = xabs*cosd(sita(k))+yabs*sind(sita(k)); % to rotate and get the new coordinates
        ynewp = yabs*cosd(sita(k))-xabs*sind(sita(k));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%The projection part%%%%%%%%%%%%%%%%%%%%%%%%%%        
        % to reshape the elevation matrix
        ynan = ynewp(IDnan); % y-coordinates of non-nan data point
        xnan = xnewp(IDnan); % x-coordinates of non-nan data point     
        DEMsubZnan = DEM.Z(IDnan); % elevation of non-nan data point
        [B,I] = sort(xnan, 'ascend'); % B is sorted x-coordinate
        DEMsubZproj = DEMsubZnan(I); % now z is also in the same order as x-coordinates       
    % to sort data into uniform-spaced intervals for area integration
        ngrid = ceil((B(2:end)-B(1))/discs); 
        divisionID = find(diff(ngrid)==1);
        dzmax = zeros(length(divisionID),1);
        dzmin = zeros(length(divisionID),1); 
        for j = 1:length(divisionID)
            if j == 1
            endID = divisionID(j);
            beginID = 1;
            dzmax(j) = max(DEMsubZproj(beginID:endID));
            dzmin(j) = min(DEMsubZproj(beginID:endID));
            else 
            endID = divisionID(j);
            beginID = divisionID(j-1);
            dzmax(j) = max(DEMsubZproj(beginID:endID));
            dzmin(j) = min(DEMsubZproj(beginID:endID));
            end        
        end
    % this is the area that are projected on the x-o-z surface
         Parea(k) = sum(dzmax-dzmin)*discs/1e+6; % km^2 

  end

% Retreat rate from the Basin projection method, [m/Myr]
Rate_shadow = Rate_convention*length(IDnan)*cs^2/1e+6./Parea; 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate +-34% uncertainty for all three rates(conventional erosion rate, two horizontal retreat rates).
% uncertainty is propogated from error of measured 10Be concentration.

% Randomly selected number
n_rand=1e4; % number of simulations
C_rand=truncnormrnd([n_rand,1],conc,conc_std,0,inf);
Pn = mean(production.neutron.Z(:),'omitnan');
Pms = mean(production.smuon.Z(:),'omitnan');
Pmf = mean(production.fmuon.Z(:),'omitnan');
IDnan =find(~isnan(production.neutron.Z));% index of non-nan data point
PnStd=sqrt((Pn*0.025).^2+(Pn*0.09).^2);
Pn_rand=truncnormrnd([n_rand,1],Pn,PnStd,0,inf); % random spallation prod rate
    
% Propagation of Uncertainties for production by muons (50%), added to scaling scheme (9%)
PmsStd=sqrt((Pms*0.5).^2+(Pms*0.09).^2);
Pms_rand=truncnormrnd([n_rand,1],Pms,PmsStd,0,inf); % random slow muon prod rate

PmfStd=sqrt((Pmf*0.5).^2+(Pmf*0.09).^2);
Pmf_rand=truncnormrnd([n_rand,1],Pmf,PmfStd,0,inf); % random fast muon prod rate

% Propagation of Uncertainties to total cosmo nuclides atoms 
M_Be_rand = (muN*Pn_rand+muMs*Pms_rand+muMf*Pmf_rand)*cs^2*length(IDnan);

% (1) calculate uncertainty for Rate_dot
dot_rate_avg = zeros(length(Rate_dot),length(C_rand));  
dot_rate_1Std_low = zeros(size(Rate_dot));
dot_rate_1Std_up = zeros(size(Rate_dot));
for i = 1:length(Rate_dot)
    
    for j=1:n_rand
        if C_rand(j)~=0
            dot_rate_avg(i,j)=(1/C_rand(j))*M_Be_rand(j)/Vrock(i)*1e+4;
            
        end
    end
    
    % generate production rate function and find upper and lower uncertainty
    rate_means=dot_rate_avg(i,:);  % simulated rates for a sample
    rate_d = Rate_dot(i);       % rate calculated from the measured concentration and given prod rates

    max_r=min(max(rate_means),1e+4); % max simulated rate
    min_r=min(rate_means);         % min simulated rate)

    x = linspace(min_r,max_r, 3000); % generate 3000 points in between min_r and max_r
    pd=fitdist(rate_means.','kernel');

    y=pdf(pd,x);

    %calculate cumulative distribution function for specified distribution
    mcdf=cdf(pd,x);
    [~,ind]=min(abs(x-rate_d)); 

    % find the one standard deviation
    mid_cdf=mcdf(1,ind); 
    low_cdf=mid_cdf-0.34; % 1 sigma 
    up_cdf=mid_cdf+0.34;  % 1 sigma

    [~,ind]=min(abs(mcdf-low_cdf));
    low_std=rate_d-x(1,ind);

    [~,ind]=min(abs(mcdf-up_cdf));
    up_std=x(1,ind)-rate_d;
    dot_rate_1Std_low(i) = low_std;
    dot_rate_1Std_up(i) = up_std;
end


% (2) calculate uncertainty for conventional erosion rate: fluxVertical
convention_rate_avg = zeros(length(Rate_convention),length(C_rand));  
convention_rate_1Std_low = zeros(size(Rate_convention));
convention_rate_1Std_up = zeros(size(Rate_convention));

for i = 1:length(Rate_convention)
    
    for j=1:n_rand
        if C_rand(j)~=0
            convention_rate_avg(i,j)=(1/C_rand(j))*M_Be_rand(j)./(length(IDnan)*cs^2)*1e+4;
          
        end
    end
    
    % generate production rate function and find upper and lower uncertainty
    VTCrate_means=convention_rate_avg(i,:);  % simulated rates for a sample
    VTCrate_d = Rate_convention(i);       % rate calculated from the measured concentration and given prod rates

    VTCmax_r=min(max(VTCrate_means),1e+4); % max simulated rate
    VTCmin_r=min(VTCrate_means);         % min simulated rate)

    x = linspace(VTCmin_r,VTCmax_r, 3000); % generate 3000 points in between min_r and max_r
    pd=fitdist(VTCrate_means.','kernel');

    y=pdf(pd,x);

    %calculate cumulative distribution function for specified distribution
    mcdf=cdf(pd,x);
    [~,ind]=min(abs(x-VTCrate_d)); 

    % find the one standard deviation
    mid_cdf=mcdf(1,ind); 
    VTClow_cdf=mid_cdf-0.34; % 1 sigma 
    VTCup_cdf=mid_cdf+0.34;  % 1 sigma

    [~,ind]=min(abs(mcdf-VTClow_cdf));
    VTClow_std=VTCrate_d-x(1,ind);

    [~,ind]=min(abs(mcdf-VTCup_cdf));
    VTCup_std=x(1,ind)-VTCrate_d;
    convention_rate_1Std_low(i) = VTClow_std;
    convention_rate_1Std_up(i) = VTCup_std;
end


% (3) calculate uncertainty for rates from Rate_shadow
shadow_rate_avg = zeros(length(Rate_shadow),length(C_rand));  
shadow_rate_1Std_low = zeros(size(Rate_shadow));
shadow_rate_1Std_up = zeros(size(Rate_shadow));
for i = 1:length(Rate_shadow)
    
    for j=1:n_rand
        if C_rand(j)~=0
            shadow_rate_avg(i,j)=(1/C_rand(j))*M_Be_rand(j)/Parea(i)/1e+6*1e+4;
            
        end
    end
    
    % generate production rate function and find upper and lower uncertainty
    rate_means=shadow_rate_avg(i,:);  % simulated rates for a sample
    rate_d = Rate_shadow(i);       % rate calculated from the measured concentration and given prod rates

    max_r=min(max(rate_means),1e+4); % max simulated rate
    min_r=min(rate_means);         % min simulated rate)

    x = linspace(min_r,max_r, 3000); % generate 3000 points in between min_r and max_r
    pd=fitdist(rate_means.','kernel');

    y=pdf(pd,x);

    %calculate cumulative distribution function for specified distribution
    mcdf=cdf(pd,x);
    [~,ind]=min(abs(x-rate_d)); 

    % find the one standard deviation
    mid_cdf=mcdf(1,ind); 
    low_cdf=mid_cdf-0.34; % 1 sigma 
    up_cdf=mid_cdf+0.34;  % 1 sigma

    [~,ind]=min(abs(mcdf-low_cdf));
    low_std=rate_d-x(1,ind);

    [~,ind]=min(abs(mcdf-up_cdf));
    up_std=x(1,ind)-rate_d;
    shadow_rate_1Std_low(i) = low_std;
    shadow_rate_1Std_up(i) = up_std;
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save data to .mat file 

filename = strcat(demtxt(1:end-4),'_retreat_rate.mat');

Rate_erosion = [Rate_convention convention_rate_1Std_low convention_rate_1Std_up];
Rate_erosion = repmat(Rate_erosion,length(Rate_dot),1);

% DATA includes 10 columns. 1-retreat direction; 2 to 4-retreat rate from
% Local scalar product method and the 34% uncertainty; 5 to 7-retreat rate
% from Basin projection method and the 34% uncertainty; 8 to
% 10-conventional erosion rate and the 34% uncertainty.
% DATA is saved into a .mat file named with the keyword in DEM tif file.
DATA = horzcat(regional_dir', Rate_dot,dot_rate_1Std_low, dot_rate_1Std_up,...
               Rate_shadow', shadow_rate_1Std_low', shadow_rate_1Std_up',...
               Rate_erosion);
           
% create a folder to store DATA
folder = 'Retreat_rate_from_10Be/';
if ~exist(folder, 'dir')
  mkdir(cd, folder);
end
save([folder,filename],'DATA')

