% this script partition the dilution of nuclides for geomorphic zones and
% then calculate horizontal retreat rates of the escarpment

% this is an updated version of undifferentiated reteat rate calculation
% from Be10 concentration. The undifferentiated retreat rate calculation
% don't consider the contribution of high land and lowland sediment to the
% nuclides

% % Date July 26, 2020 by Y.Wang
% Last modified on Aug. 11, 2021 to allow users manually input necessary
% parameters in PART 0

close all
clear

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
regional_dir = 57;

% Input 5) UTM zone numbe, For southern hemisphere, use negative values.
% e.g. Madagascar (-38), western India (43). Please refer to the World UTM
% grid zone map for your study area. 
ZONE = 43; 

% Input 6) UTM heimisphere, northern hemisphere, use 'N'; southern
% hemisphere, use 'S'.
utmhemi = 'N';

% Input 7) the average erosion rate of the plateau
edot = 10e-6; % erosion rate at plateau, m/yr


%% load DEM

DEM = GRIDobj(demtxt);
DEM = fillsinks(DEM);
FD = FLOWobj(DEM,'preprocess','carve');
crita = .1e6; % critical area of river extraction, m2
S = STREAMobj(FD,'minarea',crita/DEM.cellsize^2);
cs = DEM.cellsize;
A = flowacc(FD).*(cs^2)/1e+6; % km^2

% to define the inclination flux vector with the horizontal surface,
% 90-vertical, 0-purely horizontal
phi = 0; % the declination of the flux vector with the horizontal surface, [degree]

%% calculate production, and total nuclides atoms, and conventional erosion rate
% UTM zone of DEM, for south hemisphere, use negative values
production = CNP(DEM, ZONE); % this step will take some time, please be patient
% Calculate total Be10 nuclides atoms
rho_rock = 2.65; % density of rock, [g/cm^3]
attN = 160; %Spallation attenuation length (Braucher et al., 2011) [g/cm2]
attMs = 1500; % Slow muons attenuation length (Braucher et al., 2011) [g/cm2]
attMf = 4320;% Fast muons attenuation length (Braucher et al., 2011) [g/cm2]

muN=attN/rho_rock;
muMs=attMs/rho_rock;
muMf=attMf/rho_rock;

M_Be = muN*sum(production.neutron.Z(:),'omitnan')*cs^2+muMs*sum(production.smuon.Z(:),'omitnan')*cs^2+ ...
       muMf*sum(production.fmuon.Z(:),'omitnan')*cs^2; % total atoms production per year, 

%% manually pick escarpment edge and use them as outlets of highland rivers
figure(1)
imageschs(DEM);
hold on
plot(S, 'w-', 'LineWidth',1);

fprintf('\n #################### INSTRUCTIONS ####################\n')
fprintf('\nclick on the fiugre to choose outlets of plteau-flowing rivers, press RETURN when you finish\n');
fprintf('\noutlets should define the escarpment edge \n');
[xoutlets, youtlets] = ginputc('ShowPoints', true, 'ConnectPoints', true); 

[~,~,IXgrid] = snap2stream(S,xoutlets,youtlets);
DB = drainagebasins(FD,IXgrid);
db = logical(DB.Z);

%% divide basin into three geomorphic zones: highland, escarpment, lowland
DEMp = DEM;
DEMp.Z(~db)=nan; % high land DEM

DEMe = DEM;
DEMe.Z(db)=nan; % escarpment and lowland DEM

% define the direction of horizontal retreat rates

Vasp = regional_dir+180;
%% calculate the projected area of the escarpment for the retreat direction
DEM =DEMe;
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
IDnan =find(~isnan(DEM.Z));% index of non-nan data point
discs = DEM.cellsize; % discretize cellsize
  for k = 1:length(sita)
        xnew=xydem1; % x coordinates for the whole rectagular domain that enclose the basin
        ynew=xydem2; % y coordinates for the whole rectagular domain that enclose the basin
    % to rotate x-o-y surface with angle sita? in counterclockwise direction and use the centroid as rotate original point
        xco = mean(xnew(IDnan)); % x-coordinates of centroid
        yco = mean(ynew(IDnan)); % y-coordinates of centroid
        xabs = xnew-xco; % New coordinates in respect of the centroid
        yabs= ynew -yco;   
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
         Parea(k) = sum(dzmax-dzmin)*discs; % m^2 
  end

  

Ap = sum(db,'all')*cs^2;   
Vtotal = M_Be/conc*1e-2; % m3/yr 
Rate_shadow = (Vtotal-edot*Ap)/Parea*1e6; % m/Myr

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ***************calculate 1 sigma of uncertainty*********************

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

% (3) calculate uncertainty for rates from Rate_shadow
shadow_rate_avg = zeros(length(Rate_shadow),length(C_rand));  
shadow_Rate1Std_low = zeros(size(Rate_shadow));
shadow_Rate1Std_up = zeros(size(Rate_shadow));
for i = 1:length(Rate_shadow)
    
    for j=1:n_rand
        if C_rand(j)~=0
            shadow_rate_avg(i,j)=((1/C_rand(j))*M_Be_rand(j)*1e-2-edot*Ap)/Parea(i)*1e+6;
            
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
    shadow_Rate1Std_low(i) = low_std;
    shadow_Rate1Std_up(i) = up_std;
end


%%
% column 1-3, retreat rate from shadow method at given regional retreat
% dirction, column 4 is area of user-defined plateau in the escarpment
% basin
temp = [Rate_shadow,shadow_Rate1Std_low, shadow_Rate1Std_up,Ap/1e6];













% DEMe = crop(DEMe);
% figure
% imageschs(DEMe)
% filename = strcat(demtxt(1:length(demtxt)-4),'_plateauRemove.tif');
% 
% GRIDobj2geotiff(DEMe,filename)
% 



