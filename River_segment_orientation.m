% this script is to manually pick the characteristic rivers of a basin and
% calculate their orientations

% Author: Yanyan Wang (yanyan.wang@erdw.ethz.ch)
% 
% 

clear;close all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% addpath of the .tif file
addpath('Users/yanywang/Desktop/escarpment')

% name of the .tif file
demtxt = 'ghats_0_0.TIF';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare data

DEM = GRIDobj(demtxt);
crita =1e+6;% critical area for channel head initiation, [m^2]
DEM = fillsinks(DEM);
cs = DEM.cellsize;
DEM = fillsinks(DEM);
FD = FLOWobj(DEM); % flow direction
S = STREAMobj(FD,'minarea',crita/cs^2);
A = flowacc(FD).*(cs^2)/1e+6; % km^2

% index of non-nan data point
IDnan =find(~isnan(DEM.Z));

% find find ID of the outlet
ID = find(A.Z(:)==max(A.Z(:)));
[~,ki] = ismember(ID, S.IXgrid);% ki is the outlet ID


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Manually pick the two ends of river segments, read the instructions in
% the Command Window.

% A) to prepare figure for choose channel segments
figure
imageschs(DEM);
hold on
plot(S, 'w-', 'LineWidth',1);

% B) to manually choose characteristic channel segments that draining through the escarpment face
% B1) to choose channel heads of all characteristic channel segments
fprintf('\nclick on the fiugre to choose channel segments head,press RETURN when you finish\n');
fprintf('\nfor each channel segment, first point should be channel head\n');

[xsegh, ysegh] = ginputc('ShowPoints', true, 'ConnectPoints', true);  
plot(xsegh,ysegh,'go','MarkerFaceColor','g','MarkerSize',10); hold on
for k=1:length(xsegh)
   str = num2str(k);
   text(xsegh(k),ysegh(k), str, 'Color','r','FontSize',17 )
end
% B2) then choose the cross point of channel segments with escarpment
% top line, in the order of channel heads are chosen
fprintf('\nNow choose channel segments outlets,press RETURN when you finish\n');
[xsegot, ysegot] = ginputc('ShowPoints', true, 'ConnectPoints', true); 
hold off
close all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate channel segemnt orientation in terms of normal definition: N-0 E-90 S-180 W-270
% orientation is defined to be the dirction from channel segment outlet to
% channel segment head

Pheads = horzcat(xsegh,ysegh,zeros(size(xsegh)));
Pouts = horzcat(xsegot,ysegot,zeros(size(xsegh)));
% the user defined unit vector in north direction |north|==1
North = [0 1 0];
% channel segment orientation with respect to North, [degree]
channel_orient = zeros(length(xsegh),1);

for ti = 1:length(xsegh)
    Pend = Pheads(ti,:);
    Pstart = Pouts(ti,:); 
    % channel segments orientations with respect to normal North 
    channel_orient(ti)= vector_orientation(Pstart,Pend, North);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print orientation of rivers to screen
fprintf('\n River orientations are (degree): \n')
fprintf('%d\n',round(channel_orient'))






