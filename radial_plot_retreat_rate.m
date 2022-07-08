% this script visualize retreat rate and direction with a radial plot
% Input: 
%       .mat file from script of Retreat_rate_from_10Be_single_basin.m
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Step 1 %%%%%%%%%%%%%%%%%%%%%%%%
%% load .mat file

% add path of your folder to MATLAB search path
addpath('/Users/yanywang/Desktop/Directional_mass_flux_from_10Be/Retreat_rate_from_10Be')

% load the .mat file
load('ghats_0_0_retreat_rate.mat')

% change this number to adjust the limit of your compass radius
maxlim = 500;  

% transfer params for plot

Vasp = DATA(:,1); % retreat direction, [degree]

Rate_dot = DATA(:,2);% retreat rate from Local scalar product method, [m/Ma]

Rate_shadow = DATA(:,5);% retreat rate from Basin projection method, [m/Ma]

erosion = DATA(1,8);
erosion = round(erosion,1);% conventional erosion rate, [m/Ma]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Step 2 %%%%%%%%%%%%%%%%%%%%%%%%
%% radial plot
% (1) plot retreat rate from Basin projection method

figure(1)
% this is to create an empty compass chart
xfake = [0 maxlim 0 -maxlim];
yfake = [maxlim 0 -maxlim 0];
hfake = compass(xfake, yfake);
view([90 -90])
set(hfake,'Visible','off');
set(findall(gcf, 'String', '0'),'String', 'North'); 
set(findall(gcf, 'String', '90'),'String', 'East'); 
set(findall(gcf, 'String', '180'),'String', 'South'); 
set(findall(gcf, 'String', '270'),'String', 'West'); 
set(findall(gcf, 'String', '30', '-or','String','60','-or','String','120','-or','String','150', ...
                  '-or','String','210','-or','String','240','-or','String','300',...
                  '-or','String','330') ,'String', '  ');
set(gca, 'Color', [0.7 0.7 0.7])
fobj = findobj(gcf);
c=findall(fobj, 'Type', 'text');
set(c,'FontSize',16);
hold on

% this is to plot your data
[hh, ff] = pol2cart(Vasp*pi/180, Rate_shadow);
h=compass(hh,ff);
set(h,'LineWidth',1)
% this is to set the compass direction to be North-0, East-90, South-180,
% West-270
view([90 -90])

% to delete the arrow head that is set in compass.m 
for k = 1:length(h)
    a = get(h(k), 'xdata'); 
    b = get(h(k), 'ydata'); 
    set(h(k), 'xdata', a(1:2), 'ydata', b(1:2), 'color', 'r')
end

hold on
str= 'Basin projection method';
text(maxlim*0.95, -maxlim*1.3, str,'FontSize',12,'Color','k');

hold on
str= strcat('Erosion rate = ', num2str(erosion), ' m/Ma');
text(maxlim*0.95, maxlim*0.5, str,'FontSize',12,'Color','k');


% (2) to plot the Local scalar product method result
figure(2)
maxlim = 1.5*maxlim;
% this is to create an empty compass chart
xfake = [0 maxlim 0 -maxlim];
yfake = [maxlim 0 -maxlim 0];
hfake = compass(xfake, yfake);
view([90 -90])
set(hfake,'Visible','off');
set(findall(gcf, 'String', '0'),'String', 'North'); 
set(findall(gcf, 'String', '90'),'String', 'East'); 
set(findall(gcf, 'String', '180'),'String', 'South'); 
set(findall(gcf, 'String', '270'),'String', 'West'); 
set(findall(gcf, 'String', '30', '-or','String','60','-or','String','120','-or','String','150', ...
                  '-or','String','210','-or','String','240','-or','String','300',...
                  '-or','String','330') ,'String', '  ');
set(gca, 'Color', [0.7 0.7 0.7])
fobj = findobj(gcf);
c=findall(fobj, 'Type', 'text');
set(c,'FontSize',16);
hold on

% this is to plot your data
[hh, ff] = pol2cart(Vasp*pi/180, Rate_dot);
h=compass(hh,ff);
set(h,'LineWidth',1)
% this is to set the compass direction to be North-0, East-90, South-180,
% West-270
view([90 -90])

% to delete the arrow head that is set in compass.m 
for k = 1:length(h)
    a = get(h(k), 'xdata'); 
    b = get(h(k), 'ydata'); 
    set(h(k), 'xdata', a(1:2), 'ydata', b(1:2), 'color', 'r')
end

hold on
str= 'Local scalar product method';
text(maxlim*0.95, -maxlim*1.3, str,'FontSize',12,'Color','k');

hold on
str= strcat('Erosion rate = ', num2str(erosion), ' m/Ma');
text(maxlim*0.95, maxlim*0.5, str,'FontSize',12,'Color','k');



