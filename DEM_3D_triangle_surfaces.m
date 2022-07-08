% This function discretize 3D catchment surfaces into 3D triangulations and
% calculate the unit norm, aspect and area of these triangulations.
% Input: 
%       DEM in the form of GRIDobj
%       threshold area of triangulations, no more than area of 2 cells
% Output: a struct 
% 
% Author: Yanyan Wang (yanyan.wang@erdw.ethz.ch) 
% Date: 2 March, 2018
% 
function tri3D = DEM_3D_triangle_surfaces(varargin)
i = length(varargin);
if i==1
   DEM = varargin{1};%  
   thresArea = DEM.cellsize^2; % default value for 90m DEM
else if i==2
        DEM = varargin{1};%  
        thresArea = varargin{2};% user defined value,  
    end    
end

DEMZ = DEM.Z;
num=numel(DEMZ); % the number of all data points
IND=1:1:num;
[rows,cols] = ind2sub(DEM.size,IND);
% xy is the geographic x-y coordinates data points
xydem =  [double(rows(:)) double(cols(:)) ones(numel(rows),1)] * DEM.refmat;
xoriginal = reshape(xydem(:,1), DEM.size);% x coordinates for the whole rectagular domain that enclose the basin
yoriginal = reshape(xydem(:,2), DEM.size);% y coordinates for the whole rectagular domain that enclose the basin
zoriginal = double(DEM.Z);   
% find out non nan points and their coodinates
nanID = ~isnan(DEM.Z);
x = xoriginal(nanID);
y = yoriginal(nanID);
z = zoriginal(nanID);
    
%% discretize basin surfaces into 3D triangulations and calculate their surface norms, aspect and areas
% tri = delaunay(X,Y) creates a 2-D Delaunay triangulation. 'tri' is a matrix representing the set of triangles that make up the triangulation.
tri = delaunay(x,y);
p = [x,y,z];
% find out these triangles that 
zz = ones(size(z));
pp = [x,y,zz];
% Obtain the edges in each triangle formed by the 'delaunaytriangulation'
v1 = pp(tri(:,2),:)-pp(tri(:,1),:);
v2 = pp(tri(:,3),:)-pp(tri(:,2),:);
% Calculating the cross product of the edges in each triangle of the surface
cp = 0.5*cross(v1,v2); % vector whose value is triangulations area 
% Surface area of the entire surface is calculated as the sum of the areas of the individual triangles
Ai = sqrt(dot(cp, cp, 2));
% find the id of triangles whose projection area on the x-o-y plane is no
% more than one cell area
id = Ai<=thresArea;
ID = find(id);
% Obtain the edges in each triangle formed by the 'delaunaytriangulation'
% in 3D space
vv1 = p(tri(:,2),:)-p(tri(:,1),:);
vv2 = p(tri(:,3),:)-p(tri(:,2),:);
% find the unit norm of these 3D surfaces
FN = cross(vv1(ID,:),vv2(ID,:),2);
mag = sqrt(dot(FN,FN,2));
% unit norm of these 3D surfaces
Sunorm = ones(size(FN));
for i = 1:length(mag)
    Sunorm(i,:) = FN(i,:)/mag(i);    
end

% aspects of these 3D triangle surfaces
Sasp = cart2pol(Sunorm(:,1),Sunorm(:,2));
Sasp = mod(90+Sasp/pi*180,360);

% area of these 3D surfaces
ccp = 0.5*cross(vv1,vv2); % vector whose value is triangulations area
AAi = sqrt(dot(ccp, ccp, 2));
SA = AAi(id);


tri3D.unitnorm = Sunorm;% unit norm of these 3D surfaces
tri3D.aspect = Sasp;% aspects of these 3D triangle surfaces
tri3D.area = SA;% area of these 3D surfaces













