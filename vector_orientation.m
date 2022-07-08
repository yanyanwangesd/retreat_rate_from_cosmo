% this function calculate orientation of a 2D vector, like the east-oriented
% vector [0 1 0], with respect to a user defined North-direction vector.
% default North direction vector is [0 1 0]
% it returns the orientation of the vector, in degree
% it use the conventional definition of orientation: 0/360-North 90-East,
% 180-South 270-West
% it calcilate the orientation: first input Point--> second input point
% 

% Input: the start and end point of a vevtor p1(x1, y1,0), p2(x2,y2,0) and
%        a user defined North vector p(x0, y0, 0)
% Syntax: 
%       p = [1 0 0]; q = [1 2 0]; 
%       orient_p2q = vector_orientation(p, q); 
%     Or: 
%       p = [1 0 0]; q = [1 2 0]; North_Ref = [ -1,-2 , 0]
%       orient_p2q = vector_orientation(p, q, North_Ref);

% Author: Yanyan Wang on Jan.21, 2018

function [orient] = vector_orientation(varargin)
    Pstart = varargin{1};%  tail point of the vector
    Pend = varargin{2};%  head point of the vector
    num = length(varargin);
    switch num
        case 2
        north = [ 0 1 0]; % default North
        case 3
        north = varargin{3}; % user defined North
    end
    vector = Pend - Pstart;
    vect_scaled = vector/norm(vector);
    beta = cross(vect_scaled,north);
    temp = asin(beta(3))/pi*180; % sin(angle), asin.m returns [-pi/2 pi/2]
    beta2 = dot(vect_scaled,north);
    temp2 = acos(beta2)/pi*180; % cos(angle), acos.m returns [0 pi]       
    if beta(3)>=0
       orient = temp2;
    end
    if beta(3)<0 && beta2<0
       orient=180-temp;          
    end
    if beta(3)<0 && beta2>0
      orient=360+temp;
    end

