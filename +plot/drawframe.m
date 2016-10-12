
function drawframe(r, ORIENTATION, scale)
%DRAWFRAME  plots a graphical description of a coordinate frame
%
%	DRAWFRAME(r,A)
%	DRAWFRAME(r,A, SCALE)
%
% ORIENTATION is a rotation matrix [3x3] OR euler parameters [4x1], r is a point.
%
% adapted from:
% $Id: drawframe.m,v 1.1 2009-03-17 16:40:18 bradleyk Exp $
% Copyright (C) 2005, by Brad Kratochvil

if 2  == nargin,
    scale = 1;
end

if numel(ORIENTATION) == 9 % orientation matrix
    A = ORIENTATION;
elseif numel(ORIENTATION) == 4 % euler parameters
    A = utility.p2A(ORIENTATION); % convert to orientation matrix
else
    error('Orientation parameter is the wrong size')
end

hchek = ishold;
hold on

% scaled, orthogonal unit vectors for frame
unitX = scale*[1;0;0]; 
unitY = scale*[0;1;0];
unitZ = scale*[0;0;1];

% rotated unit vectors
x = A*unitX;
y = A*unitY;
z = A*unitZ;

if (isequal(zeros(3,1), r)) && (isequal(eye(3),A))
    % use black for the inertial frame
    quiver3(r,r,r,x,y,z,'Color',[0,0,0]);

    % label the inertial frame axes
    Atext = A+scale*0.05;
    text(Atext(1,1),Atext(2,1),Atext(3,1), 'X');
    text(Atext(1,2),Atext(2,2),Atext(3,2), 'Y');
    text(Atext(1,3),Atext(2,3),Atext(3,3), 'Z');
else
    %use rgb for body frames
    quiver3(r(1),r(2),r(3),x(1),x(2),x(3),'Color',[1,0,0]);
    quiver3(r(1),r(2),r(3),y(1),y(2),y(3),'Color',[0,1,0]);
    quiver3(r(1),r(2),r(3),z(1),z(2),z(3),'Color',[0,0,1]);
end


if hchek == 0
    hold off
end
end

