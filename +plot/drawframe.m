
function drawframe(r, ORIENTATION, scale, alternate_color)
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

if ~exist('alternate_color', 'var'),
    alternate_color = false;
end

if numel(ORIENTATION) == 9 % orientation matrix
    A = ORIENTATION;
elseif numel(ORIENTATION) == 4 % euler parameters
    A = utility.p2A(ORIENTATION); % convert to orientation matrix
else
    error('Orientation parameter is the wrong size')
end


plot3(r(1), r(2), r(3));

hchek = ishold;
hold on

if (isequal(zeros(3,1), r)) && (isequal(eye(3),A))
    % use gray for the base frame
    plot.arrow3(r, scale*A(1:3,1), 'k');
    plot.arrow3(r, scale*A(1:3,2), 'k');
    plot.arrow3(r, scale*A(1:3,3), 'k');
    
    Atext = A+scale*0.05;
    text(Atext(1,1),Atext(2,1),Atext(3,1), 'X');
    text(Atext(1,2),Atext(2,2),Atext(3,2), 'Y');
    text(Atext(1,3),Atext(2,3),Atext(3,3), 'Z');
else
    if alternate_color,
        plot.arrow3(r, scale*A(1:3,1), 'c');
        plot.arrow3(r, scale*A(1:3,2), 'm');
        plot.arrow3(r, scale*A(1:3,3), 'k');
    else
        plot.arrow3(r, scale*A(1:3,1), 'r');
        plot.arrow3(r, scale*A(1:3,2), 'g');
        plot.arrow3(r, scale*A(1:3,3), 'b');
    end
end

xlabel('x');
ylabel('y');
zlabel('z');

if hchek == 0
    hold off
end
end

