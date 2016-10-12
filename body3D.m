classdef body3D < handle
    % Filename: body.m
    % Author:   Samuel Acuña
    % Date:     11 Oct 2016
    % About:
    % a generic spatial rigid body, including all the attributes of the
    % body, as well as local points on the body.This also
    % estabilishes the local reference frame.
   properties
      ID;   % body ID number
      r; % = [x;y;z] location of body RF from GLOBAL RF
      p; % = [e0;e1;e2;e3] euler parameters of BODY RF
      pdot; % = [4x1] time derivative of euler parameters p
      m; % mass of the part
      J; % inertia tensor of the part
      point; % structure of points on the body, defined in BODY RF
   end
   properties (Dependent)
       nPoints; % number of points on the body
       A; % rotation matrix, an expression of the euler parameters
   end
   methods (Access = public)
       function body = body3D(ID,r,p,pdot) %constructor function
            body.ID = ID;
            body.r = r;
            body.p = p;
            body.pdot = pdot;
            body.point = {}; % no points defined yet
       end
        function addPoint(body,r) % add a body to the system
            % inputs:
            %    - body: body you are adding point to
            %    - r: [3x1] position of point in BODY RF
            if nargin < 2
                r = [0;0;0]; % point located at origin of body
            end
            body.point{body.nPoints+1} = r;
        end

   end
   methods % methods block with no attributes
       function nPoints = get.nPoints(body) % calculate number of bodies in system
           nPoints = length(body.point);
       end
       function A = get.A(body) % calculate rotation matrix, an expression of the euler parameters
           A = utility.p2A(body.p);
       end
   end
end 