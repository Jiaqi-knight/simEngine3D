classdef system3D < handle
    % Filename: system.m
    % Author:   Samuel Acuña
    % Date:     11 Oct 2016
    % About:
    % the dynamical system we are modeling. This is made up of the
    % collection of bodies, as well as constraints. This also
    % estabilishes the global inertial reference frame.
    
    properties
        % RF = struct('O',[0;0;0],'X',[1;0;0], 'Y',[0;1;0], 'Z',[0;0;1]); 
        r = [0;0;0]; % Global Reference Frame
        p = [1;0;0;0]; % Global Euler Parameters
        bodies; % collection of bodies in the system
    end
    
    properties (Dependent)
        nBodies; % number of bodies in the system
    end
    
    methods
        function sys = system3D() %constructor function
            sys.bodies = {};
        end
    end
    methods (Access = public)
        function addBody(sys,r,p) % add a body to the system
            % inputs:
            %    - sys: system you are adding body to
            %    - r: [3x1] position of body in GLOBAL RF
            %    - p: [4x1] euler parameters of the BODY RF
            if nargin < 3
                p = [1;0;0;0]; % no change in body orientation from GLOBAL RF
            end
            if nargin < 2
                r = [0;0;0] % body located at origin
            end
            
            ID = sys.nBodies+1; %body ID number
            sys.bodies{ID} = body(ID,r,p); % new instance of the body class
        end
        function plot(sys,frames) % plots the bodies in the system
            if nargin == 1
                frames = 0;
            end
            figure()
            plot.drawframe(sys.r,sys.p) % plot GLOBAL RF
            
            if sys.nBodies == 0; disp('No bodies in the system.'); return; end;
            
            hold on;
            for i = 1:sys.nBodies % plot bodies in system
                
                % plot marker for every body
                r = [sys.bodies{i}.r(1);sys.bodies{i}.r(2);sys.bodies{i}.r(3)];
                scatter3(r(1),r(2),r(3));
                
                if frames %optionally plot body reference frames
                    plot.drawframe(sys.bodies{i}.r,sys.bodies{i}.p)
                end
                
                % add body labels
                r = r + 0.05; %adjust so label is not right over point
                text(r(1),r(2),r(3),num2str(sys.bodies{i}.ID));
            end
            hold off;
            
        end
    end
    methods % methods block with no attributes
        function nBodies = get.nBodies(sys) % calculate number of bodies in system
            nBodies = length(sys.bodies);
        end
    end
end

