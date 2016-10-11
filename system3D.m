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
        body; % collection of bodies in the system
    end
    
    properties (Dependent)
        nBodies; % number of bodies in the system
    end
    
    methods
        function sys = system3D() %constructor function
            sys.body = {}; % no bodies defined yet
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
                r = [0;0;0]; % body located at origin
            end
            
            ID = sys.nBodies+1; %body ID number
            sys.body{ID} = body3D(ID,r,p); % new instance of the body class
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
                r = [sys.body{i}.r(1);sys.body{i}.r(2);sys.body{i}.r(3)];
                scatter3(r(1),r(2),r(3),500);
                
                if frames %optionally plot body reference frames
                    plot.drawframe(sys.body{i}.r,sys.body{i}.p)
                end
                
                % add body labels
                r_text = r + 0.05; %adjust so label is not right over body origin
                text(r_text(1),r_text(2),r_text(3),num2str(sys.body{i}.ID));
                
                % if points on body, plot them.
                if sys.body{i}.nPoints ~= 0
                    for j = 1:sys.body{i}.nPoints % plot bodies in system
                        sbar = [sys.body{i}.point{j}(1);sys.body{i}.point{j}(2);sys.body{i}.point{j}(3)]; % local position of point
                        Asbar = utility.p2A(sys.body{i}.p)*sbar; % rotated to global rf
                        rAsbar = r + Asbar; % global position of point
                        
                        plot3([r(1); rAsbar(1)],[r(2); rAsbar(2)],[r(3); rAsbar(3)],'ks-') % line from BODY RF to point
                        
                        % add point labels
                        rAsbar_text = rAsbar + 0.05; %adjust so label is not right over point
                        text(rAsbar_text(1),rAsbar_text(2),rAsbar_text(3),num2str(j));
                    end
                end
            end
            hold off;
        end
    end
    methods % methods block with no attributes
        function nBodies = get.nBodies(sys) % calculate number of bodies in system
            nBodies = length(sys.body);
        end
    end
end

