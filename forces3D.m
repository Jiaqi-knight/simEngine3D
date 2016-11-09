classdef forces3D < handle
    % Filename: forces3D.m
    % Author:   Samuel Acuña
    % Date:     31 Oct 2016
    %
    % About:
    % This class handles instances of forces and torques on a body.
    % It Uses the r-p formulation (euler parameters). Each instance of this
    % class is a force and torque coupling, since sometimes they are
    % coupled.
    % Observe that the unitDirection inputs differ between torques and
    % forces. The unitDirection for forces is experessed in the GLOBAL RF
    % and the unitDirection for torques is expressed in the BODY RF
    
    properties
        body;           % parent body3D object to which is this force/torque is applied
        forceLocation;  % vector, where the applied force is located, expressed in BODY RF
        forceMagnitude; % magnitude of force, often constant, but can be a function of time
        forceDirection; % [3x1] unit vector, the direction of the applied force, in GLOBAL RF
        torqueMagnitude; %  magnitude of torque applied, expressed in BODY RF
        torqueDirection; % [3x1] unit vector, the direction of the applied torque, in BODY RF
    end
    properties (Dependent)
        force;      % [3x1] vector, force applied at a point, expressed in GLOBAL RF
        torque;     % [3x1] vector, torque applied, expressed in BODY RF
        %force_r;    % [3x3] matrix, partial derivative of force w.r.t. body position
        %force_p;    % [3x4] matrix, partial derivative of force w.r.t. body orientation (euler parameters)
        %force_rdot; % [3x3] matrix, partial derivative of force w.r.t. body velocity
        %force_pdot; % [3x4] matrix, partial derivative of force w.r.t. body angular velocity (euler parameters)
        %torque_r;   % [3x3] matrix, partial derivative of torque w.r.t. body position
        %torque_p;   % [3x4] matrix, partial derivative of torque w.r.t. body orientation (euler parameters)
        %torque_rdot;% [3x3] matrix, partial derivative of torque w.r.t. body velocity
        %torque_pdot;% [3x4] matrix, partial derivative of torque w.r.t. body angular velocity (euler parameters)
    end
    methods
        function f = forces3D(body,forcesType,varargin) % constructor function
            f.body = body;
            % see specific forcesType functions for info on required inputs
            switch forcesType % choose type of force / torque to implement
                case 'torque'
                   f.addTorque(varargin{1},varargin{2});
                case 'force' % specify force location with body pointID
                   f.addForce(varargin{1},varargin{2},varargin{3});
                case 'force2' % specify force location as vector in BODY RF
                   f.addForce2(varargin{1},varargin{2},varargin{3});
                otherwise
                    error('Forces/Torques of this type not implemented yet.');
            end
        end
    end
    methods (Access = private)
        function addTorque(f,torqueMagnitude,unitDirection)
            % from ME751_f2016 slide 27 of lecture 10/03/16.
            % since torque acts along the whole rigid body, no need to
            % specify a point of application
            % INPUTS:
            %   torqueMagnitude   : magnitude of torque applied, expressed
            %                       in BODY RF. 
            %                       e.g.  torqueMagnitude = 5;
            %
            %                       Can also be a matlab anonymous function
            %                       of torque applied, as function of time
            %                       e.g.  torqueMagnitude = @(t) 20*sin(t)
            %
            %   unitDirection : [3x1] unit vector of the direction of the
            %                   applied torque. for example, for a torque
            %                   about the body Z axis, = [0;0;1]. This is
            %                   expressed in the BODY RF
            
            % check direction is a unit vector
            if norm(unitDirection) ~= 1
                warning('Normalized the applied torque unit direction');
                unitDirection = unitDirection/norm(unitDirection);
            end
            
            f.forceLocation = [0;0;0];
            f.forceMagnitude = 0;
            f.forceDirection = [1;0;0]; % random unit vector
            f.torqueMagnitude = torqueMagnitude;
            f.torqueDirection = unitDirection;
        end
        function addForce(f,pointID,forceMagnitude,unitDirection) % add constant active force to body
            % from ME751_f2016 slide 27 of lecture 10/03/16.
            % sometimes a force has a corresponding torque.
            % INPUTS:
            %   pointID     : ID of the body point where the force is applied
            %   forceMagnitude    : magnitude of force applied, expressed
            %                       in GLOBAL RF. 
            %                       e.g.  forceMagnitude = 5;
            %
            %                       Can also be a matlab anonymous function
            %                       of force applied, as function of time
            %                       e.g.  forceeMagnitude = @(t) 20*sin(t)
            %   unitDirection : [3x1] unit vector of the direction of the
            %                   applied force. for example, for a torque
            %                   about the body Z axis, = [0;0;1]. This is
            %                   expressed in the GLOBAL RF
            
            % check direction is a unit vector
            if norm(unitDirection) ~= 1
                warning('Normalized the applied force unit direction');
                unitDirection = unitDirection/norm(unitDirection);
            end
            
            f.forceLocation = f.body.point{pointID};
            f.forceMagnitude = forceMagnitude;
            f.forceDirection = unitDirection;
            f.torqueMagnitude = 0;
            f.torqueDirection = [1;0;0]; % random unit vector
        end
        function addForce2(f,forceLocation,forceMagnitude,unitDirection) % add constant active force to body
            % from ME751_f2016 slide 27 of lecture 10/03/16.
            % sometimes a force has a corresponding torque.
            % INPUTS:
            %   forceLocation     : [3x1] location of point where the force
            %                       is applied, in BODY RF
            %   forceMagnitude    : magnitude of force applied, expressed
            %                       in GLOBAL RF. 
            %                       e.g.  forceMagnitude = 5;
            %
            %                       Can also be a matlab anonymous function
            %                       of force applied, as function of time
            %                       e.g.  forceeMagnitude = @(t) 20*sin(t)
            %   unitDirection : [3x1] unit vector of the direction of the
            %                   applied force. for example, for a torque
            %                   about the body Z axis, = [0;0;1]. This is
            %                   expressed in the GLOBAL RF
            
            % check direction is a unit vector
            if norm(unitDirection) ~= 1
                warning('Normalized the applied force unit direction');
                unitDirection = unitDirection/norm(unitDirection);
            end
            
            f.forceLocation = forceLocation;
            f.forceMagnitude = forceMagnitude;
            f.forceDirection = unitDirection;
            f.torqueMagnitude = 0;
            f.torqueDirection = [1;0;0]; % random unit vector
        end
    end
    methods % methods block with no attributes
        function force = get.force(f) % calculate force exerted on body
            if isa(f.forceMagnitude, 'function_handle') % see if force is a function
                time = f.body.system.time;
                force = f.forceMagnitude(time)*f.forceDirection; % evaluate at system time
            else % force is constant
                force = f.forceMagnitude*f.forceDirection;
            end
        end
        function torque = get.torque(f) % calculate torque exerted on body
            if isa(f.torqueMagnitude, 'function_handle') % see if torque is a function
                time = f.body.system.time;
                torque = f.torqueMagnitude(time)*f.torqueDirection; % evaluate at system time
            else 
                if f.forceMagnitude == 0 % applied torque is constant
                    torque = f.torqueMagnitude*f.torqueDirection;
                else % torque as a result of applied force
                    torque = utility.tilde(f.forceLocation)*f.body.A'*f.force;
                end
            end
        end
%         function force_r = get.force_r(f) % partial derivative of force w.r.t. body position
%             % [3x3] matrix, partial derivative of force w.r.t. body position
%             % see ME751_f2016 slide 33 from lecture 10/17/16
%             
%             % NOTE: although force could be a function of body position (r)
%             % At present, we are only supporting force functions of time. 
%             % Thus, force_r will be zeros.
%             force_r = zeros(3,3);
%         end
%         function force_p = get.force_p(f) % partial derivative of force w.r.t. body orientation (euler parameters)
%             % [3x4] matrix, partial derivative of force w.r.t. body orientation (euler parameters)
%             % see ME751_f2016 slide 33 from lecture 10/17/16
%             
%             % NOTE: although force could be a function of body euler parameters (p)
%             % At present, we are only supporting force functions of time. 
%             % Thus, force_p will be zeros.
%             force_p = zeros(3,4);
%         end
%         function force_rdot = get.force_rdot(f) % partial derivative of force w.r.t. body velocity
%             % [3x3] matrix, partial derivative of force w.r.t. body velocity
%             % see ME751_f2016 slide 33 from lecture 10/17/16
%             
%             % NOTE: although force could be a function of body velocity (rdot)
%             % At present, we are only supporting force functions of time. 
%             % Thus, force_rdot will be zeros.
%             force_rdot = zeros(3,3);
%         end
%         function force_pdot = get.force_pdot(f) % partial derivative of force w.r.t. body angular velocity (euler parameters)
%             % [3x4] matrix, partial derivative of force w.r.t. body angular velocity (euler parameters)
%             % see ME751_f2016 slide 33 from lecture 10/17/16
%             
%             % NOTE: although force could be a function of body angular velocity (pdot)
%             % At present, we are only supporting force functions of time. 
%             % Thus, force_pdot will be zeros.
%             force_pdot = zeros(3,4);
%         end
%         function torque_r = get.torque_r(f) % partial derivative of torque w.r.t. body position
%             % [3x3] matrix, partial derivative of torque w.r.t. body position
%             % see ME751_f2016 slide 33 from lecture 10/17/16
%             
%             % NOTE: although torque could be a function of body position (r)
%             % At present, we are only supporting torque functions of time. 
%             % Thus, torque_r will be zeros.
%             torque_r = zeros(3,3);
%         end
%         function torque_p = get.torque_p(f) % partial derivative of torque w.r.t. body orientation (euler parameters)
%             % [3x4] matrix, partial derivative of torque w.r.t. body orientation (euler parameters)
%             % see ME751_f2016 slide 33 from lecture 10/17/16
%             
%             % NOTE: although torque could be a function of body euler parameters (p)
%             % At present, we are only supporting torque functions of time. 
%             % Thus, torque_p will be zeros.
%             torque_p = zeros(3,4);
%         end
%         function torque_rdot = get.torque_rdot(f) % partial derivative of torque w.r.t. body velocity
%             % [3x3] matrix, partial derivative of torque w.r.t. body velocity
%             % see ME751_f2016 slide 33 from lecture 10/17/16
%             
%             % NOTE: although torque could be a function of body velocity (rdot)
%             % At present, we are only supporting torque functions of time. 
%             % Thus, torque_rdot will be zeros.
%             torque_rdot = zeros(3,3);
%         end
%         function torque_pdot = get.torque_pdot(f) % partial derivative of torque w.r.t. body orientation (euler parameters)
%             % [3x4] matrix, partial derivative of torque w.r.t. body angular velocity (euler parameters)
%             % see ME751_f2016 slide 33 from lecture 10/17/16
%             
%             % NOTE: although torque could be a function of body angular velocity (pdot)
%             % At present, we are only supporting torque functions of time. 
%             % Thus, torque_pdot will be zeros.
%             torque_pdot = zeros(3,4);
%         end
    end
end