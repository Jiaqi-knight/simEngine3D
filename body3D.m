classdef body3D < handle
    % Filename: body3D.m
    % Author:   Samuel Acuña
    % Date:     11 Oct 2016
    %
    % About:
    % a generic spatial rigid body, including all the attributes of the
    % body, as well as local points on the body. This also
    % estabilishes the local reference frame. Assumes center of mass is the
    % center of the BODY RF.
    properties
        system;   % parent system3D object to which this body is a member. 
                  % useful because it containts system time and gravity
        ID;       % body ID number
        color;    % rgb, used for plotting
        r;        % = [x;y;z] location of body RF from GLOBAL RF        
        p;        % = [e0;e1;e2;e3] euler parameters of BODY RF
        rdot;     % = [3x1] time derivative of position r
        pdot;     % = [4x1] time derivative of euler parameters p
        rddot;    % = [3x1] time derivative of position r
        pddot;    % = [4x1] time derivative of euler parameters p
        m;        % mass of the part
        Jbar;     % inertia tensor of the part, in BODY RF
        isGround; % if body is ground
        point;    % structure of points on the body, defined in BODY RF
        forces;   % collection of forces3D objects, defining forces and torques applied to body
    end
    properties (Dependent)
        nPoints; % number of points on the body
        A; % rotation matrix, an expression of the euler parameters
        nForces; % number of forces/torques applied  
    end
    methods (Access = public)
        function body = body3D(system,ID,bodyType,r,p,rdot,pdot,rddot,pddot) %constructor function
            switch bodyType
                case 'ground'
                    isGround = 1;
                case 'free'
                    isGround = 0;
                otherwise
                    error('specified body type is not permitted');
            end
            if ~exist('r','var') || isempty(r)
                r = [0;0;0]; % body located at origin
            end
            if ~exist('p','var') || isempty(p)
                p = [1;0;0;0]; % no change in body orientation from GLOBAL RF
            end       
            if ~exist('rdot','var') || isempty(rdot) || isGround == 1
                rdot = [0;0;0]; % no velocity change in body position
            end
            if ~exist('pdot','var') || isempty(pdot) || isGround == 1
                pdot = [0;0;0;0]; % no velocity change in body orientation
            end
            if ~exist('rddot','var') || isempty(rddot) || isGround == 1
                rddot = [0;0;0]; % no velocity change in body position
            end
            if ~exist('pddot','var') || isempty(pddot) || isGround == 1
                pddot = [0;0;0;0]; % no velocity change in body orientation
            end
            
            % instantiate body properties
            body.system = system;
            body.ID = ID;
            body.color = [0.6602, 0.6602, 0.6602]; % default color is dark gray.
            body.r = r;
            body.p = p;
            body.rdot = rdot;
            body.pdot = pdot;
            body.rddot = rddot;
            body.pddot = pddot;
            body.m = 0;         % specify with setMass()
            body.Jbar = zeros(3);  % specify with setInertia()
            body.isGround = isGround;
            body.point = {}; % no points defined yet
            body.forces = {}; % no forces or torques defined yet
        end
        function addPoint(body,aBar) % add a point to the body
            % inputs:
            %    - body: body you are adding point to
            %    - aBar: [3x1] position of point in BODY RF
            if ~exist('aBar','var') || isempty(aBar)
                aBar = [0;0;0]; % point located at origin of body
            end
            body.point{body.nPoints+1} = aBar;
        end
        function setMass(body,mass) % specify body mass
            body.m = mass;
        end
        function setInertia(body,inertiaTensor) % specify body inertia tensor
            % inertia tensor specified in BODY RF
            if ~isequal(size(inertiaTensor),[3 3])
                error('body inertia must be a 3x3 matrix')
            end
            body.Jbar = inertiaTensor;
        end
        function addForceOfGravity(body) % add force of gravity
            % from ME751_f2016 slide 29 of lecture 10/03/16.
            ID = body.nForces + 1;
            forceLocation = [0;0;0]; % acts at body center of mass, assumed to be the origin
            forceMagnitude = norm(body.m*body.system.g);
            unitDirection = (body.m*body.system.g)/forceMagnitude;
            body.forces{ID} = forces3D(body,'force2',forceLocation,forceMagnitude,unitDirection); % new instance of the forces3D class
        end
        function addForces(body,varargin) % add active force / torque to body
            ID = body.nForces + 1;
            body.forces{ID} = forces3D(body,varargin{:}); % new instance of the forces3D class
        end
    end
    methods % methods block with no attributes
        function nPoints = get.nPoints(body) % calculate number of points on the body
            nPoints = length(body.point);
        end
        function A = get.A(body) % calculate rotation matrix, an expression of the euler parameters
            A = utility.p2A(body.p);
        end
        function nForces = get.nForces(body)
            % calculate the number of Force-Torque combinations applied to the body
            nForces = length(body.forces);
        end
    end
end