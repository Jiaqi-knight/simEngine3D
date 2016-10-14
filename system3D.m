classdef system3D < handle
    % Filename: system.m
    % Author:   Samuel Acuña
    % Date:     11 Oct 2016
    % About:
    % the dynamical system we are modeling. This is made up of the
    % collection of bodies, as well as constraints. This also
    % estabilishes the global inertial reference frame.
    
    properties
        r = [0;0;0]; % Global Reference Frame
        p = [1;0;0;0]; % Global Euler Parameters
        body; % collection of bodies in the system
        cons; % collection of constraints in the system
    end
    
    properties (Dependent)
        nBodies; % number of bodies in the system
        nConstraints; % number of constraints in the system
    end
    
    methods (Access = public)
        function sys = system3D() %constructor function
            sys.body = {}; % no bodies defined yet
            sys.cons = {}; % no constraints defined yet
        end
        function addBody(sys,varargin) % add a body to the system
            ID = sys.nBodies+1; %body ID number
            sys.body{ID} = body3D(ID,varargin{:}); % new instance of the body class
        end
        function addConstraint(sys,constraintName,varargin) % add constraint to the system            
            ID = sys.nConstraints+1; %constraint ID number
            switch  constraintName
                case 'dp1'
                    sys.cons{ID} = constraint.dp1(varargin{:}); % new instance of dp1 class
                case 'cd'
                    sys.cons{ID} = constraint.cd(varargin{:}); % new instance of dp1 class
                otherwise
                    error('Constraint not implemented yet.');
            end
        end
        function plot(sys,varargin) % plots the bodies in the system
            % wrapper to plot functions
            plot.plotSystem(sys,varargin{:});
        end
    end
    methods % methods block with no attributes
        function nBodies = get.nBodies(sys) % calculate number of bodies in system
            nBodies = length(sys.body);
        end
        function nConstraints = get.nConstraints(sys) % calculate number of constraints in system
            nConstraints = length(sys.cons);
        end
    end
end

