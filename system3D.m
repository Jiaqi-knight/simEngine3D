classdef system3D < handle
    % Filename: system.m
    % Author:   Samuel Acuña
    % Date:     11 Oct 2016
    %
    % About:
    % the dynamical system we are modeling. This is made up of the
    % collection of bodies, as well as constraints. This also
    % estabilishes the global inertial reference frame.
    
    properties
        r = [0;0;0]; % Global Reference Frame
        p = [1;0;0;0]; % Global Euler Parameters
        body; % collection of bodies in the system, inlcuding ground bodies
        cons; % collection of constraints in the system
    end
    
    properties (Dependent)
        nBodies;            % number of bodies in the system
        nGrounds;            % number of grounded bodies in the system
        nGenCoordinates;    % number of generalized coordinates in system
        nConstraints;       % number of constraints in the system
        nDOF;               % number of Degrees Of Freedom of the system
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
        function addConstraint(sys,constraintName,varargin) % add kinematic constraint to the system            
            ID = sys.nConstraints+1; %constraint ID number
            switch  constraintName
                case 'dp1'
                    sys.cons{ID} = constraint.dp1(varargin{:}); % new instance of constraint.dp1 class
                case 'dp2'
                    sys.cons{ID} = constraint.dp2(varargin{:}); % new instance of constraint.dp2 class
                case 'cd'
                    sys.cons{ID} = constraint.cd(varargin{:}); % new instance of constraint.cd class
                case 'd'
                    sys.cons{ID} = constraint.d(varargin{:}); % new instance of constraint.d class
                case 'p1'
                    sys.cons{ID} = constraint.p1(varargin{:}); % new instance of constraint.p1 class
                case 'p2'
                    sys.cons{ID} = constraint.p2(varargin{:}); % new instance of constraint.p2 class
                case 'sj'
                    sys.cons{ID} = constraint.sj(varargin{:}); % new instance of constraint.sj class
                case 'rj'
                    sys.cons{ID} = constraint.rj(varargin{:}); % new instance of constraint.rj class
                case 'p_norm'
                    % add p_norm constraint for every body in system (except ground)
                    for i = 1:sys.nBodies
                        if ~sys.body{i}.isGround
                            sys.cons{ID} = constraint.p_norm(sys.body{i}); % new instance of constraint.p_norm class
                            ID = ID+1;
                        end
                    end
                otherwise
                    error('Constraint not implemented yet.');
            end
        end
        function assembleConstraints(sys)
            disp('function assembleConstraints not implemented yet')
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
        function nGrounds = get.nGrounds(sys) % calculate number of grounded bodies in system
            nGrounds = 0;
            if sys.nBodies>0
                for i = 1:sys.nBodies
                    if sys.body{i}.isGround
                        nGrounds = nGrounds + 1;
                    end
                end
            end
        end
        function nGenCoordinates = get.nGenCoordinates(sys) % calculate number of generalized coordinates in system
            % 7 DOF available for each body
            nGenCoordinates = 7*(sys.nBodies-sys.nGrounds);
        end
        function nConstraints = get.nConstraints(sys) % calculate number of constraints in system
            nConstraints = length(sys.cons);
        end
        function nDOF = get.nDOF(sys) % count number of Degrees of Freedom
            % number of degrees of freedom removed by constraints
            rDOF = 0;
            if sys.nConstraints>0
                for i = 1:sys.nConstraints
                    rDOF = rDOF + sys.cons{i}.rDOF;
                end
            end
            nDOF = sys.nGenCoordinates - rDOF;
        end
    end
end

