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
        r = [0;0;0];    % Global Reference Frame
        p = [1;0;0;0];  % Global Euler Parameters
        body;           % collection of bodies in the system, inlcuding ground bodies
        cons;           % collection of constraints in the system
        phi;            % full constraint matrix
        phi_q;          % jacobian of phiF
        nu;             % RHS of velocity equation
        gammaHat;       % RHS of acceleration equation, in r-p formulation          
    end
    
    properties (Dependent)
        nBodies;            % number of bodies in the system
        nGrounds;           % number of grounded bodies in the system
        nGenCoordinates;    % number of generalized coordinates in system
        nConstraints;       % number of constraints in the system
        nConstrainedDOF;    % number of degrees of freedom that have been constrained
        nDOF;               % number of Degrees Of Freedom of the system
    end
    
    methods (Access = public)
        function sys = system3D() %constructor function
            sys.body = {}; % no bodies defined yet
            sys.cons = {}; % no constraints defined yet
            
            % empty for now, construct with assembleConstraints function:
            sys.phi = []; 
            sys.phi_q = []; 
            sys.nu = [];
            sys.gammaHat = [];
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
            sys.constructPhi();   % construct phiF matrix
            sys.constructPhi_q(); % construct phiF_q matrix
            sys.constructNu();    % construct RHS of velocity equation
            sys.constructGammaHat();%construct RHS of acceleration equation
        end
        function plot(sys,varargin) % plots the bodies in the system
            % wrapper to plot functions
            plot.plotSystem(sys,varargin{:});
        end
        function constructPhi(sys) % construct phiF matrix (full constraint matrix)
            % phi = [nConstrainedDOF x 1]
            sys.phi = zeros(sys.nConstrainedDOF,1); % initialize phi for speed
            row = 1;
            for i = 1:sys.nConstraints
                row_new = row + sys.cons{i}.rDOF - 1;
                sys.phi(row:row_new)  =  sys.cons{i}.phi; % plug in constraint phi's
                row = row_new + 1;
            end
        end
        function constructPhi_q(sys) % construct jacobian of phi
            % phi_q = partial derivative of sys.phi with respect to the
            % generalized coordinates
            % phi_q = [nConstrainedDOF x nGenCoordinates]
            %
            % For every constraint, pull phi_r and phi_p, then insert into
            % correct location for the jacobian. This location depends on
            % the number of bodies in the system. Grounded bodies not included. 
            % Takes the form:
            %   [phi_r(body1) ... phi_r(body_last) | phi_p(body1) ... phi_p(body_last)]
            
            phi_r = zeros(sys.nConstrainedDOF,3*(sys.nBodies-sys.nGrounds)); % initialize phi_r for speed
            phi_p = zeros(sys.nConstrainedDOF,4*(sys.nBodies-sys.nGrounds)); % initialize phi_p for speed
            
            % get ID's of bodies that are not grounded
            bodyID = zeros([sys.nBodies-sys.nGrounds],1); %initialize size
            j = 1; % counter 
            for i = 1:sys.nBodies
                if ~sys.body{i}.isGround
                    bodyID(j) = sys.body{i}.ID; % pull body ID
                    j = j+1; % increment counter
                end
            end
            
            % loop through constraints
            row = 1; %counter
            for i = 1:sys.nConstraints
                row_new = row + sys.cons{i}.rDOF - 1;
                if strcmp(class(sys.cons{i}),'constraint.p_norm') % if euler parameter normalization constraint
                    col = find(bodyID == sys.cons{i}.bodyi.ID); % get column position
                    phi_r(row,(3*col - 2):(3*col)) = sys.cons{i}.phi_r; % place in row
                    phi_p(row,(4*col - 3):(4*col)) = sys.cons{i}.phi_p;
                elseif ~sys.cons{i}.bodyi.isGround % if bodyi is not ground....
                    coli = find(bodyID == sys.cons{i}.bodyi.ID); % get body i column position
                    phi_r(row:row_new,(3*coli - 2):(3*coli)) = sys.cons{i}.phi_r(1:sys.cons{i}.rDOF,1:3); % place in row
                    phi_p(row:row_new,(4*coli - 3):(4*coli)) = sys.cons{i}.phi_p(1:sys.cons{i}.rDOF,1:4);
                    if ~sys.cons{i}.bodyj.isGround % and if body j is not ground
                        colj = find(bodyID == sys.cons{i}.bodyj.ID); % get body j column position
                        phi_r(row:row_new,(3*colj - 2):(3*colj)) = sys.cons{i}.phi_r(1:sys.cons{i}.rDOF,4:6); % place in row
                        phi_p(row:row_new,(4*colj - 3):(4*colj)) = sys.cons{i}.phi_p(1:sys.cons{i}.rDOF,5:8);
                    end
                elseif sys.cons{i}.bodyi.isGround % if bodyi is ground....
                    colj = find(bodyID == sys.cons{i}.bodyj.ID); % get body j column position
                    phi_r(row:row_new,(3*colj - 2):(3*colj)) = sys.cons{i}.phi_r(1:sys.cons{i}.rDOF,1:3); % place in row
                    phi_p(row:row_new,(4*colj - 3):(4*colj)) = sys.cons{i}.phi_p(1:sys.cons{i}.rDOF,1:4);
                end
                row = row_new + 1; % increment row counter
            end
            
            % create phi_q
            sys.phi_q = [phi_r phi_p];
        end
        function constructNu(sys) % construct nu, RHS Of velocity equation
            % nu = [nConstrainedDOF x 1]
            sys.nu = zeros(sys.nConstrainedDOF,1); % initialize nu for speed
            row = 1;
            for i = 1:sys.nConstraints
                row_new = row + sys.cons{i}.rDOF - 1;
                sys.phi(row:row_new)  =  sys.cons{i}.nu; % plug in constraint nu's
                row = row_new + 1;
            end
        end
        function constructGammaHat(sys) % construct gammaHat, RHS Of acceleration equation, in r-p formulation
            % gammaHat = [nConstrainedDOF x 1]
            sys.gammaHat = zeros(sys.nConstrainedDOF,1); % initialize nu for speed
            row = 1;
            for i = 1:sys.nConstraints
                row_new = row + sys.cons{i}.rDOF - 1;
                sys.phi(row:row_new)  =  sys.cons{i}.gammaHat; % plug in constraint nu's
                row = row_new + 1;
            end
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
        function nConstrainedDOF = get.nConstrainedDOF(sys) % calculate number of degrees of freedom that have been constrained
            nConstrainedDOF = sys.nGenCoordinates - sys.nDOF;
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

