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
        q;              % q = [r;p] vector, generalized coordinates of ungrounded bodies
        qdot;           % time derivaite of q vector
        qddot;          % 2nd time derivaite of q vector
        phi;            % full constraint matrix
        phi_q;          % jacobian of phiF
        nu;             % RHS of velocity equation
        gammaHat;       % RHS of acceleration equation, in r-p formulation
        time;           % current time within the system
    end
    
    properties (Dependent)
        bodyIDs;            % ID numbers of ungrounded bodies in the system
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
            sys.time = [];
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
        function assembleConstraints(sys) % construct phi matrices 
            sys.constructPhi();   % construct phi matrix
            sys.constructPhi_q(); % construct phi_q matrix
            sys.constructNu();    % construct RHS of velocity equation
            sys.constructGammaHat();%construct RHS of acceleration equation
            sys.constructQ();     % construct q (r and p) of bodies in system
        end
        function state = kinematicsAnalysis(sys,timeStart,timeEnd,timeStep) % perform kinematics analysis
            % perform kinematics analysis on the system. System must be
            % fully constrained (nDOF = 0).
            % time inputs are measured in seconds.
            % output:
            %   state = [time x [q,qdot,qddot]]
            %           e.g. for one body with 100 timesteps: [100x7x3]
            
            if ~exist('timeStep','var') || isempty(timeStep)
                timeStep = 10e-3; % default time step
            end
            
            timeGrid = timeStart:timeStep:timeEnd;
            
            state = zeros(length(timeGrid),sys.nGenCoordinates,3); % preallocate for saved state
            
            % iterate throughout the time grid            
            for iT = 1:length(timeGrid)
                t = timeGrid(iT); % current time step
                sys.setSystemTime(t); % set system time
                
                % Position Analysis
                if t ~= timeStart % except at initial conditions
                    % solve for r and p
                    tolerance = 1e-9;  % tolerance
                    maxIterations = 50; % maximum iterations
                    sys.positionAnalysis(tolerance,maxIterations)
                end
                
                % to prevent singularities, check jacobian. Later...
                
                % solve for velocities
                sys.velocityAnalysis();
                
                % solve for accelerations
                sys.accelerationAnalysis();
                
                % store system state
                state(iT,:,:,:) = sys.storeSystemState();
                
                disp(['Kinematics analysis completed for t = ' num2str(t) ' sec.']);
            end
        end
        function plot(sys,varargin) % plots the bodies in the system
            % wrapper to plot functions
            plot.plotSystem(sys,varargin{:});
        end
    end
    methods (Access = private)
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
            bodyID = sys.bodyIDs;
            
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
                sys.nu(row:row_new)  =  sys.cons{i}.nu; % plug in constraint nu's
                row = row_new + 1;
            end
        end
        function constructGammaHat(sys) % construct gammaHat, RHS Of acceleration equation, in r-p formulation
            % gammaHat = [nConstrainedDOF x 1]
            sys.gammaHat = zeros(sys.nConstrainedDOF,1); % initialize nu for speed
            row = 1;
            for i = 1:sys.nConstraints
                row_new = row + sys.cons{i}.rDOF - 1;
                sys.gammaHat(row:row_new)  =  sys.cons{i}.gammaHat; % plug in constraint gammaHats's
                row = row_new + 1;
            end
        end
        function constructQ(sys) % construct q (r and p) of bodies in system
            % get ID's of bodies that are not grounded
            bodyID = sys.bodyIDs;
            
            % initialize for speed
            r = zeros(3*length(bodyID),1);
            p = zeros(4*length(bodyID),1);
            
            % pull r and p out of ungrounded bodies
            rowr = 1; % r count
            rowp = 1; % p count
            for i = bodyID  
                r(rowr:rowr+2) = sys.body{i}.r;
                p(rowp:rowp+3) = sys.body{i}.p;
                rowr = rowr + 3;
                rowp = rowp + 4;
            end
            
            sys.q = [r;p]; % assemble q vector
        end
        function positionAnalysis(sys,tolerance,maxIterations) % position analysis of system
            % Performs kinematic position analysis using Newton-Raphson method

            guessQ = sys.q; % initial guess
            
            for i = 1:maxIterations % iterate
                sys.constructPhi();   % get constraint matrix
                sys.constructPhi_q(); % get jacobian
                correction = sys.phi_q\sys.phi; % correction 
                guessQ = guessQ - correction; % update guess
                
                % update bodies
                sys.updateSystem(guessQ,[],[]);
                
                if norm(correction) < tolerance %check for convergence
                    break;
                end
            end
            if i >= maxIterations
                error('Failed to converge position analysis, the maximum number of iterations is reached')
            end
            
        end
        function velocityAnalysis(sys) % solve for velocities of the system
            % compute velocity and derivative of euler parameters
            % for a given time (t)
            
            % calculate RHS of velocity equation
            sys.constructNu();
            
            % find jacobian of phi
            sys.constructPhi_q();
            
            % solve for velocities (qdot)
            qdot = sys.phi_q\sys.nu;
            
            % update bodies
            sys.updateSystem([],qdot,[]);
        end
        function accelerationAnalysis(sys) % solve for accelerations of the system
            % compute acceleration and 2nd derivative of euler parameters
            % for a given time (t)
            
            % calculate RHS of acceleration equation
            sys.constructGammaHat();
            
            % already find jacobian of phi
            
            % solve for accelerations (qddot)
            qddot = sys.phi_q\sys.gammaHat;
            
            % update bodies
            sys.updateSystem([],[],qddot);
        end
        function updateSystem(sys,q,qdot,qddot) %update r,p,and derivatives for each body in system
            %%% update q (r and p)
            if ~isempty(q)
                sys.q = q;
                r = q(1:3*(sys.nBodies-sys.nGrounds));
                p = q(3*(sys.nBodies-sys.nGrounds)+1:end);
                
                rowr = 1;
                rowp = 1;
                for i = sys.bodyIDs
                    sys.body{i}.r = r(rowr:rowr+2);
                    sys.body{i}.p = p(rowp:rowp+3);
                    rowr = rowr + 3;
                    rowp = rowp + 4;
                end
            end
            
            %%% update qdot (rdot and pdot)
            if ~isempty(qdot)
                sys.qdot = qdot;
                rdot = qdot(1:3*(sys.nBodies-sys.nGrounds));
                pdot = qdot(3*(sys.nBodies-sys.nGrounds)+1:end);
                
                rowr = 1;
                rowp = 1;
                for i = sys.bodyIDs
                    sys.body{i}.rdot = rdot(rowr:rowr+2);
                    sys.body{i}.pdot = pdot(rowp:rowp+3);
                    rowr = rowr + 3;
                    rowp = rowp + 4;
                end
            end
            
            %%% update qddot (rddot and pddot)
            if ~isempty(qddot)
                sys.qddot = qddot;
                rddot = qddot(1:3*(sys.nBodies-sys.nGrounds));
                pddot = qddot(3*(sys.nBodies-sys.nGrounds)+1:end);
                
                rowr = 1;
                rowp = 1;
                for i = sys.bodyIDs
                    sys.body{i}.rddot = rddot(rowr:rowr+2);
                    sys.body{i}.pddot = pddot(rowp:rowp+3);
                    rowr = rowr + 3;
                    rowp = rowp + 4;
                end
            end
        end
        function setSystemTime(sys,t) % set time for entire system
            % set system time
            sys.time = t;
            % update constraints time
            for i = 1:sys.nConstraints
                sys.cons{i}.t = t;
            end
        end
        function state = storeSystemState(sys) 
            % store position and orientation information for the current
            % time step
            state = [sys.q,sys.qdot,sys.qddot];
        end
    end
    methods % methods block with no attributes
        function bodyIDs = get.bodyIDs(sys) % calculate ID numbers of ungrounded bodies in the system
            bodyIDs = zeros([sys.nBodies-sys.nGrounds],1); %initialize size
            j = 1; % counter 
            for i = 1:sys.nBodies
                if ~sys.body{i}.isGround
                    bodyIDs(j) = sys.body{i}.ID; % pull body ID
                    j = j+1; % increment counter
                end
            end
        end
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

