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
        q;              % q = [r;p] vector, generalized coordinates of ungrounded (free) bodies
        qdot;           % time derivaite of q vector
        qddot;          % 2nd time derivaite of q vector
        phiF;           % full constraint matrix phiF = [phi; phiP]
        phiF_q;         % jacobian of phiF
        nuF;            % RHS of velocity equation
        gammaHatF;      % RHS of acceleration equation, in r-p formulation
        phi;            % set of kinematic and driving constraints
        phi_r;          % jacobian of phi, with respect to r only
        phi_p;          % jacobian of phi, with respect to p only
        phiP;           % set of Euler Parameter Normalization constraints
        M;              % M matrix for Equations of Motion (masses)
        P;              % P matrix for Equations of Motion (euler parameters)
        Jp;             % J^p matrix for Equations of Motion (inertias)
        F;              % F vector used in Equations of Motion (forces)
        TauHat;         % TauHat vector used in Equations of Motion (torques)
        lambda;         % lagrange multipliers of the system, by constraints
        rForce;         % reaction forces for each body in system
        rTorque;        % reaction torques for each body in system
        g;              % vector of acceleration due to gravity
        time;           % current time within the system
    end
    
    properties (Dependent)
        bodyIDs;            % ID numbers of ungrounded bodies in the system
        nBodies;            % number of bodies in the system
        nFreeBodies;        % number of free bodies in the system (ungrounded)
        nGroundBodies;      % number of grounded bodies in the system
        nConstraints;       % number of constraints in the system
        nForces;            % count number of forces/torques applied to bodies
        nGenCoordinates;    % number of generalized coordinates in system
        nConstrainedDOF;    % number of degrees of freedom that have been constrained
        rKDOF;              % number of degrees of freedom removed by kinematic and driving constraints
        rPDOF;              % number of degrees of freedom removed by euler parameter normalization constraints
        nDOF;               % number of Degrees Of Freedom of the system
    end
    
    methods (Access = public)
        function sys = system3D() %constructor function
            sys.body = {}; % no bodies defined yet
            sys.cons = {}; % no constraints defined yet
                        
            % empty for now, constructed later as needed
            sys.q = [];
            sys.qdot = [];
            sys.qddot = [];
            sys.phiF = [];
            sys.phiF_q = [];
            sys.nuF = [];
            sys.gammaHatF = [];
            sys.phi = []; 
            sys.phi_r = [];
            sys.phi_p = [];
            sys.phiP = [];
            sys.M = [];
            sys.P = [];
            sys.F = [];
            sys.lambda = [];
            sys.rForce = [];
            sys.rTorque= [];
            sys.setSystemTime(0);
            sys.g = [0;0;-9.81]; % default gravity is pointed -z direction
        end
        function addBody(sys,varargin) % add a body to the system
            ID = sys.nBodies+1; %body ID number
            sys.body{ID} = body3D(sys,ID,varargin{:}); % new instance of the body class
        end
        function addConstraint(sys,constraintName,varargin) % add kinematic constraint to the system            
            ID = sys.nConstraints+1; %constraint ID number
            switch  constraintName
                case 'dp1'
                    sys.cons{ID} = constraint.dp1(sys,varargin{:}); % new instance of constraint.dp1 class
                case 'dp2'
                    sys.cons{ID} = constraint.dp2(sys,varargin{:}); % new instance of constraint.dp2 class
                case 'cd'
                    sys.cons{ID} = constraint.cd(sys,varargin{:}); % new instance of constraint.cd class
                case 'd'
                    sys.cons{ID} = constraint.d(sys,varargin{:}); % new instance of constraint.d class
                case 'p1'
                    sys.cons{ID} = constraint.p1(sys,varargin{:}); % new instance of constraint.p1 class
                case 'p2'
                    sys.cons{ID} = constraint.p2(sys,varargin{:}); % new instance of constraint.p2 class
                case 'sj'
                    sys.cons{ID} = constraint.sj(sys,varargin{:}); % new instance of constraint.sj class
                case 'rj'
                    sys.cons{ID} = constraint.rj(sys,varargin{:}); % new instance of constraint.rj class
                otherwise
                    error('Constraint not implemented yet.');
            end
        end
        function assembleConstraints(sys) % construct phi matrices 
            sys.addEulerParamConstraints(); % add euler parameter normalization constraints to system
            sys.constructPhiF();   % construct phiF matrix
            sys.constructPhiF_q();   % construct phiF_q matrix
            sys.constructNuF();    % construct RHS of velocity equation
            sys.constructGammaHatF();%construct RHS of acceleration equation
            sys.constructQ();     % construct q (r and p) of bodies in system
            sys.constructPhi_q(); % construct phi_r & phi_p 
        end
        function state = kinematicsAnalysis(sys,timeStart,timeEnd,timeStep) % perform kinematics analysis
            % perform kinematics analysis on the system. System must be
            % fully constrained (nDOF = 0).
            % inputs: measured in seconds.
            %   timeStart : the starting time of simulation, usually 0
            %   timeEnd   : ending time of simulation
            %   timeStep  : step size of time during stimulation (optional)
            %
            % output:
            %   state = cell array of states of system
            %
            % to further clarify, each cell of state is a structure,
            % defining the state of system. For example, at timePoint 1:
            %   state{1}.time
            %   state{1}.q 
            %   state{1}.qdot
            %   state{1}.qddot 
            
            if sys.nDOF >0
                error('For kinematics analysis, the total number of degrees of freedom must be zero.');
            end
            if ~exist('timeStep','var') || isempty(timeStep)
                timeStep = 10^-2; % default time step
            end
            
            timeGrid = timeStart:timeStep:timeEnd;
            
            state = cell(length(timeGrid),1); % preallocate for saved state
            
            % iterate throughout the time grid            
            for iT = 1:length(timeGrid)
                t = timeGrid(iT); % current time step
                sys.setSystemTime(t); % set system time
                
                % solve for positions
                if t ~= timeStart % except at initial conditions
                    sys.positionAnalysis()
                end
               
                % solve for velocities
                sys.velocityAnalysis();
                
                % solve for accelerations
                sys.accelerationAnalysis();
                
                % store system state:
                % return position and orientation information for the
                % current time step
                state{iT}.time = sys.time;
                state{iT}.q = sys.q;
                state{iT}.qdot = sys.qdot;
                state{iT}.qddot = sys.qddot;
                
                disp(['Kinematics analysis completed for t = ' num2str(t) ' sec.']);
            end
        end
        function state = inverseDynamicsAnalysis(sys,timeStart,timeEnd,timeStep) % perform inverse Dynamics Analysis
            % inverse dynamics: specify motion of the mechanical system,
            % and we find the set of forces/torques that were actually
            % applied to the mechanical system to lead to this motion.
            % 
            % We will follow the 3 steps listed on ME751_f2016 slide 8 of
            % lecture 10/10/16
            % 
            % System must be fully constrained (nDOF = 0).
            %
            % inputs: measured in seconds.
            %   timeStart : the starting time of simulation, usually 0
            %   timeEnd   : ending time of simulation
            %   timeStep  : step size of time during stimulation (optional)
            %
            % output:
            %   state = cell array of states of system
            %
            % to further clarify, each cell of state is a structure,
            % defining the state of system. For example, at timePoint 1:
            %   state{1}.time
            %   state{1}.q 
            %   state{1}.qdot
            %   state{1}.qddot
            %   state{1}.rForce
            %   state{1}.rTorque
                        
            if sys.nDOF >0
                error('For inverse dynamics analysis, the total number of degrees of freedom must be zero.');
            end
            if ~exist('timeStep','var') || isempty(timeStep)
                timeStep = 10^-2; % default time step
            end
            
            timeGrid = timeStart:timeStep:timeEnd;
            
            state = cell(length(timeGrid),1); % preallocate for saved state
            
            
            % update constant Matrices needed:
            sys.constructMMatrix(); % update M 
            
            % iterate throughout the time grid:
            for iT = 1:length(timeGrid)
                t = timeGrid(iT); % current time step
                sys.setSystemTime(t); % set system time
                
                %%%%%%%%%%
                % STEP ONE: solve linear system for rddot and pddot
                
                % solve for positions
                if t ~= timeStart % except at initial conditions                  
                    sys.positionAnalysis()
                end
               
                % solve for velocities
                sys.velocityAnalysis();
                
                % solve for accelerations
                sys.accelerationAnalysis();
                
                %%%%%%%%%%
                % STEP TWO: Solve for the Lagrange multipliers (lambda)
                
                sys.constructPhi_q();    % update phi_r & phi_p
                sys.constructPMatrix();  % update P
                sys.constructMMatrix();  % update M
                sys.constructJpMatrix(); % update J^p         
                sys.constructFVector();  % update F
                sys.constructTauHatVector();  % update TauHat
                rddot = sys.qddot(1:3*(sys.nFreeBodies)); % pull accelerations
                pddot = sys.qddot(3*(sys.nFreeBodies)+1:end);
                
                % write the equations of motion for r-p formuation (slide 28 of
                % ME751_f2016 lecture 10/05/16) in a form for Inv. Dyn., like
                % is done on slide 8 of ME751_f2016 lecture 10/10/16.
                % Thus, the Left hand side (LHS) matrix becomes:
                %   [phi_r' zeros(3*nb,nb);
                %    phi_p'        P'      ]
                % this is also seen on slide 12 of ME751_f2016 lecture 10/10/16.
                LHS = [sys.phi_r' zeros(3*sys.nFreeBodies,sys.nFreeBodies);
                       sys.phi_p' sys.P']; % build LHS matrix
                
                % the right hand side (RHS) matrix is takes the form:
                %   -[  M*rddot - F;
                %     J^p*pddot - TauHat]
                RHS = -[sys.M*rddot - sys.F; 
                        sys.Jp*pddot - sys.TauHat]; 
                
                % solve for lagrange multipliers
                lambdaVector = LHS\RHS;
                sys.lambda  = lambdaVector(1:sys.rKDOF*sys.nFreeBodies);
                lambdaP = lambdaVector(sys.rKDOF*sys.nFreeBodies+1:end); % not sure if I need to use this later. 

                                
                %%%%%%%%%%
                % STEP THREE: recover the reaction forces and/or torques that
                % should act on each body so that the system experiences the
                % motion you prescribed.
                sys.calculateReactions(); % in sys.rForce, sys.rTorque
                
                % store system state:
                % return position and orientation information for the
                % current time step
                state{iT}.time = sys.time;
                state{iT}.q = sys.q;
                state{iT}.qdot = sys.qdot;
                state{iT}.qddot = sys.qddot;
                state{iT}.rForce = sys.rForce;
                state{iT}.rTorque = sys.rTorque;

                disp(['Inverse dynamics analysis completed for t = ' num2str(t) ' sec.']);
            end
        end
        function addGravityForces(sys) % add gravity force to body
            bodyID = sys.bodyIDs; % get ID's of bodies that are not grounded (free bodies)
            for i = bodyID
                sys.body{i}.addForceOfGravity(); % add gravity force to body
            end
        end
        function checkInitialConditions(sys) % check that initial conditions satisfy the prescribed constraints
            % (ME751_f2016, slide 40-42, lecture 10/17/16)
            % initial conditions must satisfy constraint equations, or we
            % will start off on the wrong foot and wont be able to get a
            % correct solution. We MUST start with in a healthy (consistent)
            % configuration for our solution to make sense.
            tolerance = 10^-8;
            
            % check that conditions satisfy level zero constraints
            sys.constructPhiF();   % construct phi and phiP matrix
            
            
            
        end
        function plot(sys,varargin) % plots the bodies in the system
            % wrapper to plot functions
            plot.plotSystem(sys,varargin{:});
        end
    end
    methods (Access = public)
        function addEulerParamConstraints(sys) % add euler parameter normalization constraints
            % add p_norm constraint for every body in system (except ground)
            % (remove 1 DOF for every ungrounded body in the system)
            ID = sys.nConstraints+1; %constraint ID number
            for i = 1:sys.nBodies
                if ~sys.body{i}.isGround
                    sys.cons{ID} = constraint.p_norm(sys,sys.body{i}); % new instance of constraint.p_norm class
                    ID = ID+1;
                end
            end
        end
        function constructPhi(sys) % construct phi matrix (kinematic & driving constraints)
            % phi = [(rKDOF) x 1]
            
            sys.phi = zeros(sys.rKDOF,1); % initialize phi for speed
            row = 1;
            for i = 1:(sys.nConstraints-sys.rPDOF) % only kinematic and driving constraints
                row_new = row + sys.cons{i}.rDOF - 1;
                sys.phi(row:row_new)  =  sys.cons{i}.phi; % plug in constraint phi's
                row = row_new + 1;
            end
        end
        function constructPhiP(sys) % construct phiP matrix (euler parameter constraints)
            % phiP = [rPDOF x 1]
            sys.phiP = zeros(sys.rPDOF,1); % initialize phiP for speed
            consNum = sys.nConstraints-sys.rPDOF+1;
            for i = 1:sys.rPDOF % only euler parameter constraints
                sys.phiP(i) = sys.cons{consNum}.phi; % plug in constraint phi's
                consNum = consNum+1;
            end
        end
        function constructPhiF(sys) % construct phiF matrix (full constraint matrix)
            % phiF = [nConstrainedDOF x 1]
            sys.constructPhi();
            sys.constructPhiP();
            sys.phiF = [sys.phi; sys.phiP];
        end
        function constructPhi_q(sys) % construct jacobian of phi
            % phi_q = partial derivative of sys.phi with respect to the
            % generalized coordinates
            % phi_q = [(rKDOF) x nGenCoordinates]
            % phi_q = [phi_r phi_q];
            %
            % For every constraint, pull phi_r and phi_p, then insert into
            % correct location for the jacobian. This location depends on
            % the number of bodies in the system. Grounded bodies not included.
            % Takes the form:
            %   [phi_r(body1) ... phi_r(body_last) | phi_p(body1) ... phi_p(body_last)]
            
            phi_r = zeros(sys.rKDOF,3*sys.nFreeBodies); % initialize phi_r for speed
            phi_p = zeros(sys.rKDOF,4*sys.nFreeBodies); % initialize phi_p for speed
            
            % get ID's of bodies that are not grounded (free bodies)
            bodyID = sys.bodyIDs;
            
            % loop through constraints
            row = 1; %counter
            for i = 1:(sys.nConstraints-sys.rPDOF) % only kinematic and driving constraints
                row_new = row + sys.cons{i}.rDOF - 1;
                if ~sys.cons{i}.bodyi.isGround % if bodyi is not ground....
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
            sys.phi_r = phi_r;
            sys.phi_p = phi_p;
        end
        function constructPhiF_q(sys) % construct jacobian of phiF
            % phiF_q = partial derivative of sys.phiF with respect to the
            % generalized coordinates
            % phiF_q = [nConstrainedDOF x nGenCoordinates]
            % phiF_q = [phiF_r phiF_q];
            %
            % For every constraint, pull phi_r and phi_p, then insert into
            % correct location for the jacobian. This location depends on
            % the number of bodies in the system. Grounded bodies not included. 
            % Takes the form:
            %   [phi_r(body1) ... phi_r(body_last) | phi_p(body1) ... phi_p(body_last)]
            
            phiF_r = zeros(sys.nConstrainedDOF,3*sys.nFreeBodies); % initialize phiF_r for speed
            phiF_p = zeros(sys.nConstrainedDOF,4*sys.nFreeBodies); % initialize phiF_p for speed
            
            % get ID's of bodies that are not grounded (free bodies)
            bodyID = sys.bodyIDs;
            
            % loop through constraints
            row = 1; %counter
            for i = 1:sys.nConstraints
                row_new = row + sys.cons{i}.rDOF - 1;
                if strcmp(class(sys.cons{i}),'constraint.p_norm') % if euler parameter normalization constraint
                    col = find(bodyID == sys.cons{i}.bodyi.ID); % get column position
                    phiF_r(row,(3*col - 2):(3*col)) = sys.cons{i}.phi_r; % place in row
                    phiF_p(row,(4*col - 3):(4*col)) = sys.cons{i}.phi_p;
                elseif ~sys.cons{i}.bodyi.isGround % if bodyi is not ground....
                    coli = find(bodyID == sys.cons{i}.bodyi.ID); % get body i column position
                    phiF_r(row:row_new,(3*coli - 2):(3*coli)) = sys.cons{i}.phi_r(1:sys.cons{i}.rDOF,1:3); % place in row
                    phiF_p(row:row_new,(4*coli - 3):(4*coli)) = sys.cons{i}.phi_p(1:sys.cons{i}.rDOF,1:4);
                    if ~sys.cons{i}.bodyj.isGround % and if body j is not ground
                        colj = find(bodyID == sys.cons{i}.bodyj.ID); % get body j column position
                        phiF_r(row:row_new,(3*colj - 2):(3*colj)) = sys.cons{i}.phi_r(1:sys.cons{i}.rDOF,4:6); % place in row
                        phiF_p(row:row_new,(4*colj - 3):(4*colj)) = sys.cons{i}.phi_p(1:sys.cons{i}.rDOF,5:8);
                    end
                elseif sys.cons{i}.bodyi.isGround % if bodyi is ground....
                    colj = find(bodyID == sys.cons{i}.bodyj.ID); % get body j column position
                    phiF_r(row:row_new,(3*colj - 2):(3*colj)) = sys.cons{i}.phi_r(1:sys.cons{i}.rDOF,1:3); % place in row
                    phiF_p(row:row_new,(4*colj - 3):(4*colj)) = sys.cons{i}.phi_p(1:sys.cons{i}.rDOF,1:4);
                end
                row = row_new + 1; % increment row counter
            end
            
            % create phiF_q
            sys.phiF_q = [phiF_r phiF_p];
        end
        function constructNuF(sys) % construct nuF, RHS Of velocity equation
            % nuF = [nConstrainedDOF x 1]
            sys.nuF = zeros(sys.nConstrainedDOF,1); % initialize nuF for speed
            row = 1;
            for i = 1:sys.nConstraints
                row_new = row + sys.cons{i}.rDOF - 1;
                sys.nuF(row:row_new)  =  sys.cons{i}.nu; % plug in constraint nu's
                row = row_new + 1;
            end
        end
        function constructGammaHatF(sys) % construct gammaHatF, RHS Of acceleration equation, in r-p formulation
            % gammaHatF = [nConstrainedDOF x 1]
            sys.gammaHatF = zeros(sys.nConstrainedDOF,1); % initialize nuF for speed
            row = 1;
            for i = 1:sys.nConstraints
                row_new = row + sys.cons{i}.rDOF - 1;
                sys.gammaHatF(row:row_new)  =  sys.cons{i}.gammaHat; % plug in constraint gammaHats's
                row = row_new + 1;
            end
        end
        function constructQ(sys) % construct q (r and p) of bodies in system
            % get ID's of bodies that are not grounded
            bodyID = sys.bodyIDs;
            
            % initialize for speed
            r = zeros(3*sys.nFreeBodies,1);
            p = zeros(4*sys.nFreeBodies,1);
            
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
            
            if ~exist('tolerance','var') || isempty(tolerance)
                tolerance = 1e-9;  % default tolerance
            end
            if ~exist('maxIterations','var') || isempty(maxIterations)
                maxIterations = 50; % default maximum iterations
            end
             
            guessQ = sys.q; % initial guess
            
            for i = 1:maxIterations % iterate
                sys.constructPhiF();   % get constraint matrix
                sys.constructPhiF_q(); % get jacobian
                correction = sys.phiF_q\sys.phiF; % correction 
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
            sys.constructNuF();
            
            % find jacobian of phi
            sys.constructPhiF_q();
            
            % solve for velocities (qdot)
            qdot = sys.phiF_q\sys.nuF;
            
            % update bodies
            sys.updateSystem([],qdot,[]);
        end
        function accelerationAnalysis(sys) % solve for accelerations of the system
            % compute acceleration and 2nd derivative of euler parameters
            % for a given time (t)
            
            % calculate RHS of acceleration equation
            sys.constructGammaHatF();
            
            % already found jacobian of phi
            
            % solve for accelerations (qddot)
            qddot = sys.phiF_q\sys.gammaHatF;
            
            % update bodies
            sys.updateSystem([],[],qddot);
        end
        function constructPMatrix(sys) % construct the P matrix used in EOM
            % from ME751_f2016 slide 27 of lecture 10/05/16
            sys.P = zeros(sys.nFreeBodies,4*sys.nFreeBodies); % preallocate for speed
            bodyID = sys.bodyIDs; % get ID's of bodies that are not grounded (free bodies)
            for i = 1:sys.nFreeBodies
               sys.P(i,(4*i-3):4*i) = sys.body{bodyID(i)}.p'; % plug in p vector transposed
            end  
        end
        function constructMMatrix(sys) % construct the P matrix used in EOM
            % from ME751_f2016 slide 27 of lecture 10/05/16
            sys.M = zeros(3*sys.nFreeBodies,3*sys.nFreeBodies); % preallocate for speed
            bodyID = sys.bodyIDs; % get ID's of bodies that are not grounded (free bodies)
            for i = 1:sys.nFreeBodies
               sys.M((3*i-2):3*i,(3*i-2):3*i) = sys.body{bodyID(i)}.m*eye(3); % plug in mass matrices
            end  
        end
        function constructJpMatrix(sys) % construct the J^P matrix used in EOM
            % from ME751_f2016 slide 27 of lecture 10/05/16
            sys.Jp = zeros(4*sys.nFreeBodies,4*sys.nFreeBodies); % preallocate for speed
            bodyID = sys.bodyIDs; % get ID's of bodies that are not grounded (free bodies)
            for i = 1:sys.nFreeBodies
                 J = sys.body{bodyID(i)}.Jbar; % pull inertia tensor
                 G = utility.Gmatrix(sys.body{bodyID(i)}.p); % calc G matrix
                 sys.Jp((4*i-3):4*i,(4*i-3):4*i) = 4*G'*J*G; % plug in J^p terms
            end  
        end
        function constructFVector(sys) % construct F vector used in EOM
            % from ME751_f2016 slide 26 from lecture 10/05/16:
            % F = F^m + F^a
            % from ME751_f2016 slide 29 from lecture 10/03/16:
            %   F^m : mass distributed forces (usually = mg)
            % from ME751_f2016 slide 27 from lecture 10/03/16:
            %   F^a : active forces
            
            % from ME751_f2016 slide 27 of lecture 10/05/16:
            sys.F = zeros(3*sys.nFreeBodies,1); % preallocate for speed
            
            bodyID = sys.bodyIDs; % get ID's of bodies that are not grounded (free bodies)
            for i = 1:sys.nFreeBodies % iterate through every free body
                F = zeros(3,1);
                for j = 1:sys.body{bodyID(i)}.nForces
                    F = F + sys.body{bodyID(i)}.forces{j}.force; % sum active and gravitational forces acting on that body
                end
                sys.F((3*i-2):3*i) = F;
            end
        end
        function constructTauHatVector(sys) % construct TauHat vector used in EOM
            % from ME751_f2016 slide 27 of lecture 10/05/16:
            sys.TauHat = zeros(4*sys.nFreeBodies,1); % preallocate for speed
            
            bodyID = sys.bodyIDs; % get ID's of bodies that are not grounded (free bodies)
            for i = 1:sys.nFreeBodies % iterate through every free body
                nBar = zeros(3,1);
                for j = 1:sys.body{bodyID(i)}.nForces
                    nBar = nBar + sys.body{bodyID(i)}.forces{j}.torque; % sum active and gravitational torques acting on that body
                end
                
                % TauHat constructed using formula from ME751_f2016 slide 26 from lecture 10/05/16
                J = sys.body{bodyID(i)}.Jbar; % pull inertia tensor
                G = utility.Gmatrix(sys.body{bodyID(i)}.p); % calc G matrix
                Gdot = utility.Gmatrix(sys.body{bodyID(i)}.pdot); % calc Gdot matrix
                P = sys.body{bodyID(i)}.p; % pull euler parameters
                TauHat = 2*G'*nBar + 8*Gdot'*J*Gdot*P;
                sys.TauHat((4*i-3):4*i) = TauHat;
            end
        end
        function calculateReactions(sys) % calculate reaction forces and torques
            % From ME751_f2016 slide 8 of lecture 10/10/16
            sys.rForce = cell(sys.nFreeBodies,length(sys.lambda)); % preallocate for speed
            sys.rTorque = cell(sys.nFreeBodies,length(sys.lambda)); % preallocate for speed
            bodyID = sys.bodyIDs; % get ID's of bodies that are not grounded (free bodies)
            for i = 1:sys.nFreeBodies
                % Reaction Forces
                phi_r_i = sys.phi_r(:,(3*i-2):3*i); % pull phi_r related to body i
                
                % total reaction forces for this body given by:
                % rForceTotal on body i = -phi_r_i'*sys.lambda
                % However, we want to find reactions forces for each lambda
                % we are interested in. such that:
                % rForceTotal on body i = (rForce on body i due to lamba1) + (rForce on body i due to lamba2) + ...

                for j = 1:length(sys.lambda) % look at each row of phi_r_i
                    sys.rForce{i,j} = -phi_r_i(j,:)'*sys.lambda(j); % reaction force j on body i
                end
                
                
                % Reaction Torques
                phi_p_i = sys.phi_p(:,(4*i-3):4*i); % pull phi_p related to body i
                
                % total reaction torques for this body as:
                % rTorqueTotal for body i = -phi_p_i'*sys.lambda
                % However, we want to find reactions torques for each lambda
                % we are interested in. such that:
                % rTorqueTotal on body i = (rTorque on body i due to lamba1) + (rTorque on body i due to lamba2) + ...
                % but rTorque is [4x1], need to convert to r-omega formulation
                G = utility.Gmatrix(sys.body{bodyID(i)}.p); % calc G matrix
                for j = 1:length(sys.lambda) % look at each row of phi_p_i
                    sys.rTorque{i,j} = 0.5*G*-phi_p_i(j,:)'*sys.lambda(j);
                end
            end
        end
        function updateSystem(sys,q,qdot,qddot) %update r,p,and derivatives for each body in system
            %%% update q (r and p)
            if ~isempty(q)
                sys.q = q;
                r = q(1:3*(sys.nFreeBodies));
                p = q(3*(sys.nFreeBodies)+1:end);
                
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
                rdot = qdot(1:3*(sys.nFreeBodies));
                pdot = qdot(3*(sys.nFreeBodies)+1:end);
                
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
                rddot = qddot(1:3*(sys.nFreeBodies));
                pddot = qddot(3*(sys.nFreeBodies)+1:end);
                
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
            % set system time. Bodies, constraints, and forces should be
            % able to access this time.
            sys.time = t;
        end
    end
    methods % methods block with no attributes
        function bodyIDs = get.bodyIDs(sys) % calculate ID numbers of free bodies in the system
            bodyIDs = zeros(sys.nFreeBodies,1); %initialize size
            j = 1; % counter 
            for i = 1:sys.nBodies
                if ~sys.body{i}.isGround
                    bodyIDs(j) = sys.body{i}.ID; % pull body ID
                    j = j+1; % increment counter
                end
            end
        end
        function nBodies = get.nBodies(sys) % calculate total number of bodies in system
            nBodies = length(sys.body);
        end
        function nGroundBodies = get.nGroundBodies(sys) % calculate number of grounded bodies in system
            nGroundBodies = 0;
            if sys.nBodies>0
                for i = 1:sys.nBodies
                    if sys.body{i}.isGround
                        nGroundBodies = nGroundBodies + 1;
                    end
                end
            end
        end
        function nFreeBodies = get.nFreeBodies(sys) % calculate number of free bodies in system
            nFreeBodies = sys.nBodies-sys.nGroundBodies;
        end
        function nGenCoordinates = get.nGenCoordinates(sys) % calculate number of generalized coordinates in system
            % 7 DOF available for each body
            nGenCoordinates = 7*(sys.nFreeBodies);
        end
        function nConstraints = get.nConstraints(sys) % calculate number of constraints in system
            nConstraints = length(sys.cons);
        end
        function nConstrainedDOF = get.nConstrainedDOF(sys) % calculate number of degrees of freedom that have been constrained
            nConstrainedDOF = sys.nGenCoordinates - sys.nDOF;
        end
        function rKDOF = get.rKDOF(sys) % number of degrees of freedom removed by kinematic and driving constraints
            rKDOF = 0;
            if sys.nConstraints>0
                for i = 1:sys.nConstraints
                    if ~strcmp(class(sys.cons{i}),'constraint.p_norm') % dont include euler parameter normalization constraints
                        rKDOF = rKDOF + sys.cons{i}.rDOF;
                    end
                end
            end
        end
        function rPDOF = get.rPDOF(sys) % number of degrees of freedom removed by euler parameter normalization constraints
            rPDOF = 0;
            if sys.nConstraints>0
                for i = 1:sys.nConstraints
                    if strcmp(class(sys.cons{i}),'constraint.p_norm') % only include euler parameter normalization constraints
                        rPDOF = rPDOF + sys.cons{i}.rDOF;
                    end
                end
            end
        end
        function nDOF = get.nDOF(sys) % count number of Degrees of Freedom
            nDOF = sys.nGenCoordinates - (sys.rKDOF + sys.rPDOF);
        end
        function nForces = get.nForces(sys) % count number of forces/torques applied to bodies
            bodyID = sys.bodyIDs; % get ID's of bodies that are not grounded (free bodies)
            nForces = 0;
            for i = 1:sys.nFreeBodies
                nForces = nForces + sys.body{bodyID(i)}.nForces;
            end
        end
    end
end

