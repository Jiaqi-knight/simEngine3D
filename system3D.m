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
        r_global = [0;0;0];    % Global Reference Frame
        p_global = [1;0;0;0];  % Global Euler Parameters
        body;           % collection of bodies in the system, inlcuding ground bodies
        cons;           % collection of constraints in the system
        q;              % q = [r;p] vector, generalized coordinates of ungrounded (free) bodies
        qdot;           % time derivaite of q vector
        qddot;          % 2nd time derivaite of q vector
        phi;            % set of kinematic and driving constraints
        phi_r;          % jacobian of phi, with respect to r only
        phi_p;          % jacobian of phi, with respect to p only
        phiP;           % set of Euler Parameter Normalization constraints
        phiF;           % full constraint matrix phiF = [phi; phiP]
        phiF_q;         % jacobian of phiF
        nu;             % RHS of velocity equation
        nuP;            % RHS of velocity equation, for euler parameters
        nuF;            % RHS of Full velocity equation
        gammaHat;       % RHS of acceleration equation, in r-p formulation
        gammaHatP;      % RHS of acceleration equation for euler parameter constraints
        gammaHatF;      % RHS of full acceleration equation
        M;              % M matrix for Equations of Motion (masses)
        P;              % P matrix for Equations of Motion (euler parameters)
        Jp;             % J^p matrix for Equations of Motion (inertias)
        F;              % F vector used in Equations of Motion (forces)
        TauHat;         % TauHat vector used in Equations of Motion (torques)
        lambda;         % lagrange multipliers of the system, by constraints
        lambdaP;        % lagrange multipliers of the system, from euler parameter constraints
        rForce;         % reaction forces for each body in system
        rTorque;        % reaction torques for each body in system
        g;              % vector of acceleration due to gravity
        time;           % current time within the system
    end
    
    properties (Dependent)
        r;                  % vector of positions of the bodies
        rdot;               % vector of velocities of the bodies
        rddot;              % vector of accelerations of the bodies
        p;                  % vector of euler parameters of the bodies (rotation)
        pdot;               % vector of derivative of euler parameters of the bodies (rotational velocity)
        pddot;              % vector of 2nd derivative of euler parameters of the bodies (rotational acceleration)
        bodyIDs;            % ID numbers of ungrounded bodies in the system
        nPoints;            % number of points identified in the system
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
        function assembleConstraints(sys) % add euler parameter normalization constraints and construct phi matrices 
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
                state{iT} = sys.getSystemState();
                
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
                sys.lambdaP = lambdaVector(sys.rKDOF*sys.nFreeBodies+1:end); % not sure if I need to use this later. 

                                
                %%%%%%%%%%
                % STEP THREE: recover the reaction forces and/or torques that
                % should act on each body so that the system experiences the
                % motion you prescribed.
                sys.calculateReactions(); % in sys.rForce, sys.rTorque
                
                % store system state:
                state{iT} = sys.getSystemState();
                
                disp(['Inverse dynamics analysis completed for t = ' num2str(t) ' sec.']);
            end
        end
        function state = dynamicsAnalysis(sys,timeStart,timeEnd,timeStep,newtonMethod) % perform dynamics analysis
            % perform dynamics analysis on the system. System must be not be
            % fully constrained (nDOF > 0).
            %
            % dynamics analysis peformed using the 6 steps outline on slide 43 of
            % lecture 10/17/16, ME751_f2016
            %
            % inputs: measured in seconds.
            %   timeStart : the starting time of simulation, usually 0
            %   timeEnd   : ending time of simulation
            %   timeStep  : step size of time during stimulation (optional)
            %   newtonMethod: {'NR','MN','QN'}
            %
            % output:
            %   state = cell array of states of system
            
            if sys.nDOF == 0
                warning('For  dynamics analysis, the total number of degrees of freedom should be greater than zero.');
            end
            if ~exist('timeStep','var') || isempty(timeStep)
                timeStep = 10^-2; % default time step
            end
            if ~exist('newtonMethod','var') || isempty(newtonMethod)
                newtonMethod = 'QN'; % default newton method is Quasi-Newton
            end
            
            timeGrid = timeStart:timeStep:timeEnd; % establish time grid
            state = cell(length(timeGrid),1); % preallocate for saved state
            
            
            %%%%%%%%%%%
            % check initial conditions (ME751_f2016, slide 41, lecture 10/17/16)
            iT = 1; % at timeStart
            sys.setSystemTime(timeGrid(iT)); % set system time
            sys.checkInitialConditions();
            sys.findInitialAccelerations();
            state{1} = sys.getSystemState(); % store initial system state: (iT = 1)
            disp(['Dynamics analysis completed for t = ' num2str(sys.time) ' sec.'])    
            
            % store violation of the velocity constraints (optional)
            state{1}.velocityViolation = sys.violationVelocityConstraints();
            
            %%%%%%%%%%% 
            % using initial conditions,
            % find first time step using Quasi-Newton method and BDF of order 1
            iT = 2; % first integration step
            sys.setSystemTime(timeGrid(iT)); % set system time
            sys.BDF_dynamicsAnalysis(1,newtonMethod,timeStep,state{1}); % ORDER 1 BDF, Quasi-Newton method, send initial state
            state{2} = sys.getSystemState();
            disp(['Dynamics analysis completed for t = ' num2str(sys.time) ' sec.'])
            
            % store violation of the velocity constraints (optional)
            state{2}.velocityViolation = sys.violationVelocityConstraints();
            
            % iterate throughout the time grid:
            for iT = 3:length(timeGrid)
                t = timeGrid(iT); % current time step
                sys.setSystemTime(t); % set system time
                sys.BDF_dynamicsAnalysis(2,newtonMethod,timeStep,state{iT-1},state{iT-2}); % ORDER 2 BDF, Quasi-Newton method, send past two states
                state{iT} = sys.getSystemState();
                
                disp(['Dynamics analysis completed for t = ' num2str(t) ' sec.']);
                
                % store violation of the velocity constraints (optional)
                state{iT}.velocityViolation = sys.violationVelocityConstraints();
            end
        end %dynamicsAnalysis
        function addGravityForces(sys) % add gravity force to body
            bodyID = sys.bodyIDs; % get ID's of bodies that are not grounded (free bodies)
            for i = bodyID
                sys.body{i}.addForceOfGravity(); % add gravity force to body
            end
        end
        function checkInitialConditions(sys,tolerance) % check that initial conditions satisfy the prescribed constraints
            % (ME751_f2016, slide 40-42, lecture 10/17/16)
            % initial conditions must satisfy constraint equations, or we
            % will start off on the wrong foot and wont be able to get a
            % correct solution. We MUST start with in a healthy (consistent)
            % configuration for our solution to make sense.
            
            if ~exist('tolerance','var') || isempty(tolerance)
               tolerance = 10^-5; % default tolerance
            end
            
            %%%%%%%%%%
            % check that conditions satisfy level zero constraints
            sys.constructPhiF();   % construct phi and phiP matrix
            if any(abs(sys.phiF) > tolerance)
                error(' Initial position conditions are not consistent.')
            end
            
            
            %%%%%%%%%%
            % check that conditions satisfy level one constraints
            sys.constructNuF();
            vel = sys.phi_r*sys.rdot+sys.phi_p*sys.pdot;
            if any(abs(vel-sys.nu) > tolerance)
                error(' Initial velocity conditions are not consistent.')
            end
            
            sys.constructPMatrix();
            if any(abs(sys.P*sys.pdot) > tolerance)
                error(' Initial euler parameter conditions are not consistent.')
            end
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
        function findInitialAccelerations(sys) % find initial conditions for acceleration
            % from ME751_f2016, slide 42 from lecture 10/17/2016
            
            %%%%%%
            % fill out LHS matrix
            sys.constructMMatrix();  % update M
            sys.constructJpMatrix(); % update J^p
            sys.constructPhiF();     % update phi & phiP
            sys.constructPhi_q();    % update phi_r & phi_p
            sys.constructPMatrix();  % update P
            Z1 = zeros(4*sys.nFreeBodies,3*sys.nFreeBodies); 
            Z2 = zeros(sys.nFreeBodies,3*sys.nFreeBodies); 
            Z3 = zeros(sys.nFreeBodies);
            Z4 = zeros(sys.rKDOF,sys.nFreeBodies);
            Z5 = zeros(sys.rKDOF);
            LHS = [    sys.M         Z1'     Z2'  sys.phi_r';
                          Z1     sys.Jp   sys.P'  sys.phi_p';
                          Z2      sys.P      Z3          Z4';
                   sys.phi_r  sys.phi_p      Z4          Z5];
            
            %%%%%%
            % fill out RHS matrix
            sys.constructFVector();       % update F
            sys.constructTauHatVector();  % update TauHat
            sys.constructGammaHatF();     % gammaHat and gammaHatF
            RHS = [sys.F; sys.TauHat; sys.gammaHatP; sys.gammaHat];
            
            %%%%%%
            % solve for initial accelerations
            accel = LHS\RHS;
            sys.updateSystem([],[],accel(1:7*sys.nFreeBodies)); % update qddot
            sys.lambdaP = accel(7*sys.nFreeBodies+1:7*sys.nFreeBodies+sys.rPDOF); % update lambdaP
            sys.lambda = accel(7*sys.nFreeBodies+sys.rPDOF+1:end); % update lambda
            sys.calculateReactions(); % derive sys.rForce, sys.rTorque
        end
        function QN_BDF1_check(sys,h,state1) % Quasi-Newton using BDF order 1
            
            %%%%%%%%%%%%%%%%%%
             % useful constants
            sys.constructMMatrix();  % update M
            Z1 = zeros(4*sys.nFreeBodies,3*sys.nFreeBodies); 
            Z2 = zeros(sys.nFreeBodies,3*sys.nFreeBodies); 
            Z3 = zeros(sys.nFreeBodies);
            Z4 = zeros(sys.rKDOF,sys.nFreeBodies);
            Z5 = zeros(sys.rKDOF);
            
            %%%%%%%%%%%%%%%%%%
            % Steps:
            % Find PSI
            % find Ystar
            % find G(Ystar)
            % find Ystar + delta
            % find G( Ystar + delta )
            % find (G(Ystar+delta) - g(Ystar))/delta
            % compare with original PSI
            
            %%%%%%%%%%%%%%%%%%
            % Find PSI
            % BDF constants
            beta0  =  1;
            alpha1 = -1;
            Crdot   = -alpha1*state1.rdot;
            Cpdot   = -alpha1*state1.pdot;
            Cr      = -alpha1*state1.r + beta0*h*Crdot;
            Cp      = -alpha1*state1.p + beta0*h*Cpdot;
            % initial guesses
            rddot   = state1.rddot;
            pddot   = state1.pddot;
            lambdaP = state1.lambdaP;
            lambda  = state1.lambda;
            %compute position and velocity using BDF and most recent accelerations
            r     = Cr + beta0^2*h^2*rddot;
            p     = Cp + beta0^2*h^2*pddot;
            rdot  = Crdot + beta0*h*rddot;
            pdot  = Cpdot + beta0*h*pddot;
            sys.updateSystem([r; p],[rdot; pdot],[]);
            %compute PSI (jacobian)
            sys.constructPhiF();            % update phi & phiP
            sys.constructJpMatrix();        % update J^p
            sys.constructPhi_q();           % update phi_r & phi_p
            PSI = [    sys.M         Z1'     Z2'  sys.phi_r';
                          Z1     sys.Jp   sys.P'  sys.phi_p';
                          Z2      sys.P      Z3          Z4';
                   sys.phi_r  sys.phi_p      Z4          Z5];
               
               
            %%%%%%%%%%%%%%%%%%
            % Find OG PSI
            r     = state1.r;
            p     = state1.p;
            rdot  = state1.rdot;
            pdot  = state1.pdot;
            rddot = state1.rddot;
            pddot = state1.pddot;
            sys.updateSystem([r; p],[rdot; pdot],[]);
            %compute PSI (jacobian)
            sys.constructPhiF();            % update phi & phiP
            sys.constructJpMatrix();        % update J^p
            sys.constructPhi_q();           % update phi_r & phi_p
            PSIog = [    sys.M         Z1'     Z2'  sys.phi_r';
                          Z1     sys.Jp   sys.P'  sys.phi_p';
                          Z2      sys.P      Z3          Z4';
                   sys.phi_r  sys.phi_p      Z4          Z5];
            
            %%%%%%%%%%%%%%%%%%
            % find Ystar
            r     = state1.r;
            p     = state1.p;
            rdot  = state1.rdot;
            pdot  = state1.pdot;
            rddot = state1.rddot;
            pddot = state1.pddot;
            sys.updateSystem([r; p],[rdot; pdot],[]);
            
            
            %%%%%%%%%%%%%%%%%%
            % find G(Ystar)
            sys.constructPhiF();            % update phi & phiP
            sys.constructJpMatrix();        % update J^p
            sys.constructPhi_q();           % update phi_r & phi_p
            sys.constructPMatrix();         % update P
            sys.constructFVector();         % update F
            sys.constructTauHatVector();    % update TauHat
            G = [sys.M*rddot + sys.phi_r'*lambda + -sys.F;
                sys.Jp*pddot + sys.phi_p'*lambda + sys.P'*lambdaP + -sys.TauHat;
                1/(beta0^2*h^2)*sys.phiP;
                1/(beta0^2*h^2)*sys.phi];
            
            %%%%%%%%%%%%%%%%%%
            % find Ystar + delta
            delta = 10^-3;
            r     = state1.r;
            p     = state1.p;
            rdot  = state1.rdot;
            pdot  = state1.pdot;
            rddot = state1.rddot;
            pddot = state1.pddot;
            lambda(5) = lambda(5)+delta;
            
            sys.updateSystem([r; p],[rdot; pdot],[]);
            
            
            %%%%%%%%%%%%%%%%%%
            % find G( Ystar + delta )
            sys.constructPhiF();            % update phi & phiP
            sys.constructJpMatrix();        % update J^p
            sys.constructPhi_q();           % update phi_r & phi_p
            sys.constructPMatrix();         % update P
            sys.constructFVector();         % update F
            sys.constructTauHatVector();    % update TauHat
            Gdelta = [sys.M*rddot + sys.phi_r'*lambda + -sys.F;
                sys.Jp*pddot + sys.phi_p'*lambda + sys.P'*lambdaP + -sys.TauHat;
                1/(beta0^2*h^2)*sys.phiP;
                1/(beta0^2*h^2)*sys.phi];
            
            
            %%%%%%%%%%%%%%%%%%
            % find (G(Ystar+delta) - g(Ystar))/delta
            FD = (Gdelta-G)./delta;
            
            
            %%%%%%%%%%%%%%%%%%
            % compare with original PSI
            PSIog
            PSI
            FD
            r
            
            
            
            
            
            
           

        end
        function BDF_dynamicsAnalysis(sys,orderNum,newtonMethod,h,varargin) % BDF Solution of the Dynamics Analysis Problem
            % Backwards Difference Formula, using a newton method
            % see ME751_f2016, slides 18-22, 10/14/16
            % also see ME751_f2016, slide 43, 10/17/16
            % 
            % Note: This function only solves one intergration step, at the current system time
            %
            % INPUTS:
            %     orderNum : order of BDF method {1,2,3,4,5,6}
            % newtonMethod : type of newton method used in numerical method {'NR','MN','QN'}
            %                NR : Newton-Raphson Method
            %                MN : Modified Newton Method
            %                QN : Quasi-Newton Method
            %            h : step size, in seconds
            %      varargin: a collection of states of the system, in order from
            %                most recent to least recent...
            %                varargin{1} : state1 : structure of system values at previous time step, t(n-1)
            %                varargin{2} : state2 : structure of system values at previous time step, t(n-2)
            
            % solver parameters
            maxIterations = 50; % 50;
            tolerance = 1e-2; %1e-8;
            
            switch orderNum % BDF order number
                case 1 % Order 1: (this is also known as the Backwards Euler Method)
                    % BDF constants
                    beta0  =  1;
                    alpha1 = -1;
                    state1 = varargin{1}; % system state from previous time step, t(n-1)
                    % Constant terms from slide 30-31, lecture 10/14/16
                    Crdot   = -alpha1*state1.rdot;
                    Cpdot   = -alpha1*state1.pdot;
                    Cr      = -alpha1*state1.r + beta0*h*Crdot;
                    Cp      = -alpha1*state1.p + beta0*h*Cpdot;
                
                case 2 % order 2
                    % BDF constants
                    beta0  =  2/3;
                    alpha1 = -4/3;
                    alpha2 =  1/3;
                    state1 = varargin{1}; % system state from previous time step, t(n-1)
                    state2 = varargin{2}; % system state from previous time step, t(n-2)
                    % Constant terms from slide 30-31, lecture 10/14/16
                    Crdot   = -alpha1*state1.rdot + -alpha2*state2.rdot;
                    Cpdot   = -alpha1*state1.pdot + -alpha2*state2.pdot;
                    Cr      = -alpha1*state1.r + -alpha2*state2.r + beta0*h*Crdot;
                    Cp      = -alpha1*state1.p + -alpha2*state2.p + beta0*h*Cpdot;
                otherwise
                    error(['BDF of Order ' orderNum ' not implemented yet.']);
            end
                
            
            %%%%%%%%%%%
            % STEP 0 : prime new time step
            % (follow Iteration steps from slide 43, 10/17/16)
           
            % initial guesses
            rddot   = state1.rddot;
            pddot   = state1.pddot;
            lambdaP = state1.lambdaP;
            lambda  = state1.lambda;
            

            % useful constants
            sys.constructMMatrix();  % update M
            Z1 = zeros(4*sys.nFreeBodies,3*sys.nFreeBodies); 
            Z2 = zeros(sys.nFreeBodies,3*sys.nFreeBodies); 
            Z3 = zeros(sys.nFreeBodies);
            Z4 = zeros(sys.rKDOF,sys.nFreeBodies);
            Z5 = zeros(sys.rKDOF);
            
            % Newton Iteration Method
            nu = 1; % iteration counter
            while nu < maxIterations
                %%%%%%%%%%%
                % STEP 1 : compute position and velocity using BDF and most recent accelerations
                r     = Cr + beta0^2*h^2*rddot;
                p     = Cp + beta0^2*h^2*pddot;
                rdot  = Crdot + beta0*h*rddot;
                pdot  = Cpdot + beta0*h*pddot;
                sys.updateSystem([r; p],[rdot; pdot],[rddot; pddot]);
                
                
                %%%%%%%%%%%
                % STEP 2 : compute the residual (g) in nonlinear system, which
                % is a function of most recent accelerations
                % (see slide 26, 10/17/16)
                sys.constructPhiF();            % update phi & phiP
                sys.constructJpMatrix();        % update J^p
                sys.constructPhi_q();           % update phi_r & phi_p
                sys.constructPMatrix();         % update P
                sys.constructFVector();         % update F
                sys.constructTauHatVector();    % update TauHat
                
                % descretized equations of motion
                g = [sys.M*rddot + sys.phi_r'*lambda + -sys.F;
                    sys.Jp*pddot + sys.phi_p'*lambda + sys.P'*lambdaP + -sys.TauHat;
                    1/(beta0^2*h^2)*sys.phiP;
                    1/(beta0^2*h^2)*sys.phi];

                %%%%%%%%%%%
                % STEP 3 : Solve nonlinear system to get correction deltaZ. 
                switch newtonMethod
                    % (For overview of newton methods available, see slide 26, 10/17/16)
                    case 'QN' 
                        % Use a QUASI-NEWTON iteration matrix (psi).
                        % (see slide 31, 10/17/16)
                        if nu == 1 % only compute once, on the first iteration
                            psi = [    sys.M         Z1'     Z2'  sys.phi_r';
                                          Z1     sys.Jp   sys.P'  sys.phi_p';
                                          Z2      sys.P      Z3          Z4';
                                   sys.phi_r  sys.phi_p      Z4          Z5];
                        end
                    case 'NR'
                        % Use a full NEWTON-RAPHSON iteration matrix (psi).
                        % (see slide 30, 10/17/16)
                        psi = sys.constructPsiMatrix(h,beta0,lambda,lambdaP);
                    case 'MN'
                        % Use a MODIFIED-NEWTON iteration matrix (psi).
                        % (see slide 26, 10/17/16)
                        if nu == 1 % only compute once, on the first iteration
                            psi = sys.constructPsiMatrix(h,beta0,lambda,lambdaP);
                        end
                    otherwise
                        error('Specified Newton Method not implemented yet.');
                end
                
                % calculate correction
                deltaZ = psi\-g; 
                %deltaZ(:,nu) = psi\-g; % for recording instances of deltaZ
                
                %%%%%%%%%%%
                % STEP 4 : Improve the quality of the approximate solution
                Z_old = [rddot; pddot; lambdaP; lambda];
                Z = Z_old + deltaZ;
                %Z = Z_old + deltaZ(:,nu); % for recording instances of deltaZ
                rddot = Z(1:length(rddot));
                pddot = Z(length(rddot)+1:length(rddot)+length(pddot));
                lambdaP = Z(length(rddot)+length(pddot)+1:length(rddot)+length(pddot)+length(lambdaP));
                lambda  = Z(length(rddot)+length(pddot)+length(lambdaP)+1:end);
                
                %%%%%%%%%%%
                % STEP 5 : if the norm of the correction is small enough,
                % go to step 6. Otherwise
                nu = nu + 1;
                
                if norm(deltaZ) < tolerance
                    break;
                end
            end
            if nu >= maxIterations % notify of failed convergence
                warning(['Did not converge in Quasi-Newton iteration, BDF1. Time = ' num2str(sys.time)]);
                disp('but we are going on, baby!');
            end
            
            %%%%%%%%%%%
            % STEP 6 : Accep the acclerations and lambdas computed in step
            % 4 as the solutions. Using the acclerations, find position and
            % velocities.
            r     = Cr + beta0^2*h^2*rddot;
            p     = Cp + beta0^2*h^2*pddot;
            rdot  = Crdot + beta0*h*rddot;
            pdot  = Cpdot + beta0*h*pddot;
            sys.updateSystem([r; p],[rdot; pdot],[rddot; pddot]); % final system state
            sys.lambdaP = lambdaP;
            sys.lambda = lambda;
            sys.calculateReactions(); % derive sys.rForce, sys.rTorque
            
            % At this point you just finished one integration step. Life is
            % good. Now repeat this for another time step
        end
        function violation = violationVelocityConstraints(sys)
            % evaluate violation of the velocity constraints
            % see slide 41 for lecture 10/17/16
            sys.constructNuF();
            vel = sys.phi_r*sys.rdot+sys.phi_p*sys.pdot;
            violation = vel-sys.nu;
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
                
                rowr = 1;
                rowp = 1;
                for i = sys.bodyIDs
                    sys.body{i}.r = sys.r(rowr:rowr+2);
                    sys.body{i}.p = sys.p(rowp:rowp+3);
                    rowr = rowr + 3;
                    rowp = rowp + 4;
                end
            end
            
            %%% update qdot (rdot and pdot)
            if ~isempty(qdot)
                sys.qdot = qdot;
                
                rowr = 1;
                rowp = 1;
                for i = sys.bodyIDs
                    sys.body{i}.rdot = sys.rdot(rowr:rowr+2);
                    sys.body{i}.pdot = sys.pdot(rowp:rowp+3);
                    rowr = rowr + 3;
                    rowp = rowp + 4;
                end
            end
            
            %%% update qddot (rddot and pddot)
            if ~isempty(qddot)
                sys.qddot = qddot;
                
                rowr = 1;
                rowp = 1;
                for i = sys.bodyIDs
                    sys.body{i}.rddot = sys.rddot(rowr:rowr+2);
                    sys.body{i}.pddot = sys.pddot(rowp:rowp+3);
                    rowr = rowr + 3;
                    rowp = rowp + 4;
                end
            end
        end
        function state = getSystemState(sys) % return structure of relevant state variables
            state.time = sys.time;
            state.r = sys.r;
            state.p = sys.p;
            state.rdot = sys.rdot;
            state.pdot = sys.pdot;
            state.rddot = sys.rddot;
            state.pddot = sys.pddot;
            state.lambda = sys.lambda;
            state.lambdaP = sys.lambdaP;
            state.rForce = sys.rForce;
            state.rTorque = sys.rTorque;
        end
        function setSystemTime(sys,t) % set time for entire system
            % set system time. Bodies, constraints, and forces should be
            % able to access this time.
            sys.time = t;
        end
    end
    methods (Access = private) % methods for constructing system terms 
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
        function constructNu(sys) % construct nu vector (kinematic & driving constraints)
            % nu = [(rKDOF) x 1]
            sys.nu = zeros(sys.rKDOF,1); % initialize nu for speed
            row = 1;
            for i = 1:(sys.nConstraints-sys.rPDOF) % only kinematic and driving constraints
                row_new = row + sys.cons{i}.rDOF - 1;
                sys.nu(row:row_new)  =  sys.cons{i}.nu; % plug in constraint nu's
                row = row_new + 1;
            end
        end
        function constructNuP(sys) % construct nu vector (euler parameter constraints)
            % nuP = [rPDOF x 1]
            sys.nuP = zeros(sys.rPDOF,1); % initialize phiP for speed
            consNum = sys.nConstraints-sys.rPDOF+1;
            for i = 1:sys.rPDOF % only euler parameter constraints
                sys.nuP(i) = sys.cons{consNum}.nu; % plug in constraint nu's
                consNum = consNum+1;
            end
        end
        function constructNuF(sys) % construct nuF (full velocity nu vector)
            % nuF = [nConstrainedDOF x 1]
            sys.constructNu();
            sys.constructNuP();
            sys.nuF = [sys.nu; sys.nuP];
        end
        function constructGammaHat(sys) % construct gammaHat vector (kinematic & driving constraints)
            % gammaHat = [(rKDOF) x 1]
            sys.gammaHat = zeros(sys.rKDOF,1); % initialize nu for speed
            row = 1;
            for i = 1:(sys.nConstraints-sys.rPDOF) % only kinematic and driving constraints
                row_new = row + sys.cons{i}.rDOF - 1;
                sys.gammaHat(row:row_new)  =  sys.cons{i}.gammaHat; % plug in constraint nu's
                row = row_new + 1;
            end
        end
        function constructGammaHatP(sys) % construct gammaHatP vector (euler parameter constraints)
            % gammaHatP = [rPDOF x 1]
            sys.gammaHatP = zeros(sys.rPDOF,1); % initialize gammaHatP for speed
            consNum = sys.nConstraints-sys.rPDOF+1;
            for i = 1:sys.rPDOF % only euler parameter constraints
                sys.gammaHatP(i) = sys.cons{consNum}.gammaHat; % plug in constraint gammaHat's
                consNum = consNum+1;
            end
        end
        function constructGammaHatF(sys) % construct gammaHatF vector (Full RHS Of acceleration equation)
            % gammaHatF = [nConstrainedDOF x 1]
            sys.constructGammaHat();
            sys.constructGammaHatP();
            sys.gammaHatF = [sys.gammaHat; sys.gammaHatP];
        end
        function constructQ(sys) % construct q (r and p) of bodies in system
            % get ID's of bodies that are not grounded
            bodyID = sys.bodyIDs;
            
            % initialize for speed
            r = zeros(3*sys.nFreeBodies,1);
            rdot = zeros(3*sys.nFreeBodies,1);
            rddot = zeros(3*sys.nFreeBodies,1);
            p = zeros(4*sys.nFreeBodies,1);
            pdot = zeros(4*sys.nFreeBodies,1);
            pddot = zeros(4*sys.nFreeBodies,1);
            
            % pull r,rdot,rddot and p,pdot,pddot out of ungrounded bodies
            rowr = 1; % r count
            rowp = 1; % p count
            for i = bodyID  
                r(rowr:rowr+2) = sys.body{i}.r;
                rdot(rowr:rowr+2) = sys.body{i}.rdot;
                rddot(rowr:rowr+2) = sys.body{i}.rddot;
                p(rowp:rowp+3) = sys.body{i}.p;
                pdot(rowp:rowp+3) = sys.body{i}.pdot;
                pddot(rowp:rowp+3) = sys.body{i}.pddot;
                rowr = rowr + 3;
                rowp = rowp + 4;
            end
            
            % assemble q,qdot,qddot vector
            sys.q = [r;p]; 
            sys.qdot = [rdot;pdot];
            sys.qddot = [rddot;pddot];
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
                p = sys.body{bodyID(i)}.p; % pull euler parameters
                TauHat = 2*G'*nBar + 8*Gdot'*J*Gdot*p;
                sys.TauHat((4*i-3):4*i) = TauHat;
            end
        end
        function psi = constructPsiMatrix(sys,h,beta0,lambda,lambdaP) % construct the full blown Jacobian (psi) matrix used in dynamics analysis.
            % constuct the full blown Jacobian (psi) matrix used in dynamics analysis.
            % inputs:
            %       h: step size, in seconds
            %   beta0: BDF constant
            %  lambda: current lagrange multipliers of the system
            % lambdaP: current lagrange multipliers of the system, wrt to
            % euler parameter normalization constraints
            
            % These matrices should have been updated previous
            % to the call to construct psi: 
            % sys.phi_r, sys.phi_p, sys.P, sys.F, sys.TauHat
            
            % although forces could be a function of r, rdot, p, and pdot...
            % At present, we are only supporting force functions of time. 
            % Thus, F_r, F_p, F_rdot, F_pdot will be zeros.
            F_r         = zeros(3*sys.nFreeBodies,3*sys.nFreeBodies);
            F_p         = zeros(3*sys.nFreeBodies,4*sys.nFreeBodies); 
            F_rdot      = zeros(3*sys.nFreeBodies,3*sys.nFreeBodies); 
            F_pdot      = zeros(3*sys.nFreeBodies,4*sys.nFreeBodies); 
            
            % although torques could be a function of r, rdot, p, and pdot...
            % At present, we are only supporting torque functions of time. 
            % Thus, TauHat_r, TauHat_rdot will be zeros,
            % and TauHat_p, TauHat_pdot will be simplified.
            TauHat_r    = zeros(4*sys.nFreeBodies,3*sys.nFreeBodies);
            TauHat_rdot = zeros(4*sys.nFreeBodies,3*sys.nFreeBodies);
            
            
            % build required terms: JpPddot_p, TauHat_p, TauHat_pdot
            bodyID = sys.bodyIDs; % get ID's of bodies that are not grounded (free bodies)
            JpPddot_p = zeros(4*sys.nFreeBodies,4*sys.nFreeBodies);   % preallocate for speed
            TauHat_p = zeros(4*sys.nFreeBodies,4*sys.nFreeBodies);    % preallocate for speed
            TauHat_pdot = zeros(4*sys.nFreeBodies,4*sys.nFreeBodies); % preallocate for speed
            for i = 1:sys.nFreeBodies % iterate through every free body
                %%%% find JpPddot_p (see slide 34, lecture 10/17/16)
                p = sys.body{bodyID(i)}.p;
                pddot = sys.body{bodyID(i)}.pddot;
                Jbar = sys.body{bodyID(i)}.Jbar;
                G = utility.Gmatrix(p);
                Gddot = utility.Gmatrix(pddot);
                a = Jbar*G*pddot;
                T = [0 -a'; a -utility.tilde(a)];
                JpPddot_p(4*i-3:4*i,4*i-3:4*i) = -4*G'*Jbar*Gddot + 4*T;
                
                %%%% find TauHat_p (see slide 26, lecture 10/05/16)
                pdot = sys.body{bodyID(i)}.pdot;
                Gdot = utility.Gmatrix(pdot);
                TauHat_p(4*i-3:4*i,4*i-3:4*i) = 8*Gdot'*Jbar*Gdot;  % block diagonal matrix
                
                %%%% find TauHat_pdot (see slide 26, lecture 10/05/16)
                %%%% See slide 34 from lecture 10/17/16 for similar
                %%%% derivation for finding partial of 8*Gdot'*Jbar*Gdot*p
                %%%% with respect to pdot
                a2 = Jbar*Gdot*p;
                T2 = [0 -a2'; a2 -utility.tilde(a2)];
                TauHat_pdot(4*i-3:4*i,4*i-3:4*i) = -8*Gdot'*Jbar*G + 8*T2;  % block diagonal matrix
            end
            
            % build required terms: phiRLambda_r, phiRLambda_p, phiPLambda_r, phiPLambda_p
            lambdaCounter = 1;
            phiLambda_rr = zeros(3*sys.nFreeBodies,3*sys.nFreeBodies); %preallocate for speed
            phiLambda_rp = zeros(3*sys.nFreeBodies,4*sys.nFreeBodies); %preallocate for speed
            phiLambda_pr = zeros(4*sys.nFreeBodies,3*sys.nFreeBodies); %preallocate for speed
            phiLambda_pp = zeros(4*sys.nFreeBodies,4*sys.nFreeBodies); %preallocate for speed
            for i = 1:(sys.nConstraints-sys.rPDOF) % only kinematic and driving constraints
                % need to send the right lambda to the right subconstraint
                % for calculation of the partial derivative
                nLambdas = sys.cons{i}.rDOF;
                lambdas = lambda(lambdaCounter:lambdaCounter+nLambdas-1);
                
                % calculate terms for the constraint
                phiLambda_rr_terms = sys.cons{i}.phiLambda_rr(lambdas);
                %phiLambda_rp_terms = sys.cons{i}.phiLambda_rp(lambdas);
                %phiLambda_pr_terms = sys.cons{i}.phiLambda_pr(lambdas);
                %phiLambda_pp_terms = sys.cons{i}.phiLambda_pp(lambdas);
                
                % cycle through each term in the phiLambda_terms matrices
                for j = 1:nLambdas;
                    if ~sys.cons{i}.bodyi.isGround && ~sys.cons{i}.bodyj.isGround % if bodyi and bodyj are free bodies....
                        fbNumi = find(sys.bodyIDs == sys.cons{i}.bodyi.ID); % free body number from body ID number for bodyi
                        fbNumj = find(sys.bodyIDs == sys.cons{i}.bodyj.ID); % free body number from body ID number for bodyj
                        %%%%%% phiLambda_rr term is [6x6];
                        % plug in terms for phiLambda_rrii
                        phiLambda_rr(3*fbNumi-2:3*fbNumi,3*fbNumi-2:3*fbNumi) = phiLambda_rr(3*fbNumi-2:3*fbNumi,3*fbNumi-2:3*fbNumi) + phiLambda_rr_terms(6*j-5:6*j-3,1:3);
                        % plug in terms for phiLambda_rrjj
                        phiLambda_rr(3*fbNumj-2:3*fbNumj,3*fbNumj-2:3*fbNumj) = phiLambda_rr(3*fbNumj-2:3*fbNumj,3*fbNumj-2:3*fbNumj) + phiLambda_rr_terms(6*j-2:6*j,4:6);
                        % plug in terms for phiLambda_rrij
                        phiLambda_rr(3*fbNumi-2:3*fbNumi,3*fbNumj-2:3*fbNumj) = phiLambda_rr(3*fbNumi-2:3*fbNumi,3*fbNumj-2:3*fbNumj) + phiLambda_rr_terms(6*j-5:6*j-3,4:6);
                        % plug in terms for phiLambda_rrji
                        phiLambda_rr(3*fbNumj-2:3*fbNumj,3*fbNumi-2:3*fbNumi) = phiLambda_rr(3*fbNumj-2:3*fbNumj,3*fbNumi-2:3*fbNumi) + phiLambda_rr_terms(6*j-2:6*j,1:3);
                    elseif ~sys.cons{i}.bodyi.isGround % if bodyi is free and bodyj is ground
                        fbNum = find(sys.bodyIDs == sys.cons{i}.bodyi.ID); % free body number from body ID number
                        %%%%%% phi_Lambda_rr term is [3x3];
                        phiLambda_rr(3*fbNum-2:3*fbNum,3*fbNum-2:3*fbNum) = phiLambda_rr(3*fbNum-2:3*fbNum,3*fbNum-2:3*fbNum) + phiLambda_rr_terms(3*j-2:3*j,1:3);
                    elseif sys.cons{i}.bodyi.isGround % if bodyi is ground and bodyj is free....
                        fbNum = find(sys.bodyIDs == sys.cons{i}.bodyj.ID); % free body number from body ID number
                        %%%%%% phi_Lambda_rr term is [3x3];
                        phiLambda_rr(3*fbNum-2:3*fbNum,3*fbNum-2:3*fbNum) = phiLambda_rr(3*fbNum-2:3*fbNum,3*fbNum-2:3*fbNum) + phiLambda_rr_terms(3*j-2:3*j,1:3);
                    end
                end
                
                
                % set lambdaCounter for next constraint
                lambdaCounter = lambdaCounter + nLambdas;
            end
            
            
            % NOW I JUST NEED TO DEFINE THESE TERMS... JOY!!!
            %             phiRLambda_r
            %             phiRLambda_p
            %             phiPLambda_r
            %             phiPLambda_p
            %             PLamdaP_p
            
            
            
            
            
            % create psi entries (see slide 30, 10/17/16)
            % psi11 : [3*nFreeBodies x 3*nFreeBodies]
            % psi12 : [3*nFreeBodies x 4*nFreeBodies]
            % psi21 : [4*nFreeBodies x 3*nFreeBodies]
            % psi22 : [4*nFreeBodies x 4*nFreeBodies]
            psi11 = sys.M + h^2*beta0^2*phiLambda_rr - h^2*beta0^2*F_r - h*beta0*F_rdot;
            psi12 = h^2*beta0^2*phiLambda_rp - h^2*beta0^2*F_p - h*beta0*F_pdot;
            psi21 = h^2*beta0^2*phiLambda_pr - h^2*beta0^2*TauHat_r - h*beta0*TauHat_rdot;
            psi22 = sys.Jp + h^2*beta0^2*JpPddot_p + h^2*beta0^2*PLamdaP_p + h^2*beta0^2*phiLambda_pp - h^2*beta0^2*TauHat_p - h*beta0*TauHat_pdot;
            
            % zeros matrices
            Z2 = zeros(sys.nFreeBodies,3*sys.nFreeBodies);
            Z3 = zeros(sys.nFreeBodies);
            Z4 = zeros(sys.rKDOF,sys.nFreeBodies);
            Z5 = zeros(sys.rKDOF);
            
            % construct psi, (see slide 30, 10/17/16)
            psi = [    psi11      psi12      Z2'  sys.phi_r';
                       psi21      psi22   sys.P'  sys.phi_p';
                          Z2      sys.P      Z3          Z4';
                   sys.phi_r  sys.phi_p      Z4          Z5];
            
        end
        function [F_r, F_p, F_rdot, F_pdot, TauHat_r, TauHat_p, TauHat_rdot, TauHat_pdot] = buildPartialDerivativesOfForces(sys) % partial derivatives of sys.F and sys.TauHat
            % for the current system time step, construct the partial
            % derivatives of sys.F and sys.TauHat
            % these partial derivatives are with respect to r,p,rdot,pdot
            
            % based on: 
            % ME751_f2016 slide 33 from lecture 10/17/16
            % ME751_f2016 slide 26-27 from lecture 10/05/16
            
            F_r         = zeros(3*sys.nFreeBodies,3*sys.nFreeBodies);
            F_p         = zeros(3*sys.nFreeBodies,4*sys.nFreeBodies); 
            F_rdot      = zeros(3*sys.nFreeBodies,3*sys.nFreeBodies); 
            F_pdot      = zeros(3*sys.nFreeBodies,4*sys.nFreeBodies); 
            TauHat_r    = zeros(4*sys.nFreeBodies,3*sys.nFreeBodies);
            TauHat_p    = zeros(4*sys.nFreeBodies,4*sys.nFreeBodies);
            TauHat_rdot = zeros(4*sys.nFreeBodies,3*sys.nFreeBodies);
            TauHat_pdot = zeros(4*sys.nFreeBodies,4*sys.nFreeBodies);
            
%             bodyID = sys.bodyIDs; % get ID's of bodies that are not grounded (free bodies)
%             for i = 1:sys.nFreeBodies % iterate through every free body
%                 %iF_r    = zeros(3,3); % iterative F_r
%                 %iF_p    = zeros(3,4); % iterative F_p
%                 %iF_rdot = zeros(3,3); % iterative F_r
%                 %iF_pdot = zeros(3,4); % iterative F_p
%                 
%                 for j = 1:sys.body{bodyID(i)}.nForces
%                     % sum partial derivatives acting on that body
%                     %iF_r    = iF_r + sys.body{bodyID(i)}.forces{j}.force_r; %[3x3]
%                     %iF_p    = iF_p + sys.body{bodyID(i)}.forces{j}.force_p; %[3x4]
%                     %iF_rdot = iF_rdot + sys.body{bodyID(i)}.forces{j}.force_rdot; %[3x3]
%                     %iF_pdot = iF_pdot + sys.body{bodyID(i)}.forces{j}.force_pdot; %[3x4]
%                     
%                     % sum active and gravitational torques acting
%                     nBar = nBar + sys.body{bodyID(i)}.forces{j}.torque; 
%                 end
%                 
% %                 p = sys.body{bodyID(i)}.p; % pull euler parameters
% %                 J = sys.body{bodyID(i)}.Jbar; % pull inertia tensor
% %                 G = utility.Gmatrix(p); % calc G matrix
% %                 Gdot = utility.Gmatrix(sys.body{bodyID(i)}.pdot); % calc Gdot matrix
% %                 G_r = zeros()
% %                 TauHat_p = 8*Gdot'*J*Gdot*p;
%                 
%                 % plug sums into final terms
%                 %F_r((3*i-2):3*i,(3*i-2):3*i)    = iF_r;
%                 %F_p((3*i-2):3*i,(4*i-3):4*i)    = iF_p;
%                 %F_rdot((3*i-2):3*i,(3*i-2):3*i) = iF_rdot;
%                 %F_pdot((3*i-2):3*i,(4*i-3):4*i) = iF_pdot;
%             end
        end
        function JpPddot_p = buildJpPddot_p(sys) % build partial derivative JpPddot_p
            % this is a term used in the construction of Newton Raphson psi
            % matrix. 
            bodyID = sys.bodyIDs; % get ID's of bodies that are not grounded (free bodies)
            JpPddot_p = zeros(4*sys.nFreeBodies,4*sys.nFreeBodies); % preallocate for speed
            for i = 1:sys.nFreeBodies % iterate through every free body
                %%%% find JpPddot_p (see slide 34, lecture 10/17/16)
                p = sys.body{bodyID(i)}.p;
                pddot = sys.body{bodyID(i)}.pddot;
                Jbar = sys.body{bodyID(i)}.Jbar;
                a = Jbar*utility.Gmatrix(p)*pddot;
                T = [0 -a'; a -utility.tilde(a)];
                JpPddot_p(4*i-3:4*i,4*i-3:4*i) = 4*(T-utility.Gmatrix(p)'*Jbar*utility.Gmatrix(pddot));
            end
        end
        function TauHat_p = buildTauHat_p(sys) % build partial derivative TauHat_p
            % this is a term used in the construction of Newton Raphson psi
            % matrix. 
            bodyID   = sys.bodyIDs; % get ID's of bodies that are not grounded (free bodies)
            TauHat_p = zeros(4*sys.nFreeBodies,4*sys.nFreeBodies); % preallocate for speed
            
            for i = 1:sys.nFreeBodies % iterate through every free body
                %%%% find TauHat_p (see slide 26, lecture 10/05/16)
                Jbar = sys.body{bodyID(i)}.Jbar;
                Gdot = utility.Gmatrix(sys.body{bodyID(i)}.pdot);
                TauHat_p(4*i-3:4*i,4*i-3:4*i) = 8*Gdot'*Jbar*Gdot;  % block diagonal matrix
            end
        end
    end
    methods % methods block with no attributes
        function bodyIDs = get.bodyIDs(sys) % calculate ID numbers of free bodies in the system
            bodyIDs = zeros(1,sys.nFreeBodies); %initialize size, must be row vector
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
        function nPoints = get.nPoints(sys) % calculate total number of points in the system
            nPoints = 0;
            if sys.nBodies == 0; return; end;
            for i = 1:sys.nBodies
                nPoints = nPoints + sys.body{i}.nPoints;
            end
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
        function r = get.r(sys) % vector of positions of the bodies
            r = sys.q(1:3*(sys.nFreeBodies));
        end
        function rdot = get.rdot(sys) % vector of velocities of the bodies
            rdot = sys.qdot(1:3*(sys.nFreeBodies));
        end
        function rddot = get.rddot(sys) % vector of accelerations of the bodies
            rddot = sys.qddot(1:3*(sys.nFreeBodies));
        end
        function p = get.p(sys) % vector of euler parameters of the bodies (rotation)
            p = sys.q(3*(sys.nFreeBodies)+1:end);
        end
        function pdot = get.pdot(sys) % vector of derivative of euler parameters of the bodies (rotational velocity)
            pdot = sys.qdot(3*(sys.nFreeBodies)+1:end);
        end
        function pddot = get.pddot(sys) % vector of 2nd derivative of euler parameters of the bodies (rotational acceleration)
            pddot = sys.qddot(3*(sys.nFreeBodies)+1:end);
        end
    end
end

