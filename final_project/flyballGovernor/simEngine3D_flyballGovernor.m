% Filename: simEngine3D_flyballGovernor.m
% Author:   Samuel Acuña
% Date:     14 Dec 2016
%
% About:    
% Driver file for the simEngine3D framework. Defines the system, the
% bodies in the system, and the constraints in the system.
% 
% This driver specifies the flyball governor mechanism, 
% as described in Haug's book, section 12.6 (pages 481-488)

clear; clc; close all

option = 1; % CHOOSE SIMULATION 
% option 1 = driving constraint to test steady state.
% option 2 = variable torque to test governor

% access directory with simEngine in it.
currentpath = cd('..'); cd('..');
addpath(pwd); % add grandparent folder to path
cd(currentpath); % return to initial folder

%% DEFINE SYSTEM %%

sys = system3D(); %initialize system

%% DEFINE BODIES IN THE SYSTEM %%
% with initial position and orientation estimates

% body 1, GROUND
sys.addBody('ground');

% body 2, SPINDLE
r_body = [0; 0.2; 0];
sys.addBody('free',r_body);

sys.body{2}.setMass(200); %kg
J = zeros(3); 
J(1,1) = 25;
J(2,2) = 50;
J(3,3) = 25;
sys.body{2}.setInertia(J);

% body 3, BALL 1
r_body = [-0.16*sind(45); 0.2 - 0.16*cosd(45); 0];
R = utility.R3(45*(pi/180));
p_body = utility.A2p(R);
sys.addBody('free',r_body,p_body);

sys.body{3}.setMass(1); %kg
J = zeros(3); 
J(1,1) = 0.1;
J(2,2) = 0.1;
J(3,3) = 0.1;
sys.body{3}.setInertia(J);

% body 4, BALL 2
r_body = [0.16*sind(45); 0.2 - 0.16*cosd(45); 0];
R = utility.R3(-45*(pi/180));
p_body = utility.A2p(R);
sys.addBody('free',r_body,p_body);

sys.body{4}.setMass(1); %kg
J = zeros(3); 
J(1,1) = 0.1;
J(2,2) = 0.1;
J(3,3) = 0.1;
sys.body{4}.setInertia(J);

% body 5, COLLAR
r_body = [0; 0.05; 0];
sys.addBody('free',r_body);

sys.body{5}.setMass(1); %kg
J = zeros(3); 
J(1,1) = 0.15;
J(2,2) = 0.125;
J(3,3) = 0.15;
sys.body{5}.setInertia(J);



%% DEFINE POINTS ON THE BODIES %%
% these are used to specify heads and tails of vectors

% GROUND
sys.body{1}.addPoint([     0;    0;    0]); % body 1, point 1
sys.body{1}.addPoint([  0.01;    0;    0]); % body 1, point 2
sys.body{1}.addPoint([     0;    0; 0.01]); % body 1, point 3

% SPINDLE
sys.body{2}.addPoint([     0;    0;    0]); % body 2, point 1
sys.body{2}.addPoint([     0; -0.2;    0]); % body 2, point 2
sys.body{2}.addPoint([  0.01;    0;    0]); % body 2, point 3
sys.body{2}.addPoint([     0;    0; 0.01]); % body 2, point 4

% BALL 1
sys.body{3}.addPoint([    0;    0;    0]); % body 3, point 1
sys.body{3}.addPoint([ 0.16;    0;    0]); % body 3, point 2
sys.body{3}.addPoint([    0;    0; 0.01]); % body 3, point 3
sys.body{3}.addPoint([ 0.08;    0;    0]); % body 3, point 4

% BALL 2
sys.body{4}.addPoint([    0;    0;    0]); % body 4, point 1
sys.body{4}.addPoint([-0.16;    0;    0]); % body 4, point 2
sys.body{4}.addPoint([    0;    0;-0.01]); % body 4, point 3
sys.body{4}.addPoint([-0.08;    0;    0]); % body 4, point 4

% COLLAR
sys.body{5}.addPoint([    0;    0;    0]); % body 5, point 1
sys.body{5}.addPoint([    0; 0.01;    0]); % body 5, point 2
sys.body{5}.addPoint([    0;    0; 0.01]); % body 5, point 3

%% Add TSDAs to the system

% spring, SPINDLE to COLLAR
k = 1000;
l0 = 0.15;
c = 30;
h = @(l,ldot,t)0;
sys.addTSDA(sys.body{2},1,sys.body{5},1,k,l0,c,h)

%% PLOT THE SYSTEM in 3D %%
sys.plot(1,0.05) % plot with reference frames
%view(44,52); 
view(0,90);
axis equal

%% DEFINE CONSTRAINTS AMONG THE BODIES %%

% KINEMATIC CONSTRAINTS

% distance: BALL 1 to SPINDLE 
sys.addConstraint('d',sys.body{3},4,sys.body{5},1,0.10922)

% distance: BALL 2 to SPINDLE
sys.addConstraint('d',sys.body{4},4,sys.body{5},1,0.10922)

% revolute joint: GROUND to SPINDLE 
sys.addConstraint('rj',sys.body{1},1,1,2,1,3,sys.body{2},2,1,2); 

% revolute joint: SPINDLE to BALL 1
sys.addConstraint('rj',sys.body{2},1,1,2,1,3,sys.body{3},2,1,3); 

% revolute joint: SPINDLE to BALL 2
sys.addConstraint('rj',sys.body{2},1,1,2,1,3,sys.body{4},2,1,3)

% translational joint: SPINDLE to COLLAR
sys.addConstraint('tj',sys.body{2},1,1,3,1,4,sys.body{5},1,1,3,1,2); 

if option == 1
    % DRIVING CONSTRAINTS
    % uses dp1 constraint to specify angle of spindle
    a = norm(sys.body{1}.point{3});
    b = norm(sys.body{2}.point{3});
    f = @(t) a*b*cos(11.0174*t + pi/2 + 0.0001);
    fdot = @(t) a*b*(-11.0174*sin(pi/2 + 11.0174*t + 0.0001));
    fddot = @(t) a*b*(-(11.0174)^2*cos(pi/2 + 11.0174*t + 0.0001));
    sys.addConstraint('dp1',sys.body{1},1,3,sys.body{2},1,3,f,fdot,fddot) % unit length vectors
end
% add euler parameter normalization constraints to system
sys.addEulerParamConstraints(); 


%% ASSEMBLE CONSTRAINT MATRIX 
sys.assembleConstraints();

if option == 2
    % ADD EXTERNAL TORQUE to SPINDLE
    torqueDirection = [0;1;0];
    torqueMagnitude = -25;
    sys.body{2}.addForces('piecewiseTorque',torqueMagnitude,torqueDirection,1,2,0,-25,25,-25);
    
    % ADD ENGINE COMPENSATION TORQUE to SPINDLE
    torqueMagnitude = 7500;
    sys.body{2}.addForces('specificTorque',torqueMagnitude,torqueDirection);
end
%% INITIAL CONDITIONS 

w0 = [0;11.0174;0]; % rad/2
sys.body{2}.pdot = utility.w2pdot(w0,sys.body{2}.p);

sys.setInitialVelocities(); % find initial velocities of constrained bodies
 
%sys.checkInitialConditions(1e-4);

%% ADD FORCES
sys.g = [0;-9.81;0];
sys.addGravityForces();


%% KINEMATICS ANALYSIS
% timeStart = 0; %seconds
% timeEnd =  10;
% timeStep = 10^-3;
% 
% tic
% state = sys.kinematicsAnalysis(timeStart,timeEnd,timeStep);
% timeQN = toc;
% 
% 
% return
%% DYNAMICS ANALYSIS
timeStart = 0; %seconds
timeEnd =  10;
timeStep = 10^-2;

tic
state = sys.dynamicsAnalysis(timeStart,timeEnd,timeStep,'QN');
timeQN = toc;

%save('state_flyballGovernor1.mat','state')
%save('state_flyballGovernor2.mat','state')
disp('done with dynamics analysis.')
disp(['For step-size h=' num2str(timeStep) ' seconds,'])
disp(['Time to compute Quasi-Newton Solution: ' num2str(timeQN) ' seconds'])


%% PLOT COLLAR POSITION
% pull simEngine3D data:
time = timeStart:timeStep:timeEnd;
X = zeros(length(state),1); % preallocate for speed
Y = zeros(length(state),1); % preallocate for speed
for i = 1:length(state)
    Y(i) = state{i}.r(11); % y position of collar
end

% plot collar position 
figure
fig = gcf;
fig.Color = [1 1 1]; % set background color to white

hold on
plot(time,Y,'r')
hold off
lh1 = legend('Y simEngine3D');
set(lh1,'FontSize',12)
xlabel('Time (sec)','FontSize',12)
ylabel('Collar Displacement (m)','FontSize',12)
title('FLYBALL GOVERNOR','FontSize',12)
grid on

%% ANIMATION OF SYSTEM
% play animation of the dynamics analysis
plot.animateSystem(sys,state);


