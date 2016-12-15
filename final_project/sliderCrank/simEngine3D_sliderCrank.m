% Filename: simEngine3D_sliderCrank.m
% Author:   Samuel Acuña
% Date:     17 Nov 2016
%
% About:    
% Driver file for the simEngine3D framework. Defines the system, the
% bodies in the system, and the constraints in the system.
% 
% This driver specifies the slider crank mechanism, 
% as described in Haug's book, section 12.2 (pages 455-458)

clear; clc; close all

% access directory with simEngine in it.
currentpath = cd('..'); cd('..');
addpath(pwd); % add grandparent folder to path
cd(currentpath); % return to initial folder

%% DEFINE SYSTEM %%

sys = system3D(); %initialize system

%% DEFINE BODIES IN THE SYSTEM %%
% with initial position and orientation estimates

% body 1, crank
r_body = [0; 0.1; 0.12];
R = utility.R1(pi/2); 
p_body = utility.A2p(R); % euler parameter
sys.addBody('free',r_body,p_body);

sys.body{1}.setMass(0.12); %kg
J = zeros(3); 
J(1,1) = 0.0001;
J(2,2) = 0.00001;
J(3,3) = 0.0001;
sys.body{1}.setInertia(J);


% body 2, connecting rod
r_body = [0.1; 0.05; 0.1];
p_body = [0.9;-0.21;0.4;-0.1]; % euler parameter estimate
sys.addBody('free',r_body,p_body);

sys.body{2}.setMass(0.5); %kg
J = zeros(3); 
J(1,1) = 0.0004;
J(2,2) = 0.004;
J(3,3) = 0.004;
sys.body{2}.setInertia(J);

% body 3, slider
r_body = [0.2; 0; 0];
sys.addBody('free',r_body);

sys.body{3}.setMass(2.0); %kg
J = zeros(3); 
J(1,1) = 0.0001;
J(2,2) = 0.0001;
J(3,3) = 0.0001;
sys.body{3}.setInertia(J);

% body 4, ground
sys.addBody('ground') 


%% DEFINE POINTS ON THE BODIES %%
% these are used to specify heads and tails of vectors

% crank
sys.body{1}.addPoint([    0;    0; 0]); % body 1, point 1
sys.body{1}.addPoint([ 0.01;    0; 0]); % body 1, point 2
sys.body{1}.addPoint([    0; 0.08; 0]); % body 1, point 3

% % connecting rod
sys.body{2}.addPoint([     0;    0;    0]); % body 2, point 1
sys.body{2}.addPoint([ -0.15;    0;    0]); % body 2, point 2
sys.body{2}.addPoint([  0.15;    0;    0]); % body 2, point 3
sys.body{2}.addPoint([  0.15; 0.05;    0]); % body 2, point 4


% slider
sys.body{3}.addPoint([    0;    0;    0]); % body 2, point 1
sys.body{3}.addPoint([ 0.20;    0;    0]); % body 2, point 2
sys.body{3}.addPoint([ 0.25;    0;    0]); % body 2, point 3
sys.body{3}.addPoint([    0; 0.05;    0]); % body 2, point 4
sys.body{3}.addPoint([    0;    0; 0.05]); % body 2, point 5
sys.body{3}.addPoint([ 0.05;    0;    0]); % body 2, point 6

% ground
sys.body{4}.addPoint([    0;    0;    0]); % body 4, point 1
sys.body{4}.addPoint([    0; 0.05;    0]); % body 4, point 2
sys.body{4}.addPoint([    0;    0; 0.05]); % body 4, point 3
sys.body{4}.addPoint([    0;  0.1; 0.12]); % body 4, point 4
sys.body{4}.addPoint([ 0.05;    0;    0]); % body 4, point 5


%% move connecting rod into place
crank_tip = sys.body{1}.r + utility.p2A(sys.body{1}.p)*sys.body{1}.point{3};
slider_tip = sys.body{3}.r;
r_new = (crank_tip+slider_tip)./2;

A1 = utility.R3(-(pi/2-atan2(slider_tip(1),crank_tip(2))));
%A2 = utility.R2(pi/2-atan2(slider_tip(1),crank_tip(3))); % starting value rotation
A2 = utility.R2( 0.7297276562270); % iterated until aligned perfectly
A3 = utility.R1(-0.6435011087932); % iterated until aligned perfectly
R = A1*A2*A3;
p_new = utility.A2p(R);

sys.body{2}.r = r_new;
sys.body{2}.p = p_new;

% A2 iteration
%[crank_tip,  sys.body{2}.r + R*sys.body{2}.point{2}] % to compare values
%[slider_tip, sys.body{2}.r + R*sys.body{2}.point{3}] % to compare values

% A3 iteration
%hinge = sys.body{2}.r + R*sys.body{2}.point{4} - slider_tip;
%dot(slider_tip,hinge) % to compare values

%% PLOT THE SYSTEM in 3D %%
sys.plot(1,0.05) % plot with reference frames
view(121,16);
axis equal

%% DEFINE CONSTRAINTS AMONG THE BODIES %%

% KINEMATIC CONSTRAINTS

% revolute joint: ground to crank (table 10.2.1)
sys.addConstraint('rj',sys.body{4},4,1,2,1,3,sys.body{1},1,1,2); % consistent

% spherical joint: crank to connecting rod (table 10.2.2)
sys.addConstraint('sj',sys.body{1},3,sys.body{2},2); % consistent

% revolute-cylindrical joint: connecting rod to slider (table 10.2.3)
%sys.addConstraint('rcj',sys.body{2},3,4,3,sys.body{3},2,3,1,4,1,5,2); %
%for some reason, this constraint isnt working...

% dp1 and p2, part of rcj
sys.addConstraint('dp1',sys.body{2},3,4,sys.body{3},2,3)
sys.addConstraint('p2',sys.body{3},1,4,1,5,2,sys.body{2},3)

% translational joint: slider to ground (table 10.2.4)
sys.addConstraint('tj',sys.body{3},1,1,4,1,5,sys.body{4},1,1,3,1,5); 

% distance constraint: connecting rod to slider (table 10.2.5)
sys.addConstraint('d',sys.body{2},3,sys.body{3},6,0.05); 


% DRIVING CONSTRAINTS
% uses dp1 constraint to specify angle of crank
a = norm(sys.body{1}.point{3});
b = norm(sys.body{4}.point{2});
f = @(t) a*b*cos(-2*pi*t + pi/2 + 0.0001);
fdot = @(t) a*b*(-2*pi*sin(2*pi*t - pi/2 + 0.0001));
fddot = @(t) a*b*(-4*pi^2*cos(2*pi*t - pi/2 + 0.0001));
sys.addConstraint('dp1',sys.body{1},1,3,sys.body{4},1,2,f,fdot,fddot) % unit length vectors

% add euler parameter normalization constraints to system
sys.addEulerParamConstraints(); 

%% ASSEMBLE CONSTRAINT MATRIX 
sys.assembleConstraints();

%% INITIAL CONDITIONS 

sys.setInitialVelocities(); % find initial velocities of constrained bodies
 
%sys.checkInitialConditions(1e-4);

%% ADD FORCES
sys.g = [0;-9.81;0];
sys.addGravityForces();


% %% KINEMATICS ANALYSIS
% timeStart = 0; %seconds
% timeEnd =  1;
% timeStep = 10^-3;
% 
% tic
% state = sys.kinematicsAnalysis(timeStart,timeEnd,timeStep);
% timeQN = toc;



%% DYNAMICS ANALYSIS
timeStart = 0; %seconds
timeEnd =  1;
timeStep = 10^-3;

tic
state = sys.dynamicsAnalysis(timeStart,timeEnd,timeStep,'QN');
timeQN = toc;

%save('state_sliderCrank10.mat','state')
disp('done with dynamics analysis.')
disp(['For step-size h=' num2str(timeStep) ' seconds,'])
disp(['Time to compute Quasi-Newton Solution: ' num2str(timeQN) ' seconds'])

%% ANIMATION OF SYSTEM
% play animation of the dynamics analysis
plot.animateSystem(sys,state,[121,16],10);


%% PLOT KINEMATICS
time = timeStart:timeStep:timeEnd;

% kinematics for slider BODY 3
rOprime3 = zeros(length(state),3); % preallocate for speed
rdotOprime3 = zeros(length(state),3);
rddotOprime3 = zeros(length(state),3);
for i = 1:length(state)
    rOprime3(i,1) = state{i}.r(4); % x value of rOprime3
    rOprime3(i,2) = state{i}.r(5); % y value of rOprime3
    rOprime3(i,3) = state{i}.r(6); % z value of rOprime3
    rdotOprime3(i,1) = state{i}.rdot(4); % xdot value of rOprime3
    rdotOprime3(i,2) = state{i}.rdot(5); % ydot value of rOprime3
    rdotOprime3(i,3) = state{i}.rdot(6); % zdot value of rOprime3
    rddotOprime3(i,1) = state{i}.rddot(4); % xddot value of rOprime3
    rddotOprime3(i,2) = state{i}.rddot(5); % yddot value of rOprime3
    rddotOprime3(i,3) = state{i}.rddot(6); % zddot value of rOprime3
end


% plot position, velocity, acceleration of Oprime, BODY 3
figure 
subplot(3,1,1)
hold on
plot(time,rOprime3(:,1))
%plot(time,rOprime3(:,2))
%plot(time,rOprime3(:,3))
title('Position of Slider')
xlabel('Time (sec)')
ylabel('Position (m)')
%legend('X','Y','Z')
hold off

subplot(3,1,2)
hold on
plot(time,rdotOprime3(:,1))
%plot(time,rdotOprime3(:,2))
%plot(time,rdotOprime3(:,3))
title('Velocity of Slider')
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
%legend('X','Y','Z')
hold off

subplot(3,1,3)
hold on
plot(time,rddotOprime3(:,1))
%plot(time,rddotOprime3(:,2))
%plot(time,rddotOprime3(:,3))
title('Acceleration of Slider')
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
%legend('X','Y','Z')
hold off



