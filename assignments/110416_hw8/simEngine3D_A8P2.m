% Filename: simEngine3D_A8P2.m
% Author:   Samuel Acuña
% Date:     02 Nov 2016
%
% About:    
% Driver file for the simEngine3D framework, for hw8 . Defines the system, the
% bodies in the system, and the constraints in the system.

clear; clc; close all

% access directory with simEngine in it.
currentpath = cd('..'); cd('..');
addpath(pwd); % add grandparent folder to path
cd(currentpath); % return to initial folder

%% DEFINE SYSTEM %%

sys = system3D(); %initialize system

%% DEFINE BODIES IN THE SYSTEM %%

% body 1, ground
sys.addBody('ground') % body 1, ground

% body 2, in initial configuration. All units in SI
L2 = 2; 
theta = pi/2; %initial rotation,radians
R = utility.R2(pi/2)*utility.R3(theta); % initial pendulum rotation
p = utility.A2p(R); % euler parameter for pendulum LRF rotation
sys.addBody('free',[0;L2*sin(theta);-L2*cos(theta)],p); % add body 2
sys.body{2}.color = [0.4141, 0.3516, 0.8008]; %slateblue

% set body 2 mass
rho = 7800; Length = 4; width = 0.05; 
volume = Length*width*width;
mass = rho*volume;
sys.body{2}.setMass(mass);

% set body 2 inertia
J = zeros(3); 
J(1,1) = mass/12*(width^2+width^2);
J(2,2) = mass/12*(Length^2+width^2);
J(3,3) = mass/12*(Length^2+width^2);
sys.body{2}.setInertia(J);

% body 3, in initial configuration. All units in SI
L3 = L2/2;
theta = 0; %initial rotation,radians
R = utility.R2(pi/2)*utility.R3(theta); % initial pendulum rotation
p = utility.A2p(R); % euler parameter for pendulum LRF rotation
sys.addBody('free',[0;2*L2+L3*sin(theta);-L3*cos(theta)],p); % add body 3
sys.body{3}.color = [0.6953, 0.1328, 0.1328]; %Firebrick red

% set body 3 mass
rho = 7800; Length = 2; width = 0.05; 
volume = Length*width*width;
mass = rho*volume;
sys.body{3}.setMass(mass);

% set body 3 inertia
J = zeros(3); 
J(1,1) = mass/12*(width^2+width^2);
J(2,2) = mass/12*(Length^2+width^2);
J(3,3) = mass/12*(Length^2+width^2);
sys.body{3}.setInertia(J);


%% DEFINE POINTS ON THE BODIES %%
% these are used to specify heads and tails of vectors
sys.body{1}.addPoint([0;0;0]); % body 1, point 1
sys.body{1}.addPoint([0;1;0]); % body 1, point 2
sys.body{1}.addPoint([0;0;-1]); % body 1, point 3
sys.body{2}.addPoint([0;0;0]); % body 2, point 1
sys.body{2}.addPoint([0;0;1]); % body 2, point 2
sys.body{2}.addPoint([-2;0;0]); % body 2, point 3
sys.body{2}.addPoint([2;0;0]); % body 2, point 4
sys.body{2}.addPoint([0;-1;0]); % body 2, point 5

sys.body{3}.addPoint([-1;0;0]); % body 3, point 1
sys.body{3}.addPoint([0;0;0]); % body 3, point 2
sys.body{3}.addPoint([0;0;1]); % body 3, point 3


%% PLOT THE SYSTEM in 3D %%
sys.plot(1) % plot with reference frames
%sys.plot() % plot without reference frames
view(98,12);
axis equal

%% DEFINE CONSTRAINTS AMONG THE BODIES %%

% KINEMATIC CONSTRAINTS
% revolute joint with ground body
sys.addConstraint('rj',sys.body{1},1,1,2,1,3,sys.body{2},3,1,2)

% revolute joint body 3 to body 2
sys.addConstraint('rj',sys.body{2},4,1,4,1,5,sys.body{3},1,2,3)


%% ASSEMBLE CONSTRAINT MATRIX 
sys.assembleConstraints()

%% ADD FORCES
sys.addGravityForces();

%% DYNAMICS ANALYSIS
timeStart = 0; %seconds
timeEnd =  10;
timeStep = 10^-2; %10^-3;
tic
state = sys.dynamicsAnalysis(timeStart,timeEnd,timeStep);
toc

%save('state_A8P2.mat','state')
%load('state_A7.mat')
disp('done')

%% PLOT REACTION TORQUE

time = timeStart:timeStep:timeEnd;

% kinematics for pendulum, the Oprime frame, BODY 2
rOprime2 = zeros(length(state),3); % preallocate for speed
rdotOprime2 = zeros(length(state),3);
rddotOprime2 = zeros(length(state),3);
for i = 1:length(state)
    rOprime2(i,1) = state{i}.r(1); % x value of rOprime2
    rOprime2(i,2) = state{i}.r(2); % y value of rOprime2
    rOprime2(i,3) = state{i}.r(3); % z value of rOprime2
    rdotOprime2(i,1) = state{i}.rdot(1); % xdot value of rOprime2
    rdotOprime2(i,2) = state{i}.rdot(2); % ydot value of rOprime2
    rdotOprime2(i,3) = state{i}.rdot(3); % zdot value of rOprime2
    rddotOprime2(i,1) = state{i}.rddot(1); % xddot value of rOprime2
    rddotOprime2(i,2) = state{i}.rddot(2); % yddot value of rOprime2
    rddotOprime2(i,3) = state{i}.rddot(3); % zddot value of rOprime2
end

% plot position, velocity, acceleration of Oprime, BODY 2
figure 
subplot(3,1,1)
hold on
plot(time,rOprime2(:,1))
plot(time,rOprime2(:,2))
plot(time,rOprime2(:,3))
title('Position of Body 1 (O-prime)')
xlabel('Time (sec)')
ylabel('Position (m)')
legend('X','Y','Z')
hold off

subplot(3,1,2)
hold on
plot(time,rdotOprime2(:,1))
plot(time,rdotOprime2(:,2))
plot(time,rdotOprime2(:,3))
title('Velocity of Body 1 (O-prime)')
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
legend('X','Y','Z')
hold off

subplot(3,1,3)
hold on
plot(time,rddotOprime2(:,1))
plot(time,rddotOprime2(:,2))
plot(time,rddotOprime2(:,3))
title('Acceleration of Body 1 (O-prime)')
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
legend('X','Y','Z')
hold off


%% kinematics for pendulum, the Oprime frame, BODY 3
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
plot(time,rOprime3(:,2))
plot(time,rOprime3(:,3))
title('Position of Body 2 (O-prime)')
xlabel('Time (sec)')
ylabel('Position (m)')
legend('X','Y','Z')
hold off

subplot(3,1,2)
hold on
plot(time,rdotOprime3(:,1))
plot(time,rdotOprime3(:,2))
plot(time,rdotOprime3(:,3))
title('Velocity of Body 2 (O-prime)')
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
legend('X','Y','Z')
hold off

subplot(3,1,3)
hold on
plot(time,rddotOprime3(:,1))
plot(time,rddotOprime3(:,2))
plot(time,rddotOprime3(:,3))
title('Acceleration of Body 2 (O-prime)')
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
legend('X','Y','Z')
hold off

%% plot violation of the velocity constraint
violation = zeros(length(state),1);
for i = 1:length(state)
    v = state{i}.velocityViolation;
    violation(i) = norm(v(6:10));
end

figure
plot(time,violation)
title('Norm-2 Violation of the Velocity Constraints (Rev Joint: Body 1 & Body 2)')
xlabel('Time (sec)')
ylabel('Violation (m/s)')

