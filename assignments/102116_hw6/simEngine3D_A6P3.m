
% Filename: simEngine3D_A6P3.m
% Author:   Samuel Acuña
% Date:     18 Oct 2016
%
% About:    
% Driver file for the simEngine3D framework, for hw6 . Defines the system, the
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
sys.addBody([0;0;0],[],[],[],[],[],1) 

% body 2, in initial configuration
t = 0; L = 2;
theta = pi/4*cos(2*t); %radians
R = utility.R2(pi/2)*utility.R3(theta); % initial pendulum rotation
p = utility.A2p(R); % euler parameter for pendulum LRF rotation
sys.addBody([0;L*sin(theta);-L*cos(theta)],p);

%% DEFINE POINTS ON THE BODIES %%
% these are used to specify heads and tails of vectors
sys.body{1}.addPoint([0;0;0]); % body 1, point 1
sys.body{1}.addPoint([0;1;0]); % body 1, point 2
sys.body{1}.addPoint([0;0;-1]); % body 1, point 3
sys.body{2}.addPoint([0;0;0]); % body 2, point 1
sys.body{2}.addPoint([0;0;1]); % body 2, point 2
sys.body{2}.addPoint([-2;0;0]); % body 2, point 3
sys.body{2}.addPoint([1;0;0]); % body 2, point 4

%% PLOT THE SYSTEM in 3D %%
sys.plot(1) % plot with reference frames
% sys.plot() % plot without reference frames
view(98,12);
axis equal

%% DEFINE CONSTRAINTS AMONG THE BODIES %%

% KINEMATIC CONSTRAINTS
% revolute joint with ground body
sys.addConstraint('rj',sys.body{1},1,1,2,1,3,sys.body{2},3,1,2)

% DRIVING CONSTRAINTS
% uses dp1 constraint to specify angle of pendulum
f = @(t)cos((pi*cos(2*t))/4 - pi/2);%cos((pi*cos(2*t))/4); 
fdot = @(t)((pi*sin(2*t)*sin((pi*cos(2*t))/4 - pi/2))/2);%((pi*sin(2*t)*sin((pi*cos(2*t))/4))/2);
fddot = @(t)(pi*cos(2*t)*sin((pi*cos(2*t))/4 - pi/2) - (pi^2*sin(2*t)^2*cos((pi*cos(2*t))/4 - pi/2))/4);%(pi*cos(2*t)*sin((pi*cos(2*t))/4) - (pi^2*sin(2*t)^2*cos((pi*cos(2*t))/4))/4);
sys.addConstraint('dp1',sys.body{1},1,2,sys.body{2},1,4,f,fdot,fddot,t) % unit length vectors

% EULER PARAMETER NORMALIZATION CONSTRAINTS
sys.addConstraint('p_norm')

%% ASSEMBLE CONSTRAINT MATRIX 
sys.assembleConstraints()

%% KINEMATIC ANALYSIS
timeStart = 0; %seconds
timeEnd = 10;
timeStep = 10^-2;%10^-3;

state = sys.kinematicsAnalysis(timeStart,timeEnd,timeStep);
%save('state.mat','state')
%load('state.mat')

%% PLOT KINEMATICS

time = timeStart:timeStep:timeEnd;

% origin of pendulum
rOprime(:,1)     = state(:,1,1); rOprime(:,2)     = state(:,2,1); rOprime(:,3)     = state(:,3,1);
rdotOprime(:,1)  = state(:,1,2); rdotOprime(:,2)  = state(:,2,2); rdotOprime(:,3)  = state(:,3,2);
rddotOprime(:,1) = state(:,1,3); rddotOprime(:,2) = state(:,2,3); rddotOprime(:,3) = state(:,3,3);
pOprime(:,1)     = state(:,4,1); pOprime(:,2)     = state(:,5,1); pOprime(:,3)     = state(:,6,1); pOprime(:,4)     = state(:,7,1);
pdotOprime(:,1)  = state(:,4,2); pdotOprime(:,2)  = state(:,5,2); pdotOprime(:,3)  = state(:,6,2); pdotOprime(:,4)  = state(:,7,2);
pddotOprime(:,1) = state(:,4,3); pddotOprime(:,2) = state(:,5,3); pddotOprime(:,3) = state(:,6,3); pddotOprime(:,4) = state(:,7,3);


figure
subplot(3,1,1)
hold on
plot(time,rOprime(:,1))
plot(time,rOprime(:,2))
plot(time,rOprime(:,3))
title('Position of point O-prime')
xlabel('Time (sec)')
ylabel('Position (m)')
legend('X','Y','Z')
hold off

subplot(3,1,2)
hold on
plot(time,rdotOprime(:,1))
plot(time,rdotOprime(:,2))
plot(time,rdotOprime(:,3))
title('Velocity of point O-prime')
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
legend('X','Y','Z')
hold off

subplot(3,1,3)
hold on
plot(time,rddotOprime(:,1))
plot(time,rddotOprime(:,2))
plot(time,rddotOprime(:,3))
title('Acceleration of point O-prime')
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
legend('X','Y','Z')
hold off


% point Q, at the hinge
sBar = [-2 0 0]';
rQ     = zeros(length(time),3);
rQdot  = zeros(length(time),3);
rQddot = zeros(length(time),3);
for i = 1:length(time)
    p = pOprime(i,:)';
    pdot = pdotOprime(i,:)';
    pddot = pddotOprime(i,:)';
    
    A = utility.p2A(p);
    B = utility.Bmatrix(p,sBar);
    Bdot = utility.Bmatrix(pdot,sBar);
    rQ(i,:) = rOprime(i,:)' + A*sBar;
    rQdot(i,:) = rdotOprime(i,:)' + B*pdot;
    rQddot(i,:) = rddotOprime(i,:)' + B*pddot + Bdot*pdot;
end

figure
subplot(3,1,1)
hold on
plot(time,rQ(:,1))
plot(time,rQ(:,2))
plot(time,rQ(:,3))
title('Position of point Q')
xlabel('Time (sec)')
ylabel('Position (m)')
legend('X','Y','Z')
hold off

subplot(3,1,2)
hold on
plot(time,rQdot(:,1))
plot(time,rQdot(:,2))
plot(time,rQdot(:,3))
title('Velocity of point Q')
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
legend('X','Y','Z')
hold off

subplot(3,1,3)
hold on
plot(time,rQddot(:,1))
plot(time,rQddot(:,2))
plot(time,rQddot(:,3))
title('Acceleration of point Q')
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
legend('X','Y','Z')
hold off


