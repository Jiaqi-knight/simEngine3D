
% Filename: simEngine3D_A6P3.m
% Author:   Samuel Acu�a
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
sys.addBody('ground') % body 1, ground

% body 2, in initial configuration. All units in SI
t = 0; L = 2; 
theta = pi/4*cos(2*t); %initial rotation,radians
R = utility.R2(pi/2)*utility.R3(theta); % initial pendulum rotation
p = utility.A2p(R); % euler parameter for pendulum LRF rotation
sys.addBody('free',[0;L*sin(theta);-L*cos(theta)],p); % add body 2

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


%% ASSEMBLE CONSTRAINT MATRIX 
sys.assembleConstraints()

%% ADD FORCES
sys.addGravityForces();

%% INVERSE DYNAMICS ANALYSIS
timeStart = 0; %seconds
timeEnd = 10;
timeStep = 10^-2; %10^-3;

%state = sys.inverseDynamicsAnalysis(timeStart,timeEnd,timeStep);
%save('state_A7.mat','state')
load('state_A7_quick.mat')




%% PLOT REACTION TORQUE

time = timeStart:timeStep:timeEnd;
 
% reaction torque for pendulum
torque = zeros(length(state),1); % preallocate for speed
for i = 1:length(state)
    torque(i) = state{i}.rTorque{1}(3);
end

% plot reaction torque
figure 
plot(time,torque)
title('Reaction Torque of Pendulum')
xlabel('Time (sec)')
ylabel('Torque (N-m)')
% 
% subplot(3,1,2)
% hold on
% plot(time,rdotOprime(:,1))
% plot(time,rdotOprime(:,2))
% plot(time,rdotOprime(:,3))
% title('Velocity of point O-prime')
% xlabel('Time (sec)')
% ylabel('Velocity (m/s)')
% legend('X','Y','Z')
% hold off
% 
% subplot(3,1,3)
% hold on
% plot(time,rddotOprime(:,1))
% plot(time,rddotOprime(:,2))
% plot(time,rddotOprime(:,3))
% title('Acceleration of point O-prime')
% xlabel('Time (sec)')
% ylabel('Acceleration (m/s^2)')
% legend('X','Y','Z')
% hold off
% 
% 
% % derive kinematics of point Q, at the hinge
% sBar = [-2 0 0]';
% rQ          = zeros(length(state),3);
% rQdot       = zeros(length(state),3);
% rQddot      = zeros(length(state),3);
% 
% for i = 1:length(state)
%     % get pOprime
%     p     = state{i}.q(4:7); 
%     pdot  = state{i}.qdot(4:7); 
%     pddot = state{i}.qddot(4:7); 
%     
%     % derive rQ
%     A = utility.p2A(p);
%     B = utility.Bmatrix(p,sBar);
%     Bdot = utility.Bmatrix(pdot,sBar);
%     rQ(i,:) = rOprime(i,:)' + A*sBar;
%     rQdot(i,:) = rdotOprime(i,:)' + B*pdot;
%     rQddot(i,:) = rddotOprime(i,:)' + B*pddot + Bdot*pdot;
% end
% 
% % plot position, velocity, acceleration of point Q
% figure
% subplot(3,1,1)
% hold on
% plot(time,rQ(:,1))
% plot(time,rQ(:,2))
% plot(time,rQ(:,3))
% title('Position of point Q')
% xlabel('Time (sec)')
% ylabel('Position (m)')
% legend('X','Y','Z')
% hold off
% 
% subplot(3,1,2)
% hold on
% plot(time,rQdot(:,1))
% plot(time,rQdot(:,2))
% plot(time,rQdot(:,3))
% title('Velocity of point Q')
% xlabel('Time (sec)')
% ylabel('Velocity (m/s)')
% legend('X','Y','Z')
% hold off
% 
% subplot(3,1,3)
% hold on
% plot(time,rQddot(:,1))
% plot(time,rQddot(:,2))
% plot(time,rQddot(:,3))
% title('Acceleration of point Q')
% xlabel('Time (sec)')
% ylabel('Acceleration (m/s^2)')
% legend('X','Y','Z')
% hold off
% 
% 

