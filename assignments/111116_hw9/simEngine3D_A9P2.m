% Filename: simEngine3D_A9P2.m
% Author:   Samuel Acuña
% Date:     09 Nov 2016
%
% About:    
% Driver file for the simEngine3D framework, for hw9 . Defines the system, the
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
sys.body{3}.addPoint([1;0;0]); % body 3, point 4, just for looks


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


%% ASSEMBLE CONSTRAINT MATRIX (and add euler parameter normalization constraints)
sys.assembleConstraints()

%% ADD FORCES
sys.addGravityForces();

%% DYNAMICS ANALYSIS
timeStart = 0; %seconds
timeEnd =  10;
timeStep = 10^-2; %10^-3;

% using BDF2 with newton-raphson iterations
tic
stateNR = sys.dynamicsAnalysis(timeStart,timeEnd,timeStep,'NR');
timeNR = toc;
%save('stateNR_A9P2.mat','stateNR')
disp('done with dynamics analysis.')
disp(['For step-size h=' num2str(timeStep) ' seconds,'])
disp(['Time to compute Newton-Raphson Solution: ' num2str(timeNR) ' seconds'])

% % using BDF2 with Modified-newton iterations
% tic
% stateMN = sys.dynamicsAnalysis(timeStart,timeEnd,timeStep,'MN');
% timeMN = toc;
% %save('stateMN_A9P2.mat','stateMN')
% disp('done with dynamics analysis.')
% disp(['For step-size h=' num2str(timeStep) ' seconds,'])
% disp(['Time to compute Modified-Newton Solution: ' num2str(timeMN) ' seconds'])

% % using BDF2 with quasi-newton iterations
% tic
% stateQN = sys.dynamicsAnalysis(timeStart,timeEnd,timeStep,'QN');
% timeQN = toc;
% %save('stateQN_A9P2.mat','stateQN')
% disp('done with dynamics analysis.')
% disp(['For step-size h=' num2str(timeStep) ' seconds,'])
% disp(['Time to compute Quasi-Newton Solution: ' num2str(timeQN) ' seconds'])
return

%% PLOTS - be sure to run all three newton-methods above before running the plots.

time = timeStart:timeStep:timeEnd;

% kinematics for pendulum, the Oprime frame, BODY 3
rddotOprime3_NR = zeros(length(stateNR),3); % preallocate for speed
rddotOprime3_MN = zeros(length(stateMN),3); % preallocate for speed
rddotOprime3_QN = zeros(length(stateQN),3); % preallocate for speed
for i = 1:length(stateNR)
    rddotOprime3_NR(i,1) = stateNR{i}.rddot(1); % xddot
    rddotOprime3_NR(i,2) = stateNR{i}.rddot(2); % yddot
    rddotOprime3_NR(i,3) = stateNR{i}.rddot(3); % zddot
    rddotOprime3_MN(i,1) = stateMN{i}.rddot(1);
    rddotOprime3_MN(i,2) = stateMN{i}.rddot(2); 
    rddotOprime3_MN(i,3) = stateMN{i}.rddot(3); 
    rddotOprime3_QN(i,1) = stateQN{i}.rddot(1);
    rddotOprime3_QN(i,2) = stateQN{i}.rddot(2); 
    rddotOprime3_QN(i,3) = stateQN{i}.rddot(3);
end

% differences in translational acceleration
errorAccel_QN = rddotOprime3_NR - rddotOprime3_QN;
errorAccel_MN = rddotOprime3_NR - rddotOprime3_MN;

% plot differences for MODIFIED NEWTON
figure 
fig = gcf;
fig.Color = [1 1 1]; % set background color to white
st1 = suptitle('A9P2 MODIFIED-NEWTON: acceleration error');
set(st1,'FontSize',20,'FontWeight','bold')

subplot(3,1,1)
plot(time,errorAccel_MN(:,1),'r') % xddot component
xlabel('Time (sec)')
ylabel('Difference (m/s^2)')
legend('Accel X')

subplot(3,1,2)
plot(time,errorAccel_MN(:,2),'g') % yddot component
xlabel('Time (sec)')
ylabel('Difference (m/s^2)')
legend('Accel Y')

subplot(3,1,3)
plot(time,errorAccel_MN(:,3),'b') % zddot component
xlabel('Time (sec)')
ylabel('Difference (m/s^2)')
legend('Accel Z')
hold off

% plot differences for QUASI NEWTON
figure 
fig = gcf;
fig.Color = [1 1 1]; % set background color to white
st2 = suptitle('A9P2 QUASI-NEWTON: acceleration error');
set(st2,'FontSize',20,'FontWeight','bold')

subplot(3,1,1)
plot(time,errorAccel_QN(:,1),'r') % xddot component
xlabel('Time (sec)')
ylabel('Difference (m/s^2)')
legend('Accel X')

subplot(3,1,2)
plot(time,errorAccel_QN(:,2),'g') % yddot component
xlabel('Time (sec)')
ylabel('Difference (m/s^2)')
legend('Accel Y')

subplot(3,1,3)
plot(time,errorAccel_QN(:,3),'b') % zddot component
xlabel('Time (sec)')
ylabel('Difference (m/s^2)')
legend('Accel Z')
hold off

%%%%%%% iteration plots
%% number of iterations to converge
nIterations_NR = zeros(length(stateNR),1); % preallocate for speed
nIterations_MN = zeros(length(stateMN),1); % preallocate for speed
nIterations_QN = zeros(length(stateQN),1); % preallocate for speed
for i = 1:length(stateNR)
    nIterations_NR(i) = stateNR{i}.nIterations; 
    nIterations_MN(i) = stateMN{i}.nIterations; 
    nIterations_QN(i) = stateQN{i}.nIterations; 
end

% plot number of iterations to converge
figure 
fig = gcf;
fig.Color = [1 1 1]; % set background color to white
st = suptitle('A9P2b Iterations to converge');
set(st,'FontSize',20,'FontWeight','bold')
hold on
plot(time,nIterations_NR,'gO','MarkerSize',12)
plot(time,nIterations_QN,'bx-')
plot(time,nIterations_MN,'r.-')

hold off
xlabel('Time (sec)')
ylabel('Number of iterations')
lh = legend('Newton-Raphson','Quasi-Newton','Modified-Newton')
set(lh,'FontSize',15)

%% play animation of the dynamics analysis
plot.animateSystem(sys,stateNR);


