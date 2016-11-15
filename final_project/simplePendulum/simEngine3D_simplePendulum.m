% Filename: simEngine3D_simplePendulum.m
% Author:   Samuel Acuña
% Date:     14 Nov 2016
%
% About:    
% Driver file for the simEngine3D framework. Defines the system, the
% bodies in the system, and the constraints in the system.
% 
% This driver specifies the simple pendulum benchmark model:
% http://lim.ii.udc.es/mbsbenchmark/dist/A01/A01_specification.xml
% (see ME751_f2016 slide 36 on lecture 09/28/16)

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

% although the benchmark webpage suggests that the initial configuraiton is
% [1,0,0], their benchmark data starts from [-1,0,0].
% body 2, in initial configuration. All units in SI
L = 1; %meter 
p = [0;0;0;-1]; % initial euler parameter 
sys.addBody('free',[-L;0;0],p); % add body 2
sys.body{2}.color = [0.4141, 0.3516, 0.8008]; %slateblue

% set body 2 mass
mass = 1; %kg
sys.body{2}.setMass(mass);

% set body 2 inertia
J = zeros(3); % point mass has no inertia 
sys.body{2}.setInertia(J);


%% DEFINE POINTS ON THE BODIES %%
% these are used to specify heads and tails of vectors
sys.body{1}.addPoint([0;0;0]); % body 1, point 1
sys.body{1}.addPoint([0.5;0;0]); % body 1, point 2
sys.body{1}.addPoint([0;0.5;0]); % body 1, point 3

sys.body{2}.addPoint([0;0;0]); % body 2, point 1
sys.body{2}.addPoint([-1;0;0]); % body 2, point 2
sys.body{2}.addPoint([0;0;0.5]); % body 2, point 3


%% PLOT THE SYSTEM in 3D %%
sys.plot(1) % plot with reference frames
view(0,90);
axis equal

%% DEFINE CONSTRAINTS AMONG THE BODIES %%

% KINEMATIC CONSTRAINTS
% revolute joint with ground body
sys.addConstraint('rj',sys.body{1},1,1,2,1,3,sys.body{2},2,1,3)


%% ASSEMBLE CONSTRAINT MATRIX (and add euler parameter normalization constraints)
sys.assembleConstraints()

%% ADD FORCES
sys.g = [0;-9.81;0];
sys.addGravityForces();

%% DYNAMICS ANALYSIS
timeStart = 0; %seconds
timeEnd =  10;
timeStep = 10^-2;

tic
state = sys.dynamicsAnalysis(timeStart,timeEnd,timeStep,'QN');
timeQN = toc;

%save('state_simplePendulum.mat','state')
disp('done with dynamics analysis.')
disp(['For step-size h=' num2str(timeStep) ' seconds,'])
disp(['Time to compute Quasi-Newton Solution: ' num2str(timeQN) ' seconds'])


%% COMPARE WITH BENCHMARK
% For step-size h=0.01 seconds,
% Time to compute Quasi-Newton Solution: 59.8413 seconds

% pull benchmark data
benchmark = load('benchmark_A01_solution_data.txt');
benchmark_time = benchmark(:,1);
benchmark_X = benchmark(:,2);
benchmark_Y = benchmark(:,3);

% pull simEngine3D data:
time = timeStart:timeStep:timeEnd;
X = zeros(length(state),1); % preallocate for speed
Y = zeros(length(state),1); % preallocate for speed
for i = 1:length(state)
    X(i) = state{i}.r(1); % x position of pendulum
    Y(i) = state{i}.r(2); % y position of pendulum
end

% plot position comparison
bm_index = find(benchmark_time == 10);
figure
fig = gcf;
fig.Color = [1 1 1]; % set background color to white

hold on
plot(time,X,'r')
plot(time,Y,'b--')
plot(benchmark_time(1:bm_index),benchmark_X(1:bm_index),'r^');
plot(benchmark_time(1:bm_index),benchmark_Y(1:bm_index),'bO');
hold off
lh1 = legend('X simEngine3D','Y simEngine3D','X benchmark','Y benchmark');
set(lh1,'FontSize',12)
xlabel('Time (sec)','FontSize',12)
ylabel('Displacement (m)','FontSize',12)
title('A01\_pendulum\_2D','FontSize',12)
grid on

%% CHECK DIFFERENCES IN POSITION TO BENCHMARK
error_x = X(1:50:end) - benchmark_X(1:bm_index);
error_y = Y(1:50:end) - benchmark_Y(1:bm_index);

% low precision check
lowPrecisionTol = 1.0e-1;
lpc_x = zeros(length(error_x),1);
lpc_y = zeros(length(error_y),1);
for i = 1:length(error_x)
    if abs(error_x(i)) > lowPrecisionTol
        lpc_x(i) = 1;
    end
    if abs(error_y(i)) > lowPrecisionTol
        lpc_y(i) = 1;
    end
end
nlpc_x = length(find(lpc_x == 1));
nlpc_y = length(find(lpc_y == 1));
disp(['X coordinate: There are ' num2str(nlpc_x) ' entries outside of allowable low precision error']);
disp(['Y coordinate: There are ' num2str(nlpc_y) ' entries outside of allowable low precision error']);

% high precision check
highPrecisionTol = 1.0e-3;
hpc_x = zeros(length(error_x),1);
hpc_y = zeros(length(error_y),1);
for i = 1:length(error_x)
    if abs(error_x(i)) > highPrecisionTol
        hpc_x(i) = 1;
    end
    if abs(error_y(i)) > highPrecisionTol
        hpc_y(i) = 1;
    end
end
nhpc_x = length(find(hpc_x == 1));
nhpc_y = length(find(hpc_y == 1));
disp(['X coordinate: There are ' num2str(nlpc_x) ' entries outside of allowable high precision error']);
disp(['Y coordinate: There are ' num2str(nlpc_y) ' entries outside of allowable high precision error']);




%% ANIMATION OF SYSTEM
% play animation of the dynamics analysis
plot.animateSystem(sys,state);


