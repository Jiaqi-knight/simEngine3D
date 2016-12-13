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

sys.body{1}.setMass(0.5); %kg
J = zeros(3); 
J(1,1) = 0.0004;
J(2,2) = 0.004;
J(3,3) = 0.004;
sys.body{1}.setInertia(J);

% body 3, slider
r_body = [0.2; 0; 0];
sys.addBody('free',r_body);

sys.body{1}.setMass(0.12); %kg
J = zeros(3); 
J(1,1) = 0.0001;
J(2,2) = 0.0001;
J(3,3) = 0.0001;
sys.body{1}.setInertia(J);

% body 4, ground
sys.addBody('ground') 

%% DEFINE POINTS ON THE BODIES %%
% these are used to specify heads and tails of vectors

% crank
sys.body{1}.addPoint([    0;    0; 0]); % body 1, point 1
sys.body{1}.addPoint([ 0.01;    0; 0]); % body 1, point 2
sys.body{1}.addPoint([    0; 0.08; 0]); % body 1, point 3

% connecting rod
sys.body{2}.addPoint([     0;    0;    0]); % body 2, point 1
sys.body{2}.addPoint([ -0.15;    0;    0]); % body 2, point 2
sys.body{2}.addPoint([  0.15;    0;    0]); % body 2, point 3
sys.body{2}.addPoint([  0.15; 0.01;    0]); % body 2, point 4
sys.body{2}.addPoint([  0.15;    0; 0.01]); % body 2, point 5

% slider
sys.body{3}.addPoint([    0;    0;    0]); % body 2, point 1
sys.body{3}.addPoint([ 0.20;    0;    0]); % body 2, point 2
sys.body{3}.addPoint([ 0.21;    0;    0]); % body 2, point 3
sys.body{3}.addPoint([    0; 0.01;    0]); % body 2, point 4
sys.body{3}.addPoint([    0;    0; 0.01]); % body 2, point 5
sys.body{3}.addPoint([ 0.01;    0;    0]); % body 2, point 6

% ground
sys.body{4}.addPoint([    0;    0;    0]); % body 4, point 1
sys.body{4}.addPoint([    0; 0.01;    0]); % body 4, point 2
sys.body{4}.addPoint([    0;    0; 0.01]); % body 4, point 3
sys.body{4}.addPoint([    0;  0.1; 0.12]); % body 4, point 4
sys.body{4}.addPoint([ 0.20;    0;    0]); % body 4, point 5
sys.body{4}.addPoint([ 0.21;    0;    0]); % body 4, point 6
sys.body{4}.addPoint([ 0.20; 0.01;    0]); % body 4, point 7



%% PLOT THE SYSTEM in 3D %%
sys.plot(0,0.05) % plot with reference frames
view(121,16);
axis equal

%% DEFINE CONSTRAINTS AMONG THE BODIES %%

% KINEMATIC CONSTRAINTS

% revolute joint: ground to crank (table 10.2.1)
sys.addConstraint('rj',sys.body{4},4,1,2,1,3,sys.body{1},1,1,2); % consistent

% spherical joint: crank to connecting rod (table 10.2.2)
sys.addConstraint('sj',sys.body{1},3,sys.body{2},2);

% revolute-cylindrical joint: connecting rod to slider (table 10.2.3)
sys.addConstraint('rcj',sys.body{2},3,4,3,5,3,sys.body{3},2,3,2);

% translational joint: slider to ground (table 10.2.4)
sys.addConstraint('tj',sys.body{3},1,1,4,1,5,sys.body{4},5,5,7,5,6); % consistent

% distance constraint: connecting rod to slider (table 10.2.5)
sys.addConstraint('d',sys.body{2},3,sys.body{3},6,0.01); % consistent

% add euler parameter normalization constraints to system
sys.addEulerParamConstraints(); 

%% ASSEMBLE CONSTRAINT MATRIX 

%q_assembled = sys.assemblyAnalysis(1.2)
%save('q_assembled.mat',q_assembled)
%sys.assembleConstraints(q_assembled)
return
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


