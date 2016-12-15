% Filename: simEngine3D_n4BarMechanism.m
% Author:   Samuel Acuña
% Date:     14 Nov 2016
%
% About:    
% Driver file for the simEngine3D framework. Defines the system, the
% bodies in the system, and the constraints in the system.
% 
% This driver specifies the N four-bar mechanism model:
% http://lim.ii.udc.es/mbsbenchmark/dist/A02/A02_specification.xml
% (see ME751_f2016 slide 36 on lecture 09/28/16)

clear; clc; close all

% access directory with simEngine in it.
currentpath = cd('..'); cd('..');
addpath(pwd); % add grandparent folder to path
cd(currentpath); % return to initial folder

%% DEFINE SYSTEM %%

sys = system3D(); %initialize system

%% DEFINE BODIES IN THE SYSTEM %%

N = 2; % N four-bar linkages
l_link = 1; % meters. Length of a link
m_link = 1; %kg. Mass of a link
J_link = 1/12*m_link*l_link^2; % inertia of long, slender rod.

% add grounded bodies (bodies 1:N+1)
r_link = [0;0;0];
for i = 1:N+1
    sys.addBody('ground',r_link);
    r_link = r_link + [l_link; 0;0];
end

% add driver links (bodies N+2:2N+2)
r_link = [0;l_link/2;0];
R = utility.R3(pi/2); % initially rotate to keep x axis along length of link
p = utility.A2p(R); % euler parameter
for i = 1:N+1
    sys.addBody('free',r_link,p);
    r_link =  r_link + [l_link; 0;0];
end

% add coupler links (bodies 2N+3:3N+2)
r_link = [l_link/2;l_link;0];
for i = 1:N
    sys.addBody('free',r_link);
    r_link =  r_link + [l_link; 0;0];
end

% set link masses and inertia
J = zeros(3); % all links should be rod along x axis
J(3,3) = J_link; % benchmark only has inertia term in zz
for i = sys.nGroundBodies+1:sys.nBodies
    sys.body{i}.setMass(m_link);
    sys.body{i}.setInertia(J);
end

%% DEFINE POINTS ON THE BODIES %%
% these are used to specify heads and tails of vectors

for i = 1:sys.nGroundBodies
    sys.body{i}.addPoint([0;0;0]); % body i, point 1
    sys.body{i}.addPoint([l_link/4;0;0]); % body i, point 2
    sys.body{i}.addPoint([0;l_link/4;0]); % body i, point 3
end
for i = sys.nGroundBodies+1:sys.nBodies
    sys.body{i}.addPoint([0;0;0]); % body i, point 1, origin
    sys.body{i}.addPoint([-l_link/2;0;0]); % body i, point 2, -x
    sys.body{i}.addPoint([l_link/2;0;0]); % body i, point 3, x
    sys.body{i}.addPoint([0;l_link/4;0]); % body i, point 4, +y
    sys.body{i}.addPoint([0;0;l_link/4]); % body i, point 5, +z
end


%% PLOT THE SYSTEM in 3D %%
%sys.plot(1,0.2) % plot with reference frames
sys.plot(0,0.2) % plot without reference frames
view(0,90);
axis equal

%% DEFINE CONSTRAINTS AMONG THE BODIES %%

% KINEMATIC CONSTRAINTS
% revolute joints with ground bodies and driver links
for i = 1:N+1 % ground index
    driverLink_index = N+1+i;
    sys.addConstraint('rj',sys.body{i},1,1,2,1,3,sys.body{driverLink_index},2,1,5);
end

% revolute joints with driver links and coupler link
% note: to avoid overconstraints, only between driver link and coupler link pointing to the right.
for i = 1:N 
    driverLink_index = N+1+i;
    couplerLink_index = 2*N+2+i;
    sys.addConstraint('rj',sys.body{driverLink_index},3,1,3,1,4,sys.body{couplerLink_index},2,1,5);
end

% final CD joints between driver links and coupler link pointing to the left
for i = 1:N 
    driverLink_index = N+2+i;
    couplerLink_index = 2*N+2+i;
    sys.addConstraint('cd','x',sys.body{driverLink_index},3,sys.body{couplerLink_index},3);
    sys.addConstraint('cd','y',sys.body{driverLink_index},3,sys.body{couplerLink_index},3);
end

% EULER PARAMETER NORMALIZATION CONSTRAINTS
sys.addEulerParamConstraints(); 

%% ASSEMBLE CONSTRAINT MATRIX
sys.assembleConstraints();

%% INITIAL CONDITIONS 

% from initial velocity at top left corner 
vel0 = 1; %m/s
w0 = [0;0;-1]; % rad/2
sys.body{N+2}.rdot = [vel0/2;0;0]; %initial velocity of driver link body
sys.body{N+2}.pdot = utility.w2pdot(w0,sys.body{N+2}.p);


sys.setInitialVelocities(); % find initial velocities of constrained bodies

%% ADD FORCES
sys.g = [0;-9.81;0];
sys.addGravityForces();

%% DYNAMICS ANALYSIS
timeStart = 0; %seconds
timeEnd =  10;
timeStep = 10^-3; %-3;

tic
state = sys.dynamicsAnalysis(timeStart,timeEnd,timeStep,'QN');
timeQN = toc;

%save('state_N2_h3_n4BarMechanism.mat','state')
disp('done with dynamics analysis.')
disp(['For step-size h=' num2str(timeStep) ' seconds, and N = ' num2str(N) ','])
disp(['Time to compute Quasi-Newton Solution: ' num2str(timeQN) ' seconds'])


%% ANIMATION OF SYSTEM
%play animation of the dynamics analysis
plot.animateSystem(sys,state,[0 90],10);


%% COMPARE WITH BENCHMARK


% pull benchmark data
benchmark = load('benchmark_A02_solution_data.txt');
benchmark_time = benchmark(:,1);
benchmark_X = benchmark(:,2);
benchmark_Y = benchmark(:,3);

%% pull simEngine3D data:
time = timeStart:timeStep:timeEnd;
sBar = [l_link/2;0;0]; % local position of point B0 on leftmost driver link
X_B0 = zeros(length(state),1); % preallocate for speed
Y_B0 = zeros(length(state),1); % preallocate for speed

% derive point B0
for i = 1:length(state)
    r = state{i}.r(1:3); % position of leftmost driver link
    p = state{i}.p(1:4); % euler parameter of leftmost driver link
    A = utility.p2A(p);
    rB0 = r + A*sBar; % location of point B0
    X_B0(i) = rB0(1);
    Y_B0(i) = rB0(2);
end


%% plot position comparison
bm_index = find(benchmark_time == timeEnd);
figure
fig = gcf;
fig.Color = [1 1 1]; % set background color to white


hold on
plot(time,X_B0,'g')
plot(time,Y_B0,'c--')
plot(benchmark_time(1:bm_index),benchmark_X(1:bm_index),'r');
plot(benchmark_time(1:bm_index),benchmark_Y(1:bm_index),'b--');
hold off
lh1 = legend('X_B_0 simEngine3D','Y_B_0 simEngine3D','X_B_0 benchmark','Y_B_0 benchmark');
set(lh1,'FontSize',12)
xlabel('Time (sec)','FontSize',12)
ylabel('Displacement (m)','FontSize',12)
title('A02\_four\_bar\_mechanism','FontSize',12)
grid on



%% CHECK DIFFERENCES IN POSITION TO BENCHMARK
error_x = X_B0(1:50:end) - benchmark_X(1:bm_index);
error_y = Y_B0(1:50:end) - benchmark_Y(1:bm_index);

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
disp(['Comparing ' num2str(bm_index) ' datapoints,'])
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





