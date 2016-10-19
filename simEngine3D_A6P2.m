
% Filename: simEngine3D_A6P2.m
% Author:   Samuel Acuña
% Date:     18 Oct 2016
%
% About:    
% Driver file for the simEngine3D framework, for hw6 . Defines the system, the
% bodies in the system, and the constraints in the system.

clear; clc; close all

%% DEFINE SYSTEM %%

sys = system3D(); %initialize system

%% DEFINE BODIES IN THE SYSTEM %%
sys.addBody([0;0;0],[],[],[],[],[],1) % body 1, ground

t = 0;
L = 2;
theta = pi/4*cos(2*t); %radians
R = utility.R2(pi/2)*utility.R3(theta); % initial pendulum rotation
p = utility.A2p(R); % euler parameter for pendulum LRF rotation

sys.addBody([0;L*sin(theta);-L*cos(theta)],p) % body 2

%% DEFINE POINTS ON THE BODIES %%
% these are used to specify heads and tails of vectors
sys.body{1}.addPoint([0;0;0]); % body 1, point 1
sys.body{1}.addPoint([0;1;0]); % body 1, point 2
sys.body{1}.addPoint([0;0;1]); % body 1, point 3
sys.body{2}.addPoint([0;0;0]); % body 2, point 1
sys.body{2}.addPoint([0;0;1]); % body 2, point 2
sys.body{2}.addPoint([-2;0;0]); % body 2, point 3


%% PLOT THE SYSTEM in 3D %%
sys.plot(1) % plot with reference frames
% sys.plot() % plot without reference frames
view(98,12);
axis equal

%% DEFINE CONSTRAINTS AMONG THE BODIES %%

% revolute joint with ground body
sys.addConstraint('rj',sys.body{1},1,1,2,1,3,sys.body{2},3,1,2)

% add driving constraint
% theta = pi/4*cos(2*t); %radians

% add euler parameter normalization constraints
sys.addConstraint('p_norm')


%% ASSEMBLE CONSTRAINT MATRIX 

% next steps:
% jacobian of p_norm constraint
% driving constraints
% generate full constraint matrix
% display full jacobian, nu, gamma

%% DISPLAY CONSTRAINT PROPERTIES
disp('For t=0,')
disp('phi^F:')
sys.cons{1}.phi

disp(' ')
disp('nu:')
sys.cons{1}.nu

disp(' ')
disp('gamma:')
sys.cons{1}.nu

