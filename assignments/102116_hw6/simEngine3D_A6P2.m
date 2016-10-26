
% Filename: simEngine3D_A6P2.m
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
sys.addBody('ground') % body 1, ground

t = 0;
L = 2;
theta = pi/4*cos(2*t); %radians
R = utility.R2(pi/2)*utility.R3(theta); % initial pendulum rotation
p = utility.A2p(R); % euler parameter for pendulum LRF rotation

sys.addBody('free',[0;L*sin(theta);-L*cos(theta)],p) % body 2

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

% revolute joint with ground body
sys.addConstraint('rj',sys.body{1},1,1,2,1,3,sys.body{2},3,1,2)

% add driving constraint
t = 0;
f = cos((pi*cos(2*t))/4); 
fdot = (pi*sin(2*t)*sin((pi*cos(2*t))/4))/2;
fddot = pi*cos(2*t)*sin((pi*cos(2*t))/4) - (pi^2*sin(2*t)^2*cos((pi*cos(2*t))/4))/4;

sys.addConstraint('dp1',sys.body{1},1,3,sys.body{2},1,4,f,fdot,fddot) % unit length vectors

% add euler parameter normalization constraints
sys.addConstraint('p_norm')


%% ASSEMBLE CONSTRAINT MATRIX 
sys.assembleConstraints()

%% DISPLAY CONSTRAINT PROPERTIES
disp('For t=0,')
disp('FULL CONSTRAINT MATRIX:')
sys.phi

disp('JACOBIAN OF FULL CONSTRAINT MATRIX:')
sys.phi_q

disp(' ')
disp('NU:')
sys.nu

disp(' ')
disp('GAMMA:')
sys.gammaHat

