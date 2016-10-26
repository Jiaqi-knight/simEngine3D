
% Filename: simEngine3D_A6P1.m
% Author:   Samuel Acuña
% Date:     14 Oct 2016
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
p = utility.A2p(utility.R1(45)); % euler parameter for 45 deg rot about X axis
sys.addBody('free',[1;1;1],p) % body 1
sys.addBody('free',[1;1;0],p) % body 2
sys.addBody('ground',[0;1;0],p) % body 3, ground

%% DEFINE POINTS ON THE BODIES %%
sys.body{1}.addPoint([1;0;0]); % body 1, point 1
sys.body{1}.addPoint([0;0;0]); % body 1, point 2
sys.body{2}.addPoint([0;1;0]); % body 2, point 1
sys.body{2}.addPoint([0;0;0]); % body 2, point 2
sys.body{3}.addPoint([0;0;0]); % body 3, point 1
sys.body{3}.addPoint([0;0;1]); % body 3, point 2

%% PLOT THE SYSTEM in 3D %%
sys.plot(1) % plot with reference frames
% sys.plot()
view(52,28);

%% DEFINE CONSTRAINTS AMONG THE BODIES %%

sys.addConstraint('d',sys.body{2},2,sys.body{1},2,1)

PiQj = utility.dij(sys.body{1},2,sys.body{3},2);
C = norm(PiQj) % distance between points Pi and Qj
sys.addConstraint('d',sys.body{1},2,sys.body{3},2,C)

sys.addConstraint('dp2',sys.body{1},2,1,2,sys.body{2},2)
sys.addConstraint('dp2',sys.body{2},2,1,2,sys.body{3},1)


% DISPLAY CONSTRAINT PROPERTIES
sys.cons{1} % d, 2 bodies
sys.cons{2} % d, body with ground
sys.cons{3} % dp2, 2 bodies
sys.cons{4} % dp2, body with ground
