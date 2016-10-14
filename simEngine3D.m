
% Filename: simEngine3D.m
% Author:   Samuel Acuña
% Date:     11 Oct 2016
% About:    
% Driver file for the simEngine3D framework. Defines the system, the
% bodies in the system, and the constraints and forces in the system. Then
% calls for the solution of the equations of motion for the system.

clear; clc; close all

%% DEFINE SYSTEM %%

sys = system3D(); %initialize system

p = utility.A2p(utility.R1(45)); % euler parameter for 45 deg rot about X axis

%% DEFINE BODIES IN THE SYSTEM %%
sys.addBody([1;1;1],p) % body 1
sys.addBody([1;1;0],p) % body 2
sys.addBody([0;1;0],p,[],[],[],1) % body 3, ground

%% DEFINE POINTS ON THE BODIES %%
sys.body{1}.addPoint([1;0;0]); % body 1, point 1
sys.body{1}.addPoint([0;0;0]); % body 1, point 2
sys.body{2}.addPoint([0;1;0]); % body 2, point 1
sys.body{2}.addPoint([0;0;0]); % body 2, point 2
sys.body{3}.addPoint([0;0;0]); % body 3, point 1
sys.body{3}.addPoint([0;0;1]); % body 3, point 2

%% PLOT THE SYSTEM %%
sys.plot(1) % plot with reference frames
% sys.plot()

%% DEFINE CONSTRAINTS AMONG THE BODIES %%

sys.addConstraint('dp1',sys.body{1},1,2,sys.body{2},1,2)
sys.addConstraint('dp1',sys.body{2},1,2,sys.body{3},1,2)
sys.addConstraint('cd','x',sys.body{2},1,sys.body{1},2)

sys.cons{1} % dp1, 2 bodies
sys.cons{2} % dp1, body with ground
sys.cons{3} % cd 2 bodies
