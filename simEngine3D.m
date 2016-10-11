
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



%% DEFINE BODIES IN THE SYSTEM %%
sys.addBody([10;0;0]) % body 1
sys.addBody([5;20;7]) % body 2

%% DEFINE POINTS ON THE BODIES %%
sys.body{1}.addPoint([2;2;2])

%% PLOT THE SYSTEM %%
% sys.plot(1) % plot with reference frames
sys.plot()

%% DEFINE CONSTRAINTS AMONG THE BODIES %%

