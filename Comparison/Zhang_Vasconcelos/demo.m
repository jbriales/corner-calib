%% Demo script: 
%  To use lccMinSol and lccZhang see the structure of the variables 
%  'Planes' and 'Points'. The output calibrations are T and Tz.
%
%  Type 'help lccMinSol' and 'help lccZhang' for details

clear all;
close all;
clc;

addpath(genpath('./MinimalSolution'));
addpath(genpath('./Zhang'));
addpath('./Optimization');
addpath('./Rotation');

% load input data
load('./DemoData.mat');

%% Minimal Solution

disp('Method using minimal solution');
[T, r, inl, Ti, ri] = lccMinSol(Planes, Points);
fprintf('\n\n');

%% Zhang

disp('Method using Zhang algorithm');
[Tz, rz, inlz, Tiz, riz] = lccZhang(Planes, Points);
fprintf('\n\n');

