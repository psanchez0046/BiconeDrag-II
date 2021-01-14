% This is an example of a script that calls postprocessingBiconeCPC.m
close all
clear, clc
% Geometry parameters
h = 0.022;% distance between interface and cup bottom [m]
R1 = 0.0340;% bicone radius [m]
R = 0.04;% cup radius [m]
% Parameters of the rheometer dynamics
inertia = 0.0000242019;% system (rotor + bicone) inertia [Kg·m^2]
b = 3.2e-8;% frictional torque coefficient [Kg·m^2·rad/s]
% Mesh parameters
N = 1000;% Subintervals in r direction
M = N/2;% Subintervals in z direction
% Subphase physical parameters
rho_bulk = 1000;% subphase density [Kg/m^3]
eta_bulk = 1e-3;% complex subphase viscosity [Pa·s]
% Iterative scheme parameters
iteMax = 100;% maximum number of iterations
tolMin = 1e-5;% threshold tolerance
% Input/output data
colIndexAR = 2;% ordinal number of the data of the column that contains the modulus of the amplitude ratio
colIndexDelta = 3;% ordinal number of the data of the column that contains the modulus of the amplitude ratio
colIndexFreq = 1;% ordinal number of the data of the column that contains the modulus of the amplitude ratio
inputFilepath = pwd;% input filepath
outputFilepath = pwd;% output filepath
% Execute postprocessingBiconeCPC.m with the specified input data
[GData,etasData,bouData,ARcalcData,deltaARcalcData,iterationsTimesData,iterationsData,timeElapsedTotal] = postprocessingBiconeCPC(h,R1,R,inertia,b,N,M,rho_bulk,eta_bulk,iteMax,tolMin,colIndexAR,colIndexDelta,colIndexFreq,inputFilepath,outputFilepath);
