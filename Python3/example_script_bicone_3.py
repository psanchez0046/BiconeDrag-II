# This is an example of a script that calls postprocessingBiconeCPC
import os
import bicfuns_3

# Geometry parameters
h = 0.022# distance between interface and cup bottom [m]
R1 = 0.0340# bicone radius [m]
R = 0.04# cup radius [m]
# Parameters of the dynamics rheometer 
inertia = 0.0000242019# system (rotor + bicone) inertia [Kg·m^2]
b = 3.2e-8# frictional torque coefficient [Kg·m^2·rad/s]
# Mesh parameters
n = 500# Subintervals in r direction
m = 250# Subintervals in z direction
# Subphase physical parameters
rho_bulk = 1000# subphase density [Kg/m^3]
eta_bulk = 1e-3# subphase viscosity [Pa·s]
# Iterative scheme parameters
iteMax = 100# maximum number of iterations
tolMin = 0.00001# threshold tolerance
# Input/output data
colIndexAR = 1# ordinal number of the data of the column that contains the modulus of the amplitude ratio
colIndexDelta = 2# ordinal number of the data of the column that contains the modulus of the amplitude ratio
colIndexFreq = 0# ordinal number of the data of the column that contains the modulus of the amplitude ratio
inputFilepath=os.getcwd()# input filepath
outputFilepath=os.getcwd()# output filepath
# Execute postprocessingBiconeCPC with the specified input data
GData,etasData,bouData,ARcalcData,deltaARcalcData,iterationsTimesData,iterationsData,timeElapsedTotal=bicfuns_3.postprocessingBiconeCPC(h,R1,R,inertia,b,n,m,rho_bulk,eta_bulk,iteMax,tolMin,colIndexAR,colIndexDelta,colIndexFreq,inputFilepath,outputFilepath)
