function [ g ] = solve_NS_bicono( Re, Bo, N, M, R1_adim, delta_z )
% This function solves the Navier-Stokes equations with no-slip
% and Boussinesq_Scriven boundary conditions with second order centered 
% finite differences for the Bicone bob-cylindrical cup configuration

% INPUTS:
% Re: Reynolds number
% Bo: Boussinesq number
% N: subintervals on r
% M: subintervals on z
% R1_adim: bicone radius (nondimensional)
% delta_z: mesh spacing in z (nondimensional)

% OUTPUT:
% g: velocity field (nondimensional)

% # of Subintervals on the bicone lower surface
Nb=floor(N*R1_adim);
% Dimensions of the square coefficients matrix A 
dim=(N+1)*(M+1);
%% ACUMULATING INDEX AND VALUES
% tic
%bicone nodes
rows = 1:Nb+1;
cols = 1:Nb+1;
coefs = ones(1, Nb+1);
%simmetry axis and wall nodes
rows = [rows (1:M)*(N+1) (1:M-1)*(N+1)+1];
cols = [cols (1:M)*(N+1) (1:M-1)*(N+1)+1];
coefs = [coefs ones(1, 2*M-1)];
%ground nodes
rows = [rows (1:N+1)+M*(N+1)];
cols = [cols (1:N+1)+M*(N+1)];
coefs = [coefs ones(1, N+1)];

% interface nodes (BS + NS)
rows = [rows (Nb+2:N) (Nb+2:N) (Nb+2:N) (Nb+2:N)];
cols = [cols (Nb+2:N)-1 (Nb+2:N) (Nb+2:N)+1 (Nb+2:N)+(N+1)];
coefs = [coefs N*N*(1+2*Bo*(1/delta_z))*(1-(1./(2*((Nb+2:N)-1)))) ....
    -1i*Re - N*N*(1+2*Bo*(1/delta_z))*(2+(1./(((Nb+2:N)-1).*((Nb+2:N)-1)))) - 2*(1/delta_z)*(1/delta_z) ....
    N*N*(1+2*Bo*(1/delta_z))*(1+(1./(2*((Nb+2:N)-1)))) ....
    repmat(2*(1/delta_z)*(1/delta_z), [1 N+1-(Nb+2)])];

% internal nodes
[P, Q] = ndgrid(2:M, 2:N);
rows = [rows reshape(((P-1)*(N+1)+Q)', 1, []) ....
    reshape(((P-1)*(N+1)+Q)', 1, []) ....
    reshape(((P-1)*(N+1)+Q)', 1, []) ....
    reshape(((P-1)*(N+1)+Q)', 1, []) ....
    reshape(((P-1)*(N+1)+Q)', 1, [])];
cols = [cols reshape(((P-2)*(N+1)+Q)', 1, []) ....
    reshape(((P-1)*(N+1)+Q-1)', 1, []) ....
    reshape(((P-1)*(N+1)+Q)', 1, []) ....
    reshape(((P-1)*(N+1)+Q+1)', 1, []) ....
    reshape((P*(N+1)+Q)', 1, [])];
coefs = [coefs repmat((1/delta_z)*(1/delta_z), [1 (N-1)*(M-1)]) ....
    repmat(N*N*(1-(1./(2*((2:N)-1)))), [1 (M-1)]) ....
    repmat(-1i*Re-N*N*(2+(1./(((2:N)-1).*((2:N)-1))))-2*(1/delta_z)*(1/delta_z), [1 (M-1)]) ....
    repmat(N*N*(1+(1./(2*((2:N)-1)))), [1 (M-1)]) ....
    repmat((1/delta_z)*(1/delta_z), [1 (N-1)*(M-1)])];
% toc
%% FILL A
% tic
A = sparse(rows, cols, coefs);
% toc
%% FILL b with no-slip boundary conditions over the bicone lower surface
b = sparse(1:Nb+1, ones(1, Nb+1), (0:Nb)*(1/(R1_adim*N)), dim, 1);
%% Solving the linear system of equations
% tic
g = A\b;
% toc
end