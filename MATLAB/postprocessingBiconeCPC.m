function [GData,etasData,bouData,ARcalcData,deltaARcalcData,iterationsTimesData,iterationsData,timeElapsedTotal]=postprocessingBiconeCPC(h,R1,R,inertia,b,N,M,rho_bulk,eta_bulk,iteMax,tolMin,colIndexAR,colIndexDelta,colIndexFreq,inputFilepath,outputFilepath)
timerTotal = tic;
% This code calculates the rheological properties from amplitude ratio
% and phase lag data of a rotational rheometer with a bicone fixture

% INPUTS (all units are in the International System)
% %Geometry parameters:
% h: distance between interface and cup bottom [m]
% R1: bicone radius [m]
% R: cup radius [m]
% % Parameters of the Rheometer dynamics:
% inertia: system (rotor + bicone) inertia [Kg·m^2]
% b: frictional torque coefficient [Kg·m^2·s]
% %Mesh parameters:
% N: Subintervals in r direction
% M: Subintervals in z direction
% % Subphase physical parameters:
% rho_bulk: subphase density [Kg/m^3]
% eta_bulk: subphase viscosity [Pa·s]
% %Iterative scheme parameters:
% iteMax: maximum number of iterations
% tolMin: threshold tolerance
% % Input data
% colIndexAR: ordinal number of the data of the column that contains the modulus of the amplitude ratio
% colIndexDelta: ordinal number of the data of the column that contains the modulus of the amplitude ratio
% colIndexFreq: ordinal number of the data of the column that contains the modulus of the amplitude ratio
% inputFilepath: input filepath
% outputFilepath: output filepath

% OUTPUTS(optional)
% timeElapsedTotal: total execution time
% iterationsTimesData: time elapsed in processing each line of data
% iterationsData: number of iterations to process each line of data
% bouData: Boussinesq number for each line of each experimental data file
% etasData: Surface viscoelasticity for each line of data
% gData: Complex dynamic surface moduli for each line of data

h_adim = h/R;
% delta_r=1/N;% mesh spacing in r
delta_z = h_adim/M;% mesh spacing in z
R1_adim = R1/R;

expFilenames = GetFilenames('_exp.txt',inputFilepath);

% Initializing optional output data
iterationsTimesData = cell([1 length(expFilenames)]);
iterationsData = cell([1 length(expFilenames)]);
bouData = cell([1 length(expFilenames)]);
etasData = cell([1 length(expFilenames)]);
GData = cell([1 length(expFilenames)]);
ARcalcData = cell([1 length(expFilenames)]);
deltaARcalcData = cell([1 length(expFilenames)]);

% Prompting the number of input data files to process
fprintf('Experiment files: %s \n\n', num2str(length(expFilenames)))

% for-end looping on each input data file (*_exp.txt file)
for ite = 1:length(expFilenames)
    fprintf('Analyzing file %s...\n', num2str(ite))
    expFileData = importdata(char(expFilenames(ite)));
    % Prompting the number of lines to process in the file
    [filasexpFiledata, ~] = size(expFileData);
    fprintf('Lines: %s\n\n', num2str(filasexpFiledata))    
    % for-end looping on each line under analysis
    timeElapsedIT = zeros(filasexpFiledata, 1);
    Bou_final = zeros(filasexpFiledata, 1);
    lambda_final = zeros(filasexpFiledata, 1);
    ARcalc_final = zeros(filasexpFiledata, 1);
    errorAR_final = zeros(filasexpFiledata, 1);
    delta_AR_final = zeros(filasexpFiledata, 1);
    frec = expFileData(:, colIndexFreq);
    for lin = 1:filasexpFiledata        
        omegarad = 2*pi*frec(lin);
        % Displaying the number of the line under analysis
        fprintf('Analyzing line %s from file %s...\n', num2str(lin), num2str(ite))
        Re = rho_bulk*omegarad*R*R/eta_bulk;% Reynolds number
        ARexp = expFileData(lin,colIndexAR)*(cos(expFileData(lin,colIndexDelta))+1i*sin(expFileData(lin,colIndexDelta)));
        % Intializing variables...
        lambda = 1;
        Bou = [];
        ARcalc = [];
        errorAR = [];
        Bou(lambda) = 0;
        timerVal = tic;
        % Solving the Navier-Stokes equation
        g = solve_NS_bicono(Re, Bou(lambda), N, M, R1_adim, delta_z);
        Nb = floor(N*R1_adim);
        % Calculating subphase drag integral by the compound trapezium rule
%         integral=((M*R*R*R)/(2*h*N*N*N))*(sum(((2:Nb)-1)'.*((2:Nb)-1)'*2.*(g((2:Nb))-g((2:Nb)+(N+1))))+Nb*Nb*(g(Nb+1)-g((Nb+1)+(N+1))));
        integral = ((M*R*R*R)/(h*N*N*N))*(trapz(((2:Nb+1)-1)'.*((2:Nb+1)-1)'.*(g((2:Nb+1))-g((2:Nb+1)+(N+1)))));
        C = 1i*omegarad*2*pi*R1*eta_bulk;
        Tsub = -C*integral;% Subphase drag
        Tsurf = C*R1*R*Bou(lambda)*(R1*(N/R)*(g(Nb+2)-g(Nb+1))-1);% Surface drag
        ARcalc(lambda) = -Tsub - Tsurf - inertia*omegarad*omegarad + 1i*b*omegarad;
        % Calculating first value of the relative error 
        errorAR(lambda) = abs((ARcalc(lambda))/(ARexp)-1);
        % while loop performing the successive iterations
        while lambda<iteMax && errorAR(lambda)>tolMin 
            lambda = lambda+1;
            Bou(lambda) = (1i*omegarad*2*pi*R1*eta_bulk*integral-ARexp-inertia*omegarad*omegarad+1i*b*omegarad)/....
                (1i*omegarad*2*pi*R1*R1*R*eta_bulk*(R1*(N/R)*(g(Nb+2)-g(Nb+1))-1));       
            % Solving the Navier-Stokes equation
            g = solve_NS_bicono(Re, Bou(lambda), N, M, R1_adim, delta_z);
            Nb = floor(N*R1_adim);
            % Calculating subphase drag integral by the compound trapezium rule
%             integral=((M*R*R*R)/(2*h*N*N*N))*(sum(((2:Nb)-1)'.*((2:Nb)-1)'*2.*(g((2:Nb))-g((2:Nb)+(N+1))))+Nb*Nb*(g(Nb+1)-g((Nb+1)+(N+1))));
            integral = ((M*R*R*R)/(h*N*N*N))*(trapz(((2:Nb+1)-1)'.*((2:Nb+1)-1)'.*(g((2:Nb+1))-g((2:Nb+1)+(N+1)))));
            Tsub = -C*integral;% Subphase drag
            Tsurf = C*R1*R*Bou(lambda)*(R1*(N/R)*(g(Nb+2)-g(Nb+1))-1);% Surface drag
            ARcalc(lambda) = -Tsub - Tsurf - inertia*omegarad*omegarad + 1i*b*omegarad;
            % Calculating successive values of the relative error
            errorAR(lambda) = abs((ARcalc(lambda))/(ARexp)-1);            
        end  
        % Exit of the iterative process
        if lambda == iteMax               
            fprintf('Limit reached\n') % Maximum # of iterations reached
        else
            fprintf('OK convergence at %s iterations\n', num2str(lambda)) % Iterations have converged!!!
        end    
        timeElapsedIT(lin) = toc(timerVal);
        % Displaying the time used in the iterative process for each line
        fprintf('Iterative process time = %s s\n\n', num2str(timeElapsedIT(lin)))
        Bou_final(lin) = Bou(lambda);
        ARcalc_final(lin) = ARcalc(lambda);
        delta_AR_final(lin) = atan(imag(ARcalc(lambda))/real(ARcalc(lambda)));
        if delta_AR_final(lin) < 0
            delta_AR_final(lin) = delta_AR_final(lin)+pi;
        end    
        errorAR_final(lin) = errorAR(lambda);
        lambda_final(lin) = lambda;   
    end
    % Calculating variables depending on the Boussinesq number
    eta_s_final = Bou_final*R*eta_bulk;% converged viscoelasticity
    G_complex = 1i*R*eta_bulk*2*pi*frec.*Bou_final;% converged dynamic surface moduli    
    % Exporting results to the output data file
    results = [frec real(G_complex) imag(G_complex) real(eta_s_final) imag(eta_s_final) real(Bou_final) imag(Bou_final) abs(ARcalc_final) delta_AR_final timeElapsedIT lambda_final];
    if ispc
        dlmwrite(strcat(outputFilepath,'\',strrep(char(expFilenames(ite)),'exp','out')),results,'precision','%.15g','delimiter','\t')
    else
        dlmwrite(strcat(outputFilepath,'/',strrep(char(expFilenames(ite)),'exp','out')),results,'precision','%.15g','delimiter','\t')
    end    
    % Optional output data
    iterationsTimesData{ite} = timeElapsedIT;
    iterationsData{ite} = lambda_final;
    bouData{ite} = Bou_final;
    etasData{ite} = eta_s_final;
    GData{ite} = G_complex;
    ARcalcData{ite} = abs(ARcalc_final);
    deltaARcalcData{ite} = delta_AR_final;
end
timeElapsedTotal = toc(timerTotal);
fprintf('Total postprocessing program time = %s s\n', num2str(timeElapsedTotal))
end