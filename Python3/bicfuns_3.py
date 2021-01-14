import numpy as np
import glob
import sys
import time 
from scipy import sparse
import scipy.sparse.linalg
# Defining timer depending on the operative system
if sys.platform == 'win32':
    # on Windows, the best timer is time.clock
    timer = time.clock
else:
    # on most other platforms the best timer is time.time
    timer = time.time      
    
def solve_NS_bicono(Re, Bo, n, m, R1_adim, delta_z):
    # This function solves the Navier-Stokes equations with no-slip
    # and Boussinesq_Scriven boundary conditions with second order centered 
    # finite differences for the Bicone bob-cylindrical cup configuration
    
    # INPUTS:
    # Re: Reynolds number
    # Bo: Boussinesq number
    # N: subintervals on r
    # M: subintervals on z
    # R1_adim: bicone radius (nondimensional)
    # delta_z: mesh spacing in z (nondimensional)
	
    # OUTPUT:
    # g: velocity field (nondimensional)
	
    # # of Subintervals on the bicone lower surface
    nb = np.floor(n*R1_adim)
    nb = int(nb)
    # Dimensions of the square coefficients matrix A 
    dim = (m+1)*(n+1)
    #ACUMULATING INDEX AND VALUES FOR A
    # bicone nodes
    rows = np.array(np.arange(nb+1, dtype='int'))
    cols = np.array(np.arange(nb+1, dtype='int'))
    coefs = np.ones(nb+1, dtype='complex')
    # simmetry axis and wall nodes
    rows = np.append(rows, [n*np.arange(1, m+1, dtype='int')+np.arange(m), n*np.arange(1, m+1, dtype='int')+np.arange(m)+1])
    cols = np.append(cols, [n*np.arange(1, m+1, dtype='int')+np.arange(m), n*np.arange(1, m+1, dtype='int')+np.arange(m)+1])
    coefs = np.append(coefs, np.ones(2*m, dtype='complex'))
    # ground nodes
    rows = np.append(rows, [np.arange(1, n+1, dtype='int')+m*(n+1)])
    cols = np.append(cols, [np.arange(1, n+1, dtype='int')+m*(n+1)])
    coefs = np.append(coefs, np.ones(n, dtype='complex'))
    # interface nodes (Boussinesq_Scriven condition)
    rows = (np.append(rows, [np.arange(nb+1,n, dtype='int'), 
                             np.arange(nb+1,n, dtype='int'), 
                             np.arange(nb+1,n, dtype='int'), 
                             np.arange(nb+1,n, dtype='int')]))
    cols = (np.append(cols, [np.arange(nb+1,n, dtype='int')-1, 
                             np.arange(nb+1,n, dtype='int'), 
                             np.arange(nb+1,n, dtype='int')+1, 
                             np.arange(nb+1,n, dtype='int')+(n+1)]))    
    # coefs = (np.append(coefs, [Bo*n*n*(1-(1/(2*np.arange(nb+1,n, dtype='complex')))), 
    #                             -Bo*n*n*(2+(1/(np.arange(nb+1,n, dtype='complex')*np.arange(nb+1,n, dtype='complex'))))-(1/delta_z),
    #                             Bo*n*n*(1+(1/(2*np.arange(nb+1,n, dtype='complex')))),
    #                             np.repeat(1/delta_z, np.arange(nb+1,n).size)]))
    coefs = (np.append(coefs, [n*n*(1+2*Bo*(1/delta_z))*(1-(1/(2*np.arange(nb+1,n, dtype='complex')))), 
                                -1j*Re - n*n*(1+2*Bo*(1/delta_z))*(2+(1/(np.arange(nb+1,n, dtype='complex')*np.arange(nb+1,n, dtype='complex')))) - 2*(1/delta_z)*(1/delta_z),
                                n*n*(1+2*Bo*(1/delta_z))*(1+(1/(2*np.arange(nb+1,n, dtype='complex')))),
                                np.repeat(2*(1/delta_z)*(1/delta_z), np.arange(nb+1,n).size)]))
    
    #internal nodes
    P,Q = np.mgrid[1:m, 1:n]
    rows = (np.append(rows, [(P*(n+1)+Q).reshape(1, -1),
                             (P*(n+1)+Q).reshape(1, -1),
                             (P*(n+1)+Q).reshape(1, -1),
                             (P*(n+1)+Q).reshape(1, -1),
                             (P*(n+1)+Q).reshape(1, -1)]))
    cols = (np.append(cols, [((P-1)*(n+1)+Q).reshape(1, -1),
                             (P*(n+1)+Q-1).reshape(1, -1),
                             (P*(n+1)+Q).reshape(1, -1),
                             (P*(n+1)+Q+1).reshape(1, -1),
                             ((P+1)*(n+1)+Q).reshape(1, -1)]))    
    coefs = (np.append(coefs, [np.array([np.repeat((1/delta_z)*(1/delta_z), n-1)]*(m-1), dtype='complex').reshape(1, -1),
                               np.array([n*n*(1-(1/(2*np.arange(1,n))))]*(m-1), dtype='complex').reshape(1, -1),                               
                               np.array([-1j*Re-n*n*(2+(1/(np.arange(1,n)*np.arange(1,n))))-2*(1/delta_z)*(1/delta_z)]*(m-1), dtype='complex').reshape(1, -1),                               
                               np.array([n*n*(1+(1/(2*np.arange(1,n))))]*(m-1), dtype='complex').reshape(1, -1),                               
                               np.array([np.repeat((1/delta_z)*(1/delta_z), n-1)]*(m-1), dtype='complex').reshape(1, -1)]))
    # FILL A
    A = sparse.csc_matrix((coefs, (rows, cols)), shape=(dim, dim))
    #ACUMULATING INDEX AND VALUES FOR b
    bcoefs = np.arange(0, nb+1, dtype='float64')*(1/(R1_adim*n))
    brows = np.arange(0, nb+1, dtype='int')
    bcols = np.zeros(nb+1, dtype='int')
    # FILL b with no-slip boundary conditions over the bicone lower surface
    b = sparse.csc_matrix((bcoefs, (brows, bcols)), shape=(dim, 1))
	# Solving the linear system of equations
    g = sparse.linalg.spsolve(A, b)
    return g

def postprocessingBiconeCPC(h, R1, R, inertia, b, n, m, rho_bulk, eta_bulk, iteMax, tolMin, colIndexAR, colIndexDelta, colIndexFreq, inputFilepath, outputFilepath):
	# This code calculates the rheological properties from amplitude ratio
	# and phase lag data of a rotational rheometer with a bicone fixture
    
    # INPUTS (all units are in the International System)
    # # Geometry parameters:
    # h: distance between interface and cup bottom [m]
    # R1: bicone radius [m]
    # R: cup radius [m]
    # # Rheometer parameters :
    # inertia: system (rotor + bicone) inertia [Kg·m^2]
    # b: frictional torque coefficient [Kg·m^2·s]
    # # Mesh parameters:
    # N: Subintervals in r direction
    # M: Subintervals in z direction
    # # Subphase physical parameters:
    # rho_bulk: subphase density [Kg/m^3]
    # eta_bulk: subphase viscosity [Pa·s]
    # # Iterative scheme parameters:
    # iteMax: maximum number of iterations
    # tolMin: threshold tolerance
    # # Input data
    # colIndexAR: ordinal number of the data of the column that contains the modulus of the amplitude ratio
    # colIndexDelta: ordinal number of the data of the column that contains the modulus of the amplitude ratio
    # colIndexFrec: ordinal number of the data of the column that contains the modulus of the amplitude ratio
    # inputDataFilepath: input filepath
    # inputDataFilepath: output filepath
    
    # OUTPUTS(optional)
	# gData: Complex dynamic surface moduli for each line of data
	# etasData: Surface viscoelasticity for each line of data
	# bouData: Boussinesq number for each line of each experimental data file
	# ARcalcData: modulus of the amplitude ratio for each line of data
	# deltaARcalcData: phase of the amplitude ratio for each line pf data
	# iterationsTimesData: time elapsed in processing each line of data
    # iterationsData: number of iterations to process each line of data
	# timeElapsedTotal: total execution time
	
    timerTotal = timer()

    h_adim = h/R
    # delta_r=1/N #mesh spacing in r
    delta_z = h_adim/m# mesh spacing in z
    R1_adim = R1/R
    
    expFilenames = []
    if sys.platform == 'win32':
        for name  in glob.glob(inputFilepath + '\\*_exp.txt'):
            expFilenames.append(name)
    else:
        for name  in glob.glob(inputFilepath + '/*_exp.txt'):
            expFilenames.append(name)
        
    # Initializing optional output data
    iterationsTimesData = []
    iterationsData = []
    bouData = []
    etasData = []
    GData = []
    ARcalcData = []
    deltaARcalcData = []    
	# Displaying the number of input data files to process
    print('Experiment files: {} \n'.format(len(expFilenames)))
	
    # for-end looping on each input data file (*_exp.txt file)
    for ite in range(len(expFilenames)):
        print('Analyzing file {}...\n'.format(ite+1))
        if sys.platform == 'win32':
            expFileData = np.genfromtxt(expFilenames[ite].split('\\')[-1], delimiter='\t', dtype=None)
        else:
            expFileData = np.genfromtxt(expFilenames[ite].split('/')[-1], delimiter='\t', dtype=None)
            
        if len(expFileData.shape)==1:
            expFileData=np.array([expFileData])    
        # Displaying the number of lines to process in the file
        filasexpFiledata = int(expFileData.shape[0])
        
        print('Lines: {} \n'.format(filasexpFiledata))
        expFileData = expFileData.T
        # for-end looping on each line under analysis
        timeElapsedIT = np.zeros(filasexpFiledata, 'float64')
        Bou_final = np.zeros(filasexpFiledata, 'complex')
        lambda_final = np.zeros(filasexpFiledata, 'int')
        ARcalc_final = np.zeros(filasexpFiledata, 'complex')
        errorAR_final = np.zeros(filasexpFiledata, 'float64')
        delta_AR_final = np.zeros(filasexpFiledata, 'float64')
        
        for lin in range(filasexpFiledata):
            omegarad = 2*np.pi*expFileData[colIndexFreq][lin]
            # Displaying the number of the line under analysis
            print('Analyzing line {} from file {}...\n'.format(lin+1, ite+1))
            Re = rho_bulk*omegarad*R*R/eta_bulk# Reynolds number
            
            ARexp = expFileData[colIndexAR][lin]*(np.cos(expFileData[colIndexDelta][lin])+1j*np.sin(expFileData[colIndexDelta][lin]))
            ARclean = 1j*((np.pi*omegarad*R1*R1*R1*R1*eta_bulk)/(2*h))-inertia*omegarad*omegarad
            # Intializing variables...
            lmbd = 1
            Bou = []
            ARcalc = []
            errorAR = []
            
            Bou.append(((R-R1)/(1j*2*np.pi*omegarad*eta_bulk*R1*R1*R*R))*(ARexp-ARclean))
            timerVal = timer()
            # Solving the Navier-Stokes equation
            g = solve_NS_bicono(Re, Bou[-1], n, m, R1_adim, delta_z)
            nb = np.floor(n*R1_adim)
            nb = int(nb)
            # Calculating drag integral by the compound trapezium rule
            integral = ((m*R*R*R)/(2*h*n*n*n))*(np.sum(np.arange(1,nb)*np.arange(1,nb)*2*(g[1:nb]-g[1+(n+1):nb+(n+1)]))+nb*nb*(g[nb]-g[nb+(n+1)]))
            C = 1j*omegarad*2*np.pi*R1*eta_bulk
            Tsub = -C*integral# Subphase drag
            Tsurf = C*R1*R*Bou[-1]*(R1*(n/R)*(g[nb+1]-g[nb])-1)# Surface drag
            ARcalc.append(-Tsub-Tsurf-inertia*omegarad*omegarad+1j*b*omegarad)
            # Calculating first value of the relative error 
            errorAR.append(np.absolute((ARcalc[-1])/(ARexp)-1))
            # while loop performing the successive iterations     
            while lmbd < iteMax and errorAR[-1] > tolMin: 
                lmbd += 1
                Bou.append((1j*omegarad*2*np.pi*R1*eta_bulk*integral-ARexp-inertia*omegarad*omegarad+1j*b*omegarad)/(1j*omegarad*2*np.pi*R1*R1*R*eta_bulk*(R1*(n/R)*(g[nb+1]-g[nb])-1)))      
                # Solving the Navier-Stokes equation
                g = solve_NS_bicono( Re, Bou[-1], n, m, R1_adim, delta_z )
                # Calculating subphase drag integral by the compound trapezium rule
                integral = ((m*R*R*R)/(2*h*n*n*n))*(np.sum(np.arange(1,nb)*np.arange(1,nb)*2*(g[1:nb]-g[1+(n+1):nb+(n+1)]))+nb*nb*(g[nb]-g[nb+(n+1)]))
                Tsub = -C*integral# Subphase drag
                Tsurf = C*R1*R*Bou[-1]*(R1*(n/R)*(g[nb+1]-g[nb])-1)# Surface drag
                ARcalc.append(-Tsub-Tsurf-inertia*omegarad*omegarad+1j*b*omegarad)
                # Calculating successive values of the relative error
                errorAR.append(np.absolute((ARcalc[-1])/(ARexp)-1))            
            # Exit of the iterative process
            if lmbd == iteMax:
                print('Limit reached \n')# Maximum # of iterations reached
            else:
                print('OK convergence at {} iterations\n'.format(lmbd))# Iterations have converged!!!
                
            timeElapsedIT[lin] = timer() - timerVal
			# Displaying the time used in the iterative process for each line
            print('Iterative process time = {} s\n'.format(timeElapsedIT[lin]))
            Bou_final[lin] = Bou[-1]
            ARcalc_final[lin] = ARcalc[-1]
            delta_AR_final[lin] = np.arctan(np.imag(ARcalc[-1])/np.real(ARcalc[-1]))
            if delta_AR_final[lin] < 0:
                delta_AR_final[lin] = delta_AR_final[lin] + np.pi
                
            errorAR_final[lin] = errorAR[-1]
            lambda_final[lin] = lmbd;   
        # Calculating variables depending on the Boussinesq number
        eta_s_final = Bou_final*R*eta_bulk# converged viscoelasticity
        G_complex = 1j*2*np.pi*expFileData[colIndexFreq]*R*eta_bulk*Bou_final# converged dynamic surface moduli 
        
        # Exporting results to the output data file
        results = (np.array([expFileData[colIndexFreq], np.real(G_complex), np.imag(G_complex), np.real(eta_s_final),
                             np.imag(eta_s_final), np.real(Bou_final), np.imag(Bou_final), np.absolute(ARcalc_final),
                             delta_AR_final, timeElapsedIT, lambda_final]))
        results = results.T
        if sys.platform == 'win32':
            np.savetxt(outputFilepath+'\\'+expFilenames[ite].split('\\')[-1].replace('exp','out'), results, fmt='%.14f', delimiter='\t')
        else:
            np.savetxt(outputFilepath+'/'+expFilenames[ite].split('/')[-1].replace('exp','out'), results, fmt='%.14f', delimiter='\t')
        # Optional output data
        iterationsTimesData.append(timeElapsedIT)
        iterationsData.append(lambda_final)
        bouData.append(Bou_final)
        etasData.append(eta_s_final)
        GData.append(G_complex)
        ARcalcData.append(np.absolute(ARcalc_final))
        deltaARcalcData.append(delta_AR_final)
    timeElapsedTotal = timer() - timerTotal
    print('Total postprocessing program time = {} s\n'.format(timeElapsedTotal))
    return GData,etasData,bouData,ARcalcData,deltaARcalcData,iterationsTimesData,iterationsData,timeElapsedTotal