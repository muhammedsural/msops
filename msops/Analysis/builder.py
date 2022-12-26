import openseespy.opensees as ops
import math
from ..Units.Unit import Unit
import matplotlib.pyplot as plt
import opsvis as opsv
import numpy as np
import logging
from datetime import datetime
import sys
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
from scipy.integrate import cumtrapz
pd.reset_option('all')
import warnings
warnings.simplefilter(action='ignore', category=ResourceWarning)

 
class opsbuild:
    _version_ = "0.4.0"
    def __init__(self) -> None:
        pass
    # Define some local handy procs
    def addMaterial(MaterialType,InputArray):
        """
            MaterialType: Openseespy material tag
            InputArray  : Openseespy material input in array
        """
        global matTag
        matTag += 1
        ops.uniaxialMaterial(MaterialType,matTag,*InputArray)
        return matTag
    
    def modelbuild(ndm=2,ndf=3):
        """
        Modeli siler ve ndm=2 ndf=3 çerçeve model oluşturur.

        """
        # clear memory of past model definitions
        # ------------------------------------------------------------------------
        ops.wipe('all')                                    
        
        # Create ModelBuilder (with two-dimensions and 3 DOF/node)
        # ------------------------------------------------------------------------
        ops.model('basic', '-ndm', ndm, '-ndf', ndf)
    
    def calcnodalmass(node_dict,column_dict,beam_dict,colLenDict,beamLenDict,Acol,Abeam,distLoadbeam,nodalLoadcol=0,concDensity=23.53):
        """
        Elemanların zati ağırlıklarını düğüm noktalarına dağıtıp sözlük olarak değerleri döndürür ve total kütleyi verir
        Inputs :
                node_dict   : Düğüm noktalarının bilgisinin tutulduğu sözlük
                column_dict : Kolon düğüm noktalarının tutulduğu sözlük
                beam_dict   : Kiriş düğüm noktalarının bilgisinin tutulduğu sözlük
                colLenDict  : Kolon boylarının tutulduğu sözlük
                beamLenDict : Kiriş boylarının tutulduğu sözlük
                Acol        : Kolon alanı
                Abeam       : Kiriş alanı
                distLoadbeam: Kiriş üzerindeki yayılı yük modal analiz için
                nodalLoadcol: Kolon üzerindeki tekil yük modal analiz için default 0 kN
                concDensity : Beton yoğunluğu dedault 23.53 kN/m3

        Outputs:
                nodal_mass_dict
        """
        nodal_mass_dict ={}
        totalWeight = 0

        for node in node_dict.keys():
            nodal_mass_dict[node] = 0
            for colId in column_dict.keys():
                if node in column_dict[colId]:
                    wnode = round(concDensity*colLenDict[colId]*Acol/2,3)+ nodalLoadcol
                    nodal_mass_dict[node] += wnode/Unit.g
                    totalWeight += wnode
            for beamId in beam_dict.keys():
                if node in beam_dict[beamId]:
                    wnode = round(concDensity*beamLenDict[beamId]*Abeam/2,3) + distLoadbeam/2
                    nodal_mass_dict[node] += wnode/Unit.g
                    totalWeight += wnode
        print(f" total zati kütle  = {round(totalWeight,3)}")
        return nodal_mass_dict
    
    def createnodes(node_dict,nodal_mass_dict=None):
        """
        node_dict       : Düğüm noktalarının sözlüğü
        nodal_mass_dict : Düğüm noktalarının zati kütlesinin tutulduğu sözlük
        Elemanların zati ağırlıkları düğüm noktalarına eşit olarak paylaştırılır.
        """
        for node in node_dict.keys():
            if nodal_mass_dict is not None:
                ops.node(node,node_dict[node][0],node_dict[node][1],'-mass',nodal_mass_dict[node], nodal_mass_dict[node], 0.0)
            else:
                ops.node(node,node_dict[node][0],node_dict[node][1])

    def rayleigh(numMode=2,xDamp = 0.05,betaKcurr = 0.0,betaKinit = 0.0):
        """
        xDamp     = Default %5 damping ratio  percentage of critical damping
        alphaM    = M-prop. damping; D=alphaM*M
        betaKcurr = K-proportional damping; +beatKcurr*Kcurrent
        betaKinit = initial-stiffness proportional damping  +beatKinit*Kinit
        """
        # Opensees örneklerinde rayleigh damping kullanımı
        #Lambda = ops.eigen('-fullGenLapack',numMode) #eigenvalue mode 2
        #w1 = math.pow(Lambda[0],0.5)
        #w2 = math.pow(Lambda[1],0.5)
        # calculate damping parameters
        #zeta = 0.02                                       # percentage of critical damping
        #n = 10.0                                          # stiffness multiplier for rotational spring
        #a0  = xDamp*2.0*w1*w2/(w1 + w2)                   # mass damping coefficient based on first and second modes
        #a1  = xDamp*2.0/(w1 + w2)                         # stiffness damping coefficient based on first and second modes
        #a1_mod = a1*(1.0+n)/n                             # modified stiffness damping coefficient used for n modified elements. See Zareian & Medina 2010.
        #ops.rayleigh(a0 0.0 0.0 0.0)                      # assign mass proportional damping to structure (only assigns to nodes with mass)
        
        #Lambda = ops.eigen('-fullGenLapack',numMode) #eigenvalue mode 2
        
        w = ops.eigen(numMode)
        T1 = 2*np.pi/w[0]**0.5
        print(f'Fundamental period, T = {T1} sec')
        
        zi = xDamp # Mode 1
        zj = 0.02 # Mode 2
        wi = math.pow(w[0],0.5)
        wj = math.pow(w[1],0.5)
        A = np.zeros((2,2))
        A[0,0] = 1/wi; A[0,1] = wi
        A[1,0] = 1/wj; A[1,1] = wj
        b = np.zeros(2)
        b[0] = zi
        b[1] = zj
        x = np.linalg.solve(0.5*A,b)
        ops.rayleigh(x[0],betaKcurr,betaKinit,x[1])
        
        print(f" mass damping coefficient {round(x[0],4)} stiffness damping coefficient {round(x[1],4)} betaKcurr = {betaKcurr} betaKinit = {betaKinit}")
        """
        #method 1
        #w1 = math.pow(Lambda[0],0.5)
        #w2 = math.pow(Lambda[1],0.5)
        #a0  = xDamp*2.0*w1*w2/(w1 + w2)                           # mass damping coefficient based on first and second modes
        #a1  = xDamp*2.0/(w1 + w2)                                 # stiffness damping coefficient based on first and second modes
        #print(f"w1 = {w1} w2 = {w2}")
        #print(f" mass damping coefficient {round(a0,4)} stiffness damping coefficient {round(a1,4)} betaKcurr = {betaKcurr} betaKinit = {betaKinit}")
        #ops.rayleigh(round(a0,4),betaKcurr,betaKinit,round(a1,4))
        
        #method 2
        # define DAMPING--------------------------------------------------------------------------------------
        # apply Rayleigh DAMPING from $xDamp
        # D=$alphaM*M + $betaKcurr*Kcurrent + $betaKcomm*KlastCommit + $beatKinit*$Kinitial
        #alphaM = 0.0				# M-prop. damping; D = alphaM*M	
        #Omega = math.pow(Lambda[0],0.5)
        #betaKcomm = 2*(0.02/Omega)
        #print(f"mode 1 için Lambda = {Lambda} Omega = {Omega} betaKcomm = {betaKcomm}")
        #ops.rayleigh(alphaM,betaKcurr, betaKinit, betaKcomm) # RAYLEIGH damping
        """
    
    def BuildRCrectSection (SecID, HSec, BSec, coverH, coverB, coreID, coverID, steelID,
                        numBarsTop, barAreaTop, numBarsBot, barAreaBot, numBarsIntTot, barAreaInt, 
                        nfCoreY=20, nfCoreZ=20, nfCoverY=20, nfCoverZ=20, pflag=0):
        """
        Build fiber rectangular RC section, 1 steel layer top, 1 bot, 1 skin, confined core
        Define a procedure which generates a rectangular reinforced concrete section
        with one layer of steel at the top & bottom, skin reinforcement and a 
        confined core.
                by: Volkan Ozsarac, 2021
                    adapted from Silvia Mazzoni, 2006
        
        Formal arguments
            SecID:         tag for the section that is generated by this procedure
            HSec:          depth of section, along local-y axis
            BSec:          width of section, along local-z axis
            cH:            distance from section boundary to neutral axis of reinforcement
            cB:            distance from section boundary to side of reinforcement
            coreID:        material tag for the core patch
            coverID:       material tag for the cover patches
            steelID:       material tag for the reinforcing steel
            numBarsTop:    number of reinforcing bars in the top layer
            numBarsBot:    number of reinforcing bars in the bottom layer
            numBarsIntTot: TOTAL number of reinforcing bars on the intermediate layers, symmetric about z axis and 2 bars per layer-- needs to be an even integer
            barAreaTop:    cross-sectional area of each reinforcing bar in top layer
            barAreaBot:    cross-sectional area of each reinforcing bar in bottom layer
            barAreaInt:    cross-sectional area of each reinforcing bar in intermediate layer 
            nfCoreY:       number of fibers in the core patch in the y direction
            nfCoreZ:       number of fibers in the core patch in the z direction
            nfCoverY:      number of fibers in the cover patches with long sides in the y direction
            nfCoverZ:      number of fibers in the cover patches with long sides in the z direction
            pflag:         optional flag to plot fiber sections
            
                                y
                                ^
                                |     
                    ----------------------- --        |
                    |   o      o      o   | -- coverH |
                    |                     |           |
                    |   o             o   |           |
            z <---  |          +          |       HSec|
                    |   o             o   |           |
                    |                     |           |
                    |   o o o o o o o o   | --coverH  |
                    ----------------------- --        |
                    |--------Bsec---------|
                    |---|   coverB   |----|
        
                            y
                            ^
                            |    
                    ----------------------
                    |\      cover       /|
                    | \------Top-------/ |
                    |c|                |c|
                    |o|                |o|
            z <-----|v|      core      |v|  HSec
                    |e|                |e|
                    |r|                |r|
                    | /-------Bot------\ |
                    |/       cover      \|
                    ----------------------
                            Bsec
            
        
        Notes
            The core concrete ends at the NA of the reinforcement
            The center of the section is at (0,0) in the local axis system
        """
        # Define some parameters
        # ------------------------------------------------------------------------
        coverY = HSec/2.0                 # The distance from the section z-axis to the edge of the cover concrete -- outer edge of cover concrete
        coverZ = BSec/2.0                 # The distance from the section y-axis to the edge of the cover concrete -- outer edge of cover concrete
        coreY = coverY-coverH             # The distance from the section z-axis to the edge of the core concrete --  edge of the core concrete/inner edge of cover concrete
        coreZ = coverZ-coverB             # The distance from the section y-axis to the edge of the core concrete --  edge of the core concrete/inner edge of cover concrete
        numBarsInt = int(numBarsIntTot/2) # number of intermediate bars per side
        IntYStart = -coreY + 2*coreY/(numBarsInt+1) # Bottom (y-axis) coordinate for intermediate bars (assuming that they are uniformly distributed)
        IntYEnd = coreY - 2*coreY/(numBarsInt+1)    # Top (y-axix) coordinate for for intermediate bars (assuming that they are uniformly distributed)
        
        # Define the fiber section
        # ------------------------------------------------------------------------
        ops.section('Fiber', SecID)
        # Define the core patch
        #   patch('quad', matTag, numSubdivIJ, numSubdivJK,    *crdsI    ,     *crdsJ    ,    *crdsK    ,   *crdsL    )
        ops.patch('quad', coreID, nfCoreZ    , nfCoreY    , -coreY, coreZ, -coreY, -coreZ, coreY, -coreZ, coreY, coreZ)

        # Define the four cover patches
        ops.patch('quad', coverID, 2, nfCoverY, -coverY, coverZ, -coreY,   coreZ,   coreY,   coreZ,  coverY, coverZ)
        ops.patch('quad', coverID, 2, nfCoverY, -coreY, -coreZ,  -coverY, -coverZ,  coverY, -coverZ, coreY,  -coreZ)
        ops.patch('quad', coverID, nfCoverZ, 2, -coverY, coverZ, -coverY, -coverZ, -coreY,  -coreZ, -coreY,   coreZ)
        ops.patch('quad', coverID, nfCoverZ, 2,  coreY,  coreZ,   coreY,  -coreZ,   coverY, -coverZ, coverY, coverZ)    

        # Define reinforcing layers
        if numBarsInt != 0:
            ops.layer('straight', steelID, numBarsInt ,barAreaInt, IntYStart,  coreZ,  IntYEnd,  coreZ) # intermediate skin reinf. +z
            ops.layer('straight', steelID, numBarsInt ,barAreaInt, IntYStart,  -coreZ,  IntYEnd,  -coreZ) # intermediate skin reinf. -z
        ops.layer('straight', steelID, numBarsTop ,barAreaTop,  coreY,  coreZ,  coreY, -coreZ) # top layer reinfocement
        ops.layer('straight', steelID, numBarsBot ,barAreaBot, -coreY,  coreZ, -coreY, -coreZ) # bottom layer reinforcement
        
        # Plot the fiber section data
        # ------------------------------------------------------------------------
        
        fiber_sec = [['section', 'Fiber', SecID],
                        ['patch', 'quad', coreID, nfCoreZ, nfCoreY, -coreY, coreZ, -coreY, -coreZ, coreY, -coreZ, coreY, coreZ],
                        ['patch', 'quad', coverID, 2, nfCoverY, -coverY, coverZ, -coreY,   coreZ,   coreY,   coreZ,  coverY, coverZ],
                        ['patch', 'quad', coverID, 2, nfCoverY, -coreY, -coreZ,  -coverY, -coverZ,  coverY, -coverZ, coreY,  -coreZ],
                        ['patch', 'quad', coverID, nfCoverZ, 2, -coverY, coverZ, -coverY, -coverZ, -coreY,  -coreZ, -coreY,   coreZ],
                        ['patch', 'quad', coverID, nfCoverZ, 2,  coreY,  coreZ,   coreY,  -coreZ,   coverY, -coverZ, coverY, coverZ],
                        ['layer', 'straight', steelID ,numBarsInt ,barAreaInt, IntYStart,  coreZ,  IntYEnd,  coreZ],
                        ['layer', 'straight', steelID ,numBarsInt ,barAreaInt, IntYStart,  -coreZ,  IntYEnd,  -coreZ],
                        ['layer', 'straight', steelID ,numBarsTop ,barAreaTop,  coreY,  coreZ,  coreY, -coreZ],
                        ['layer', 'straight', steelID ,numBarsBot ,barAreaBot, -coreY,  coreZ, -coreY, -coreZ]
                        ]
        return fiber_sec
        
    def monotonic_test(dref, nSteps, Material,mass=1*Unit.tonne,name='core',plot=True):
        """
        Run monotonic loading test on a uniaxial material
        Args:
            Material:    List containg, Material properties Such as Material = ['Steel01', F_y, k_el, r_post] # Bilinear with kinematic hardening
            dref:        Reference displacement to which cycles are run
            nSteps:      Number of displacement increments
            mass :       mass of nodes default = 1
        Returns:
            Deformation: List containing deformation of material thoughout the analysis
            Forces:      List containing forces in material throughout the analysis
        """
        
        # Wipe any existing model
        # ------------------------------------------------------------------------
        ops.wipe()

        # Create ModelBuilder (with 1-dimension and 1 DOF/node)
        # ------------------------------------------------------------------------
        ops.model('basic', '-ndm', 1, '-ndf', 1)

        # Define nodes
        # ------------------------------------------------------------------------
        node1 = 1 # Tag for node 1 (fixed node)
        node2 = 2 # Tag for node 2 (free node)
        coord1 = 0.0 # 1 dimensional coordinate for node 1
        coord2 = 0.0 # 1 dimensional coordinate for node 2
        ops.node(node1, coord1)
        ops.node(node2, coord2)

        # Define single-point constraints
        # ------------------------------------------------------------------------
        ops.fix(node1, 1) # Fix node 1, 
        ops.fix(node2, 0) # release node 2 (this is optional, by default it is unrestrained)

        # Define the nodal mass
        # ------------------------------------------------------------------------
        
        ops.mass(node2, mass)

        # Define materials
        # ------------------------------------------------------------------------
        # spring
        spring_tag = 1
        spring_type = Material[0]
        spring_props = Material[1:]
        ops.uniaxialMaterial(spring_type, spring_tag, *spring_props)
        
        # low stiffness material (for numerical purposes)
        dummy_tag = 2 # tag for low stiffness material
        ops.uniaxialMaterial('Elastic', dummy_tag, 1e-5)
        
        # Define elements
        # ------------------------------------------------------------------------
        ops.element('zeroLength', spring_tag, node1, node2, "-mat", spring_tag, "-dir", 1)
        ops.element('zeroLength', dummy_tag, node1, node2, "-mat", dummy_tag, "-dir", 1)

        # Define the load pattern
        # ------------------------------------------------------------------------
        ops.timeSeries('Linear', 1)
        ops.pattern('Plain', 1, 1)
        ops.load(2, 1)

        # Set analysis settings
        # ------------------------------------------------------------------------
        # Wipe any previous analysis object
        ops.wipeAnalysis()

        # Convergence Test -- determines when convergence has been achieved.
        tol = 1.0e-8  # Set the tolerance (default)
        iterMax = 50  # Set the max bumber of iterations (default)
        pFlag = 0     # Optional print flag (default is 0). Valid options: 0-5
        nType = 2     # optional type of norm (default is 2). Valid options: 0-2
        ops.test('NormDispIncr', tol, iterMax, pFlag, nType)

        # SolutionAlgorithm -- determines the sequence of steps taken to solve the non-linear equation at the current time step
        ops.algorithm('Newton', '-initial')

        # DOF_Numberer -- determines the mapping between equation numbers and degrees-of-freedom
        ops.numberer('RCM')

        # SystemOfEqn/Solver -- within the solution algorithm, it specifies how to store and solve the system of equations in the analysis
        ops.system('BandGeneral')

        # Constraints handler: determines how the constraint equations are enforced in the analysis -- how it handles the boundary conditions/imposed displacements
        ops.constraints('Transformation')
        
        # Integrator -- determines the predictive step for time t+dt
        dU = dref/nSteps # displacement increment
        ops.integrator('DisplacementControl', 2, 1, dU)

        # AnalysisType -- defines what type of analysis is to be performed ('Static', 'Transient' etc.)
        ops.analysis('Static')

        # Initialize some parameters
        # ------------------------------------------------------------------------
        Force = [0]
        Deformation = [0]
        
        # Perform step by step analysis
        # ------------------------------------------------------------------------
        for l in range(0,nSteps,1):
            ok = ops.analyze(1)
            Force.append(ops.basicForce(spring_tag)[0])
            Deformation.append(ops.basicDeformation(spring_tag)[0])
            if ok !=0:
                print("DispControl Analysis is FAILED")
                print("Analysis failed at nSteps: %s" %(l))
                print('-------------------------------------------------------------------------')
                break
        outputs = [Deformation,Force]
        if plot:
            # Plot the results
            # ------------------------------------------------------------------------
            fig, ax = plt.subplots(figsize=(20,5))
            fig.subplots_adjust(bottom=0.15, left=0.2)
            ax.grid(True)
            ax.plot(outputs[0], outputs[1], label=f"{Material[0]} {name}")
            ax.axhline(0, color='black', lw=2)
            ax.set_xlabel('Deformation [m]')
            ax.set_ylabel('Force')
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
            plt.show()
        return outputs

    def cyclic_test(dref, numCycles, nSteps, Material,mass=1*Unit.tonne,name='core',plot=True):
        """
        Run cyclic loading test on a uniaxial material
        Args:
            Material:    List containg, Material properties
            dref:        Reference displacement to which cycles are run.
            numCycles:   No. of cycles. Valid options either 1,2,3,4,5,6
            nSteps:      Number of displacement increments per load reversal
        Returns:
            Deformation: List containing deformation of material thoughout the analysis
            Forces:      List containing forces in material throughout the analysis
        """
        
        # Wipe any existing model
        # ------------------------------------------------------------------------
        ops.wipe()

        # Create ModelBuilder (with 1-dimension and 1 DOF/node)
        # ------------------------------------------------------------------------
        ops.model('basic', '-ndm', 1, '-ndf', 1)

        # Define nodes
        # ------------------------------------------------------------------------
        node1 = 1 # Tag for node 1 (fixed node)
        node2 = 2 # Tag for node 2 (free node)
        coord1 = 0.0 # 1 dimensional coordinate for node 1
        coord2 = 0.0 # 1 dimensional coordinate for node 2
        ops.node(node1, coord1)
        ops.node(node2, coord2)

        # Define single-point constraints
        # ------------------------------------------------------------------------
        ops.fix(node1, 1) # Fix node 1, 
        ops.fix(node2, 0) # release node 2 (this is optional, by default it is unrestrained)

        # Define the nodal mass
        # ------------------------------------------------------------------------
        ops.mass(node2, mass)

        # Define materials
        # ------------------------------------------------------------------------
        # spring
        spring_tag = 1
        spring_type = Material[0]
        spring_props = Material[1:]
        ops.uniaxialMaterial(spring_type, spring_tag, *spring_props)
        
        # low stiffness material (for numerical purposes)
        dummy_tag = 2 # tag for low stiffness material
        ops.uniaxialMaterial('Elastic', dummy_tag, 1e-5)
        
        # Define elements
        # ------------------------------------------------------------------------
        ops.element('zeroLength', spring_tag, node1, node2, "-mat", spring_tag, "-dir", 1)
        ops.element('zeroLength', dummy_tag, node1, node2, "-mat", dummy_tag, "-dir", 1)

        # Define the load pattern
        # ------------------------------------------------------------------------
        ops.timeSeries('Linear', 1)
        ops.pattern('Plain', 1, 1)
        ops.load(2, 1)

        # Initialize some parameters
        # ------------------------------------------------------------------------
        Force = [0]
        Deformation = [0]
        
        # Set analysis settings
        # ------------------------------------------------------------------------
        # Wipe any previous analysis object
        ops.wipeAnalysis()

        # Convergence Test -- determines when convergence has been achieved.
        tol = 1.0e-8  # Set the tolerance (default)
        iterMax = 50  # Set the max bumber of iterations (default)
        pFlag = 0     # Optional print flag (default is 0). Valid options: 0-5
        nType = 2     # optional type of norm (default is 2). Valid options: 0-2
        ops.test('NormDispIncr', tol, iterMax, pFlag, nType)

        # SolutionAlgorithm -- determines the sequence of steps taken to solve the non-linear equation at the current time step
        ops.algorithm('Newton', '-initial')

        # DOF_Numberer -- determines the mapping between equation numbers and degrees-of-freedom
        ops.numberer('RCM')

        # SystemOfEqn/Solver -- within the solution algorithm, it specifies how to store and solve the system of equations in the analysis
        ops.system('BandGeneral')

        # Constraints handler: determines how the constraint equations are enforced in the analysis -- how it handles the boundary conditions/imposed displacements
        ops.constraints('Transformation')
        
        # Create the list of displacements
        if numCycles == 1:
            dispList = [dref, -2*dref, dref]
            dispNoMax = 3
        elif numCycles == 2:
            dispList = [dref, -2*dref, dref,
                        dref, -2*dref, dref]
            dispNoMax = 6
        elif numCycles == 3:
            dispList = [dref, -2*dref, dref,
                        dref, -2*dref, dref,
                        dref, -2*dref, dref]
            dispNoMax = 9
        elif numCycles == 4:
            dispList = [dref, -2*dref, dref,
                        dref, -2*dref, dref,
                        dref, -2*dref, dref,
                        dref, -2*dref, dref]
            dispNoMax = 12
        elif numCycles == 5:
            dispList = [dref, -2*dref, dref,
                        dref, -2*dref, dref,
                        dref, -2*dref, dref,
                        dref, -2*dref, dref,
                        dref, -2*dref, dref]
            dispNoMax = 15
        elif numCycles == 6:
            dispList = [dref, -2*dref, dref,
                        dref, -2*dref, dref,
                        dref, -2*dref, dref,
                        dref, -2*dref, dref,
                        dref, -2*dref, dref,
                        dref, -2*dref, dref]
            dispNoMax = 18
        else:
            print("ERROR: Value for numCycles not a valid choice. Choose between 1 and 6")
            print('-------------------------------------------------------------------------') 
        
        # Iterate for each load reversal
        for d in range(1,dispNoMax+1,1):
            
            # Integrator -- determines the predictive step for time t+dt
            dU = dispList[d-1]/nSteps # displacement increment
            ops.integrator('DisplacementControl', 2, 1, dU)

            # AnalysisType -- defines what type of analysis is to be performed ('Static', 'Transient' etc.)
            ops.analysis('Static')
            
            # Perform step by step analysis
            # ------------------------------------------------------------------------
            for l in range(0,nSteps,1):
                ok = ops.analyze(1)
                Force.append(ops.basicForce(spring_tag)[0])
                Deformation.append(ops.basicDeformation(spring_tag)[0])
                if ok !=0:
                    print("DispControl Analysis is FAILED")
                    print("Analysis failed at cycle: %s and nSteps: %s" %(d,l))
                    print('-------------------------------------------------------------------------')
                    break
        outputs = [Deformation,Force]
        if plot:
            # Plot the results
            # ------------------------------------------------------------------------
            fig, ax = plt.subplots(figsize=(15,10))
            fig.subplots_adjust(bottom=0.15, left=0.2)
            ax.grid(True)
            ax.plot(outputs[0], outputs[1], label=f"{Material[0]}{name}")
            ax.axhline(0, color='black', lw=2)
            ax.set_xlabel('Deformation [m]')
            ax.set_ylabel('Force')
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
            
        return outputs
    
    def PMinteraction(secTag, kappa,HSec,epsc=-0.0035):
        
        # Define two nodes at (0,0)
        ops.node(1, 0.0, 0.0)
        ops.node(2, 0.0, 0.0)

        # Fix all degrees of freedom except axial and bending
        ops.fix(1, 1, 1, 1)
        ops.fix(2, 0, 1, 0)
        
        # Define element
        #                               tag   ndI  ndJ  secTag
        ops.element('zeroLengthSection',  1,   1,   2,  secTag)

        ops.timeSeries('Constant',2)
        ops.pattern('Plain',2,2)
        ops.sp(2,1,-epsc + kappa*HSec/2) # Strain at the centroid
        ops.sp(2,3,kappa) # Curvature
        
        
        
        ops.constraints('Transformation')
        ops.system('UmfPack')

        ops.analysis('Static') 
        ops.analyze(1)

        ops.reactions()

        P = -ops.nodeReaction(2,1)
        M = ops.nodeReaction(2,3)
        
        return P,M

    def MomentCurvature(secTag : int, axialLoad : float, maxK : int,HSec : int,BSec : int, numIncr : int =100) -> list:
        """
        Moment curvature analysis
        
        Input
            secTag          :  Kesit id
            axial load      :  Kesite etkitilen eksenel yük
            maxK            :  Hedef dönme
            HSec            :  depth of section, along local-y axis
            BSec            :  width of section, along local-z axis
            numIncr         :  Number of analysis increments
        
        Output
            Load -> List
            Curvature -> List
        """
        y1 = HSec / 2
        z1 = BSec / 2


        # Define two nodes at (0,0)
        ops.node(1, 0.0, 0.0)
        ops.node(2, 0.0, 0.0)

        # Fix all degrees of freedom except axial and bending
        ops.fix(1, 1, 1, 1)
        ops.fix(2, 0, 1, 0)
        
        # Define element
        #                               tag   ndI  ndJ  secTag
        ops.element('zeroLengthSection',  1,   1,   2,  secTag)

        # Save element forces and nodal displacements, base reactions
        # ------------------------------------------------------------------------
        #ops.reactions('-dynamic', '-rayleigh') # Must call this command before using nodeReaction() command.
        
        # Define constant axial load
        ops.timeSeries('Constant', 1)
        ops.pattern('Plain', 1, 1)
        ops.load(2, axialLoad, 0.0, 0.0)

        
        # Define analysis parameters
        ops.integrator('LoadControl', 0.0)
        ops.system('SparseGeneral', '-piv')
        ops.test('NormUnbalance', 1e-9, 10)
        ops.numberer('Plain')
        ops.constraints('Plain')
        ops.algorithm('Newton')
        ops.analysis('Static')

        # Do one analysis for constant axial load
        ops.analyze(1)

        # Define reference moment
        ops.timeSeries('Linear', 2)
        ops.pattern('Plain',2, 2)
        ops.load(2, 0.0, 0.0, 1.0)

        # Compute curvature increment
        dK = maxK / numIncr
        # Use displacement control at node 2 for section analysis
        #   integrator('DisplacementControl', nodeTag, dof, incr, numIter=1, dUmin=incr, dUmax=incr)
        ops.integrator('DisplacementControl',    2   ,  3 ,  dK ,    1     ,    dK     ,    dK)
        step = 0
        ok = 0
        #loadf = 1.0        # This feature of disabling the possibility of having a negative loading has been included
        LoadFactor = [0] 
        Curvature = [0]
        
        while step <= numIncr and ok == 0:
            step +=1
            ok = ops.analyze(1)
            if ok == 0:
                loadf = ops.getTime()                                # get load factor
                LoadFactor.append(ops.getTime())                     # append loadfactor
                Curvature.append(HSec*ops.nodeDisp(2, 3)) # append moment
            else:
                print(f"analiz converge etmedi {step}/{numIncr}")
                break
        return LoadFactor,Curvature

    def modal_analys(numMode):
        """This function use only openseespy3.4.0.2"""
        lambdas = ops.eigen('-fullGenLapack', numMode) # eigenvalue analysis 
        ops.modalProperties("-print", "-file", "ModalReport.txt", "-unorm") # perform modal analysis, and print results
        for i in range(1,numMode+1):
            opsv.plot_mode_shape(i, 20) 
            
    def modal_analys2(numEigen, pflag=1, outname=None):
        """
        Details
        -------
            This script will return the modal properties of an OpenSeespy model.
        
        Information
        -----------
            Author: Volkan Ozsarac, Earthquake Engineering PhD Candidate
            Affiliation: University School for Advanced Studies IUSS Pavia
            e-mail: volkanozsarac@iusspavia.it
        
        References
        ----------
            Chopra, A.K. 2012. Dynamics of Structures: Theory and 
            Applications to Earthquake Engineering, Prentice Hall.

        Notes
        -----
            Total (activated) mass is obtained by summing the masses assigned to the
            unrestrained degrees of freedoms (DOFs). Thus, it should not be confused 
            with total mass assigned to all DOFs. Influence vectors for rotational 
            excitation are not correct at the moment, this addition remains as future work.
            Which reference point to use is not clear for rotational excitations.
            SAP2000 and Seismostruct use different reference points.
            
        Parameters
        ----------
        numEigen : int
            Number of eigenvalues to calculate.
        pflag    : int (1 or 0)
            flag to print output information on screen
        outname  : str, optional (The default is None)
            if not None and pFlag==1, the modal properties for the 
            first numEigen modes will be writtend into outname.csv.

        Returns
        -------
        T        : numpy.ndarray
            Period array for the first numEigen modes.
        Mratios  : dictionary
            Effective modal mass participation ratios for the first numEigen modes.
        Mfactors : dictionary
            Modal particpation factors for the first numEigen modes.
        Mtots    : dictionary
            Total activated masses.

        """

        #ops.wipeAnalysis()
        ops.numberer("Plain")
        ops.system('FullGeneral')
        ops.algorithm('Linear')
        ops.analysis('Transient')

        # Extract the Mass Matrix
        # Note that this is not the global mass matrix, but unrestrained part (Muu)
        ops.integrator('GimmeMCK',1.0,0.0,0.0)
        ops.analyze(1,0.0) 
        # Number of equations in the model
        N = ops.systemSize()         # Has to be done after analyze
        Mmatrix = ops.printA('-ret') # Or use op.printA('-file','M.out')
        Mmatrix = np.array(Mmatrix) # Convert the list to an array
        Mmatrix.shape = (N,N)       # Make the array an NxN matrix
        print( '\n************************************************************', \
            '\nExtracting the mass matrix, ignore the warnings...')
            
        # Determine maximum number of DOFs/node used in the system
        NDF = 0
        for node in ops.getNodeTags():
            temp = len(ops.nodeDOFs(node))
            if temp > NDF: NDF = temp

        DOFs = []       # List containing indices of unrestrained DOFs
        used = {}       # Dictionary with nodes and associated unrestrained DOFs
        ldict = {}      # Dictionary containing influence vectors
        Mratios = {}    # Dictionary containing effective modal masses ratios
        Mfactors = {}   # Dictionary containing modal participation factors
        for i in range(1,NDF+1):
            ldict[i] = np.zeros([N,1])
            Mratios[i] = np.zeros(numEigen)
            Mfactors[i] = np.zeros(numEigen)
            
        # Create the influence vectors, and get the unrestrained DOFs assigned to the nodes
        # TODO -1: The influence vectors are not correct in case of rotational excitations
        # One typical approach is to use center of mass on plane
        idx = 0                                     # Counter for unrestrained DOFs
        for node in ops.getNodeTags():               # Start iterating over each node
            used[node] = []                         # Unrestrain local DOF ids
            ndof = len(ops.nodeDOFs(node))           # Total number of DOFs assigned
            for j in range(ndof):                   # Iterate over each DOF
                temp = ops.nodeDOFs(node)[j]         # Get the global DOF id (-1 if restrained)
                if temp not in DOFs and temp >= 0:  # Check if this DOF is unrestrained and is not known before
                    DOFs.append(temp)               # Save the global id of DOF 
                    used[node].append(j+1)          # Save the local id of DOF 
                    ldict[j+1][idx,0] = 1           # Influence vectors for horizontal and vertical excitations
                    idx += 1                        # Increase the counter

        # This does not seem necessary when numberer is "Plain"
        # But lets reorganize the mass matrix anyway
        Mmatrix = Mmatrix[DOFs,:][:,DOFs]  

        # Calculate the total masses assigned to the unrestrained DOFs
        Mtots = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0}
        for i in range(1,NDF+1):
            Mtots[i] = (ldict[i].T@Mmatrix@ldict[i])[0,0]

        # Perform eigenvalue analysis
        ops.wipeAnalysis()
        listSolvers = ['-genBandArpack','-fullGenLapack','-symmBandLapack']
        ok = 1  
        for s in listSolvers:
            print("Using %s as solver..." % s[1:])
            try:
                eigenValues = ops.eigen(s,numEigen)
                catchOK = 0
                ok = 0
            except: 
                catchOK = 1
            
            if catchOK==0:
                for i in range(numEigen):
                    if eigenValues[i] < 0: 
                        ok = 1
                if ok==0: 
                    print('Eigenvalue analysis is completed.')
                    break
        if ok!=0:
            print("Error on eigenvalue something is wrong...")
            sys.exit()
        else:  
            Lambda = np.asarray(eigenValues)
            Omega = Lambda**0.5
            T = 2*np.pi/Omega
            frq = 1/T

        # Note: influence factors for rotational excitation is wrong! 
        # Obtain modal properties
        for mode in range(1,numEigen+1):
            idx = 0
            phi = np.zeros([N,1]) # Eigen vector
            for node in used:
                for dof in used[node]:
                    phi[idx,0]=ops.nodeEigenvector(node,mode,dof)
                    idx += 1
                    
            phi = phi/(phi.T@Mmatrix@phi)**0.5  # Normalize the eigen vector by modal mass
            Mn = phi.T@Mmatrix@phi              # Modal mass (should always be equal to 1)

            for j in range(1,NDF+1):
                if Mtots[j] != 0:                              # Check if any mass is assigned
                    Ln = phi.T@Mmatrix@ldict[j]                # Modal excitation factor
                    Mnstar = (Ln**2/Mn)[0,0]                   # Effective modal mass
                    Mfactors[j][mode-1] = Ln/Mn                # Modal participation factor
                    Mratios[j][mode-1] = (Mnstar/Mtots[j]*100) # Effective modal mass participation ratio [%]
        
        for j in range(1,7):
            try: 
                Mratios[j]
            except: 
                Mratios[j]  = np.zeros(numEigen)
                Mfactors[j] = np.zeros(numEigen)

        # TODO-1: Results are not correct for rotational excitation cases, for now ignore those.
        del Mratios[6], Mratios[5], Mratios[4]
        del Mfactors[6], Mfactors[5], Mfactors[4]

        # Calculate cumulative modal mass participation ratio
        sM1 = np.cumsum(Mratios[1]); sM2 = np.cumsum(Mratios[2]); sM3 = np.cumsum(Mratios[3])  
        
        # Print modal analysis results
        if pflag == 1:
            arguments = []
            arguments.append('Modal Periods and Frequencies')
            arguments.append('%4s|%8s|%10s|%12s|%12s' \
                % ('Mode', 'T [sec]','f [Hz]','\u03C9 [rad/sec]', '\u03BB [rad\u00b2/sec\u00b2]'))
            for mode in range(numEigen):      
                arguments.append('%4s|%8s|%10s|%12s|%12s' \
                    % ("{:.0f}".format(mode+1), "{:.4f}".format(T[mode]), "{:.3f}".format(frq[mode]), \
                        "{:.2f}".format(Omega[mode]), "{:.2f}".format(Lambda[mode])))
            arguments.append('Total Activated Masses')
            arguments.append('%8s|%8s|%8s' \
                % ('M\u2081','M\u2082','M\u2083'))
            arguments.append('%8s|%8s|%8s' \
                    % ( "{:.2f}".format(Mtots[1]), "{:.2f}".format(Mtots[2]), "{:.2f}".format(Mtots[3])))
            arguments.append('Modal Mass Participation Factors') 
            arguments.append('%4s|%7s|%7s|%7s' \
                % ('Mode','\u0393\u2081','\u0393\u2082','\u0393\u2083') )             
            for mode in range(numEigen):
                arguments.append('%4s|%7s|%7s|%7s' % ("{:.0f}".format(mode+1), \
                    "{:.3f}".format(Mfactors[1][mode]), "{:.3f}".format(Mfactors[2][mode]), "{:.3f}".format(Mfactors[3][mode])))  
            arguments.append('Effective Modal Mass Participation Ratios [%]') 
            arguments.append('%4s|%7s|%7s|%7s' \
                % ('Mode','U\u2081','U\u2082','U\u2083') )              
            for mode in range(numEigen):
                arguments.append('%4s|%7s|%7s|%7s' % ("{:.0f}".format(mode+1), \
                    "{:.3f}".format(Mratios[1][mode]), "{:.3f}".format(Mratios[2][mode]), "{:.3f}".format(Mratios[3][mode])))  
            arguments.append('Cumulative Effective Modal Mass Participation Ratios [%]') 
            arguments.append('%4s|%7s|%7s|%7s' \
                % ('Mode','\u2211U\u2081','\u2211U\u2082','\u2211U\u2083') )              
            for mode in range(numEigen):
                arguments.append('%4s|%7s|%7s|%7s' % ("{:.0f}".format(mode+1), \
                    "{:.3f}".format(sM1[mode]), "{:.3f}".format(sM2[mode]), "{:.3f}".format(sM3[mode])))  

            # To the screen
            arguments = '\n'.join(arguments); print(arguments)
    
            # To the .csv file
            if outname != None:
                with open(outname+'.csv','w', encoding='utf-32') as f:
                    f.write(arguments)

        return T, Mratios, Mfactors, Mtots
    
    def do_gravity(nSteps = 500,tol = 1.0e-8,iterMax = 10,pFlag = 0,nType = 2):
            """
            Procedure to carry out gravity-analysis: load-controlled static analysis
            Args:
                nSteps:       Number of steps to run analysis
                pflag:        Optional flag to print base reactions

            # Convergence Test -- determines when convergence has been achieved.
                tol = 1.0e-8  # Set the tolerance (default)
                iterMax = 50  # Set the max bumber of iterations (default)
                pFlag = 0     # Optional print flag (default is 0). Valid options: 0-5
                nType = 2     # optional type of norm (default is 2). Valid options: 0-2

            Returns:
                eleForces:    Dictionary containing element joint forces
                nodalDisps:   Dictionary containing nodal displacements

            """
            # Set analysis parameters
            # ------------------------------------------------------------------------
            # Wipe any previous analysis object
            ops.wipeAnalysis()
            
            ops.test('NormDispIncr', tol, iterMax, pFlag, nType)
            
            # SolutionAlgorithm -- determines the sequence of steps taken to solve the non-linear equation at the current time step
            ops.algorithm('Newton', '-initial')
            
            # DOF_Numberer -- determines the mapping between equation numbers and degrees-of-freedom
            ops.numberer('RCM') #AMD , Plain , RCM
            
            # SystemOfEqn/Solver -- within the solution algorithm, it specifies how to store and solve the system of equations in the analysis
            ops.system('BandGeneral') #FullGeneral ,BandGeneral
            
            # Constraints handler: determines how the constraint equations are enforced in the analysis -- how it handles the boundary conditions/imposed displacements
            ops.constraints('Transformation')
            
            # Integrator -- determines the predictive step for time t+dt
            dLambda = 1/nSteps # the load factor increment
            ops.integrator('LoadControl', dLambda)
            
            # AnalysisType -- defines what type of analysis is to be performed ('Static', 'Transient' etc.)
            ops.analysis('Static')
            
            # Perform the analysis
            # ------------------------------------------------------------------------
            ops.analyze(nSteps)
            
            
            # Maintain constant gravity loads and reset time to zero
            # ------------------------------------------------------------------------
            ops.loadConst('-time', 0.0)
            
            # Save element forces and nodal displacements, base reactions
            # ------------------------------------------------------------------------
            
            
            BaseReactions = 0
            for node in ops.getNodeTags():
                forces = np.array(ops.nodeReaction(node))
                BaseReactions += forces
        
            eleForces = {}
            for ele in ops.getEleTags():
                forces = ops.eleForce(ele)
                eleForces[ele] = np.array(forces)
            
            nodalDisps = {}
            for node in ops.getNodeTags():
                disps = ops.nodeDisp(node)
                nodalDisps[node] = np.array(disps)
                
            eleSection_top = {}
            eleSection_bot = {}
            eleSection_steel = {}
            
            for ele in ops.getEleTags():
                # Compression fiber
                fiber_stressStrain_top   = ops.eleResponse(ele, 'section', 'fiber', str(0.30/2),str(0), str(3), 'stressStrain' )
                # Tension fiber
                fiber_stressStrain_bot   = ops.eleResponse(ele, 'section', 'fiber',str(-0.30/2),str(0), str(3), 'stressStrain' )
                # Tension fiber - steel
                fiber_stressStrain_steel = ops.eleResponse(ele, 'section', 'fiber',str(-0.30/2),str(0),str(2) , 'stressStrain' )
                
                if len(eleSection_top.keys()) !=  len(ops.getEleTags()):
                    eleSection_top[ele]   = np.array(fiber_stressStrain_top) 
                    eleSection_bot[ele]   = np.array(fiber_stressStrain_bot)
                    eleSection_steel[ele] = np.array(fiber_stressStrain_steel)
                    continue
                else:
                    eleSection_top[ele]   = np.append(eleSection_top[ele],fiber_stressStrain_top)
                    eleSection_bot[ele]   = np.append(eleSection_bot[ele],fiber_stressStrain_bot)
                    eleSection_steel[ele] = np.append(eleSection_steel[ele],fiber_stressStrain_steel)
                    
            eleSection = [eleSection_top,eleSection_bot,eleSection_steel]
            return eleForces, nodalDisps, BaseReactions,eleSection

    def do_nspa(dmax, ctrlNode, ctrlDOF, nSteps,tol = 1.0e-8,iterMax = 300,pFlag = 0,nType = 2):
        """
            Procedure to carry out a non-cylic pushover of a model
            Args:
                dmax            :     Maximum displacement to run analysis
                ctrlNode        :     Node to control with the displacement integrator
                ctrlDOF         :     DOF the loading is applied
                nSteps          :     Number of steps
                tol             :     Set the tolerance (default)
                iterMax         :     Set the max bumber of iterations (default)
                pFlag           :     Optional print flag (default is 0). Valid options: 0-5
                nType           :     Optional type of norm (default is 0). Valid options: 0-2
                capacityCurve   :     Plot Capacity curve default is True
                outFiber        :     record fiber response and plotting for selected one element [visualazation,eleid,Hsec,cover,idCoverMat,idCoreMat,idSteel] default:[True,1,0.30,0.05,1,2,3]
                outNodalDisp    :       False
                animotions      :     False
                spyMatrix       :     False

            Returns:
                LoadFactor:   List containing load factors used throughout the analysis
                DispCtrlNode: List containing displacement of control node throughout the analysis  
                sectionStrainStress
                nodalDisps

        """
        # Set analysis parameters
        # ------------------------------------------------------------------------
        # Wipe any previous analysis object
        ops.wipeAnalysis()
        
        # Set the gravity loads to be constant & reset the time in the domain
        #ops.loadConst('-time', 0.0)
        # Save element forces and nodal displacements, base reactions
        # ------------------------------------------------------------------------
        #ops.reactions() # Must call this command before using nodeReaction() command.

        # Convergence Test -- determines when convergence has been achieved.
        ops.test('NormDispIncr', tol, iterMax, pFlag, nType)
        
        # SolutionAlgorithm -- determines the sequence of steps taken to solve the non-linear equation at the current time step
        ops.algorithm('Newton','-initial')
        
        # DOF_Numberer -- determines the mapping between equation numbers and degrees-of-freedom
        ops.numberer('RCM')
        
        # SystemOfEqn/Solver -- within the solution algorithm, it specifies how to store and solve the system of equations in the analysis
        ops.system('BandGeneral')
        
        # Constraints handler: determines how the constraint equations are enforced in the analysis -- how it handles the boundary conditions/imposed displacements
        ops.constraints('Transformation')
        
        # Integrator -- determines the predictive step for time t+dt
        dU = dmax/nSteps # Displacement increment 
        ops.integrator('DisplacementControl', ctrlNode, ctrlDOF, dU)
        
        # AnalysisType -- defines what type of analysis is to be performed ('Static', 'Transient' etc.)
        ops.analysis('Static')

        # Initialize some parameters
        # ------------------------------------------------------------------------
        ok = 0.0           # analysis result
        step = 1           # current step number
        loadf = 1.0        # This feature of disabling the possibility of having a negative loading has been included
        LoadFactor = [0]   # List containing load factors used throughout the analysis
        DispCtrlNode = [0] # List containing displacement of control node throughout the analysis  

        el_tags = ops.getEleTags()
        nels = len(el_tags)
        timeV = np.zeros(nSteps)
        eleDeform = np.zeros((nSteps, nels, 6))

        # Perform the analysis
        # ------------------------------------------------------------------------
        while step <= nSteps and ok == 0 and loadf > 0:
            
            ok = ops.analyze(1) # Run a step of the analysis
            
            if ok != 0: # If the analysis fails, we can change the strategy to achieve convergence
                print("~~~ Analysis did not converge at step: %d/%d, trying KrylovNewton ..." % (step,nSteps))
                ops.algorithm('KrylovNewton')
                ok = ops.analyze(1)
                if ok == 0: # If the analysis works, we can go back to previous strategy
                    print("~~~ That worked, back to regular Newton ...")
                    ops.algorithm('Newton', '-initial')

            if ok == 0: # If the analysis successful save some stuff
                print("~~~ Analysis converge at step: %d/%d"% (step,nSteps))
                loadf = ops.getTime()                                # get load factor
                LoadFactor.append(loadf)                     # append loadfactor
                DispCtrlNode.append(ops.nodeDisp(ctrlNode, ctrlDOF)) # append displacement
                maxBaseShear = max(LoadFactor)
                if maxBaseShear > 1.5*loadf:
                    ok = -3
                    print("Structure is probably collapse.. analysis stopped...")
                    break
                
            step += 1           # Update the current step number
                
        # Print the final status of the analysis
        # ------------------------------------------------------------------------
        if ok != 0:
            Analysis =  "Displacement Control Analysis is FAILED to converge at step: %d/%d" % (step,nSteps)
        else:
            Analysis = "Displacement Control Analysis is SUCCESSFUL"
        if loadf <= 0:
            Analysis = "Stopped because of Load factor below zero: %.f" % loadf
        print('------------------------------------------------------------------------')
        print(Analysis)

        return LoadFactor, DispCtrlNode     

    def define_timeseries(record,dt,tsTag,IDloadTag = 400,GMfactSeries=9.81,GMdirection=1,GMFactPattern=1):
        """
        Zaman serisi tanımlar
 
        record       : ivme serisi   ; Time series record file name or signals list
        dt           : Kayıt zaman aralığı
        tsTag        : Zaman serisinin tag'i 
        IDloadTag    : Pattern Tag default 400,
        GMfact       : Zaman serisi faktör çarpanı default 1 bu değer fonksiyonda g ile çarpılıp zaman serisi oluşturulur; Groundmotion factor (A factor multiply load factors)
        GMdirection  : Zaman serisi doğrultusu default 1
        """
        
        if type(record) == str:
            print("dosya yolu")
            ops.timeSeries('Path', tsTag  , '-dt', dt, '-filePath', record, '-factor', GMfactSeries)
        elif type(record) == list:
            print("liste")
            ops.timeSeries('Path', tsTag, '-dt', dt, '-values', *record, '-factor', GMfactSeries) # time series object
        
        ops.pattern('UniformExcitation', IDloadTag,GMdirection, '-accel', tsTag,'-factor',GMFactPattern)# pattern object

    def do_nrha(tNode, bNode, Dt, Tmax, Dc):
        """
        Function to perform Nonlinear Response History Analysis
        Args:
            tNode:    top nodes for the drift calculation
            bNode:    bottom nodes for the drift calculation
            Dt:       Analysis time step
            Tmax:     Length of the record (including padding of 0's)
            Dc:       Drift capacity for pier drift (%)

        Returns:
            mdrft:    Peak Interstorey Drifts [%]
            
        Note:
            About Newmark Integrator;
            gamma = 1/2, beta = 1/4 --> Average Acceleration Method; Unconditionally stable
            gamma = 1/2, beta = 1/6 --> Linear Acceleration Method; Conditionally stable: Dt / T > 0.551   
        """
        # Set analysis settings
        # ------------------------------------------------------------------------
        # Wipe any previous analysis object
        ops.wipeAnalysis()
        
        # Convergence Test -- determines when convergence has been achieved.
        tol = 1.0e-8  # Set the tolerance (default)
        iterMax = 50  # Set the max bumber of iterations (default)
        pFlag = 0     # Optional print flag (default is 0). Valid options: 0-5
        nType = 2     # optional type of norm (default is 2). Valid options: 0-2
        ops.test('NormDispIncr', tol, iterMax, pFlag, nType)
        
        # SolutionAlgorithm -- determines the sequence of steps taken to solve the non-linear equation at the current time step
        ops.algorithm('Newton', '-initial')
        
        # DOF_Numberer -- determines the mapping between equation numbers and degrees-of-freedom
        ops.numberer('RCM')
        
        # SystemOfEqn/Solver -- within the solution algorithm, it specifies how to store and solve the system of equations in the analysis
        ops.system('BandGeneral')
        
        # Constraints handler: determines how the constraint equations are enforced in the analysis -- how it handles the boundary conditions/imposed displacements
        ops.constraints('Transformation')
        
        # Integrator -- determines the predictive step for time t+dt
        gamma = 0.5   # Set Newmark gamma coefficient
        beta = 0.25   # Set Newmark beta coefficient
        ops.integrator('Newmark', gamma, beta)
        
        # AnalysisType -- defines what type of analysis is to be performed ('Static', 'Transient' etc.)
        ops.analysis('Transient')
        
        # Initialize some parameters
        # ------------------------------------------------------------------------
        cIndex = 0         # Initially define the control index (-1 for non-converged, 0 for stable, 1 for global collapse)
        controlTime = 0.0  # Start the controlTime
        ok = 0             # Set the convergence to 0 (initially converged)
        mflr = 0           # Set the initial pier collapse location
        h = []             # storey heights
        mdrft = []         # the interstorey drift values
        
        for i in range(len(tNode)):
            # Find the coordinates of the nodes in Global Y (2)
            top2 = ops.nodeCoord(tNode[i], 2)
            bot2 = ops.nodeCoord(bNode[i], 2)
            dist = top2 - bot2
            h.append(dist)     # Current pier height
            mdrft.append(0.0)  # We will populate the lists with zeros initially
            if dist == 0: print("WARNING: Zerolength found in drift check")
                
        # Perform the analysis
        # ------------------------------------------------------------------------
        while cIndex == 0 and controlTime <= Tmax and ok == 0:
            ok = ops.analyze(1, Dt)      # Run a step of the analysis
            controlTime = ops.getTime()  # Update the control time

            # If the analysis fails, we can change some stuff to achieve convergence
            if ok != 0:
                print("~~~ Analysis did not converge at %.2f, trying KrylovNewton ..." % controlTime)
                ops.algorithm('KrylovNewton')
                ok = ops.analyze(1)
                if ok == 0:
                    print("~~~ That worked .. back to regular Newton")
                    ops.algorithm('Newton', '-initial')
                else: # Bye bye...  Failed to converge, exit the analysis.
                    cIndex = -1
                    
            # If the analysis is successful some results can be stored, 
            if ok == 0:        
                # Let's get the peak interstorey drift ratios at each floor (mdrft).
                for i in range(len(tNode)):
                    tNode_disp = ops.nodeDisp(tNode[i], 1)              # Current top node disp in dof 1
                    bNode_disp = ops.nodeDisp(bNode[i], 1)              # Current bottom node disp in dof 1
                    cHt = h[i]                                          # Current interstorey height
                    cdrft = 100.0 * abs(tNode_disp - bNode_disp) / cHt  # Current interstorey drift in dof 1 [%]
                    if cdrft >= mdrft[i]: mdrft[i] = cdrft              # Update the interstorey drift in dof 1 [%]
                
                # Stop the analysis if the interstorey drift limit is exceeded 
                if any(i >= Dc for i in mdrft): 
                    cIndex = 1 # Set the state of the model to local collapse (=1)

        # Print the final status of the analysis
        # ------------------------------------------------------------------------
        if cIndex == -1:
            Analysis = "Analysis is FAILED to converge at %.3f of %.3f" % (controlTime, Tmax)
        if cIndex == 0:
            text = ["\nInterstorey drift: %.2f%% at floor %d" % (mdrft[i],i+1) for i in range(len(mdrft))]
            Analysis = ''.join(['Analysis is SUCCESSFULLY completed']+text)
        if cIndex == 1:
            Analysis = "Analysis is STOPPED, peak interstorey drift ratio, %d%%, is exceeded, global COLLAPSE is observed" % Dc
        print('------------------------------------------------------------------------')
        print(Analysis)
        
        return mdrft
   
    def analysis_define(solver=0,Tol=1e-8,maxNumIter=300,pFlag=0,nType=2):
        """
            solver= default 0 other options:
                                            0:'BandGen',
                                            1:'BandSPD',
                                            2:'ProfileSPD',
                                            3:'SuperLU',
                                            4:'UmfPack',
                                            5:'FullGeneral',
                                            6:'SparseSYM',
            Tol=1e-8,
            maxNumIter=300,
            pFlag=0,
            nType=2
        """
        
        system = {
            0:'BandGen',
            1:'BandSPD',
            2:'ProfileSPD',
            3:'SuperLU',
            4:'UmfPack',
            5:'FullGeneral',
            6:'SparseSYM',
        }
        # Convergence Test -- determines when convergence has been achieved.
        ops.test('NormDispIncr', Tol, maxNumIter, pFlag, nType)
        #ops.test('EnergyIncr', Tol, maxNumIter,pFlag,nType)

        # SolutionAlgorithm -- determines the sequence of steps taken to solve the non-linear equation at the current time step
        ops.algorithm('Newton', '-initial')
        #ops.algorithm('ModifiedNewton',secantStiff,'-initial')

        # DOF_Numberer -- determines the mapping between equation numbers and degrees-of-freedom
        ops.numberer('RCM')

        # SystemOfEqn/Solver -- within the solution algorithm, it specifies how to store and solve the system of equations in the analysis
        ops.system(system[solver])
        #ops.system('UmfPack')

        # Constraints handler: determines how the constraint equations are enforced in the analysis -- how it handles the boundary conditions/imposed displacements
        ops.constraints('Transformation')

        # Integrator -- determines the predictive step for time t+dt
        gamma = 0.5   # Set Newmark gamma coefficient
        beta = 0.25   # Set Newmark beta coefficient
        ops.integrator('Newmark', gamma, beta)
        
        #ops.integrator('HHT', alpha, gamma=1.5-alpha, beta=(2-alpha)^2/4)

        # AnalysisType -- defines what type of analysis is to be performed ('Static', 'Transient' etc.)
        ops.analysis('Transient')

    def run_timehistory(columndict,DtAnalysis   =0.01,TmaxAnalysis=10,
                        outEleForces = False,
                        outNodalDisp = True,
                        outFiber     = True,
                        animotions   = False,
                        outSection   = True,
                        fiberData    = False
                        ):
        """
        DtAnalysis=0.01,TmaxAnalysis=10,outEleForces =False,outNodalDisp=True,
        #outFiber : [record,Hsec,cover,idCoverMat,idCoreMat,idSteel,lpl] default:[True,0.30,0.05,1,2,3]
        
        """
        # Konfigürasyon
        an = datetime.now()
        year = an.year
        month = an.month
        day = an.day
        logging.basicConfig(
        filename=f"{day}_{month}_{year}.txt",
        format="%(asctime)s - %(levelname)s - %(message)s ",
        filemode="a",
        level=logging.DEBUG)
        # Logger'ımızı bir değişkene atayalım
        logger = logging.getLogger()
        
        DtAnalysis = DtAnalysis*Unit.sec
        TmaxAnalysis = TmaxAnalysis*Unit.sec
        

        tCurrent = ops.getTime()

        algorithm = {1:'KrylovNewton', 2: 'SecantNewton' , 3:'ModifiedNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS',
                     7: 'Broyden', 8: 'NewtonLineSearch'}
        
        time = [tCurrent]
        
        
        
        
        outFiberData  = {}
        nodalOutput   = {}
        for node in ops.getNodeTags():
            nodalOutput[node] = {
                                    "nodeX_disp"       :[0],
                                    "nodeY_disp"       :[0],
                                    "node_rotation"    :[0]
                                }
        
        n_steps = int((TmaxAnalysis)/DtAnalysis)
        el_tags = ops.getEleTags()
        nels = len(el_tags)
        Eds = np.zeros((n_steps+100, nels, 6))
        step = 1           # current step number

                
        MomentRotation    = pd.DataFrame(columns=["Eletags","iMoment","iRotation","jMoment","jRotation"],index=([i for i in range(1,n_steps*max(ops.getEleTags()))]))
        NodalDisplacement = pd.DataFrame(columns=["Nodetags","NodeDispX","NodeDispY","NodeRotation"],index=([i for i in range(1,n_steps*max(ops.getNodeTags()))]))
        ElementForces     = pd.DataFrame(columns=["Eletags","Ni","Nj","Ti","Tj","Mi","Mj"],index=([i for i in range(1,n_steps*max(ops.getEleTags()))]))
        FiberStressStrain = pd.DataFrame(columns=["Eletags","TopCoverStress","TopCoverStrain","TopCoreStress","TopCoreStrain","TopSteelStress","TopSteelStrain","BotCoverStress","BotCoverStrain","BotCoreStress","BotCoreStrain","BotSteelStress","BotSteelStrain"])
        FiberData         = pd.DataFrame()

        ok = 0
        indexForce = 1
        indexDisp = 1
        indexMomrot = 1
        indexFibers = 1

        while ok == 0 and tCurrent <= TmaxAnalysis:
            ok = ops.analyze(1,DtAnalysis)
            if ok != 0:
                logger.warning(f"analiz, {tCurrent}/{TmaxAnalysis} adiminda converge edemedi newton algoritmasi degistiriliyor...")
                for j in algorithm:
                    if j < 4:
                        ops.algorithm(algorithm[j], '-initial')
                    else:
                        ops.algorithm(algorithm[j])
                    logger.info(f"~~~ {algorithm[j]} algoritmasi deneniyor ...")
                    ok = ops.analyze(1,DtAnalysis)
                    if ok == 0:
                        logger.info(f"~~~ {algorithm[j]} algoritmasi ile cozuldu, Yeniden Newton algoritmasina donuluyor ...")
                        ops.algorithm('Newton','-initial')
                        break
            if ok == 0 :
                #print(f"Analiz converge etti {tCurrent}/{TmaxAnalysis}")
                logger.info(f"~~~ analiz adimi Newton algoritmasiyla basari ile cozuldu {tCurrent}/{TmaxAnalysis} ...")
                tCurrent = ops.getTime()                
                time.append(tCurrent)
                step += 1           # Update the current step number
                
                
                #Dataframe yapısına çevrildi 03.12.2022
                if outEleForces:
                    for ele in ops.getEleTags():
                        forces = ops.eleForce(ele,-1) #eleman kuvvetlerini atıyorum
                        ElementForces.loc[indexForce] = [int(ele),forces[0],forces[1],forces[2],forces[3],forces[4],forces[5]]
                        indexForce += 1

                #Dataframe yapısına çevrildi 03.12.2022       
                if outNodalDisp:
                    for node in ops.getNodeTags():
                        disps = ops.nodeDisp(node,-1)
                        NodalDisplacement.loc[indexDisp] = [int(node),disps[0],disps[1],disps[2]]
                        indexDisp += 1
                
                #            [0    , 1  ,  2  ,    3     ,    4    ,   5   , 6]    
                #outFiber : [record,Hsec,cover,idCoverMat,idCoreMat,idSteel,lpl]
                #  0 ,1 , 2  , 3  ,  4 ,  5  ,     6    ,   7     ,   8   ,  9
                #  iN,jN, L  ,Bsec,Hsec,cover,idCoverMat,idCoreMat,idSteel, lpl
                # [ 1, 5, 3.0, 0.7, 0.3,0.025,     1    ,   3     ,   4   , 0.15]
                #Dataframe yapısına çevrildi 03.12.2022
                if outFiber:
                    for ele in columndict.keys():
                        fiber_stressStrain_top_cover        = ops.eleResponse(int(ele), 'sectionX', columndict[ele][9],  'fiber',  columndict[ele][4]/2                  ,    0  , columndict[ele][6]  , 'stressStrain' )
                        fiber_stressStrain_top_core         = ops.eleResponse(int(ele), 'sectionX', columndict[ele][9],  'fiber', (columndict[ele][4]-columndict[ele][5])/2,    0  , columndict[ele][7]  , 'stressStrain' )
                        # Tension fiberint
                        fiber_stressStrain_bot_cover        = ops.eleResponse(int(ele), 'sectionX', columndict[ele][9],  'fiber', -columndict[ele][4]/2                  ,    0  , columndict[ele][6]  , 'stressStrain' )
                        fiber_stressStrain_bot_core         = ops.eleResponse(int(ele), 'sectionX', columndict[ele][9],  'fiber',-(columndict[ele][4]-columndict[ele][5])/2,    0  , columndict[ele][7]  , 'stressStrain' )
                        # Tension fiber - steelint
                        fiber_stressStrain_steel_top        = ops.eleResponse(int(ele), 'sectionX', columndict[ele][9],  'fiber', columndict[ele][4]/2                   ,    0  , columndict[ele][8]  , 'stressStrain' )
                        fiber_stressStrain_steel_bot        = ops.eleResponse(int(ele), 'sectionX', columndict[ele][9],  'fiber',-columndict[ele][4]/2                   ,    0  , columndict[ele][8]  , 'stressStrain' )
                        
                        recordDatas = [int(ele),
                        fiber_stressStrain_top_cover[0],fiber_stressStrain_top_cover[1],
                        fiber_stressStrain_top_core[0],fiber_stressStrain_top_core[1],
                        fiber_stressStrain_steel_top[0],fiber_stressStrain_steel_top[1],
                        fiber_stressStrain_bot_cover[0],fiber_stressStrain_bot_cover[1],
                        fiber_stressStrain_bot_core[0],fiber_stressStrain_bot_core[1],
                        fiber_stressStrain_steel_bot[0],fiber_stressStrain_steel_bot[1]]
                        FiberStressStrain.loc[indexFibers] = recordDatas
                        indexFibers += 1

                        
                        """
                        if outFiber[0]:
                            for ele in ops.getEleTags():
                                if ele not in sectionOutput.keys():
                                    fiberOutput[ele]={  "top_cover":{"stress":[],"strain":[]},
                                                        "top_core" :{"stress":[],"strain":[]},
                                                        "bot_cover":{"stress":[],"strain":[]},
                                                        "bot_core" :{"stress":[],"strain":[]},
                                                        "steel_top":{"stress":[],"strain":[]},
                                                        "steel_bot":{"stress":[],"strain":[]},
                                                    }
                                #            [0    , 1  ,  2  ,    3     ,    4    ,   5   , 6]    
                                #outFiber : [record,Hsec,cover,idCoverMat,idCoreMat,idSteel,lpl]
                                # Compression fiber                                               eleman lokasyonu,  'fiber',        y=HSec/2            ,   z=0 ,   matTag 
                                # ilgili eleman uzunluğunda verilen plastik mafsal uzunluğuna en yakın fiber kesitteki 6 noktada verilen malzeme türünde okuma yapılır
                                fiber_stressStrain_top_cover        = ops.eleResponse(ele, 'sectionX', outFiber[6],  'fiber', outFiber[2]/2              ,    0  , outFiber[3]  , 'stressStrain' )
                                fiber_stressStrain_top_core         = ops.eleResponse(ele, 'sectionX', outFiber[6],  'fiber', (outFiber[2]-outFiber[3])/2,    0  , outFiber[4]  , 'stressStrain' )
                                # Tension fiber
                                fiber_stressStrain_bot_cover        = ops.eleResponse(ele, 'sectionX', outFiber[6],  'fiber',-outFiber[2]/2              ,    0  , outFiber[3]  , 'stressStrain' )
                                fiber_stressStrain_bot_core         = ops.eleResponse(ele, 'sectionX', outFiber[6],  'fiber',-(outFiber[2]-outFiber[3])/2,    0  , outFiber[4]  , 'stressStrain' )
                                # Tension fiber - steel
                                fiber_stressStrain_steel_top        = ops.eleResponse(ele, 'sectionX', outFiber[6],  'fiber',outFiber[2]/2               ,    0  , outFiber[5]  , 'stressStrain' )
                                fiber_stressStrain_steel_bot        = ops.eleResponse(ele, 'sectionX', outFiber[6],  'fiber',-outFiber[2]/2              ,    0  , outFiber[5]  , 'stressStrain' )
                                
                                
                                if len(fiber_stressStrain_top_cover)>1:
                                    fiberOutput[ele]["top_cover"]["stress"]   .append(fiber_stressStrain_top_cover[0])
                                    fiberOutput[ele]["top_cover"]["strain"]   .append(fiber_stressStrain_top_cover[1])
                                
                                if len(fiber_stressStrain_top_core)>1:
                                    fiberOutput[ele]["top_core"]["stress"]    .append(fiber_stressStrain_top_core[0])
                                    fiberOutput[ele]["top_core"]["strain"]    .append(fiber_stressStrain_top_core[1])
                                
                                if len(fiber_stressStrain_bot_cover)>1:
                                    fiberOutput[ele]["bot_cover"]["stress"]   .append(fiber_stressStrain_bot_cover[0])
                                    fiberOutput[ele]["bot_cover"]["strain"]   .append(fiber_stressStrain_bot_cover[1])
                                    
                                if len(fiber_stressStrain_bot_core)>1:    
                                    fiberOutput[ele]["bot_core"] ["stress"]   .append(fiber_stressStrain_bot_core [0])
                                    fiberOutput[ele]["bot_core"] ["strain"]   .append(fiber_stressStrain_bot_core [1])
                                
                                if len(fiber_stressStrain_steel_top)>1:
                                    fiberOutput[ele]["steel_top"]["stress"]   .append(fiber_stressStrain_steel_top[0])
                                    fiberOutput[ele]["steel_top"]["strain"]   .append(fiber_stressStrain_steel_top[1])
                                
                                if len(fiber_stressStrain_steel_bot)>1:
                                    fiberOutput[ele]["steel_bot"]["stress"]   .append(fiber_stressStrain_steel_bot[0])
                                    fiberOutput[ele]["steel_bot"]["strain"]   .append(fiber_stressStrain_steel_bot[1])
                        """
                                                                            
                    if animotions:
                        for el_i, ele_tag in enumerate(el_tags):
                            nd1, nd2 = ops.eleNodes(ele_tag)
                            Eds[step, el_i, :] = [ops.nodeDisp(nd1)[0],
                                                    ops.nodeDisp(nd1)[1],
                                                    ops.nodeDisp(nd1)[2],
                                                    ops.nodeDisp(nd2)[0],
                                                    ops.nodeDisp(nd2)[1],
                                                    ops.nodeDisp(nd2)[2]]
                    
                    #Dataframe yapısına çevrildi 03.12.2022
                    if outSection:
                        for ele in ops.getEleTags():
                            basicForces        = ops.eleResponse(ele,'basicForce')                   #[axial , iMoment,jMoment] 
                            basicDeforms       = ops.eleResponse(ele,'basicDeformation')             #[pinch , irotation,jrotation] 
                            MomentRotation.loc[indexMomrot] = [int(ele),basicForces[1],basicDeforms[1],basicForces[2],basicDeforms[2]]
                            indexMomrot += 1
                            
                        """
                        for ele in ops.getEleTags():

                            if ele not in sectionOutput.keys():
                                sectionOutput[ele]={
                                                    #"force"         : [],
                                                    "ibasicForce"    : [],
                                                    "jbasicForce"    : [],
                                                    "itotalRot"      : [],
                                                    #"iplasticRot"   : [],
                                                    #"ielasticRot"   : [],
                                                    "jtotalRot"      : [],
                                                    #"jplasticRot"   : [],
                                                    #"jelasticRot"   : [],
                                                    #"intgLocation"  : ops.eleResponse(ele,"integrationPoints"),
                                                    #"intgForces"    : [],
                                                    #"intgCurvature" : []
                                                }
                                
                            #intgpnt = len(ops.sectionLocation(ele))
                            basicForces        = ops.eleResponse(ele,'basicForce')                   #[axial , iMoment,jMoment] 
                            basicDeforms       = ops.eleResponse(ele,'basicDeformation')             #[pinch , irotation,jrotation] 
                            #plasticDeforms     = ops.eleResponse(ele,'plasticDeformation')           #[pinch , irotation,jrotation] 
                            #elasticDeforms     = [(i-z) for i,z in zip(basicDeforms,plasticDeforms)] #[pinch , irotation,jrotation]
                            #forces             = ops.eleResponse(ele,"force")                        #[iaxial,ishear,imoment,jaxial,jshear,jmoment]
                            
                            ibasicForces       = basicForces[1]   #i ucunun momenti tutuluyor her zaman adımında
                            jbasicForces       = basicForces[2]   #j ucunun momenti tutuluyor her zaman adımında
                            
                            ibasicDeforms      = basicDeforms[1]    #i ucunun dönmesi tutuluyor her zaman adımında
                            #iplasticDeforms    = plasticDeforms[1]  #i ucunun dönmesi tutuluyor her zaman adımında
                            #ielasticDeforms    = elasticDeforms[1]  #i ucunun dönmesi tutuluyor her zaman adımında
                            
                            jbasicDeforms      = basicDeforms[2]  #i ucunun dönmesi tutuluyor her zaman adımında
                            #jplasticDeforms    = plasticDeforms[2]  #i ucunun dönmesi tutuluyor her zaman adımında
                            #jelasticDeforms    = elasticDeforms[2]  #i ucunun dönmesi tutuluyor her zaman adımında
                            
                            #sectionForces      = ops.eleResponse(ele,"section",1,"force")       #integrasyon noktalarındaki kuvvetler tutuluyor
                            #sectionCurvatures  = ops.eleResponse(ele,"section",1,"deformation") #integrasyon noktalarındaki eğilmeler tutuluyor

                            #sectionOutput[ele]["force"]        .append(forces)        #eleman iç kuvvetleri tutuluyor
                            sectionOutput[ele]["ibasicForce"]  .append(ibasicForces)
                            
                            sectionOutput[ele]["jbasicForce"]  .append(jbasicForces)
                            sectionOutput[ele]["itotalRot"]    .append(ibasicDeforms)
                            #sectionOutput[ele]["iplasticRot"]  .append(iplasticDeforms)
                            #sectionOutput[ele]["ielasticRot"]  .append(ielasticDeforms)
                            sectionOutput[ele]["jtotalRot"]    .append(jbasicDeforms)

                            #sectionOutput[ele]["jplasticRot"]  .append(jplasticDeforms)
                            #sectionOutput[ele]["jelasticRot"]  .append(jelasticDeforms)
                            #sectionOutput[ele]["intgForces"]   .append(sectionForces)
                            #sectionOutput[ele]["intgCurvature"].append(sectionCurvatures)

                            #Bilgilendirme için normalde kullanılmıyor...
                            #axial_moment_i     = ops.eleResponse(ele,'section',2,'force')      # Column section forces, axial and moment,    node i return list
                            #axial_curvature_i  = ops.eleResponse(ele,'section',2,'deformation')# section deformations,  axial(ezilme) and curvature(eğrilik), node i return list
                            #axial_moment_j     = ops.eleResponse(ele,'section',5,'force')      # Column section forces, axial and moment,    node j return list
                            #axial_curvature_j  = ops.eleResponse(ele,'section',5,'deformation')# section deformations,  axial and curvature, node j return list
                            #sectionOutput[ele]["axial_moments_i"]   .append(axial_moment_i)
                            #sectionOutput[ele]["axial_curvatures_i"].append(axial_curvature_i)
                            #sectionOutput[ele]["axial_moments_j"]   .append(axial_moment_j)
                            #sectionOutput[ele]["axial_curvatures_j"].append(axial_curvature_j)
                            #sectionOutput[ele]["section_rotations_i"].append( ops.sectionDeformation(ele,2,2))
                            #sectionOutput[ele]["section_rotations_j"].append(ops.sectionDeformation(ele,5,2))
                            #sectionOutput[ele]["section_moments_i"]  .append(ops.sectionForce(ele,2,2))
                            #sectionOutput[ele]["section_moments_j"]  .append( ops.sectionForce(ele,5,2))"""  

                if fiberData:
                    for ele in ops.getEleTags():
                        if ele not in outFiberData.keys():
                            outFiberData[ele]={"yCoord":[],
                                           "zCoord":[],
                                           "sigma" :[],
                                           "eps"   :[]}
                        data = ops.eleResponse(ele,'section',1,'fiberData')
                        Ndata = len(data)
                        Nfibers = int(Ndata/5)
                        y = data[0:Ndata:5]
                        z = data[1:Ndata:5]
                        # A = data[2:Ndata:5] # If you want fiber areas
                        sig = data[3:Ndata:5]
                        eps = data[4:Ndata:5]
                        
                        if len(outFiberData[ele]["yCoord"]) == 0:
                            outFiberData[ele]["yCoord"].append(y)
                        if len(outFiberData[ele]["zCoord"]) == 0:
                            outFiberData[ele]["zCoord"].append(z)
                        outFiberData[ele]["sigma"] .append(sig)
                        outFiberData[ele]["eps"]   .append(eps)
            else:
                break
        # Print the final status of the analysis
            # ------------------------------------------------------------------------
        if ok != 0:
            logger.warning(f"~~~ Deplasman kontrol analizi hata verdi hata veren analiz adimi:{tCurrent}/{TmaxAnalysis}")
        else:
            logger.info("Deplasman kontrol analizi basarili...")
        
       
        return ElementForces,NodalDisplacement,MomentRotation,FiberStressStrain,Eds,time
    
    def IDA(DtAnalysis=0.01,TmaxAnalysis=10,patternTag=1 ,incrGMFactor=0.05,MaxGMFactor=2.0):
        # Analysis duration and time step
        Nsteps = int(TmaxAnalysis/DtAnalysis)

        # Arrays for plotting
        Uplot = []
        gmPlot = []

        # Maximum and increment in ground motion factor
        

        gmFact = 0.0
        while gmFact < MaxGMFactor:
            
            gmFact += incrGMFactor
            gmPlot.append(gmFact)
            
            ops.remove('loadPattern',patternTag)
            
            ops.reset()
            
            # Redefine ground motion with new factor
           #ops.pattern('UniformExcitation', IDloadTag,GMdirection, '-accel', tsTag)# pattern object
            ops.pattern('UniformExcitation',patternTag,     1     ,'-accel' ,   1,'-factor',gmFact)
            
            # Perform analysis and record maximum displacement
            Umax = 0
            for i in range(Nsteps):
                ok = ops.analyze(1,DtAnalysis)
                if ok != 0:
                    break
                U = ops.nodeDisp(max(ops.getNodeTags()),1)
                if abs(U) > Umax:
                    Umax = abs(U)
                    
            Uplot.append(Umax)
            print(gmFact)
        fig, ax = plt.subplots(figsize=(10,10),constrained_layout=True, sharey=True)
        plt.rcParams.update({'font.size': 16})
        ax.grid()
        ax.plot(Uplot,gmPlot)
        plt.xlim(0)
        plt.ylim(0)
        ax.axhline(0, color='black', lw=2)
        ax.axvline(0, color='black', lw=2)
        ax.set(xlabel="Roof Horizontal disp (m) ", ylabel="Scale Factor",title="IDA")
    
    def allMatrix():
        
        import scipy.linalg as slin
        NDF = 3
        ops.wipeAnalysis()
        ops.system('FullGeneral')
        ops.analysis('Transient')
        
        # Mass
        ops.integrator('GimmeMCK',1.0,0.0,0.0)
        ops.analyze(1,0.0)
        
        # Number of equations in the model
        N = ops.systemSize() # Has to be done after analyze
        
        M = ops.printA('-ret') # Or use ops.printA('-file','M.out')
        M = np.array(M) # Convert the list to an array
        M.shape = (N,N) # Make the array an NxN matrix
        
        # Stiffness
        ops.integrator('GimmeMCK',0.0,0.0,1.0)
        ops.analyze(1,0.0)
        K = ops.printA('-ret')
        K = np.array(K)
        K.shape = (N,N)
        
        # Damping
        ops.integrator('GimmeMCK',0.0,1.0,0.0)
        ops.analyze(1,0.0)
        C = ops.printA('-ret')
        C = np.array(C)
        C.shape = (N,N)


        # Determine number of DOFs with mass
        massDOFs = []
        for nd in ops.getNodeTags():
            for j in range(NDF): # NDF is number of DOFs/node
                if ops.nodeMass(nd,j+1) > 0.0:
                    massDOFs.append(ops.nodeDOFs(nd)[j])
        
        # Number of DOFs with mass
        Nmass = len(massDOFs)
        
        # DOFs without mass
        masslessDOFs = np.setdiff1d(range(N),massDOFs)
        Nmassless = len(masslessDOFs)
        
        # Form matrices for D*x = -lam*B*x
        B = np.zeros((2*Nmass,2*Nmass)) # = [ 0 M; M C]
        D = np.zeros((2*Nmass,2*Nmass)) # = [-M 0; 0 K]
        
        # Mass
        B[:Nmass,:][:,Nmass:2*Nmass] =  M[massDOFs,:][:,massDOFs]
        B[Nmass:2*Nmass,:][:,:Nmass] =  M[massDOFs,:][:,massDOFs]
        D[:Nmass,:][:,:Nmass]        = -M[massDOFs,:][:,massDOFs]
        
        # Damping
        B[Nmass:2*Nmass,:][:,Nmass:2*Nmass] = C[massDOFs,:][:,massDOFs]
        
        # Static condensation
        Kmm = K[massDOFs,:][:,massDOFs];     Kmn = K[massDOFs,:][:,masslessDOFs]
        Knm = K[masslessDOFs,:][:,massDOFs]; Knn = K[masslessDOFs,:][:,masslessDOFs]
        # Kc = Kmm - Kmn*inv(Knn)*Knm
        if Nmassless > 0:
            Kc = Kmm - np.dot(Kmn,np.linalg.solve(Knn,Knm))
        else:
            Kc = K
        
        # Stiffness at DOFs with mass
        D[Nmass:2*Nmass,:][:,Nmass:2*Nmass] = Kc
        
        # State space eigenvalue analysis
        [lam,x] = slin.eig(D,-B)
        
        return M,K,C,Nmass,
        

