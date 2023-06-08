#%%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#This files defines the simulation of hydraulically actuated flexible structure.

# Author        : Qasim Khadim
# Contact       : qasim.khadim@outlook.com,qkhadim22 (Github)
# Dated         : 02-05-2023
# Organization  : University of Oulu in the collaboration of LUT University and University of Innsbruck.

# Copyright     :
#%%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

import sys
sys.exudynFast = True #this variable is used to signal to load the fast exudyn module

import exudyn as exu
from exudyn.itemInterface import *
from exudyn.utilities import *
from exudyn.FEM import *

import math as mt
from math import sin, cos, sqrt, pi, tanh
import time

SC  = exu.SystemContainer()
mbs = SC.AddSystem()
#%%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# physical parameters
g           = [0, -9.81, 0]  # Gravity
tEnd        = 24  # simulation time
h           = 1e-3  # step size
nSteps      = int(tEnd/h)+2

# Loading Graphics of bodies
fileName1       = 'Graphics_Exudyn/Pillar.stl'
fileName4       = 'Graphics_Exudyn/Bracket1.stl'
fileName5       = 'Graphics_Exudyn/Bracket2.stl'

#fileNameT       = 'TiltBoomANSYS/TiltBoom' #for load/save of FEM data

feL             = FEMinterface()
feT             = FEMinterface()

#Ground body
oGround         = mbs.AddObject(ObjectGround())
markerGround    = mbs.AddMarker(MarkerBodyRigid(bodyNumber=oGround, localPosition=[0, 0, 0]))

# Pillar definition in Exudyn
L1              = 0.365    # Length in x-direction
H1              = 1.4769      # Height in y-direction
W1              = 0.25    # Width in z-direction
bodyDim1        = [L1, H1, W1]  # body dimensions
pMid1           = np.array([-0.017403, 0.577291, 0])  # center of mass, body0
PillarP         = np.array([0, 0, 0])
iCube1          = RigidBodyInertia(mass=93.26, com=pMid1,
                                   inertiaTensor=np.array([[16.358844,-1.27808, 1.7e-5],
                                                           [-1.27808, 0.612552, -5.9e-5],
                                                           [1.7e-5,  -5.9e-5  , 16.534255]]),
                                                               inertiaTensorAtCOM=True)

graphicsBody1   = GraphicsDataFromSTLfile(fileName1, color4black,verbose=False, invertNormals=True,invertTriangles=True)
graphicsBody1   = AddEdgesAndSmoothenNormals(graphicsBody1, edgeAngle=0.25*pi,addEdges=True, smoothNormals=True)
graphicsCOM1    = GraphicsDataBasis(origin=iCube1.com, length=2*W1)

# Definintion of pillar as body in Exudyn and node n1
[n1, b1]        = AddRigidBody(mainSys=mbs,
                             inertia=iCube1,  # includes COM
                             nodeType=exu.NodeType.RotationEulerParameters,
                             position=PillarP,
                             rotationMatrix=np.diag([1, 1, 1]),
                             gravity=g,
                             graphicsDataList=[graphicsCOM1, graphicsBody1])

#Pillar
Marker3         = mbs.AddMarker(MarkerBodyRigid(bodyNumber=b1, localPosition=[0, 0, 0]))                     #With Ground
Marker4         = mbs.AddMarker(MarkerBodyRigid(bodyNumber=b1, localPosition=[-0.09, 1.4261, 0]))            #Lift Boom
Marker5         = mbs.AddMarker(MarkerBodyRigid(bodyNumber=b1, localPosition=[0.17, 0.386113249, 0]))        # Cylinder 1 position

# Fixed joint between Pillar and Ground
mbs.AddObject(GenericJoint(markerNumbers=[markerGround, Marker3],constrainedAxes=[1, 1, 1, 1, 1, 1],
visualization=VObjectJointGeneric(axesRadius=0.2*W1,axesLength=1.4*W1)))


def FlexibleLiftBoom(FineMesh, nModes, loadFromSavedNPY, ComputeModes, FreeModes, HCB, ModeAnimation, StaticCase, Hydraulics,
                   useFriction, Visualization,Plotting):
    
    global Fc, Fs, sig2, vs, Av0, Av1
    
    if FineMesh:
        
        fileNameL       = 'LiftBoom/ABAQUS/liftboom-free-050623' #To load fine mesh data from Abaqus
        
        if not loadFromSavedNPY: 
                start_time                      = time.time()
                nodes                           = feL.ImportFromAbaqusInputFile(fileNameL+'.inp', typeName='Part', name='P000524_A_1-Nostopuomi_v2_fem_st')
                feL.ReadMassMatrixFromAbaqus(fileName=fileNameL + '_MASS2.mtx')             #Load mass matrix
                feL.ReadStiffnessMatrixFromAbaqus(fileName=fileNameL + '_STIF2.mtx')        #Load stiffness matrix
                feL.SaveToFile(fileNameL)
                print("--- saving LiftBoom FEM Abaqus data took: %s seconds ---" % (time.time() - start_time)) 
           
                if ComputeModes:
                    
                     from Models.ComputeModes import LiftBoomModes
                     LiftBoomModes(feL, FineMesh, nModes, FreeModes, HCB, ModeAnimation)
                                        
        else:       
                print('importing Abaqus FEM data structure of Lift Boom...')
                start_time = time.time()
                feL.LoadFromFile(fileNameL)
                cpuTime = time.time() - start_time
                print("--- importing FEM data took: %s seconds ---" % (cpuTime))
                
                # For computing modes
                
                if ComputeModes:
                    
                     from Models.ComputeModes import LiftBoomModes
                     LiftBoomModes(feL, FineMesh, nModes, FreeModes, HCB, ModeAnimation)
                     
                     
    else:
        
        fileNameL       = 'LiftBoom/ANSYS/LiftBoom' #To load fine mesh data from AANSYS
        
        if not loadFromSavedNPY: 
                start_time                          = time.time()
                
                inputFileNameStiffnessMatrixL       = 'LiftBoom/ANSYS/'+'stiffnessMatrix.txt'
                inputFileNameMassMatrixL            = 'LiftBoom/ANSYS/'+'massMatrix.txt'
                inputFileNameNodalCoordinatesL      = 'LiftBoom/ANSYS/'+'nodalCoordinates.txt'
                inputFileNameElementsL              = 'LiftBoom/ANSYS/'+'elements.txt'
                inputFileNameNodalMappingVectorL    = 'LiftBoom/ANSYS/'+'nodalMappingVector.txt'                
                feL.ReadStiffnessMatrixFromAnsys(inputFileNameStiffnessMatrixL, inputFileNameNodalMappingVectorL, verbose=True)
                feL.ReadMassMatrixFromAnsys(inputFileNameMassMatrixL, inputFileNameNodalMappingVectorL, verbose=True)
                feL.ReadNodalCoordinatesFromAnsys(inputFileNameNodalCoordinatesL, verbose=True)
                feL.ReadElementsFromAnsys(inputFileNameElementsL, verbose=True)                
                feL.SaveToFile(fileNameL)
                print("--- saving LiftBoom FEM ANSYS data took: %s seconds ---" % (time.time() - start_time)) 
                
                if ComputeModes:
                    
                     from Models.ComputeModes import LiftBoomModes
                     LiftBoomModes(feL, FineMesh, nModes, FreeModes, HCB, ModeAnimation)
        else:
                print('importing ANSYS FEM data structure of Lift Boom...')      
                start_time = time.time()
                feL.LoadFromFile(fileNameL)                
                cpuTime = time.time() - start_time
                print("--- importing FEM data took: %s seconds ---" % (cpuTime))
        
                if ComputeModes:
                    
                     from Models.ComputeModes import LiftBoomModes
                     LiftBoomModes(feL, FineMesh, nModes, FreeModes, HCB, ModeAnimation)
    
    # Boundary condition at pillar
    p2                  = [0, 0,-10e-2]
    p1                  = [0, 0, 10e-2]
    radius1             = 2.5e-002
    nodeListJoint1      = feL.GetNodesOnCylinder(p1, p2, radius1, tolerance=1e-2) 
    pJoint1             = feL.GetNodePositionsMean(nodeListJoint1)
    nodeListJoint1Len   = len(nodeListJoint1)
    noodeWeightsJoint1  = [1/nodeListJoint1Len]*nodeListJoint1Len
           
    # Boundary condition at Piston 1
    p4                  = [0.3025,-0.1049,-10e-2]
    p3                  = [0.3025,-0.1049, 10e-2]
    radius2             = 3.6e-002
    nodeListPist1       = feL.GetNodesOnCylinder(p3, p4, radius2, tolerance=1e-2)  
    pJoint2             = feL.GetNodePositionsMean(nodeListPist1)
    nodeListPist1Len    = len(nodeListPist1)
    noodeWeightsPist1   = [1/nodeListPist1Len]*nodeListPist1Len
    
    if not Hydraulics:
        # STEP 2: Craig-Bampton Modes
        boundaryList    = [nodeListJoint1]
        start_time      = time.time()
        feL.ComputeHurtyCraigBamptonModes(boundaryNodesList=boundaryList, nEigenModes=nModes, 
                                                    useSparseSolver=True,computationMode = HCBstaticModeSelection.RBE2)
        print("Hurty-Craig Bampton modes... ")
        print("eigen freq.=", feL.GetEigenFrequenciesHz())
        print("HCB modes needed %.3f seconds" % (time.time() - start_time))  
    
        colLift = color4blue
           
        LiftBoom            = ObjectFFRFreducedOrderInterface(feL)
        
        LiftBoomFFRF        = LiftBoom.AddObjectFFRFreducedOrder(mbs, positionRef=np.array([-0.09, 1.4261, 0]), 
                                              initialVelocity=[0,0,0], 
                                              initialAngularVelocity=[0,0,0],
                                              rotationMatrixRef  = RotationMatrixZ(mt.radians(0)),
                                              gravity=g,
                                              color=colLift,)
        
        Marker7             = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber=LiftBoomFFRF['oFFRFreducedOrder'],
                                              meshNodeNumbers=np.array(nodeListJoint1), #these are the meshNodeNumbers
                                              weightingFactors=noodeWeightsJoint1))
        if StaticCase:
            # Add Fixed Joint
            mbs.AddObject(GenericJoint(markerNumbers=[Marker4, Marker7],constrainedAxes=[1, 1, 1, 1, 1, 1],
                           visualization=VObjectJointGeneric(axesRadius=0.18*0.263342,axesLength=1.1*0.263342)))
            
        else:
            # Add Revolute Joint
            mbs.AddObject(GenericJoint(markerNumbers=[Marker4, Marker7],constrainedAxes=[1,1,1,1,1,0],
                         visualization=VObjectJointGeneric(axesRadius=0.18*0.263342,axesLength=1.1*0.263342)))   
        
        
        # Joint 3
        p10                 = [2.89,0.0246,-7.4e-2]
        p9                  = [2.89,0.0246, 7.4e-2]
        pdef                = [2.89,0.0246, 0]
        radius5             = 5.2e-002
        nodeListJoint3      = feL.GetNodesOnCylinder(p9, p10, radius5, tolerance=1e-2)  
        pJoint5             = feL.GetNodePositionsMean(nodeListJoint3)
        nodeListJoint3Len   = len(nodeListJoint3)
        noodeWeightsJoint3  = [1/nodeListJoint3Len]*nodeListJoint3Len
        
        Marker11            = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber= LiftBoomFFRF['oFFRFreducedOrder'], 
                                                                 meshNodeNumbers=np.array(nodeListJoint3),
                                                                 #referencePosition=pdef,
                                                                 #useAlternativeApproach=altApproach,
                                                                weightingFactors=noodeWeightsJoint3))

           
        # Add Sensor for deflection
        Deflecion          = mbs.AddSensor(SensorMarker(markerNumber=Marker11, storeInternal=True, 
                                                     fileName='ExData/LiftFlexible/SimpleDef.txt', outputVariableType=exu.OutputVariableType.Displacement))
                      
    else:
        # STEP 2: Craig-Bampton Modes
        boundaryList        = [nodeListJoint1, nodeListPist1]
        start_time          = time.time()
        feL.ComputeHurtyCraigBamptonModes(boundaryNodesList=boundaryList, nEigenModes=nModes, 
                                                    useSparseSolver=True,computationMode = HCBstaticModeSelection.RBE2)
        print("Hurty-Craig Bampton modes... ")
        print("eigen freq.=", feL.GetEigenFrequenciesHz())
        print("HCB modes needed %.3f seconds" % (time.time() - start_time))  
    
        colLift = color4blue
           
        LiftBoom            = ObjectFFRFreducedOrderInterface(feL)
        
        LiftBoomFFRF        = LiftBoom.AddObjectFFRFreducedOrder(mbs, positionRef=np.array([-0.09, 1.4261, 0]), 
                                              initialVelocity=[0,0,0], 
                                              initialAngularVelocity=[0,0,0],
                                              rotationMatrixRef  = RotationMatrixZ(mt.radians(-27.227343749651500)),
                                              gravity=g,
                                             color=colLift,)
        
        Marker7             = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber=LiftBoomFFRF['oFFRFreducedOrder'],
                                              meshNodeNumbers=np.array(nodeListJoint1), #these are the meshNodeNumbers
                                              weightingFactors=noodeWeightsJoint1))
        Marker8             = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber=LiftBoomFFRF['oFFRFreducedOrder'],
                                              meshNodeNumbers=np.array(nodeListPist1), #these are the meshNodeNumbers
                                              weightingFactors=noodeWeightsPist1))
        
        #Revolute Joint
        mbs.AddObject(GenericJoint(markerNumbers=[Marker4, Marker7],constrainedAxes=[1,1,1,1,1,0],
                           visualization=VObjectJointGeneric(axesRadius=0.18*0.263342,axesLength=1.1*0.263342)))
        
        colCyl              = color4orange
        colPis              = color4grey 
        
        import scipy.io
        
        Data             = scipy.io.loadmat('ExData/Lift@1DOF')
        DataM            =  Data['Xt']
        DataU            =  Data['U1']
        tspan            = Data['t']
        Datas            =  Data['s1']
        Datads           =  Data['ds1']
        DataFm           =  Data['F_mu_1']
        sM               =  Datas[0, :]
        dsM              =  Datads[0, :]
        p1M              =  DataM[2, :]
        p2M              =  DataM[3, :]
        aM               =  DataM[0, :]
        dAM              =  DataM[1, :]
        Fm               =  DataFm[0,:]
        U                =  DataU[0, :]
        pP               =  100e5
        L_Cyl            =  820e-3                        # Cylinder length
        D_Cyl            =  100e-3                        # Cylinder dia
        A_1              = (pi/4)*(D_Cyl)**2             # Area of cylinder side
        L_Pis1           = 535e-3                        # Piston length, also equals to stroke length
        d_pis            = 56e-3                         # Piston dia
        A_2              = A_1-(pi/4)*(d_pis)**2         # Area on piston-rod side

        # Hydraulic Parameters
        pT               = 1e5                           # Tank pressure
        Qn               = 2.138e-8                      # Nominal flow rate of valve at 18 l/min under
        d_1              = 12.7e-3                       # Dia of volume 1
        V1               = (pi/4)*(d_1)**2*1.5             # Volume V1 = V0
        V2               = (pi/4)*(d_1)**2*1.5             # Volume V2 = V1
        A                = [A_1, A_2]
        dampingHA        = 0
        Bh               = 700e6
        Bc               = 2.1000e+11
        Bo               = 1650e6
        Fc               = 210
        Fs               = 300
        sig2             = 330
        vs               = 5e-3
        
        #ODE1 for pressures:
        nODE1            = mbs.AddNode(NodeGenericODE1(referenceCoordinates=[0,0],
                                    initialCoordinates=[2.390224339049967e+06,
                                                        2.087875598556137e+06], #initialize with 20 bar
                                    numberOfODE1Coordinates=2))
        
        
        def UFfrictionSpringDamper(mbs, t, itemIndex, u, v, k, d, f0):

            return   4*(Fc*tanh(4*(abs(v    )/vs))+(Fs-Fc)*((abs(v    )/vs)/((1/4)*(abs(v    )/vs)**2+3/4)**2))*np.sign(v )+sig2*v    *tanh(4)


        if useFriction:    
            print('USE friction')
            oFriction       = mbs.AddObject(ObjectConnectorSpringDamper(markerNumbers=[Marker5, Marker8], referenceLength=0.001,stiffness=0,
                                               damping=0, force=0, velocityOffset = 0., activeConnector = True,
                                                springForceUserFunction=UFfrictionSpringDamper,
                                                  visualization=VSpringDamper(show=False) ))
            
        oHA1                = mbs.AddObject(HydraulicActuatorSimple(name='LiftCylinder', markerNumbers=[ Marker5, Marker8], 
                                            nodeNumbers=[nODE1], offsetLength=L_Pis1, strokeLength=L_Cyl, chamberCrossSection0=A[0], 
                                            chamberCrossSection1=A[1], hoseVolume0=V1, hoseVolume1=V2, valveOpening0=0, 
                                            valveOpening1=0, actuatorDamping=1e5*0, oilBulkModulus=Bo, cylinderBulkModulus=Bc, 
                                            hoseBulkModulus=Bh, nominalFlow=Qn, systemPressure=pP, tankPressure=pT, 
                                            useChamberVolumeChange=True, activeConnector=True, 
                                            visualization={'show': True, 'cylinderRadius': 50e-3, 'rodRadius': 28e-3, 
                                                            'pistonRadius': 0.04, 'pistonLength': 0.001, 'rodMountRadius': 0.0, 
                                                            'baseMountRadius': 20.0e-3, 'baseMountLength': 20.0e-3, 'colorCylinder': color4orange,
                                                            'colorPiston': color4grey}))
            
        def PreStepUserFunction(mbs, t):
                        
            Av0 = U[mt.trunc(t/h)]
            Av1 = -Av0
            
            mbs.SetObjectParameter(oHA1, "valveOpening0", Av0)
            mbs.SetObjectParameter(oHA1, "valveOpening1", Av1)

            return True
        
        mbs.SetPreStepUserFunction(PreStepUserFunction) 
        
        sForce          = mbs.AddSensor(SensorObject(objectNumber=oHA1, storeInternal=True, 
                                                     fileName='ExData/LiftFlexible/sForce.txt', outputVariableType=exu.OutputVariableType.Force))
        sDistance       = mbs.AddSensor(SensorObject(objectNumber=oHA1, storeInternal=True, 
                                                     fileName='ExData/LiftFlexible/sDistance.txt', outputVariableType=exu.OutputVariableType.Distance))
        sDisplacement   = mbs.AddSensor(SensorObject(objectNumber=oHA1, storeInternal=True,
                                                     fileName='ExData/LiftFlexible/sDisplacement.txt', outputVariableType=exu.OutputVariableType.Displacement))
        sVelocity       = mbs.AddSensor(SensorObject(objectNumber=oHA1, storeInternal=True,
                                                     fileName='ExData/LiftFlexible/sVelocity.txt', outputVariableType=exu.OutputVariableType.Velocity))
        sPressures      = mbs.AddSensor(SensorNode(nodeNumber=nODE1, storeInternal=True, 
                                                   fileName='ExData/LiftFlexible/sPressures.txt',outputVariableType=exu.OutputVariableType.Coordinates))   
    
    mbs.Assemble()
    simulationSettings = exu.SimulationSettings()

    if not StaticCase:
        simulationSettings.timeIntegration.numberOfSteps            = nSteps
        simulationSettings.timeIntegration.endTime                  = tEnd
        simulationSettings.timeIntegration.verboseMode              = 1
        
        #simulationSettings.timeIntegration.simulateInRealtime       = True
        simulationSettings.solutionSettings.solutionWritePeriod     = 50*h
        #simulationSettings.solutionSettings.sensorsWritePeriod      = h
        simulationSettings.linearSolverType                         = exu.LinearSolverType.EigenSparse
        simulationSettings.timeIntegration.newton.useModifiedNewton = True
        simulationSettings.timeIntegration.stepInformation         += 8
        simulationSettings.timeIntegration.verboseModeFile          = 0
        simulationSettings.solutionSettings.outputPrecision         =6
        simulationSettings.solutionSettings.exportVelocities = False 
        simulationSettings.timeIntegration.generalizedAlpha.useNewmark = True
        #simulationSettings.displayComputationTime                   = True
        simulationSettings.displayStatistics                        = True 
        simulationSettings.timeIntegration.generalizedAlpha.spectralRadius=0.2            
        
        exu.SolveDynamic(mbs, simulationSettings=simulationSettings,
                 solverType=exu.DynamicSolverType.TrapezoidalIndex2)
       
       # mbs.SolveDynamic(simulationSettings)  
    else:
        simulationSettings.solutionSettings.solutionWritePeriod = 2e-2  #output interval general
        simulationSettings.solutionSettings.sensorsWritePeriod = 1e-1  #output interval of sensors
        simulationSettings.timeIntegration.numberOfSteps = int(tEnd/h) #must be integer
        simulationSettings.timeIntegration.endTime = tEnd
        simulationSettings.solutionSettings.coordinatesSolutionFileName = "LiftFlexible/staticSolution.txt"
        simulationSettings.solutionSettings.appendToFile = False
        simulationSettings.staticSolver.newton.numericalDifferentiation.relativeEpsilon = 1e-4
        simulationSettings.staticSolver.newton.relativeTolerance = 1e-6
        simulationSettings.staticSolver.newton.absoluteTolerance = 1e-1
        simulationSettings.staticSolver.verboseMode = 2
        simulationSettings.timeIntegration.newton.useModifiedNewton = True
        simulationSettings.linearSolverType = exu.LinearSolverType.EigenSparse
        simulationSettings.staticSolver.numberOfLoadSteps = 5
        
        simulationSettings.staticSolver.adaptiveStep = True
        
        mbs.SolveStatic(simulationSettings)
    
    if Visualization:
        SC.visualizationSettings.window.renderWindowSize            = [1600, 1200]        
        SC.visualizationSettings.openGL.multiSampling               = 4        
        SC.visualizationSettings.openGL.lineWidth                   = 3  
        SC.visualizationSettings.general.autoFitScene               = False      
        SC.visualizationSettings.nodes.drawNodesAsPoint             = False        
        SC.visualizationSettings.nodes.showBasis                    = True  
        
        from exudyn.interactive import SolutionViewer
        SolutionViewer(mbs) 
    
    #%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    if Plotting:
        import matplotlib.pyplot as plt
        
        if not Hydraulics:
            TipDef  = mbs.GetSensorStoredData(Deflecion)
            
            fig, ax = plt.subplots()
            plt.plot(TipDef[:,0], 1000*TipDef[:,2], label='Exudyn', color='blue')
            plt.xlim(0, max(TipDef[:,0]))
            plt.ylim(-1,1)
            plt.grid(True)
            ax.set_xlabel('Time, s')
            ax.set_ylabel('Deflection, mm')
            plt.legend()
            plt.savefig('ExData/LiftFlexible/SimpleDeflection.png', dpi=300)
            plt.show()
        
        else: 
            #PLOTTING Results
            dE      = mbs.GetSensorStoredData(sDisplacement)
            sE      = mbs.GetSensorStoredData(sDistance)
            dsE     = mbs.GetSensorStoredData(sVelocity)
            dpE     = mbs.GetSensorStoredData(sPressures)
            dFE     = mbs.GetSensorStoredData(sForce)  
    
            vAct    = np.zeros((len(sE),1))
            FE      = np.zeros((len(sE),1))
      
            for i in range(len(sE)):
                LHact   = sE[i,1]
                uSD     = dE[i,1:4]
                vSD     = dsE[i,1:4]
                vAct[i]=( vSD@uSD)/LHact
                FE[i]  =UFfrictionSpringDamper(0,0,0, 0, vAct[i], 0, 0, 0)
    
    
            fig, ax = plt.subplots()
            plt.plot(vAct, FE, label='FF-NN', color='blue')
            plt.plot(dsM,Fm, label='MATLAB', color='red')
            plt.xlim(-0.1, 0.1)
            plt.ylim(-300*4,4*300)
            plt.grid(True)
            ax.set_xlabel('Velocity, m/s')
            ax.set_ylabel('$F_{\mu}$, N')
            plt.legend()
            plt.savefig('ExData/LiftFlexible/Friction.png', dpi=300)
            plt.show()
    
    
            fig, ax = plt.subplots()
            plt.plot(sE[:, 0], 1000*sE[:, 1], label='EXUDYN', color='blue')
            plt.plot(tspan, 1000*sM, label='MATLAB', color='red')
            plt.xlim(0, tEnd)
            plt.ylim(750, 1500)
            plt.grid(True)
            # Add a title and axis labels
            ax.set_xlabel('Time, s')
            ax.set_ylabel('Actuator position, mm')
            plt.legend()
            plt.savefig('ExData/LiftFlexible/ActuatorPosition.png', dpi=300)
            plt.show()


            fig, ax = plt.subplots()
            plt.plot(sE[:, 0], dsE[:, 2], label='EXUDYN', color='blue')
            plt.plot(tspan, dsM, label='MATLAB', color='red')
            plt.xlim(0, tEnd)
            plt.ylim(-0.1, 0.1)
            plt.grid(True)
            # Add a title and axis labels
            ax.set_xlabel('Time, s')
            ax.set_ylabel('Actuator velocity, m/s')
            plt.legend()
            plt.savefig('ExData/LiftFlexible/Actuatorvelocity.png', dpi=300)
            plt.show()
        
            #np.mean(dpE[:, 1])/1e6


            fig, ax = plt.subplots()
            plt.plot(sE[:, 0], dpE[:, 1], label='EXUDYN', color='blue')
            plt.plot(tspan, p1M, label='MATLAB', color='red')
            plt.xlim(0, tEnd)
            plt.ylim(-20e5, 100e5)
            plt.grid(True)
            # Add a title and axis labels
            ax.set_xlabel('Time, s')
            ax.set_ylabel('Pressure $p_1$, Pa')
            plt.legend()
            plt.savefig('ExData/LiftFlexible/p_1.png', dpi=300)
            plt.show()



            fig, ax = plt.subplots()
            plt.plot(sE[:, 0], dpE[:, 2], label='EXUDYN', color='blue')
            plt.plot(tspan, p2M, label='MATLAB', color='red')
            plt.xlim(0, tEnd)
            plt.grid(True)
            # Add a title and axis labels
            plt.ylim(-20e5, 125e5)
            ax.set_ylabel('Pressure $p_2$, Pa')
            plt.legend()
            plt.savefig('ExData/LiftFlexible/p_2.png', dpi=300)
            plt.show()
        
    return




def FlexiblePatuCrane(FineMesh, nModes, loadFromSavedNPY, ComputeModes, FreeModes, HCB, ModeAnimation, StaticCase, Hydraulics,
                   useFriction, Visualization,Plotting):
    
    global Fc, Fs, sig2, vs, Av0, Av1
    
    if FineMesh:
        
        fileNameL       = 'LiftBoom/ABAQUS/liftboom-free-050623' #To load fine mesh data from Abaqus
        fileNameT       = 'LiftBoom/ABAQUS/TiltBoom' #To load fine mesh data from AANSYS

        if not loadFromSavedNPY: 
                start_time                      = time.time()
                nodes                           = feL.ImportFromAbaqusInputFile(fileNameL+'.inp', typeName='Part', name='P000524_A_1-Nostopuomi_v2_fem_st')
                feL.ReadMassMatrixFromAbaqus(fileName=fileNameL + '_MASS2.mtx')             #Load mass matrix
                feL.ReadStiffnessMatrixFromAbaqus(fileName=fileNameL + '_STIF2.mtx')        #Load stiffness matrix
                feL.SaveToFile(fileNameL)
                print("--- saving LiftBoom FEM Abaqus data took: %s seconds ---" % (time.time() - start_time)) 
           
                if ComputeModes:
                    
                     from Models.ComputeModes import LiftBoomModes
                     LiftBoomModes(feL, FineMesh, nModes, FreeModes, HCB, ModeAnimation)
                                        
        else:       
                print('importing Abaqus FEM data structure of Lift Boom...')
                start_time = time.time()
                feL.LoadFromFile(fileNameL)
                cpuTime = time.time() - start_time
                print("--- importing FEM data took: %s seconds ---" % (cpuTime))
                
                # For computing modes
                
                if ComputeModes:
                    
                     from Models.ComputeModes import LiftBoomModes
                     LiftBoomModes(feL, FineMesh, nModes, FreeModes, HCB, ModeAnimation)
                     
                     
    else:
        
        fileNameL       = 'LiftBoom/ANSYS/LiftBoom' #To load fine mesh data from AANSYS
        fileNameT       = 'TiltBoom/ANSYS/TiltBoom' #To load fine mesh data from AANSYS
        
        if not loadFromSavedNPY: 
                start_time                          = time.time()
                
                # Lift Boom
                inputFileNameStiffnessMatrixL       = 'LiftBoom/ANSYS/'+'stiffnessMatrix.txt'
                inputFileNameMassMatrixL            = 'LiftBoom/ANSYS/'+'massMatrix.txt'
                inputFileNameNodalCoordinatesL      = 'LiftBoom/ANSYS/'+'nodalCoordinates.txt'
                inputFileNameElementsL              = 'LiftBoom/ANSYS/'+'elements.txt'
                inputFileNameNodalMappingVectorL    = 'LiftBoom/ANSYS/'+'nodalMappingVector.txt'                
                feL.ReadStiffnessMatrixFromAnsys(inputFileNameStiffnessMatrixL, inputFileNameNodalMappingVectorL, verbose=True)
                feL.ReadMassMatrixFromAnsys(inputFileNameMassMatrixL, inputFileNameNodalMappingVectorL, verbose=True)
                feL.ReadNodalCoordinatesFromAnsys(inputFileNameNodalCoordinatesL, verbose=True)
                feL.ReadElementsFromAnsys(inputFileNameElementsL, verbose=True)                
                feL.SaveToFile(fileNameL)
             
                # Tilt Boom
                inputFileNameStiffnessMatrixT       = 'TiltBoom/ANSYS/'+'stiffnessMatrix.txt'
                inputFileNameMassMatrixT            = 'TiltBoom/ANSYS/'+'massMatrix.txt'
                inputFileNameNodalCoordinatesT      = 'TiltBoom/ANSYS/'+'NodalCoordinates.txt'
                inputFileNameElementsT              = 'TiltBoom/ANSYS/'+'Elements.txt'
                inputFileNameNodalMappingVectorT    = 'TiltBoom/ANSYS/'+'nodalMappingVector.txt'                
                feT.ReadStiffnessMatrixFromAnsys(inputFileNameStiffnessMatrixT, inputFileNameNodalMappingVectorT, verbose=True)
                feT.ReadMassMatrixFromAnsys(inputFileNameMassMatrixT, inputFileNameNodalMappingVectorT, verbose=True)
                feT.ReadNodalCoordinatesFromAnsys(inputFileNameNodalCoordinatesT, verbose=True)
                feT.ReadElementsFromAnsys(inputFileNameElementsT, verbose=True)                
                feT.SaveToFile(fileNameT)
                                
                print("--- saving LiftBoom FEM ANSYS data took: %s seconds ---" % (time.time() - start_time)) 
                
                if ComputeModes:
                    
                     from Models.ComputeModes import LiftBoomModes
                     from Models.ComputeModes import TiltBoomModes
                     
                     LiftBoomModes(feL, FineMesh, nModes, FreeModes, HCB, ModeAnimation)
                     TiltBoomModes(feT, FineMesh=False, nModes=4, FreeModes=True, HCB=False, ModeAnimation=True)
        else:
                print('importing ANSYS FEM data structure...')      
                start_time = time.time()
                feL.LoadFromFile(fileNameL)  
                feT.LoadFromFile(fileNameT)               
                cpuTime = time.time() - start_time
                print("--- importing FEM data took: %s seconds ---" % (cpuTime))
        
                if ComputeModes:
                    
                     from Models.ComputeModes import LiftBoomModes
                     from Models.ComputeModes import TiltBoomModes
                     
                     LiftBoomModes(feL, FineMesh, nModes, FreeModes, HCB, ModeAnimation)
                     TiltBoomModes(feT, FineMesh, nModes, FreeModes, HCB, ModeAnimation)
    
    # Boundary condition at pillar
    p2                  = [0, 0,-10e-2]
    p1                  = [0, 0, 10e-2]
    radius1             = 2.5e-002
    nodeListJoint1      = feL.GetNodesOnCylinder(p1, p2, radius1, tolerance=1e-2) 
    pJoint1             = feL.GetNodePositionsMean(nodeListJoint1)
    nodeListJoint1Len   = len(nodeListJoint1)
    noodeWeightsJoint1  = [1/nodeListJoint1Len]*nodeListJoint1Len
           
    # Boundary condition at Piston 1
    p4                  = [0.3025,-0.1049,-10e-2]
    p3                  = [0.3025,-0.1049, 10e-2]
    radius2             = 3.6e-002
    nodeListPist1       = feL.GetNodesOnCylinder(p3, p4, radius2, tolerance=1e-2)  
    pJoint2             = feL.GetNodePositionsMean(nodeListPist1)
    nodeListPist1Len    = len(nodeListPist1)
    noodeWeightsPist1   = [1/nodeListPist1Len]*nodeListPist1Len
           
           
    # Boundary condition at cylinder 1
    p6                  = [1.265,0.2080,-8.02e-2]
    p5                  = [1.265,0.2080, 8.02e-2]
    radius3             = 3.2e-002
    nodeListCyl2        = feL.GetNodesOnCylinder(p5, p6, radius3, tolerance=1e-2)  
    pJoint3             = feL.GetNodePositionsMean(nodeListCyl2)
    nodeListCyl2Len     = len(nodeListCyl2)
    noodeWeightsCyl2    = [1/nodeListCyl2Len]*nodeListCyl2Len       
           
    # Boundary condition at Joint 2
    p8                  = [2.69,0.0066,-7.4e-2]
    p7                  = [2.69,0.0066, 7.4e-2]
    radius4             = 3.7e-002
    nodeListJoint2      = feL.GetNodesOnCylinder(p7, p8, radius4, tolerance=1e-2)  
    pJoint4             = feL.GetNodePositionsMean(nodeListJoint2)
    nodeListJoint2Len   = len(nodeListJoint2)
    noodeWeightsJoint2  = [1/nodeListJoint2Len]*nodeListJoint2Len
           
    # Joint 3
    p10                 = [2.89,0.0246,-7.4e-2]
    p9                  = [2.89,0.0246, 7.4e-2]
    radius5             = 5.2e-002
    nodeListJoint3      = feL.GetNodesOnCylinder(p9, p10, radius5, tolerance=1e-2)  
    pJoint5             = feL.GetNodePositionsMean(nodeListJoint3)
    nodeListJoint3Len   = len(nodeListJoint3)
    noodeWeightsJoint3  = [1/nodeListJoint3Len]*nodeListJoint3Len
    
    # Boundary condition at pillar
    p12                 = [0, 0,-9.63e-2]
    p11                 = [0, 0, 9.63e-2]
    radius6             = 4.8e-002
    nodeListJoint1T     = feT.GetNodesOnCylinder(p11, p12, radius6, tolerance=1e-2) 
    pJoint1T            = feT.GetNodePositionsMean(nodeListJoint1T)
    nodeListJoint1TLen  = len(nodeListJoint1T)
    noodeWeightsJoint1T = [1/nodeListJoint1TLen]*nodeListJoint1TLen
           
    # Boundary condition at Piston 1
    p14                 = [9.5e-2,0.24,-7.15e-2]
    p13                 = [9.5e-2,0.24, 7.15e-2]
    radius7             = 2.5e-002
    nodeListPist1T      = feT.GetNodesOnCylinder(p13, p14, radius7, tolerance=1e-2)  
    pJoint2T            = feT.GetNodePositionsMean(nodeListPist1T)
    nodeListPist1TLen   = len(nodeListPist1T)
    noodeWeightsPist1T  = [1/nodeListPist1TLen]*nodeListPist1TLen
                  
    # Boundary condition at cylinder 1
    # p18                 = [1.84,0.343,-4.80e-2]
    # p18s                = [1.84,0.339,-4.37e-2]
    # p17                 = [1.84,0.189,-4.80e-2]
    # p17s                = [1.84,0.193,-4.37e-2]
    # p16                 = [1.84,0.189, 4.80e-2]
    # p16s                = [1.84,0.193, 4.37e-2]
    # p15                 = [1.84,0.343, 4.80e-2]
    # p15s                = [1.84,0.339, 4.37e-2]
           
    # radius8             = 3.2e-002
    # nodeListCyl2T       = feT.GetNodesOn(p15, p16, radius8, tolerance=1e-2)  
    # pJoint3T            = feT.GetNodePositionsMean(nodeListCyl2T)
    # nodeListCyl2TLen    = len(nodeListCyl2T)
    # noodeWeightsCyl2T   = [1/nodeListCyl2TLen]*nodeListCyl2TLen       
     
    if not Hydraulics:
        # STEP 2: Craig-Bampton Modes
        print("Compute Craig-Bampton modes... ")
        boundaryListL   = [nodeListJoint1, nodeListJoint2,nodeListJoint3]
        boundaryListT  = [nodeListJoint1T, nodeListPist1T]
        
        start_time      = time.time()
        feL.ComputeHurtyCraigBamptonModes(boundaryNodesList=boundaryListL, nEigenModes=nModes, 
                                                    useSparseSolver=True,computationMode = HCBstaticModeSelection.RBE2)
        
        feT.ComputeHurtyCraigBamptonModes(boundaryNodesList=boundaryListT, nEigenModes=nModes, 
                                                    useSparseSolver=True,computationMode = HCBstaticModeSelection.RBE2)        
        
        print("Hurty-Craig Bampton modes... ")
        print("eigen freq. Lift Boom=", feL.GetEigenFrequenciesHz())
        print("eigen freq. Tilt Boom=", feT.GetEigenFrequenciesHz())
        print("HCB modes needed %.3f seconds" % (time.time() - start_time))  
    
        colLift = color4blue
           
        LiftBoom            = ObjectFFRFreducedOrderInterface(feL)
        TiltBoom            = ObjectFFRFreducedOrderInterface(feT)
        
        LiftBoomFFRF        = LiftBoom.AddObjectFFRFreducedOrder(mbs, positionRef=np.array([-0.09, 1.4261, 0]), 
                                              initialVelocity=[0,0,0], 
                                              initialAngularVelocity=[0,0,0],
                                              rotationMatrixRef  = RotationMatrixZ(mt.radians(0)),
                                              gravity=g,
                                              color=colLift,)
        
        TiltBoomFFRF        = TiltBoom.AddObjectFFRFreducedOrder(mbs, positionRef=np.array([-0.09, 1.4261, 0])+ 
                                                                 np.array([2.879420180699481+27e-3, -0.040690041435711005+8.3e-2, 0]), 
                                              initialVelocity=[0,0,0], 
                                              initialAngularVelocity=[0,0,0],
                                              rotationMatrixRef  = RotationMatrixZ(mt.radians(0)),
                                              gravity=g,
                                              color=colLift,)
        
        # 4th Body: Bracket 1
        L4              = 0.557227    # Length in x-direction
        H4              = 0.1425      # Height in y-direction
        W4              = 0.15        # Width in z-direction
        Bracket1L       = np.array([-0.09, 1.4261, 0]) + np.array([2.689524056550459, 0.00827551426741, 0])
        pMid4           = np.array([-0.257068, 0.004000, 0])  # center of mass, body0
        iCube4          = RigidBodyInertia(mass=11.524039, com=pMid4,
                                            inertiaTensor=np.array([[0.333066, 0.017355, 0],
                                                                    [0.017355, 0.081849, 0],
                                                                    [0,              0, 0.268644]]),
                                                                     inertiaTensorAtCOM=True)

        graphicsBody4   = GraphicsDataFromSTLfile(fileName4, color4blue,verbose=False, invertNormals=True,invertTriangles=True)
        graphicsBody4   = AddEdgesAndSmoothenNormals(graphicsBody4, edgeAngle=0.25*pi,addEdges=True, smoothNormals=True)
        graphicsCOM4    = GraphicsDataBasis(origin=iCube4.com, length=2*W4)
        [n4, b4]        = AddRigidBody(mainSys=mbs,inertia=iCube4,  # includes COM
                                        nodeType=exu.NodeType.RotationEulerParameters,
                                        position=Bracket1L,  # pMid2
                                        rotationMatrix=RotationMatrixZ(mt.radians(-33.195078970135967)),
                                        gravity=g,graphicsDataList=[graphicsCOM4, graphicsBody4])
        Marker15        = mbs.AddMarker(MarkerBodyRigid(bodyNumber=b4, localPosition=[0, 0, 0]))                        #With LIft Boom 
        Marker16        = mbs.AddMarker(MarkerBodyRigid(bodyNumber=b4, localPosition=[-0.457227402, 0.0425, 0]))  

       # 5th Body: Bracket 2
        L5              = 0.569009       # Length in x-direction
        H5              = 0.078827       # Height in y-direction
        W5              = 0.15           # Width in z-direction
        pMid5           = np.array([-0.212792, 0, 0])  # center of mass, body0
        Bracket1B       = np.array([-0.09, 1.4261, 0]) + np.array([2.329958953740425, 0.29034355756469298, 0])  #0.285710892999728-1.8*0.0425, -0.356968041652145+0.0525
        iCube5          = RigidBodyInertia(mass=7.900191, com=pMid5,
                                             inertiaTensor=np.array([[0.052095, 0, 0],
                                                                     [0,  0.260808, 0],
                                                                     [0,              0,  0.216772]]),
                                                                     inertiaTensorAtCOM=True)
        
        graphicsBody5   = GraphicsDataFromSTLfile(fileName5, color4blue,verbose=False, invertNormals=True,invertTriangles=True)
        graphicsBody5   = AddEdgesAndSmoothenNormals(graphicsBody5, edgeAngle=0.25*pi,addEdges=True, smoothNormals=True)
        graphicsCOM5    = GraphicsDataBasis(origin=iCube5.com, length=2*W5)
        [n5, b5]        = AddRigidBody(mainSys=mbs,inertia=iCube5,  # includes COM
                                         nodeType=exu.NodeType.RotationEulerParameters,
                                         position=Bracket1B,  # pMid2
                                         rotationMatrix=RotationMatrixZ(mt.radians(-181.7373) ),   #-5
                                         gravity=g,graphicsDataList=[graphicsCOM5, graphicsBody5])
        
        Marker18        = mbs.AddMarker(MarkerBodyRigid(bodyNumber=b5, localPosition=[0, 0, 0]))                        #With LIft Boom 
        Marker19        = mbs.AddMarker(MarkerBodyRigid(bodyNumber=b5, localPosition=[-0.475, 0, 0]))                   #With LIft Boom,-0.475 
                
        
        Marker7             = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber=LiftBoomFFRF['oFFRFreducedOrder'],
                                              meshNodeNumbers=np.array(nodeListJoint1), #these are the meshNodeNumbers
                                              weightingFactors=noodeWeightsJoint1))
        
        Marker10            = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber=LiftBoomFFRF['oFFRFreducedOrder'],
                                              meshNodeNumbers=np.array(nodeListJoint2), #these are the meshNodeNumbers
                                              weightingFactors=noodeWeightsJoint2))
        
        Marker11            = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber=LiftBoomFFRF['oFFRFreducedOrder'],
                                              meshNodeNumbers=np.array(nodeListJoint3), #these are the meshNodeNumbers
                                              weightingFactors=noodeWeightsJoint3))

        
        Marker13            = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber=TiltBoomFFRF['oFFRFreducedOrder'], 
                                                                    meshNodeNumbers=np.array(nodeListJoint1T), #these are the meshNodeNumbers
                                                                    weightingFactors=noodeWeightsJoint1T)) 

        Marker14            = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber=TiltBoomFFRF['oFFRFreducedOrder'], 
                                                                    meshNodeNumbers=np.array(nodeListPist1T), #these are the meshNodeNumbers
                                                                    weightingFactors=noodeWeightsPist1T))                        #With LIft Boom 
       
       
        if StaticCase:
            # Add Fixed Joint
            mbs.AddObject(GenericJoint(markerNumbers=[Marker4, Marker7],constrainedAxes=[1, 1, 1, 1, 1, 1],
                           visualization=VObjectJointGeneric(axesRadius=0.18*0.263342,axesLength=1.1*0.263342)))
            # Add Fixed Joint
            mbs.AddObject(GenericJoint(markerNumbers=[Marker11, Marker13],constrainedAxes=[1,1,1,1,1,1],
                         visualization=VObjectJointGeneric(axesRadius=0.18*0.263342,axesLength=1.1*0.263342)))  
            
            #Add fixed join between LiftBoom and Bracket 1
            mbs.AddObject(GenericJoint(markerNumbers=[Marker10, Marker15],constrainedAxes=[1,1,1,1,1,1],
                             visualization=VObjectJointGeneric(axesRadius=0.32*W4,axesLength=0.96*W4))) 
        
            #Add fixed joint between Bracket 1 and Bracket 2
            mbs.AddObject(GenericJoint(markerNumbers=[Marker16, Marker18],constrainedAxes=[1,1,1,1,1,1],
                             visualization=VObjectJointGeneric(axesRadius=0.23*W5,axesLength=1.0*W5)))
               
        else:
            # Add Revolute Joint
            mbs.AddObject(GenericJoint(markerNumbers=[Marker4, Marker7],constrainedAxes=[1,1,1,1,1,0],
                         visualization=VObjectJointGeneric(axesRadius=0.22*0.220,axesLength=0.72*0.220)))   
            
            # Add Revolute Joint
            mbs.AddObject(GenericJoint(markerNumbers=[Marker11, Marker13],constrainedAxes=[1,1,1,1,1,0],
                         visualization=VObjectJointGeneric(axesRadius=0.18*0.263342,axesLength=1.1*0.263342))) 
            
            # Revolute joint between LiftBoom and Bracket 1
            mbs.AddObject(GenericJoint(markerNumbers=[Marker10, Marker15],constrainedAxes=[1,1,1,1,1,0],
                             visualization=VObjectJointGeneric(axesRadius=0.32*W4,axesLength=0.96*W4))) 
        
            # Revolute joint between Bracket 1 and Bracket 2
            mbs.AddObject(GenericJoint(markerNumbers=[Marker16, Marker18],constrainedAxes=[1,1,1,1,1,0],
                             visualization=VObjectJointGeneric(axesRadius=0.23*W5,axesLength=1.0*W5)))
        
            # Revolute joint between Bracket 2 and TiltBoom
            mbs.AddObject(GenericJoint(markerNumbers=[Marker19, Marker14],constrainedAxes=[1,1,0,0,0,0],
                             visualization=VObjectJointGeneric(axesRadius=0.23*W5,axesLength=1.0*W5)))    
        
        
        # Joint 3
        p10                 = [2.89,0.0246,-7.4e-2]
        p9                  = [2.89,0.0246, 7.4e-2]
        pdef                = [2.89,0.0246, 0]
        radius5             = 5.2e-002
        nodeListJoint3      = feL.GetNodesOnCylinder(p9, p10, radius5, tolerance=1e-2)  
        pJoint5             = feL.GetNodePositionsMean(nodeListJoint3)
        nodeListJoint3Len   = len(nodeListJoint3)
        noodeWeightsJoint3  = [1/nodeListJoint3Len]*nodeListJoint3Len
        
        Marker11            = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber= LiftBoomFFRF['oFFRFreducedOrder'], 
                                                                 meshNodeNumbers=np.array(nodeListJoint3),
                                                                 #referencePosition=pdef,
                                                                 #useAlternativeApproach=altApproach,
                                                                weightingFactors=noodeWeightsJoint3))

           
        # Add Sensor for deflection
        Deflecion          = mbs.AddSensor(SensorMarker(markerNumber=Marker11, storeInternal=True, 
                                                     fileName='ExData/PatuFlexible/SimpleDef.txt', outputVariableType=exu.OutputVariableType.Displacement))
                      
    else:
        
        print("Compute Craig-Bampton modes... ")
        boundaryListL   = [nodeListJoint1, nodeListPist1, nodeListCyl2, nodeListJoint2,nodeListJoint3]
        boundaryListT  = [nodeListJoint1T, nodeListPist1T]
        
        start_time      = time.time()
        feL.ComputeHurtyCraigBamptonModes(boundaryNodesList=boundaryListL, nEigenModes=nModes, 
                                                    useSparseSolver=True,computationMode = HCBstaticModeSelection.RBE2)
        
        feT.ComputeHurtyCraigBamptonModes(boundaryNodesList=boundaryListT, nEigenModes=nModes, 
                                                    useSparseSolver=True,computationMode = HCBstaticModeSelection.RBE2)        
        
        print("Hurty-Craig Bampton modes... ")
        print("eigen freq. Lift Boom=", feL.GetEigenFrequenciesHz())
        print("eigen freq. Tilt Boom=", feT.GetEigenFrequenciesHz())
        print("HCB modes needed %.3f seconds" % (time.time() - start_time))  
    
        colLift = color4blue
           
        LiftBoom            = ObjectFFRFreducedOrderInterface(feL)
        TiltBoom            = ObjectFFRFreducedOrderInterface(feT)
        
        LiftBoomFFRF        = LiftBoom.AddObjectFFRFreducedOrder(mbs, positionRef=np.array([-0.09, 1.4261, 0]), 
                                              initialVelocity=[0,0,0], 
                                              initialAngularVelocity=[0,0,0],
                                              rotationMatrixRef  = RotationMatrixZ(mt.radians(0)),
                                              gravity=g,
                                              color=colLift,)
        
        TiltBoomFFRF        = TiltBoom.AddObjectFFRFreducedOrder(mbs, positionRef=np.array([-0.09, 1.4261, 0])+ 
                                                                 np.array([2.879420180699481+27e-3, -0.040690041435711005+8.3e-2,0]), 
                                              initialVelocity=[0,0,0], 
                                              initialAngularVelocity=[0,0,0],
                                              rotationMatrixRef  = RotationMatrixZ(mt.radians(0)),
                                              gravity=g,
                                              color=colLift,)
        
        # 4th Body: Bracket 1
        L4              = 0.557227    # Length in x-direction
        H4              = 0.1425      # Height in y-direction
        W4              = 0.15        # Width in z-direction
        Bracket1L       = np.array([-0.09, 1.4261, 0]) + np.array([2.689524056550459, 0.00827551426741, 0])
        pMid4           = np.array([-0.257068, 0.004000, 0])  # center of mass, body0
        iCube4          = RigidBodyInertia(mass=11.524039, com=pMid4,
                                            inertiaTensor=np.array([[0.333066, 0.017355, 0],
                                                                    [0.017355, 0.081849, 0],
                                                                    [0,              0, 0.268644]]),
                                                                     inertiaTensorAtCOM=True)

        graphicsBody4   = GraphicsDataFromSTLfile(fileName4, color4blue,verbose=False, invertNormals=True,invertTriangles=True)
        graphicsBody4   = AddEdgesAndSmoothenNormals(graphicsBody4, edgeAngle=0.25*pi,addEdges=True, smoothNormals=True)
        graphicsCOM4    = GraphicsDataBasis(origin=iCube4.com, length=2*W4)
        [n4, b4]        = AddRigidBody(mainSys=mbs,inertia=iCube4,  # includes COM
                                        nodeType=exu.NodeType.RotationEulerParameters,
                                        position=Bracket1L,  # pMid2
                                        rotationMatrix=RotationMatrixZ(mt.radians(-33.195078970135967)),
                                        gravity=g,graphicsDataList=[graphicsCOM4, graphicsBody4])
        Marker15        = mbs.AddMarker(MarkerBodyRigid(bodyNumber=b4, localPosition=[0, 0, 0]))                        #With LIft Boom 
        Marker16        = mbs.AddMarker(MarkerBodyRigid(bodyNumber=b4, localPosition=[-0.457227402, 0.0425, 0]))  

       # 5th Body: Bracket 2
        L5              = 0.569009       # Length in x-direction
        H5              = 0.078827       # Height in y-direction
        W5              = 0.15           # Width in z-direction
        pMid5           = np.array([-0.212792, 0, 0])  # center of mass, body0
        Bracket1B       = np.array([-0.09, 1.4261, 0]) + np.array([2.329958953740425, 0.29034355756469298, 0])  #0.285710892999728-1.8*0.0425, -0.356968041652145+0.0525
        iCube5          = RigidBodyInertia(mass=7.900191, com=pMid5,
                                             inertiaTensor=np.array([[0.052095, 0, 0],
                                                                     [0,  0.260808, 0],
                                                                     [0,              0,  0.216772]]),
                                                                     inertiaTensorAtCOM=True)
        
        graphicsBody5   = GraphicsDataFromSTLfile(fileName5, color4blue,verbose=False, invertNormals=True,invertTriangles=True)
        graphicsBody5   = AddEdgesAndSmoothenNormals(graphicsBody5, edgeAngle=0.25*pi,addEdges=True, smoothNormals=True)
        graphicsCOM5    = GraphicsDataBasis(origin=iCube5.com, length=2*W5)
        [n5, b5]        = AddRigidBody(mainSys=mbs,inertia=iCube5,  # includes COM
                                         nodeType=exu.NodeType.RotationEulerParameters,
                                         position=Bracket1B,  # pMid2
                                         rotationMatrix=RotationMatrixZ(mt.radians(-181.7373) ),   #-5
                                         gravity=g,graphicsDataList=[graphicsCOM5, graphicsBody5])
        
        Marker18        = mbs.AddMarker(MarkerBodyRigid(bodyNumber=b5, localPosition=[0, 0, 0]))                        #With LIft Boom 
        Marker19        = mbs.AddMarker(MarkerBodyRigid(bodyNumber=b5, localPosition=[-0.475, 0, 0]))                   #With LIft Boom,-0.475 
                
        
        Marker7             = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber=LiftBoomFFRF['oFFRFreducedOrder'],
                                              meshNodeNumbers=np.array(nodeListJoint1), #these are the meshNodeNumbers
                                              weightingFactors=noodeWeightsJoint1))
        Marker8         = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber=LiftBoomFFRF['oFFRFreducedOrder'],
                                              meshNodeNumbers=np.array(nodeListCyl2), #these are the meshNodeNumbers
                                              weightingFactors=noodeWeightsCyl2))           #With Cylinder 1
        
        Marker9         = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber=LiftBoomFFRF['oFFRFreducedOrder'],
                                              meshNodeNumbers=np.array(nodeListJoint2), #these are the meshNodeNumbers
                                              weightingFactors=noodeWeightsJoint2))        #With Cylinder 2
        
        
        Marker10            = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber=LiftBoomFFRF['oFFRFreducedOrder'],
                                              meshNodeNumbers=np.array(nodeListJoint2), #these are the meshNodeNumbers
                                              weightingFactors=noodeWeightsJoint2))
        
        Marker11            = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber=LiftBoomFFRF['oFFRFreducedOrder'],
                                              meshNodeNumbers=np.array(nodeListJoint3), #these are the meshNodeNumbers
                                              weightingFactors=noodeWeightsJoint3))

        
        Marker13            = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber=TiltBoomFFRF['oFFRFreducedOrder'], 
                                                                    meshNodeNumbers=np.array(nodeListJoint1T), #these are the meshNodeNumbers
                                                                    weightingFactors=noodeWeightsJoint1T)) 

        Marker14            = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber=TiltBoomFFRF['oFFRFreducedOrder'], 
                                                                    meshNodeNumbers=np.array(nodeListPist1T), #these are the meshNodeNumbers
                                                                    weightingFactors=noodeWeightsPist1T))                        #With LIft Boom
        
        # Add Revolute Joint
        mbs.AddObject(GenericJoint(markerNumbers=[Marker4, Marker7],constrainedAxes=[1,1,1,1,1,0],
                         visualization=VObjectJointGeneric(axesRadius=0.22*0.220,axesLength=0.72*0.220)))   
            
        # Add Revolute Joint
        mbs.AddObject(GenericJoint(markerNumbers=[Marker11, Marker13],constrainedAxes=[1,1,1,1,1,0],
                         visualization=VObjectJointGeneric(axesRadius=0.18*0.263342,axesLength=1.1*0.263342))) 
            
        # Revolute joint between LiftBoom and Bracket 1
        mbs.AddObject(GenericJoint(markerNumbers=[Marker10, Marker15],constrainedAxes=[1,1,1,1,1,0],
                             visualization=VObjectJointGeneric(axesRadius=0.32*W4,axesLength=0.96*W4))) 
        
        # Revolute joint between Bracket 1 and Bracket 2
        mbs.AddObject(GenericJoint(markerNumbers=[Marker16, Marker18],constrainedAxes=[1,1,1,1,1,0],
                             visualization=VObjectJointGeneric(axesRadius=0.23*W5,axesLength=1.0*W5)))
        
        # Revolute joint between Bracket 2 and TiltBoom
        mbs.AddObject(GenericJoint(markerNumbers=[Marker19, Marker14],constrainedAxes=[1,1,0,0,0,0],
                             visualization=VObjectJointGeneric(axesRadius=0.23*W5,axesLength=1.0*W5))) 
        
        colCyl              = color4orange
        colPis              = color4grey 
        
        import scipy.io
        
        Data            = scipy.io.loadmat('ExData/PATU@2DOF')
        DataM           = Data['Xt']
        DataU1          = Data['U1']
        DataU2          = Data['U2']
        tspan           = Data['t']
        sM1             = Data['s1']
        sM2             = Data['s2']
        a1M             = DataM[0, :]
        a2M             = DataM[1, :]
        da1M            = DataM[2, :]
        da2M            = DataM[3, :]       
        p1M             = DataM[4, :]
        p2M             = DataM[5, :]
        p3M             = DataM[6, :]
        p4M             = DataM[7, :]
        U1              = DataU1[0, :]
        U2              = DataU2[0, :]
        pS              = 80e5
        pT              = 1e5                           # Tank pressure
        Qn1             = (24/60000)/(sqrt(35e5)*9.9)                      # Nominal flow rate of valve at 18 l/min under
        Qn2             = (54/60000)/(sqrt(35e5)*9.9)                      # Nominal flow rate of valve at 18 l/min under
        
        # Cylinder and piston parameters
        L_Cyl1          = 820e-3                            # Cylinder length
        D_Cyl1          = 100e-3                            # Cylinder dia
        A_1             = (pi/4)*(D_Cyl1)**2                # Area of cylinder side
        L_Pis1          = 535e-3                            # Piston length, also equals to stroke length
        d_pis1          = 56e-3                             # Piston dia
        A_2             = A_1-(pi/4)*(d_pis1)**2            # Area on piston-rod side
        L_Cyl2          = 1050e-3                           # Cylinder length
        L_Pis2          = 780e-3                            # Piston length, also equals to stroke length
        d_1             = 12.7e-3                         # Dia of volume 1
        V1              = (pi/4)*(d_1)**2*1.5             # Volume V1 = V0
        V2              = (pi/4)*(d_1)**2*1.5             # Volume V2 = V1
        A               = [A_1, A_2]
        Bh              = 700e6
        Bc              = 2.1000e+11
        Bo              = 1650e6
        Fc              = 210
        Fs              = 300
        sig2            = 330
        vs              = 5e-3
        
        #add hydraulics actuator:
        colCyl              = color4orange
        colPis              = color4grey 
        LH1                 = L_Cyl1                        #zero length of actuator
        LH2                 = L_Cyl2                        #zero length of actuator
        
        
        nODE1               = mbs.AddNode(NodeGenericODE1(referenceCoordinates=[0, 0],
                                    initialCoordinates=[5.346679683966570e+06,
                                                        3.100522576711962e+06],  # initialize with 20 bar
                                                        numberOfODE1Coordinates=2))

        nODE2               = mbs.AddNode(NodeGenericODE1(referenceCoordinates=[0, 0],
                                    initialCoordinates=[4.806409161717769e+06,
                                                        7.714843744901530e+06],  # initialize with 20 bar
                                                        numberOfODE1Coordinates=2))
        def UFfrictionSpringDamper1(mbs, t, itemIndex, u, v, k, d, f0):

            return   4*(Fc*tanh(4*(abs(v    )/vs))+(Fs-Fc)*((abs(v    )/vs)/((1/4)*(abs(v    )/vs)**2+3/4)**2))*np.sign(v )+sig2*v    *tanh(4)

        
        def UFfrictionSpringDamper2(mbs, t, itemIndex, u, v, k, d, f0):

            return   (Fc*tanh(4*(abs(v    )/vs))+(Fs-Fc)*((abs(v    )/vs)/((1/4)*(abs(v    )/vs)**2+3/4)**2))*np.sign(v )+sig2*v    *tanh(4)
        
        if useFriction:
            oFriction1       = mbs.AddObject(ObjectConnectorSpringDamper(markerNumbers=[Marker5, Marker8], referenceLength=0.001,stiffness=0,
                                                damping=0, force=0, velocityOffset = 0., activeConnector = True,
                                                springForceUserFunction=UFfrictionSpringDamper1,
                                                  visualization=VSpringDamper(show=False) ))
            
            oFriction2       = mbs.AddObject(ObjectConnectorSpringDamper(markerNumbers=[Marker9, Marker18], referenceLength=0.001,stiffness=0,
                                                 damping=0, force=0, velocityOffset = 0., activeConnector = True,
                                                 springForceUserFunction=UFfrictionSpringDamper2,
                                                   visualization=VSpringDamper(show=False) ))


        oHA1                = mbs.AddObject(HydraulicActuatorSimple(name='LiftCylinder', markerNumbers=[ Marker5, Marker8], 
                                            nodeNumbers=[nODE1], offsetLength=L_Pis1, strokeLength=LH1, chamberCrossSection0=A[0], 
                                            chamberCrossSection1=A[1], hoseVolume0=V1, hoseVolume1=V2, valveOpening0=0, 
                                            valveOpening1=0, actuatorDamping=1e5*0, oilBulkModulus=Bo, cylinderBulkModulus=Bc, 
                                            hoseBulkModulus=Bh, nominalFlow=Qn1, systemPressure=pS, tankPressure=pT, 
                                            useChamberVolumeChange=True, activeConnector=True, 
                                            visualization={'show': True, 'cylinderRadius': 50e-3, 'rodRadius': 28e-3, 
                                                            'pistonRadius': 0.04, 'pistonLength': 0.001, 'rodMountRadius': 0.0, 
                                                            'baseMountRadius': 20.0e-3, 'baseMountLength': 20.0e-3, 'colorCylinder': color4orange,
                                                            'colorPiston': color4grey}))
        
        oHA2 = mbs.AddObject(HydraulicActuatorSimple(name='TiltCylinder', markerNumbers=[Marker9, Marker18], 
                                             nodeNumbers=[nODE2], offsetLength=L_Pis2, strokeLength=LH2, chamberCrossSection0=A[0], 
                                             chamberCrossSection1=A[1], hoseVolume0=V1, hoseVolume1=V2, valveOpening0=0, 
                                             valveOpening1=0, actuatorDamping=1e5*0, oilBulkModulus=Bo, cylinderBulkModulus=Bc, 
                                             hoseBulkModulus=Bh, nominalFlow=Qn2, systemPressure=pS, tankPressure=pT, 
                                             useChamberVolumeChange=True, activeConnector=True, 
                                             visualization={'show': True, 'cylinderRadius': 50e-3, 'rodRadius': 28e-3, 
                                                             'pistonRadius': 0.04, 'pistonLength': 0.001, 'rodMountRadius': 0.0, 
                                                             'baseMountRadius': 0.0, 'baseMountLength': 0.0, 'colorCylinder': color4orange,
                                                             'colorPiston': color4grey}))
        
        def PreStepUserFunction(mbs, t):
            
            Av0 = U1[mt.trunc(t/h)]
            Av1 = -Av0
            Av2 = U2[mt.trunc(t/h)]
            Av3 = -Av2

            mbs.SetObjectParameter(oHA1, "valveOpening0", Av0)
            mbs.SetObjectParameter(oHA1, "valveOpening1", Av1)
            mbs.SetObjectParameter(oHA2, "valveOpening0", Av2)
            mbs.SetObjectParameter(oHA2, "valveOpening1", Av3)
            return True

        mbs.SetPreStepUserFunction(PreStepUserFunction)
        
       # Add sensors in the system
        sForce1         = mbs.AddSensor(SensorObject(objectNumber=oHA1, storeInternal=True, 
                                                     fileName='ExData/PatuFlexible/sForce1.txt', outputVariableType=exu.OutputVariableType.Force))
        sForce2         = mbs.AddSensor(SensorObject(objectNumber=oHA2, storeInternal=True, 
                                                     fileName='ExData/PatuFlexible/sForce2.txt', outputVariableType=exu.OutputVariableType.Force))
        sDistance1      = mbs.AddSensor(SensorObject(objectNumber=oHA1, storeInternal=True, 
                                                     fileName='ExData/PatuFlexible/sDistance1.txt', outputVariableType=exu.OutputVariableType.Distance))
        sDistance2      = mbs.AddSensor(SensorObject(objectNumber=oHA2, storeInternal=True, 
                                                     fileName='ExData/PatuFlexible/sDistance2.txt', outputVariableType=exu.OutputVariableType.Distance))
        sVelocity1      = mbs.AddSensor(SensorObject(objectNumber=oHA1, storeInternal=True,
                                                     fileName='ExData/PatuRigid/sVelocity1.txt', outputVariableType=exu.OutputVariableType.VelocityLocal))
        sVelocity2      = mbs.AddSensor(SensorObject(objectNumber=oHA2, storeInternal=True,
                                                     fileName='ExData/PatuFlexible/sVelocity2.txt', outputVariableType=exu.OutputVariableType.VelocityLocal))  
        sPressures1     = mbs.AddSensor(SensorNode(nodeNumber=nODE1, storeInternal=True, 
                                                   fileName='ExData/PatuFlexible/sPressures1.txt',outputVariableType=exu.OutputVariableType.Coordinates))
        sPressures2     = mbs.AddSensor(SensorNode(nodeNumber=nODE2, storeInternal=True, 
                                                   fileName='ExData/PatuFlexible/sPressures2.txt',outputVariableType=exu.OutputVariableType.Coordinates))
    mbs.Assemble()
    simulationSettings = exu.SimulationSettings()

    if not StaticCase:
        simulationSettings.timeIntegration.numberOfSteps            = nSteps
        simulationSettings.timeIntegration.endTime                  = tEnd
        simulationSettings.timeIntegration.verboseMode              = 1
        #simulationSettings.timeIntegration.simulateInRealtime       = True
        simulationSettings.solutionSettings.solutionWritePeriod     = h
        simulationSettings.solutionSettings.sensorsWritePeriod      = h
        simulationSettings.linearSolverType                         = exu.LinearSolverType.EigenSparse
        simulationSettings.timeIntegration.newton.useModifiedNewton = True
        simulationSettings.timeIntegration.stepInformation         += 8
        #simulationSettings.displayComputationTime                   = True
        simulationSettings.displayStatistics                        = True             
        
        exu.SolveDynamic(mbs, simulationSettings=simulationSettings,
                 solverType=exu.DynamicSolverType.TrapezoidalIndex2)
       
       # mbs.SolveDynamic(simulationSettings)  
    else:
        simulationSettings.solutionSettings.solutionWritePeriod = 2e-2  #output interval general
        simulationSettings.solutionSettings.sensorsWritePeriod = 1e-1  #output interval of sensors
        simulationSettings.timeIntegration.numberOfSteps = int(tEnd/h) #must be integer
        simulationSettings.timeIntegration.endTime = tEnd
        simulationSettings.solutionSettings.coordinatesSolutionFileName = "staticSolution.txt"
        simulationSettings.solutionSettings.appendToFile = False
        simulationSettings.staticSolver.newton.numericalDifferentiation.relativeEpsilon = 1e-4
        simulationSettings.staticSolver.newton.relativeTolerance = 1e-6
        simulationSettings.staticSolver.newton.absoluteTolerance = 1e-1
        simulationSettings.staticSolver.verboseMode = 2
        simulationSettings.timeIntegration.newton.useModifiedNewton = True
        simulationSettings.linearSolverType = exu.LinearSolverType.EigenSparse
        simulationSettings.staticSolver.numberOfLoadSteps = 5
        
        simulationSettings.staticSolver.adaptiveStep = True
        
        mbs.SolveStatic(simulationSettings)
    
    if Visualization:
        SC.visualizationSettings.window.renderWindowSize            = [1600, 1200]        
        SC.visualizationSettings.openGL.multiSampling               = 4        
        SC.visualizationSettings.openGL.lineWidth                   = 3  
        SC.visualizationSettings.general.autoFitScene               = False      
        SC.visualizationSettings.nodes.drawNodesAsPoint             = False        
        SC.visualizationSettings.nodes.showBasis                    = True  
        
        from exudyn.interactive import SolutionViewer
        SolutionViewer(mbs) 
    
    #%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    if Plotting:
        import matplotlib.pyplot as plt
        
        if not Hydraulics:
            TipDef  = mbs.GetSensorStoredData(Deflecion)
            
            fig, ax = plt.subplots()
            plt.plot(TipDef[:,0], 1000*TipDef[:,2], label='Exudyn', color='blue')
            plt.xlim(0, max(TipDef[:,0]))
            plt.ylim(-1,1)
            plt.grid(True)
            ax.set_xlabel('Time, s')
            ax.set_ylabel('Deflection, mm')
            plt.legend()
            plt.savefig('ExData/SimpleDeflection.png', dpi=300)
            plt.show()
        
        else: 
            #PLOTTING Results
            dE      = mbs.GetSensorStoredData(sDisplacement)
            sE      = mbs.GetSensorStoredData(sDistance)
            dsE     = mbs.GetSensorStoredData(sVelocity)
            dpE     = mbs.GetSensorStoredData(sPressures)
            dFE     = mbs.GetSensorStoredData(sForce)  
    
            vAct    = np.zeros((len(sE),1))
            FE      = np.zeros((len(sE),1))
      
            for i in range(len(sE)):
                LHact   = sE[i,1]
                uSD     = dE[i,1:4]
                vSD     = dsE[i,1:4]
                vAct[i]=( vSD@uSD)/LHact
                FE[i]  =UFfrictionSpringDamper(0,0,0, 0, vAct[i], 0, 0, 0)
    
    
            fig, ax = plt.subplots()
            plt.plot(vAct, FE, label='FF-NN', color='blue')
            plt.plot(dsM,Fm, label='MATLAB', color='red')
            plt.xlim(-0.1, 0.1)
            plt.ylim(-300*4,4*300)
            plt.grid(True)
            ax.set_xlabel('Velocity, m/s')
            ax.set_ylabel('$F_{\mu}$, N')
            plt.legend()
            plt.savefig('ExData/Friction.png', dpi=300)
            plt.show()
    
    
            fig, ax = plt.subplots()
            plt.plot(sE[:, 0], 1000*sE[:, 1], label='EXUDYN', color='blue')
            plt.plot(tspan, 1000*sM, label='MATLAB', color='red')
            plt.xlim(0, tEnd)
            plt.ylim(750, 1500)
            plt.grid(True)
            # Add a title and axis labels
            ax.set_xlabel('Time, s')
            ax.set_ylabel('Actuator position, mm')
            plt.legend()
            plt.savefig('ExData/ActuatorPosition.png', dpi=300)
            plt.show()


            fig, ax = plt.subplots()
            plt.plot(sE[:, 0], dsE[:, 2], label='EXUDYN', color='blue')
            plt.plot(tspan, dsM, label='MATLAB', color='red')
            plt.xlim(0, tEnd)
            plt.ylim(-0.1, 0.1)
            plt.grid(True)
            # Add a title and axis labels
            ax.set_xlabel('Time, s')
            ax.set_ylabel('Actuator velocity, m/s')
            plt.legend()
            plt.savefig('ExData/Actuatorvelocity.png', dpi=300)
            plt.show()
        
            #np.mean(dpE[:, 1])/1e6


            fig, ax = plt.subplots()
            plt.plot(sE[:, 0], dpE[:, 1], label='EXUDYN', color='blue')
            plt.plot(tspan, p1M, label='MATLAB', color='red')
            plt.xlim(0, tEnd)
            plt.ylim(-20e5, 100e5)
            plt.grid(True)
            # Add a title and axis labels
            ax.set_xlabel('Time, s')
            ax.set_ylabel('Pressure $p_1$, Pa')
            plt.legend()
            plt.savefig('ExData/p_1.png', dpi=300)
            plt.show()



            fig, ax = plt.subplots()
            plt.plot(sE[:, 0], dpE[:, 2], label='EXUDYN', color='blue')
            plt.plot(tspan, p2M, label='MATLAB', color='red')
            plt.xlim(0, tEnd)
            plt.grid(True)
            # Add a title and axis labels
            plt.ylim(-20e5, 125e5)
            ax.set_ylabel('Pressure $p_2$, Pa')
            plt.legend()
            plt.savefig('ExData/p_2.png', dpi=300)
            plt.show()
        
    return