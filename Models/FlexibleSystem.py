#%%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#This files defines the simulation of hydraulically actuated flexible structure.

# Author        : Qasim Khadim
# Contact       : qasim.khadim@outlook.com,qkhadim22 (Github)
# Dated         : 02-05-2023
# Organization  : University of Oulu in the collaboration of LUT University and University of Innsbruck.

# Copyright     :
#%%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

import sys
sys.exudynFast = True

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
#feT             = FEMinterface()

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


def FFRFHydraulics(FineMesh, nModes, loadFromSavedNPY, ComputeModes, FreeModes, HCB, ModeAnimation, Hydraulics,
                   useFriction, Visualization,Plotting):
    
    global Fc, Fs, sig2, vs, Av0, Av1
    
    if FineMesh:
        
        fileNameL       = 'LiftBoom/ABAQUS/LiftBoom-freebody' #To load fine mesh data from Abaqus
        
        if not loadFromSavedNPY: 
                start_time                      = time.time()
                nodes                           = feL.ImportFromAbaqusInputFile(fileNameL+'.inp', typeName='Part', name='P000524_A_1-Nostopuomi_v2')
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
    
    # STEP 2: Craig-Bampton Modes
    boundaryList        = [nodeListJoint1, nodeListPist1]
    start_time          = time.time()
    feL.ComputeHurtyCraigBamptonModes(boundaryNodesList=boundaryList, nEigenModes=nModes, 
                                                    useSparseSolver=True,computationMode = HCBstaticModeSelection.RBE2)
    print("Hurty-Craig Bampton modes... ")
    print("eigen freq.=", feL.GetEigenFrequenciesHz())
    print("HCB modes needed %.3f seconds" % (time.time() - start_time))  
           
    LiftBoom            = ObjectFFRFreducedOrderInterface(feL)
    LiftBoomFFRF        = LiftBoom.AddObjectFFRFreducedOrder(mbs, positionRef=np.array([-0.09, 1.4261, 0]), 
                                              initialVelocity=[0,0,0], 
                                              initialAngularVelocity=[0,0,0],
                                              rotationMatrixRef  = RotationMatrixZ(mt.radians(-27.227343749651500)),
                                              gravity=g,
                                              color=[0.1,0.9,0.1,1.],)
    
    Marker7             = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber=LiftBoomFFRF['oFFRFreducedOrder'],
                                              meshNodeNumbers=np.array(nodeListJoint1), #these are the meshNodeNumbers
                                              weightingFactors=noodeWeightsJoint1))
    
    Marker8             = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber=LiftBoomFFRF['oFFRFreducedOrder'],
                                              meshNodeNumbers=np.array(nodeListPist1), #these are the meshNodeNumbers
                                              weightingFactors=noodeWeightsPist1))
    
    if not Hydraulics:
        
        #Revolute Joint
        mbs.AddObject(GenericJoint(markerNumbers=[Marker4, Marker7],constrainedAxes=[1,1,1,1,1,0],
                           visualization=VObjectJointGeneric(axesRadius=0.18*0.263342,axesLength=1.1*0.263342)))
    else:
        
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
                                                     fileName='ExData/sForce.txt', outputVariableType=exu.OutputVariableType.Force))
        sDistance       = mbs.AddSensor(SensorObject(objectNumber=oHA1, storeInternal=True, 
                                                     fileName='ExData/sDistance.txt', outputVariableType=exu.OutputVariableType.Distance))
        sDisplacement   = mbs.AddSensor(SensorObject(objectNumber=oHA1, storeInternal=True,
                                                     fileName='ExData/sDisplacement.txt', outputVariableType=exu.OutputVariableType.Displacement))
        sVelocity       = mbs.AddSensor(SensorObject(objectNumber=oHA1, storeInternal=True,
                                                     fileName='ExData/sVelocity.txt', outputVariableType=exu.OutputVariableType.Velocity))
        sPressures      = mbs.AddSensor(SensorNode(nodeNumber=nODE1, storeInternal=True, 
                                                   fileName='ExData/sPressures.txt',outputVariableType=exu.OutputVariableType.Coordinates))   
    
        
    mbs.Assemble()
    simulationSettings = exu.SimulationSettings()
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

    SC.visualizationSettings.window.renderWindowSize            = [1600, 1200]        
    SC.visualizationSettings.openGL.multiSampling               = 4        
    SC.visualizationSettings.openGL.lineWidth                   = 3  
    SC.visualizationSettings.general.autoFitScene               = False      
    SC.visualizationSettings.nodes.drawNodesAsPoint             = False        
    SC.visualizationSettings.nodes.showBasis                    = True     

    exu.SolveDynamic(mbs, simulationSettings=simulationSettings,
                 solverType=exu.DynamicSolverType.TrapezoidalIndex2)
    
    if Visualization:
        
        from exudyn.interactive import SolutionViewer
        SolutionViewer(mbs) 
    
    #%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    if Plotting:
        import matplotlib.pyplot as plt
        
        #PLOTTING Results
        dE  = mbs.GetSensorStoredData(sDisplacement)
        sE  = mbs.GetSensorStoredData(sDistance)
        dsE = mbs.GetSensorStoredData(sVelocity)
        dpE = mbs.GetSensorStoredData(sPressures)
        dFE = mbs.GetSensorStoredData(sForce)  
    
        vAct = np.zeros((len(sE),1))
        FE= np.zeros((len(sE),1))
      
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