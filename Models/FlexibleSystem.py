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
g           = [0, -9.8066, 0]  # Gravity
tEnd        = 30  # simulation time
h           = 1e-3  # step size
nSteps      = int(tEnd/h)+2
colLift     = color4blue

# Loading Graphics of bodies
fileName1       = 'Graphics_Exudyn/LowQualityPATU/Pillar.stl'
fileName4       = 'Graphics_Exudyn/HighQualityPATU/Bracket1.stl'
fileName5       = 'Graphics_Exudyn/HighQualityPATU/Bracket2.stl'

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
LiftP           = np.array([-0.09, 1.4261, 0])
TiltL           = np.array([2.879420180699481-5e-3, -0.040690041435711005, 0])  #np.array([2.879420180699481-5e-3, -0.040690041435711005+6.2e-2, 0])
Bracket1        = np.array([2.689524056550459, -0.066027551426741+15e-3, 0])
Bracket2        = np.array([-0.382, 0.28, 0]) 
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


def FlexiblePatuCrane(FineMesh, nModes, loadFromSavedNPY, ComputeModes, FreeModes, HCB, ModeAnimation, StaticCase, Hydraulics,
                   useFriction, Visualization,Plotting):
    
    if FineMesh:
        fileNameL       = 'LiftBoom/ABAQUS/liftboom-free-050623'        #To load fine mesh data from ABAQUS
        fileNameT       = 'TiltBoom/ABAQUS/tiltboom-120623'             #To load fine mesh data from ABAQUS
        
        if not loadFromSavedNPY: 
                start_time                      = time.time()
                nodesL                          = feL.ImportFromAbaqusInputFile(fileNameL+'.inp', typeName='Part', name='P000524_A_1-Nostopuomi_v2_fem_st')
                feL.ReadMassMatrixFromAbaqus(fileName=fileNameL + '_MASS2.mtx')             #Load mass matrix
                feL.ReadStiffnessMatrixFromAbaqus(fileName=fileNameL + '_STIF2.mtx')        #Load stiffness matrix
                feL.SaveToFile(fileNameL)
                
                nodesT                          = feT.ImportFromAbaqusInputFile(fileNameT+'.inp', typeName='Part', name='P000516_A_1-Taittopuomi_fem')
                feT.ReadMassMatrixFromAbaqus(fileName=fileNameT + '_MASS1.mtx')             #Load mass matrix
                feT.ReadStiffnessMatrixFromAbaqus(fileName=fileNameT + '_STIF1.mtx')        #Load stiffness matrix
                feT.SaveToFile(fileNameT)
                
                print("--- saving LiftBoom and TiltBoom FEM Abaqus data took: %s seconds ---" % (time.time() - start_time)) 
           
                if ComputeModes:
                    
                     from Models.ComputeModes import LiftBoomModes
                     from Models.ComputeModes import TiltBoomModes
                     
                     LiftBoomModes(feL, FineMesh, nModes, FreeModes, HCB, ModeAnimation)
                     TiltBoomModes(feT, FineMesh=False, nModes=4, FreeModes=True, HCB=False, ModeAnimation=True)
                                        
        else:       
                print('importing Abaqus FEM data structure of Lift Boom...')
                start_time = time.time()
                feL.LoadFromFile(fileNameL)
                feT.LoadFromFile(fileNameT)
                cpuTime = time.time() - start_time
                print("--- importing FEM data took: %s seconds ---" % (cpuTime))
                
                # For computing modes
                
                if ComputeModes:
                    
                     from Models.ComputeModes import LiftBoomModes
                     from Models.ComputeModes import TiltBoomModes
                     
                     #LiftBoomModes(feL, FineMesh, nModes, FreeModes, HCB, ModeAnimation)
                     TiltBoomModes(feT, FineMesh=True, nModes=2, FreeModes=True, HCB=False, ModeAnimation=True)
        
              # Boundary condition at pillar
        p2                  = [0, 0,-10e-2]
        p1                  = [0, 0, 10e-2]
        radius1             = 2.5e-002
        nodeListJoint1      = feL.GetNodesOnCylinder(p1, p2, radius1, tolerance=1e-2) 
        pJoint1             = feL.GetNodePositionsMean(nodeListJoint1)
        nodeListJoint1Len   = len(nodeListJoint1)
        noodeWeightsJoint1  = [1/nodeListJoint1Len]*nodeListJoint1Len 
        
        # Boundary condition at Joint 2
        p8                  = [2.69+15e-3,2.38e-02,-7.4e-2]
        p7                  = [2.69+15e-3,2.38e-02, 7.4e-2]
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
        p12                 = [9.92e-15, -0.85e-3,-9.63e-2]
        p11                 = [9.92e-15, -0.85e-3, 9.63e-2]
        radius6             = 4.82e-002
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
        
        print("Compute Craig-Bampton modes... ")
        boundaryListL   = [nodeListJoint1, nodeListJoint2,  nodeListJoint3] 
        boundaryListT  = [nodeListJoint1T]
         
        start_time      = time.time()
        feL.ComputeHurtyCraigBamptonModes(boundaryNodesList=boundaryListL, nEigenModes=nModes, 
                                                    useSparseSolver=True,computationMode = HCBstaticModeSelection.RBE2)
        
        feT.ComputeHurtyCraigBamptonModes(boundaryNodesList=boundaryListT, nEigenModes=nModes, 
                                                    useSparseSolver=True,computationMode = HCBstaticModeSelection.RBE2) 
        
        print("Hurty-Craig Bampton modes... ")
        print("eigen freq. Lift Boom=", feL.GetEigenFrequenciesHz())
        print("eigen freq. Tilt Boom=", feT.GetEigenFrequenciesHz())
        print("HCB modes needed %.3f seconds" % (time.time() - start_time)) 
        
        LiftBoom            = ObjectFFRFreducedOrderInterface(feL)
        TiltBoom            = ObjectFFRFreducedOrderInterface(feT)
        
        LiftBoomFFRF        = LiftBoom.AddObjectFFRFreducedOrder(mbs, positionRef=LiftP, 
                                             initialVelocity=[0,0,0], 
                                              initialAngularVelocity=[0,0,0],
                                              rotationMatrixRef  = RotationMatrixZ(mt.radians(-1.227343749651500)),
                                              gravity=g,
                                              color=colLift,)
        
        TiltBoomFFRF        = TiltBoom.AddObjectFFRFreducedOrder(mbs, positionRef=LiftP+TiltL, #2.879420180699481+27e-3, -0.040690041435711005+8.3e-2, 0
                                              initialVelocity=[0,0,0], 
                                              initialAngularVelocity=[0,0,0],
                                              rotationMatrixRef  = RotationMatrixZ(mt.radians(15.82031486865900)),
                                              gravity=g,
                                              color=colLift,)
        
        # 4th Body: Bracket 1
        L4              = 0.557227    # Length in x-direction
        H4              = 0.1425      # Height in y-direction
        W4              = 0.15        # Width in z-direction
        Bracket1L       = LiftP + Bracket1
        pMid4           = np.array([0.004000, 0.257068, 0])  # center of mass, body0,0.004000,-0.257068
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
                                         rotationMatrix=RotationMatrixZ(mt.radians(58.803311277817050)),
                                         gravity=g,graphicsDataList=[graphicsCOM4, graphicsBody4])
        
        # 5th Body: Bracket 2
        L5              = 0.569009       # Length in x-direction
        H5              = 0.078827       # Height in y-direction
        W5              = 0.15           # Width in z-direction
        pMid5           = np.array([0.212792, 0, 0])  # center of mass, body0
        Bracket1B       = LiftP + Bracket1+ Bracket2  #0.285710892999728-1.8*0.0425, -0.356968041652145+0.0525
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
                                          rotationMatrix=RotationMatrixZ(mt.radians(-3.151443859554505) ),   #-5
                                          gravity=g,graphicsDataList=[graphicsCOM5, graphicsBody5])
        
        Marker7         = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber=LiftBoomFFRF['oFFRFreducedOrder'],
                                             meshNodeNumbers=np.array(nodeListJoint1), #these are the meshNodeNumbers
                                              weightingFactors=noodeWeightsJoint1))
        
        Marker10        = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber=LiftBoomFFRF['oFFRFreducedOrder'],
                                              meshNodeNumbers=np.array(nodeListJoint2), #these are the meshNodeNumbers
                                              weightingFactors=noodeWeightsJoint2))      
        
        Marker11        = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber=LiftBoomFFRF['oFFRFreducedOrder'],
                                            meshNodeNumbers=np.array(nodeListJoint3), #these are the meshNodeNumbers
                                            weightingFactors=noodeWeightsJoint3))
        
        Marker13        = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber=TiltBoomFFRF['oFFRFreducedOrder'], 
                                                                    meshNodeNumbers=np.array(nodeListJoint1T), #these are the meshNodeNumbers
                                                                    weightingFactors=noodeWeightsJoint1T))
        
        Marker14        = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber=TiltBoomFFRF['oFFRFreducedOrder'], 
                                                                    meshNodeNumbers=np.array(nodeListPist1T), #these are the meshNodeNumbers
                                                                    weightingFactors=noodeWeightsPist1T))      
        
        Marker15        = mbs.AddMarker(MarkerBodyRigid(bodyNumber=b4, localPosition=[0, 0, 0]))                        #With LIft Boom 
        Marker16        = mbs.AddMarker(MarkerBodyRigid(bodyNumber=b4, localPosition=[0.0425, 0.472227402, 0]))  
        Marker18        = mbs.AddMarker(MarkerBodyRigid(bodyNumber=b5, localPosition=[0, 0, 0]))                        #With LIft Boom 
        Marker19        = mbs.AddMarker(MarkerBodyRigid(bodyNumber=b5, localPosition=[0.475, 0, 0]))                   #With LIft Boom,-0.475 
        
        # Add joints
        
        #Add Revolute Joint btw Pillar and LiftBoom
        mbs.AddObject(GenericJoint(markerNumbers=[Marker4, Marker7],constrainedAxes=[1,1,1,1,1,0],
                         visualization=VObjectJointGeneric(axesRadius=radius1,axesLength=0.18)))   
        
        #Add Revolute Joint btw LiftBoom and TiltBoom
        mbs.AddObject(GenericJoint(markerNumbers=[Marker11, Marker13],constrainedAxes=[1,1,1,1,1,0],
                         visualization=VObjectJointGeneric(axesRadius=radius6,axesLength=0.16)))     
 
        #Add Revolute Joint btw LiftBoom and Bracket 1
        mbs.AddObject(GenericJoint(markerNumbers=[Marker10, Marker15],constrainedAxes=[1,1,1,1,1,0],
                         visualization=VObjectJointGeneric(axesRadius=radius6,axesLength=0.20)))   
        
        #Add Revolute Joint btw Bracket 1 and Bracket 2
        mbs.AddObject(GenericJoint(markerNumbers=[Marker16, Marker18],constrainedAxes=[1,1,1,1,1,0],
                         visualization=VObjectJointGeneric(axesRadius=0.28*W5,axesLength=2.0*W5)))  
        
        # Revolute joint between Bracket 2 and TiltBoom
        
        mbs.AddObject(GenericJoint(markerNumbers=[Marker19, Marker14],constrainedAxes=[1,1,0,0,0,0],
                             visualization=VObjectJointGeneric(axesRadius=radius7,axesLength=2.0*W5)))   
        
    else:    
        fileNameL       = 'LiftBoom/ANSYS/LiftBoom'                     #To load fine mesh data from AANSYS
        fileNameT       = 'TiltBoom/ANSYS/TiltBoom'                     #To load fine mesh data from AANSYS
        
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
                     TiltBoomModes(feT, FineMesh=False, nModes=2, FreeModes=True, HCB=False, ModeAnimation=True)
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
        
        # Boundary condition at Joint 2
        p8                  = [2.69+15e-3,2.38e-02,-7.4e-2]
        p7                  = [2.69+15e-3,2.38e-02, 7.4e-2]
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
        p12                 = [9.92e-15, -0.85e-3,-9.63e-2]
        p11                 = [9.92e-15, -0.85e-3, 9.63e-2]
        radius6             = 4.82e-002
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
        
        print("Compute Craig-Bampton modes... ")
        boundaryListL   = [nodeListJoint1, nodeListJoint2,  nodeListJoint3] 
        boundaryListT  = [nodeListJoint1T]
         
        start_time      = time.time()
        feL.ComputeHurtyCraigBamptonModes(boundaryNodesList=boundaryListL, nEigenModes=nModes, 
                                                    useSparseSolver=True,computationMode = HCBstaticModeSelection.RBE2)
        
        feT.ComputeHurtyCraigBamptonModes(boundaryNodesList=boundaryListT, nEigenModes=nModes, 
                                                    useSparseSolver=True,computationMode = HCBstaticModeSelection.RBE2) 
        
        print("Hurty-Craig Bampton modes... ")
        print("eigen freq. Lift Boom=", feL.GetEigenFrequenciesHz())
        print("eigen freq. Tilt Boom=", feT.GetEigenFrequenciesHz())
        print("HCB modes needed %.3f seconds" % (time.time() - start_time)) 
        
        LiftBoom            = ObjectFFRFreducedOrderInterface(feL)
        TiltBoom            = ObjectFFRFreducedOrderInterface(feT)
        
        LiftBoomFFRF        = LiftBoom.AddObjectFFRFreducedOrder(mbs, positionRef=LiftP, 
                                             initialVelocity=[0,0,0], 
                                              initialAngularVelocity=[0,0,0],
                                              rotationMatrixRef  = RotationMatrixZ(mt.radians(-1.227343749651500)),
                                              gravity=g,
                                              color=colLift,)
        
        TiltBoomFFRF        = TiltBoom.AddObjectFFRFreducedOrder(mbs, positionRef=LiftP+TiltL, #2.879420180699481+27e-3, -0.040690041435711005+8.3e-2, 0
                                              initialVelocity=[0,0,0], 
                                              initialAngularVelocity=[0,0,0],
                                              rotationMatrixRef  = RotationMatrixZ(mt.radians(-1.582031486865900)),
                                              gravity=g,
                                              color=colLift,)
        
        # 4th Body: Bracket 1
        L4              = 0.557227    # Length in x-direction
        H4              = 0.1425      # Height in y-direction
        W4              = 0.15        # Width in z-direction
        Bracket1L       = LiftP + Bracket1
        pMid4           = np.array([0.004000, 0.257068, 0])  # center of mass, body0,0.004000,-0.257068
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
                                         rotationMatrix=RotationMatrixZ(mt.radians(58.803311277817050)),
                                         gravity=g,graphicsDataList=[graphicsCOM4, graphicsBody4])
        
        # 5th Body: Bracket 2
        L5              = 0.569009       # Length in x-direction
        H5              = 0.078827       # Height in y-direction
        W5              = 0.15           # Width in z-direction
        pMid5           = np.array([0.212792, 0, 0])  # center of mass, body0
        Bracket1B       = LiftP + Bracket1+ Bracket2  #0.285710892999728-1.8*0.0425, -0.356968041652145+0.0525
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
                                          rotationMatrix=RotationMatrixZ(mt.radians(-3.151443859554505) ),   #-5
                                          gravity=g,graphicsDataList=[graphicsCOM5, graphicsBody5])
        
        Marker7         = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber=LiftBoomFFRF['oFFRFreducedOrder'],
                                             meshNodeNumbers=np.array(nodeListJoint1), #these are the meshNodeNumbers
                                              weightingFactors=noodeWeightsJoint1))
        
        Marker10        = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber=LiftBoomFFRF['oFFRFreducedOrder'],
                                              meshNodeNumbers=np.array(nodeListJoint2), #these are the meshNodeNumbers
                                              weightingFactors=noodeWeightsJoint2))      
        
        Marker11        = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber=LiftBoomFFRF['oFFRFreducedOrder'],
                                            meshNodeNumbers=np.array(nodeListJoint3), #these are the meshNodeNumbers
                                            weightingFactors=noodeWeightsJoint3))
        
        Marker13        = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber=TiltBoomFFRF['oFFRFreducedOrder'], 
                                                                    meshNodeNumbers=np.array(nodeListJoint1T), #these are the meshNodeNumbers
                                                                    weightingFactors=noodeWeightsJoint1T))
        
        Marker14        = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber=TiltBoomFFRF['oFFRFreducedOrder'], 
                                                                    meshNodeNumbers=np.array(nodeListPist1T), #these are the meshNodeNumbers
                                                                    weightingFactors=noodeWeightsPist1T))      
        
        Marker15        = mbs.AddMarker(MarkerBodyRigid(bodyNumber=b4, localPosition=[0, 0, 0]))                        #With LIft Boom 
        Marker16        = mbs.AddMarker(MarkerBodyRigid(bodyNumber=b4, localPosition=[0.0425, 0.472227402, 0]))  
        Marker18        = mbs.AddMarker(MarkerBodyRigid(bodyNumber=b5, localPosition=[0, 0, 0]))                        #With LIft Boom 
        Marker19        = mbs.AddMarker(MarkerBodyRigid(bodyNumber=b5, localPosition=[0.475, 0, 0]))                   #With LIft Boom,-0.475 
        
        # Add joints
        
        #Add Revolute Joint btw Pillar and LiftBoom
        mbs.AddObject(GenericJoint(markerNumbers=[Marker4, Marker7],constrainedAxes=[1,1,1,1,1,0],
                         visualization=VObjectJointGeneric(axesRadius=radius1,axesLength=0.18)))   
        
        #Add Revolute Joint btw LiftBoom and TiltBoom
        mbs.AddObject(GenericJoint(markerNumbers=[Marker11, Marker13],constrainedAxes=[1,1,1,1,1,0],
                         visualization=VObjectJointGeneric(axesRadius=radius6,axesLength=0.16)))     
 
        #Add Revolute Joint btw LiftBoom and Bracket 1
        mbs.AddObject(GenericJoint(markerNumbers=[Marker10, Marker15],constrainedAxes=[1,1,1,1,1,0],
                         visualization=VObjectJointGeneric(axesRadius=radius6,axesLength=0.20)))   
        
        #Add Revolute Joint btw Bracket 1 and Bracket 2
        mbs.AddObject(GenericJoint(markerNumbers=[Marker16, Marker18],constrainedAxes=[1,1,1,1,1,0],
                         visualization=VObjectJointGeneric(axesRadius=0.28*W5,axesLength=2.0*W5)))  
        
        # Revolute joint between Bracket 2 and TiltBoom
        
        mbs.AddObject(GenericJoint(markerNumbers=[Marker19, Marker14],constrainedAxes=[1,1,0,0,0,0],
                             visualization=VObjectJointGeneric(axesRadius=radius7,axesLength=2.0*W5)))    
        
    if  Hydraulics:
        
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
        
        Marker8             = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber=LiftBoomFFRF['oFFRFreducedOrder'],
                                              meshNodeNumbers=np.array(nodeListPist1), #these are the meshNodeNumbers
                                              weightingFactors=noodeWeightsPist1))           #With Cylinder 1
        
        Marker9             = mbs.AddMarker(MarkerSuperElementRigid(bodyNumber=LiftBoomFFRF['oFFRFreducedOrder'],
                                              meshNodeNumbers=np.array(nodeListCyl2), #these are the meshNodeNumbers
                                              weightingFactors=noodeWeightsCyl2)) 
        
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
        pS              = 100e5
        pT              = 1e5                           # Tank pressure
        Qn1             = (18/60000)/(sqrt(20e5)*9.9)                      # Nominal flow rate of valve at 18 l/min under
        Qn2             = (22/60000)/(sqrt(20e5)*9.9)                      # Nominal flow rate of valve at 18 l/min under
        
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
                                    initialCoordinates=[6.204205450539635e+06,
                                                        3.100522576711962e+06],  # initialize with 20 bar
                                                        numberOfODE1Coordinates=2))

        nODE2               = mbs.AddNode(NodeGenericODE1(referenceCoordinates=[0, 0],
                                    initialCoordinates=[4.806409161717769e+06,
                                                        6.754810665699534e+06],  # initialize with 20 bar
                                                        numberOfODE1Coordinates=2))
        def UFfrictionSpringDamper1(mbs, t, itemIndex, u, v, k, d, f0):

            return   6*(Fc*tanh(4*(abs(v    )/vs))+(Fs-Fc)*((abs(v    )/vs)/((1/4)*(abs(v    )/vs)**2+3/4)**2))*np.sign(v )+sig2*v    *tanh(4)

        
        def UFfrictionSpringDamper2(mbs, t, itemIndex, u, v, k, d, f0):

            return   4*(Fc*tanh(4*(abs(v    )/vs))+(Fs-Fc)*((abs(v    )/vs)/((1/4)*(abs(v    )/vs)**2+3/4)**2))*np.sign(v )+sig2*v    *tanh(4)
        
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
    
    mbs.Assemble()
    simulationSettings = exu.SimulationSettings() 
    
    if not StaticCase:
        simulationSettings.timeIntegration.numberOfSteps            = nSteps
        simulationSettings.timeIntegration.endTime                  = tEnd
        simulationSettings.timeIntegration.verboseModeFile          = 0
        simulationSettings.timeIntegration.verboseMode              = 1
        simulationSettings.solutionSettings.recordImagesInterval    = -1  
        simulationSettings.solutionSettings.solutionWritePeriod     = 0.005  # store every 5 ms
        simulationSettings.timeIntegration.newton.useModifiedNewton = True
        simulationSettings.linearSolverType                         = exu.LinearSolverType.EigenSparse
        simulationSettings.linearSolverSettings.ignoreSingularJacobian=True
        simulationSettings.timeIntegration.stepInformation         += 8
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
   
    return