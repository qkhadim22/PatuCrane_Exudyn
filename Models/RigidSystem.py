#%%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#This files defines the simulation of hydraulically actuated rigid structure.

#%%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
import sys
sys.exudynFast = True

import exudyn as exu
import numpy as np
import math as mt
import scipy.io
import matplotlib.pyplot as plt

from exudyn.utilities import * #includes graphics and rigid body utilities
from math import sin, cos, sqrt, pi, tanh


SC              = exu.SystemContainer()
mbs             = SC.AddSystem()
#%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
g               = [0,-9.81,0] # Gravity
tEnd            = 20 #simulation time
nRigidBodyNodes = 4
h               = 1e-3 #step size
nSteps          = int(tEnd/h)+2

ListBodies      = ['Pillar', 'LiftBoom', 'TiltBoom+ExtensionBoom', 
                                            'Bracket 1', 'Bracket 2']  

# Loading Graphics of bodies
fileName1       = 'Graphics_Exudyn/Pillar.stl'
fileName2       = 'Graphics_Exudyn/LiftBoom.stl'
fileName3       = 'Graphics_Exudyn/TiltBoom+ExtensionBoom.stl'
fileName4       = 'Graphics_Exudyn/Bracket1.stl'
fileName5       = 'Graphics_Exudyn/Bracket2.stl'
fileName6       = 'Graphics_Exudyn/TiltBoom+ExtensionBoom_MC.stl'


#Ground body
oGround         = mbs.AddObject(ObjectGround())
Marker1         = mbs.AddMarker(MarkerBodyRigid(bodyNumber=oGround, 
                                                localPosition=[0, 0, 0]))
# First Body: Pillar
L1              = 0.365                                 # Length in x-direction
H1              = 1.4769                                # Height in y-direction
W1              = 0.25                                  # Width in z-direction
pMid1           = np.array([-0.017403, 0.577291, 0])    # Center of mass
PillarP         = np.array([0, 0, 0])
iCube1          = RigidBodyInertia(mass=93.26, com=pMid1,
                                   inertiaTensor=np.array([[16.358844,-1.27808, 1.7e-5],
                                                           [-1.27808, 0.612552, -5.9e-5],
                                                           [1.7e-5,  -5.9e-5  , 16.534255]]),
                                                               inertiaTensorAtCOM=True)

graphicsBody1   = GraphicsDataFromSTLfile(fileName1, color4black,verbose=False, 
                                                  invertNormals=True,invertTriangles=True)
graphicsBody1   = AddEdgesAndSmoothenNormals(graphicsBody1, edgeAngle=0.25*pi,addEdges=True, 
                                                     smoothNormals=True)
graphicsCOM1    = GraphicsDataBasis(origin=iCube1.com, length=2*W1)

# Definintion of pillar as body in Exudyn and node n1
[n1, b1]        = AddRigidBody(mainSys=mbs,
                             inertia=iCube1,  # includes COM
                             nodeType=exu.NodeType.RotationEulerParameters,
                             position=PillarP,
                             rotationMatrix=np.diag([1, 1, 1]),
                             gravity=g,
                             graphicsDataList=[graphicsCOM1, graphicsBody1])

Marker3         = mbs.AddMarker(MarkerBodyRigid(bodyNumber=b1, localPosition=[0, 0, 0]))                     #With Ground
Marker4         = mbs.AddMarker(MarkerBodyRigid(bodyNumber=b1, localPosition=[-0.09, 1.4261, 0]))            #Lift Boom
Marker5         = mbs.AddMarker(MarkerBodyRigid(bodyNumber=b1, localPosition=[0.17, 0.386113249, 0]))        # Cylinder 1 position
       
#Fixed joint btw Ground and Pillar
mbs.AddObject(GenericJoint(markerNumbers=[Marker1, Marker3],constrainedAxes=[1, 1, 1,1,1,1],
                                visualization=VObjectJointGeneric(axesRadius=0.2*W1,axesLength=1.4*W1)))


def RigidMultibodyHydraulics(RedundantCoordinates, Hydraulics, useFriction, Plotting):
    
    global Fc, Fs, sig2, vs, Av0, Av1, Av2, Av3

    
    if RedundantCoordinates:
        #%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                        # BODIES IN PATU CRANE
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Second Body: LiftBoom
        L2              = 3.01055           # Length in x-direction
        H2              = 0.45574           # Height in y-direction
        W2              = 0.263342          # Width in z-direction
        LiftP           = np.array([-0.09, 1.4261, 0])
        pMid2           = np.array([1.229248, 0.055596, 0])  # center of mass
        iCube2          = RigidBodyInertia(mass=143.66, com=pMid2,
                                        inertiaTensor=np.array([[1.055433, 1.442440,  -0.000003],
                                                                [ 1.442440,  66.577004, 0],
                                                                [ -0.000003,              0  ,  67.053707]]),
                                        inertiaTensorAtCOM=True)
        
        graphicsBody2   = GraphicsDataFromSTLfile(fileName2, color4blue,
                                        verbose=False, invertNormals=True,
                                        invertTriangles=True)
        graphicsBody2   = AddEdgesAndSmoothenNormals(graphicsBody2, edgeAngle=0.25*pi,
                                            addEdges=True, smoothNormals=True)
        graphicsCOM2    = GraphicsDataBasis(origin=iCube2.com, length=2*W2)
        [n2, b2]        = AddRigidBody(mainSys=mbs,
                            inertia=iCube2,  # includes COM
                            nodeType=exu.NodeType.RotationEulerParameters,
                            position=LiftP,  # pMid2
                            rotationMatrix=RotationMatrixZ(mt.radians(-1.227343749651500)),
                            gravity=g,
                            graphicsDataList=[graphicsCOM2, graphicsBody2])
        
        Marker7         = mbs.AddMarker(MarkerBodyRigid(bodyNumber=b2, localPosition=[0, 0, 0]))                       #With Pillar    
        Marker8         = mbs.AddMarker(MarkerBodyRigid(bodyNumber=b2, localPosition=[0.3025, -0.105, 0]))             #With Cylinder 1
        Marker9         = mbs.AddMarker(MarkerBodyRigid(bodyNumber=b2, localPosition=[1.263, 0.206702194, 0]))         #With Cylinder 2
        Marker10        = mbs.AddMarker(MarkerBodyRigid(bodyNumber=b2, localPosition=[2.69, -0.006592554, 0]))         #With Bracket 1  
        Marker11        = mbs.AddMarker(MarkerBodyRigid(bodyNumber=b2, localPosition=[2.881080943, 0.021592554, 0]))   #With Tilt Boom

        # Third Body: TiltBoom+ExtensionBoom   
        L3              = 2.580         # Length in x-direction
        H3              = 0.419         # Height in y-direction
        W3              = 0.220         # Width in z-direction
        TiltL           = LiftP + np.array([2.879420180699481, -0.040690041435711005, 0])
        pMid3           = np.array([ 0.659935,  0.251085, 0])  # center of mass
        iCube3          = RigidBodyInertia(mass=141.942729+ 15.928340, com=pMid3,
                                            inertiaTensor=np.array([[1.055433, 1.442440,  -0.000003],
                                                                    [1.442440,  66.577004,    0],
                                                                    [ -0.000003, 0,        67.053707]]),
                                                                     inertiaTensorAtCOM=True)
        graphicsBody3   = GraphicsDataFromSTLfile(fileName3, color4blue,verbose=False, invertNormals=True,invertTriangles=True)
        graphicsBody3   = AddEdgesAndSmoothenNormals(graphicsBody3, edgeAngle=0.25*pi,addEdges=True, smoothNormals=True)
        graphicsCOM3    = GraphicsDataBasis(origin=iCube3.com, length=2*W3)
        [n3, b3]        = AddRigidBody(mainSys=mbs,
                               inertia=iCube3,  # includes COM
                               nodeType=exu.NodeType.RotationEulerParameters,
                               position=TiltL,  # pMid2
                               rotationMatrix=RotationMatrixZ(mt.radians(-1.582031486865900)),
                               gravity=g,
                               graphicsDataList=[graphicsCOM3, graphicsBody3])
        Marker13        = mbs.AddMarker(MarkerBodyRigid(bodyNumber=b3, localPosition=[0, 0, 0]))                        #With LIft Boom 
        Marker14        = mbs.AddMarker(MarkerBodyRigid(bodyNumber=b3, localPosition=[-0.095, 0.24043237, 0]))                        #With LIft Boom 
       
        # 4th Body: Bracket 1
        L4              = 0.557227    # Length in x-direction
        H4              = 0.1425      # Height in y-direction
        W4              = 0.15        # Width in z-direction
        Bracket1L       = LiftP + np.array([2.689524056550459, -0.066027551426741, 0])
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
        Bracket1B       = LiftP + np.array([2.329958953740425, 0.22034355756469298, 0])  #0.285710892999728-1.8*0.0425, -0.356968041652145+0.0525
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
 
        #%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                                # JOINT DEFINITION# 
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
           
        #Revolute joint between Pillar and Lift boom
        mbs.AddObject(GenericJoint(markerNumbers=[Marker4, Marker7],constrainedAxes=[1,1,1,1,1,0],
                                visualization=VObjectJointGeneric(axesRadius=0.18*W2,axesLength=1.1*W2)))
        
        # Revolute joint between Lift Boom and Tilt Boom
        mbs.AddObject(GenericJoint(markerNumbers=[Marker11, Marker13],constrainedAxes=[1,1,1,1,1,0],
                             visualization=VObjectJointGeneric(axesRadius=0.22*W3,axesLength=0.72*W3)))
        
        # Revolute joint between LiftBoom and Bracket 1
        mbs.AddObject(GenericJoint(markerNumbers=[Marker10, Marker15],constrainedAxes=[1,1,1,1,1,0],
                             visualization=VObjectJointGeneric(axesRadius=0.32*W4,axesLength=0.96*W4))) 
        
        # Revolute joint between Bracket 1 and Bracket 2
        mbs.AddObject(GenericJoint(markerNumbers=[Marker16, Marker18],constrainedAxes=[1,1,1,1,1,0],
                             visualization=VObjectJointGeneric(axesRadius=0.23*W5,axesLength=1.0*W5)))
        
        # Revolute joint between Bracket 2 and TiltBoom
        mbs.AddObject(GenericJoint(markerNumbers=[Marker19, Marker14],constrainedAxes=[1,1,0,0,0,0],
                             visualization=VObjectJointGeneric(axesRadius=0.23*W5,axesLength=1.0*W5)))        
        
        # Add sensors
        Angle1          = mbs.AddSensor(SensorBody(bodyNumber=b2, localPosition=[0,0,0.0],
                                        fileName='ExData/LiftAngle1.txt', storeInternal=True,outputVariableType = exu.OutputVariableType.Rotation))
        Angle2          = mbs.AddSensor(SensorBody(bodyNumber=b3, localPosition=[0,0,0.0],
                                        fileName='ExData/TiltAngle2.txt', storeInternal=True,outputVariableType = exu.OutputVariableType.Rotation))
        Angle1_t        = mbs.AddSensor(SensorBody(bodyNumber=b2, localPosition=[0,0,0.0],
                                        fileName='ExData/LiftAngle1_t.txt', storeInternal=True,outputVariableType = exu.OutputVariableType.AngularVelocityLocal))
        Angle2_t        = mbs.AddSensor(SensorBody(bodyNumber=b3, localPosition=[0,0,0.0],
                                        fileName='ExData/TiltAngle2_t.txt', storeInternal=True,outputVariableType = exu.OutputVariableType.AngularVelocityLocal))
    
    else:
        
        # Second Body: LiftBoom
        W2              = 0.263342          # Width in z-direction
        pMid2           = np.array([1.229248, 0.055596, 0])  # center of mass
        iCube2          = RigidBodyInertia(mass=143.66, com=pMid2,
                                        inertiaTensor=np.array([[1.055433, 1.442440,  -0.000003],
                                                                [ 1.442440,  66.577004, 0],
                                                                [ -0.000003,              0  ,  67.053707]]),
                                        inertiaTensorAtCOM=True)
        
        graphicsBody2   = GraphicsDataFromSTLfile(fileName2, color4blue,
                                        verbose=False, invertNormals=True,
                                        invertTriangles=True)
        graphicsBody2   = AddEdgesAndSmoothenNormals(graphicsBody2, edgeAngle=0.25*pi,
                                            addEdges=True, smoothNormals=True)
        graphicsCOM2    = GraphicsDataBasis(origin=iCube2.com, length=2*W2)  
        
       # Third Body: TiltBoom+ExtensionBoom   
        W3              = 0.220         # Width in z-direction
        pMid3           = np.array([ 0.754935,  0.010653, 0])  # center of mass
        iCube3          = RigidBodyInertia(mass=141.942729+ 15.928340, com=pMid3,
                                            inertiaTensor=np.array([[1.055433, 1.442440,  -0.000003],
                                                                    [1.442440,  66.577004,    0],
                                                                    [ -0.000003, 0,        67.053707]]),
                                                                     inertiaTensorAtCOM=True)
        graphicsBody3   = GraphicsDataFromSTLfile(fileName6, color4blue,verbose=False, invertNormals=True,invertTriangles=True)
        graphicsBody3   = AddEdgesAndSmoothenNormals(graphicsBody3, edgeAngle=0.25*pi,addEdges=True, smoothNormals=True)
        graphicsCOM3    = GraphicsDataBasis(origin=iCube3.com, length=2*W3)
   
        # 4th Body: Bracket 1
        W4              = 0.15        # Width in z-direction
        pMid4           = np.array([-0.257068, 0.004000, 0])  # center of mass, body0
        iCube4          = RigidBodyInertia(mass=11.524039, com=pMid4,
                                            inertiaTensor=np.array([[0.333066, 0.017355, 0],
                                                                    [0.017355, 0.081849, 0],
                                                                    [0,              0, 0.268644]]),
                                                                     inertiaTensorAtCOM=True)

        graphicsBody4   = GraphicsDataFromSTLfile(fileName4, color4blue,verbose=False, invertNormals=True,invertTriangles=True)
        graphicsBody4   = AddEdgesAndSmoothenNormals(graphicsBody4, edgeAngle=0.25*pi,addEdges=True, smoothNormals=True)
        graphicsCOM4    = GraphicsDataBasis(origin=iCube4.com, length=2*W4) 
            
       # 5th Body: Bracket 2
        W5              = 0.15           # Width in z-direction
        pMid5           = np.array([-0.212792, 0, 0])  # center of mass, body0
        iCube5          = RigidBodyInertia(mass=7.900191, com=pMid5,
                                             inertiaTensor=np.array([[0.052095, 0, 0],
                                                                     [0,  0.260808, 0],
                                                                     [0,              0,  0.216772]]),
                                                                     inertiaTensorAtCOM=True)
        
        graphicsBody5   = GraphicsDataFromSTLfile(fileName5, color4blue,verbose=False, invertNormals=True,invertTriangles=True)
        graphicsBody5   = AddEdgesAndSmoothenNormals(graphicsBody5, edgeAngle=0.25*pi,addEdges=True, smoothNormals=True)
        graphicsCOM5    = GraphicsDataBasis(origin=iCube5.com, length=2*W5)
                        
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        nGeneric        = mbs.AddNode(NodeGenericODE2(referenceCoordinates=[0.]*nRigidBodyNodes,initialCoordinates=[0.]*nRigidBodyNodes,
                                       initialCoordinates_t=[0.]*nRigidBodyNodes,numberOfODE2Coordinates=nRigidBodyNodes))
        inertiaList     = [iCube2, iCube4, iCube5, iCube3]
        
        graphicsList    = [[graphicsCOM2, graphicsBody2], [graphicsCOM4, graphicsBody4],
                           [graphicsCOM5, graphicsBody5], [graphicsCOM3, graphicsBody3]] 
        
        #create KinematicTree
        JointPos        = [[-0.09, 1.4261, 0], [2.689524056550459, -6e-3, 0 ],
                           [-0.459565102810030, 0.040343557564693, 0],
                           [-0.474565102810030, 0, 0] ]

        
        jointTypes      = [exu.JointType.RevoluteZ, exu.JointType.RevoluteZ,
                           exu.JointType.RevoluteZ,exu.JointType.RevoluteZ] 
        BodiesMasses    = []
        BodiesCOMs      = exu.Vector3DList() 
        BodiesInertias  = exu.Matrix3DList()
        jointTrans      = exu.Matrix3DList()
        jointOffsets    = exu.Vector3DList()
        A               = np.eye(3)
        
        for i in range(nRigidBodyNodes):    
            inertia          = inertiaList[i]
            BodiesMasses    += [inertia.Mass()]
            BodiesCOMs.Append(inertia.COM())
            BodiesInertias.Append(inertia.InertiaCOM())
            
            if i == 0:
                A = RotationMatrixZ(mt.radians(-1.227343749651500))
            if i == 1:
                A = RotationMatrixZ(mt.radians(-33.195078970135967))
            if i == 2:
                A = RotationMatrixZ(mt.radians(-147.7373))  
            if i == 3:
                A = RotationMatrixZ(mt.radians(-181.082031486865900+1.9))                             
                    
            jointTrans.Append(A)
            jointOffsets.Append(JointPos[i])
        
        
        oKT             = mbs.AddObject(ObjectKinematicTree(nodeNumber=nGeneric, jointTypes=jointTypes, linkParents=np.arange(nRigidBodyNodes)-1,
                                  jointTransformations=jointTrans, jointOffsets=jointOffsets, linkInertiasCOM=BodiesInertias,
                                  linkCOMs=BodiesCOMs, linkMasses=BodiesMasses, gravity=g,
                                  visualization=VObjectKinematicTree(graphicsDataList = graphicsList)))

        # Adding Markers in the bodies
        Marker8         = mbs.AddMarker(MarkerKinematicTreeRigid(objectNumber=oKT, linkNumber=0,localPosition=[0.3025, -0.105, 0]))
        Marker9         = mbs.AddMarker(MarkerKinematicTreeRigid(objectNumber=oKT, linkNumber=0,localPosition=[1.263, 0.206702194, 0]))
        Marker14        = mbs.AddMarker(MarkerKinematicTreeRigid(objectNumber=oKT, linkNumber=0,localPosition=[2.879420180699481, -0.040690041435711005, 0]))        
        Marker18        = mbs.AddMarker(MarkerKinematicTreeRigid(objectNumber=oKT, linkNumber=2,localPosition=[0, 0, 0]))                        #With LIft Boom 
        Marker19        = mbs.AddMarker(MarkerKinematicTreeRigid(objectNumber=oKT, linkNumber=3,localPosition=[0.095, -0.24043237, 0]))                   #With LIft Boom,-0.475                                                                    
 
        # Revolute joint between Bracket 2 and TiltBoom
        mbs.AddObject(GenericJoint(markerNumbers=[Marker19, Marker14],constrainedAxes=[1,1,0,0,0,0],
                             visualization=VObjectJointGeneric(axesRadius=0.23*W5,axesLength=1.0*W5)))        
            
    if Hydraulics:
        
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
                                                     fileName='ExData/sForce1.txt', outputVariableType=exu.OutputVariableType.Force))
        sForce2         = mbs.AddSensor(SensorObject(objectNumber=oHA2, storeInternal=True, 
                                                     fileName='ExData/sForce2.txt', outputVariableType=exu.OutputVariableType.Force))
        sDistance1      = mbs.AddSensor(SensorObject(objectNumber=oHA1, storeInternal=True, 
                                                     fileName='ExData/sDistance1.txt', outputVariableType=exu.OutputVariableType.Distance))
        sDistance2      = mbs.AddSensor(SensorObject(objectNumber=oHA2, storeInternal=True, 
                                                     fileName='ExData/sDistance2.txt', outputVariableType=exu.OutputVariableType.Distance))
        sVelocity1      = mbs.AddSensor(SensorObject(objectNumber=oHA1, storeInternal=True,
                                                     fileName='ExData/sVelocity1.txt', outputVariableType=exu.OutputVariableType.VelocityLocal))
        sVelocity2      = mbs.AddSensor(SensorObject(objectNumber=oHA2, storeInternal=True,
                                                     fileName='ExData/sVelocity2.txt', outputVariableType=exu.OutputVariableType.VelocityLocal))  
        sPressures1     = mbs.AddSensor(SensorNode(nodeNumber=nODE1, storeInternal=True, 
                                                   fileName='ExData/sPressures1.txt',outputVariableType=exu.OutputVariableType.Coordinates))
        sPressures2     = mbs.AddSensor(SensorNode(nodeNumber=nODE2, storeInternal=True, 
                                                   fileName='ExData/sPressures2.txt',outputVariableType=exu.OutputVariableType.Coordinates))
        if RedundantCoordinates:
            Angle1h         = mbs.AddSensor(SensorBody(bodyNumber=b2, localPosition=[0,0,0.0],
                                        fileName='ExData/Angle1.txt', storeInternal=True,outputVariableType = exu.OutputVariableType.Rotation))
            Angle2h         = mbs.AddSensor(SensorBody(bodyNumber=b3, localPosition=[0,0,0.0],
                                        fileName='ExData/Angle2.txt', storeInternal=True,outputVariableType = exu.OutputVariableType.Rotation))
            Angle1_th       = mbs.AddSensor(SensorBody(bodyNumber=b2, localPosition=[0,0,0.0],
                                        fileName='ExData/Angle1_t.txt', storeInternal=True,outputVariableType = exu.OutputVariableType.AngularVelocityLocal))
            Angle2_th       = mbs.AddSensor(SensorBody(bodyNumber=b3, localPosition=[0,0,0.0],
                                        fileName='ExData/Angle2_t.txt', storeInternal=True,outputVariableType = exu.OutputVariableType.AngularVelocityLocal))
        else:
            Angle1h         = mbs.AddSensor(SensorKinematicTree(objectNumber=oKT, linkNumber=0, localPosition=[0,0,0.0],
                                        storeInternal=True,outputVariableType = exu.OutputVariableType.Rotation))
            Angle2h         = mbs.AddSensor(SensorKinematicTree(objectNumber=oKT, linkNumber=3, localPosition=[0.095, -0.24043237,0],
                                          storeInternal=True,outputVariableType = exu.OutputVariableType.Rotation))
            Angle1_th       = mbs.AddSensor(SensorKinematicTree(objectNumber=oKT, linkNumber=0, localPosition=[0,0,0.0],
                                         storeInternal=True,outputVariableType = exu.OutputVariableType.AngularVelocityLocal))
            Angle2_th       = mbs.AddSensor(SensorKinematicTree(objectNumber=oKT, linkNumber=3, localPosition=[0.095, -0.24043237,0],
                                           storeInternal=True,outputVariableType = exu.OutputVariableType.AngularVelocityLocal))
    #%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    mbs.Assemble()
    
    #Simulation settings
    simulationSettings                                                  = exu.SimulationSettings()
    simulationSettings.timeIntegration.numberOfSteps                    = nSteps
    simulationSettings.timeIntegration.endTime                          = tEnd
    simulationSettings.timeIntegration.verboseMode                      = 1
    # simulationSettings.timeIntegration.simulateInRealtime     = True
    simulationSettings.solutionSettings.sensorsWritePeriod              = h
    simulationSettings.linearSolverType                                 = exu.LinearSolverType.EigenSparse
    simulationSettings.timeIntegration.newton.useModifiedNewton         = True
    simulationSettings.timeIntegration.stepInformation                  += 8
    # simulationSettings.displayComputationTime                 = True
    simulationSettings.displayStatistics                                = True
    simulationSettings.timeIntegration.newton.relativeTolerance         = 1e-6
    simulationSettings.timeIntegration.generalizedAlpha.spectralRadius  = 0.7
    simulationSettings.timeIntegration.relativeTolerance                = 1e-6
    
    #Visualization
    SC.visualizationSettings.window.renderWindowSize                    = [1600, 1200]
    SC.visualizationSettings.openGL.multiSampling                       = 4
    SC.visualizationSettings.openGL.lineWidth                           = 3
    SC.visualizationSettings.general.autoFitScene                       = False
    SC.visualizationSettings.nodes.drawNodesAsPoint                     = False
    SC.visualizationSettings.nodes.showBasis                            = True
    
    exu.SolveDynamic(mbs, simulationSettings=simulationSettings,solverType=exu.DynamicSolverType.TrapezoidalIndex2)
    
    #from exudyn.interactive import SolutionViewer
    #SolutionViewer(mbs)
    
    if Plotting:
        
        if not Hydraulics:
            Angle1   = mbs.GetSensorStoredData(Angle1)
            Angle2   = mbs.GetSensorStoredData(Angle2)
            Angle1_t = mbs.GetSensorStoredData(Angle1_t)
            Angle2_t = mbs.GetSensorStoredData(Angle2_t)
            
        
        # Angle1 of Liftboom
            fig, ax = plt.subplots()
            plt.plot(Angle1[:, 0], np.rad2deg(Angle1[:, 3]), label='EXUDYN', color='blue')
            plt.xlim(0, tEnd)
            plt.grid(True)
            # Add a title and axis labels
            plt.ylim(-360, 360)
            #ax.set_ylabel('Angle ${\Theta}_{1_z}$, deg')
            plt.legend()
            plt.savefig('ExData/LiftTheta1.png', dpi=300)
            plt.show()
 
            # Angle of Tiltboom
            fig, ax = plt.subplots()
            plt.plot(Angle2[:, 0], np.rad2deg(Angle2[:, 3]), label='EXUDYN', color='blue')
            plt.xlim(0, tEnd)
            plt.grid(True)
            # Add a title and axis labels
            plt.ylim(-360, 360)
            #ax.set_ylabel('Angle ${\Theta}_{2_z}$, deg')
            plt.legend()
            plt.savefig('ExData/TiltTheta2.png', dpi=300)
            plt.show()       

            # Angular velocity of Liftboom
            fig, ax = plt.subplots()
            plt.plot(Angle1_t[:, 0], np.rad2deg(Angle1_t[:, 3]), label='EXUDYN', color='blue')
            plt.xlim(0, tEnd)
            plt.grid(True)
            # Add a title and axis labels
            plt.ylim(-360, 360)
            #ax.set_ylabel('Angle ${\Dot{\Theta}}_{1_z}$, deg')
            plt.legend()
            plt.savefig('ExData/LiftdTheta1.png', dpi=300)
            plt.show()
 
            # Angular velocity of Tiltboom
            fig, ax = plt.subplots()
            plt.plot(Angle2_t[:, 0], np.rad2deg(Angle2_t[:, 3]), label='EXUDYN', color='blue')
            plt.xlim(0, tEnd)
            plt.grid(True)
            # Add a title and axis labels
            plt.ylim(-360, 360)
            #ax.set_ylabel('Angle ${\Dot{\Theta}}_{2_z}$, deg')
            plt.legend()
            plt.savefig('ExData/TiltdTheta2.png', dpi=300)
            plt.show()  
        
        else:
            
            sDistance1   = mbs.GetSensorStoredData(sDistance1)
            sDistance2   = mbs.GetSensorStoredData(sDistance2)
            sPressures1  = mbs.GetSensorStoredData(sPressures1)
            sPressures2  = mbs.GetSensorStoredData(sPressures2)
            Angle1       = mbs.GetSensorStoredData(Angle1h)
            Angle2       = mbs.GetSensorStoredData(Angle2h)
            Angle1_t     = mbs.GetSensorStoredData(Angle1_th)
            Angle2_t     = mbs.GetSensorStoredData(Angle2_th)
        
            # Angle1 of Liftboom
            fig, ax = plt.subplots()
            plt.plot(Angle1[:, 0], np.rad2deg(Angle1[:, 3]), label='EXUDYN', color='blue')
            plt.plot(tspan, np.rad2deg(a1M), label='MATLAB', color='red')
            plt.xlim(0, tEnd)
            plt.grid(True)
            # Add a title and axis labels
            plt.ylim(-25, 25)
            #ax.set_ylabel('Angle ${\Theta}_{1_z}$, deg')
            plt.legend()
            plt.savefig('ExData/Theta1.png', dpi=300)
            plt.show()
 
            # Angle of Tiltboom
            fig, ax = plt.subplots()
            plt.plot(Angle2[:, 0], np.rad2deg(Angle2[:, 3]), label='EXUDYN', color='blue')
            plt.plot(tspan, np.rad2deg(a1M)+np.rad2deg(a2M), label='MATLAB', color='red')
            plt.xlim(0, tEnd)
            plt.grid(True)
            # Add a title and axis labels
            plt.ylim(-60, 30)
            #ax.set_ylabel('Angle ${\Theta}_{2_z}$, deg')
            plt.legend()
            plt.savefig('ExData/Theta2.png', dpi=300)
            plt.show()       

            # Angular velocity of Liftboom
            fig, ax = plt.subplots()
            plt.plot(Angle1_t[:, 0], np.rad2deg(Angle1_t[:, 3]), label='EXUDYN', color='blue')
            plt.plot(tspan, np.rad2deg(da1M), label='MATLAB', color='red')
            plt.xlim(0, tEnd)
            plt.grid(True)
            # Add a title and axis labels
            plt.ylim(-12, 12)
            #ax.set_ylabel('Angle ${\Dot{\Theta}}_{1_z}$, deg')
            plt.legend()
            plt.savefig('ExData/dTheta1.png', dpi=300)
            plt.show()
 
            # Angular velocity of Tiltboom
            fig, ax = plt.subplots()
            plt.plot(Angle2_t[:, 0], np.rad2deg(Angle2_t[:, 3]), label='EXUDYN', color='blue')
            plt.plot(tspan, np.rad2deg(da1M)+np.rad2deg(da2M), label='MATLAB', color='red')
            plt.xlim(0, tEnd)
            plt.grid(True)
            # Add a title and axis labels
            plt.ylim(-35, 35)
            #ax.set_ylabel('Angle ${\Dot{\Theta}}_{2_z}$, deg')
            plt.legend()
            plt.savefig('ExData/dTheta2.png', dpi=300)
            
            fig, ax = plt.subplots()
            plt.plot(sDistance1[:, 0], 1000*sDistance1[:, 1], label='EXUDYN', color='blue')
            plt.plot(tspan, 1000*sM1, label='MATLAB', color='red')
            plt.xlim(0, tEnd)
            plt.ylim(850, 1050)
            plt.grid(True)
            # Add a title and axis labels
            ax.set_xlabel('Time, s')
            ax.set_ylabel('Actuator 1, mm')
            plt.legend()
            plt.savefig('ExData/ActuatorPosition1.png', dpi=300)
            plt.show() 

            fig, ax = plt.subplots()
            plt.plot(sDistance1[:, 0], 1000*sDistance2[:, 1], label='EXUDYN', color='blue')
            plt.plot(tspan, 1000*sM2, label='MATLAB', color='red')
            plt.xlim(0, tEnd)
            plt.ylim(1000, 1400)
            plt.grid(True)
            # Add a title and axis labels
            ax.set_xlabel('Time, s')
            ax.set_ylabel('Actuator 2, mm')
            plt.legend()
            plt.savefig('ExData/ActuatorPosition2.png', dpi=300)
            plt.show()                  
   

            fig, ax = plt.subplots()
            plt.plot(sPressures1[:, 0], sPressures1[:, 1], label='EXUDYN', color='blue')
            plt.plot(tspan, p1M, label='MATLAB', color='red')
            plt.xlim(0, tEnd)
            plt.ylim(0, 100e5)
            plt.grid(True)
            # Add a title and axis labels
            ax.set_xlabel('Time, s')
            ax.set_ylabel('Pressure $p_1$, Pa')
            plt.legend()
            plt.savefig('ExData/Pressure1.png', dpi=300)
            plt.show()

            fig, ax = plt.subplots()
            plt.plot(sPressures1[:, 0], sPressures1[:, 2], label='EXUDYN', color='blue')
            plt.plot(tspan, p2M, label='MATLAB', color='red')
            plt.xlim(0, tEnd)
            plt.ylim(0, 70e5)
            plt.grid(True)
            # Add a title and axis labels
            ax.set_xlabel('Time, s')
            ax.set_ylabel('Pressure $p_2$, Pa')
            plt.legend()
            plt.savefig('ExData/Pressure2.png', dpi=300)
            plt.show()
            
            fig, ax = plt.subplots()
            plt.plot(sPressures2[:, 0], sPressures2[:, 1], label='EXUDYN', color='blue')
            plt.plot(tspan, p3M, label='MATLAB', color='red')
            plt.xlim(0, tEnd)
            plt.ylim(0, 60e5)
            plt.grid(True)
            # Add a title and axis labels
            ax.set_xlabel('Time, s')
            ax.set_ylabel('Pressure $p_3$, Pa')
            plt.legend()
            plt.savefig('ExData/Pressure3.png', dpi=300)
            plt.show()
            
            fig, ax = plt.subplots()
            plt.plot(sPressures2[:, 0], sPressures2[:, 2], label='EXUDYN', color='blue')
            plt.plot(tspan, p4M, label='MATLAB', color='red')
            plt.xlim(0, tEnd)
            plt.ylim(0, 80e5)
            plt.grid(True)
            # Add a title and axis labels
            ax.set_xlabel('Time, s')
            ax.set_ylabel('Pressure $p_4$, Pa')
            plt.legend()
            plt.savefig('ExData/Pressure4.png', dpi=300)
            plt.show()
            
            
            
    return