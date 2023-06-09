import exudyn as exu
from exudyn.itemInterface import *
from exudyn.utilities import *
from exudyn.FEM import *

import math as mt
from math import sin, cos, sqrt, pi, tanh
import time

SC  = exu.SystemContainer()
mbs = SC.AddSystem()


def LiftBoomModes(feL, FineMesh , nModes, FreeModes, HCB, ModeAnimation):
    
    if FineMesh:
        
        if FreeModes:
            
            print("Compute Free modes... ")
            feL.ComputeEigenmodes(nModes, excludeRigidBodyModes = 6, useSparseSolver = True)
            print("eigen freq.=", feL.GetEigenFrequenciesHz())
            
        if HCB:
           
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
           
           
           print("Compute Craig-Bampton modes... ")
           #boundaryList        = [nodeListJoint1, nodeListPist1, nodeListCyl2, nodeListJoint2, nodeListJoint3]
           boundaryList         = [nodeListJoint1, nodeListPist1]
           
           start_time           = time.time()
           feL.ComputeHurtyCraigBamptonModes(boundaryNodesList=boundaryList, nEigenModes=nModes, 
                                                    useSparseSolver=True,computationMode = HCBstaticModeSelection.RBE2)
                                                       
           print("Hurty-Craig Bampton modes... ")
           print("eigen freq.=", feL.GetEigenFrequenciesHz())
           print("HCB modes needed %.3f seconds" % (time.time() - start_time))  
            
            
        if ModeAnimation:
            from exudyn.interactive import AnimateModes
                
            cms     = ObjectFFRFreducedOrderInterface(feL)
                
            objFFRF = cms.AddObjectFFRFreducedOrder(mbs, positionRef=[0,0,0],initialVelocity=[0,0,0],initialAngularVelocity=[0,0,0],
                                                        gravity=[0,-9.81,0],color=[0.1,0.9,0.1,1.],)             
            
            mbs.Assemble()
                
            SC.visualizationSettings.nodes.show = False
            SC.visualizationSettings.openGL.showFaceEdges = True
            SC.visualizationSettings.openGL.multiSampling=4   
            SC.visualizationSettings.window.renderWindowSize = [1600,1080]
            SC.visualizationSettings.contour.outputVariable = exu.OutputVariableType.DisplacementLocal
            SC.visualizationSettings.contour.outputVariableComponent = 0 #component
                
            SC.visualizationSettings.general.autoFitScene = False #otherwise, model may be difficult to be moved
            nodeNumber = objFFRF['nGenericODE2'] #this is the node with the generalized coordinates
            AnimateModes(SC, mbs, nodeNumber)
            import sys
            sys.exit()   
                
    else:
        
        if  FreeModes:
            
            print("Compute Free modes... ")
            feL.ComputeEigenmodes(nModes, excludeRigidBodyModes = 6, useSparseSolver = True)
            print("eigen freq.=", feL.GetEigenFrequenciesHz())  
            
        if HCB:
           
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
           
           
           print("Compute Craig-Bampton modes... ")
           #boundaryList        = [nodeListJoint1, nodeListPist1, nodeListCyl2, nodeListJoint2, nodeListJoint3]
           boundaryList        = [nodeListJoint1, nodeListPist1]
           
           
           start_time           = time.time()
           feL.ComputeHurtyCraigBamptonModes(boundaryNodesList=boundaryList, nEigenModes=nModes, 
                                                    useSparseSolver=True,computationMode = HCBstaticModeSelection.RBE2)
                                                       
           print("Hurty-Craig Bampton modes... ")
           print("eigen freq.=", feL.GetEigenFrequenciesHz())
           print("HCB modes needed %.3f seconds" % (time.time() - start_time))  
           
           
        if ModeAnimation:
                from exudyn.interactive import AnimateModes
                
                cms     = ObjectFFRFreducedOrderInterface(feL)
                
                objFFRF = cms.AddObjectFFRFreducedOrder(mbs, positionRef=[0,0,0],initialVelocity=[0,0,0],initialAngularVelocity=[0,0,0],
                                                        gravity=[0,-9.81,0],color=[0.1,0.9,0.1,1.],)             
            
                mbs.Assemble()
                
                SC.visualizationSettings.nodes.show = False
                SC.visualizationSettings.openGL.showFaceEdges = True
                SC.visualizationSettings.openGL.multiSampling=4   
                SC.visualizationSettings.window.renderWindowSize = [1600,1080]
                SC.visualizationSettings.contour.outputVariable = exu.OutputVariableType.DisplacementLocal
                SC.visualizationSettings.contour.outputVariableComponent = 0 #component
                
                SC.visualizationSettings.general.autoFitScene = False #otherwise, model may be difficult to be moved
                nodeNumber = objFFRF['nGenericODE2'] #this is the node with the generalized coordinates
                AnimateModes(SC, mbs, nodeNumber)
                import sys
                sys.exit()  
                
    
    return




def TiltBoomModes(feT, FineMesh , nModes, FreeModes, HCB, ModeAnimation):
    
    if FineMesh:
        
        if FreeModes:
            
            print("Compute Free modes... ")
            feT.ComputeEigenmodes(nModes, excludeRigidBodyModes = 6, useSparseSolver = True)
            print("eigen freq.=", feT.GetEigenFrequenciesHz())
            
        if HCB:
           
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
          
           print("Compute Craig-Bampton modes... ")
           boundaryList         = [nodeListJoint1T, nodeListPist1T]
           
           start_time           = time.time()
           feT.ComputeHurtyCraigBamptonModes(boundaryNodesList=boundaryList, nEigenModes=nModes, 
                                                    useSparseSolver=True,computationMode = HCBstaticModeSelection.RBE2)
                                                       
           print("Hurty-Craig Bampton modes... ")
           print("eigen freq.=", feT.GetEigenFrequenciesHz())
           print("HCB modes needed %.3f seconds" % (time.time() - start_time))  
            
            
        if ModeAnimation:
            from exudyn.interactive import AnimateModes
                
            cms     = ObjectFFRFreducedOrderInterface(feT)
                
            objFFRF = cms.AddObjectFFRFreducedOrder(mbs, positionRef=[0,0,0],initialVelocity=[0,0,0],initialAngularVelocity=[0,0,0],
                                                        gravity=[0,-9.81,0],color=[0.1,0.9,0.1,1.],)             
            
            mbs.Assemble()
                
            SC.visualizationSettings.nodes.show = False
            SC.visualizationSettings.openGL.showFaceEdges = True
            SC.visualizationSettings.openGL.multiSampling=4   
            SC.visualizationSettings.window.renderWindowSize = [1600,1080]
            SC.visualizationSettings.contour.outputVariable = exu.OutputVariableType.DisplacementLocal
            SC.visualizationSettings.contour.outputVariableComponent = 0 #component
                
            SC.visualizationSettings.general.autoFitScene = False #otherwise, model may be difficult to be moved
            nodeNumber = objFFRF['nGenericODE2'] #this is the node with the generalized coordinates
            AnimateModes(SC, mbs, nodeNumber)
            import sys
            sys.exit()   
                
    else:
        
        if  FreeModes:
            
            print("Compute Free modes... ")
            feT.ComputeEigenmodes(nModes, excludeRigidBodyModes = 6, useSparseSolver = True)
            print("eigen freq.=", feT.GetEigenFrequenciesHz())  
            
        if HCB:
           
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
           #p18                 = [1.84,0.343,-4.80e-2]
           #p18s                = [1.84,0.339,-4.37e-2]
           #p17                 = [1.84,0.189,-4.80e-2]
           #p17s                = [1.84,0.193,-4.37e-2]
           #p16                 = [1.84,0.189, 4.80e-2]
           #p16s                = [1.84,0.193, 4.37e-2]
           #p15                 = [1.84,0.343, 4.80e-2]
           #p15s                = [1.84,0.339, 4.37e-2]
           
           #radius8             = 3.2e-002
           #nodeListCyl2T       = feT.GetNodesOn(p15, p16, radius8, tolerance=1e-2)  
           #pJoint3T            = feT.GetNodePositionsMean(nodeListCyl2T)
           #nodeListCyl2TLen    = len(nodeListCyl2T)
           #noodeWeightsCyl2T   = [1/nodeListCyl2TLen]*nodeListCyl2TLen       
          
           print("Compute Craig-Bampton modes... ")
           boundaryList         = [nodeListJoint1T, nodeListPist1T]
           
           
           start_time           = time.time()
           feT.ComputeHurtyCraigBamptonModes(boundaryNodesList=boundaryList, nEigenModes=nModes, 
                                                    useSparseSolver=True,computationMode = HCBstaticModeSelection.RBE2)
                                                       
           print("Hurty-Craig Bampton modes... ")
           print("eigen freq.=", feT.GetEigenFrequenciesHz())
           print("HCB modes needed %.3f seconds" % (time.time() - start_time))  
           
           
        if ModeAnimation:
                from exudyn.interactive import AnimateModes
                
                cms     = ObjectFFRFreducedOrderInterface(feT)
                
                objFFRF = cms.AddObjectFFRFreducedOrder(mbs, positionRef=[0,0,0],initialVelocity=[0,0,0],initialAngularVelocity=[0,0,0],
                                                        gravity=[0,-9.81,0],color=[0.1,0.9,0.1,1.],)             
            
                mbs.Assemble()
                
                SC.visualizationSettings.nodes.show = False
                SC.visualizationSettings.openGL.showFaceEdges = True
                SC.visualizationSettings.openGL.multiSampling=4   
                SC.visualizationSettings.window.renderWindowSize = [1600,1080]
                SC.visualizationSettings.contour.outputVariable = exu.OutputVariableType.DisplacementLocal
                SC.visualizationSettings.contour.outputVariableComponent = 0 #component
                
                SC.visualizationSettings.general.autoFitScene = False #otherwise, model may be difficult to be moved
                nodeNumber = objFFRF['nGenericODE2'] #this is the node with the generalized coordinates
                AnimateModes(SC, mbs, nodeNumber)
                import sys
                sys.exit()  
                
    
    return
