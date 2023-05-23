#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# This is an EXUDYN example
#
# Details:  Stiff flyball governor (iftomm benchmark problem) with kinematic tree
#           Ref.: https://www.iftomm-multibody.org/benchmark/problem/Stiff_flyball_governor/
#
# Author:   Johannes Gerstmayr
# Date:     2022-8-22
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

import exudyn as exu
from exudyn.itemInterface import *
from exudyn.utilities import *
from exudyn.graphicsDataUtilities import *

#from modelUnitTests import ExudynTestStructure, exudynTestGlobals
import numpy as np
from numpy import linalg as LA

SC          = exu.SystemContainer()
mbs         = SC.AddSystem()

#%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
g               = [0,-9.81,0] # Gravity
tEnd            = 20 #simulation time
h               = 1e-3 #step size
nSteps          = int(tEnd/h)+2

# Loading Graphics of bodies
fileName1       = 'Graphics_Exudyn/Pillar.stl'
fileName2       = 'Graphics_Exudyn/LiftBoom.stl'
fileName3       = 'Graphics_Exudyn/TiltBoom+ExtensionBoom.stl'
fileName4       = 'Graphics_Exudyn/Bracket1.stl'
fileName5       = 'Graphics_Exudyn/Bracket2.stl'


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

graphicsBody1   = GraphicsDataFromSTLfile(fileName1, color4black,verbose=False, 
                                                  invertNormals=True,invertTriangles=True)
graphicsBody1   = AddEdgesAndSmoothenNormals(graphicsBody1, edgeAngle=0.25*pi,addEdges=True, 
                                                     smoothNormals=True)
graphicsCOM1    = GraphicsDataBasis(origin=pMid1, length=2*W1)


nRigidBodyNodes = 1
graphicsList    = [[graphicsCOM1, graphicsBody1]]

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


inertiaLinkCOM = np.array([[16.358844,-1.27808, 1.7e-5],
                           [-1.27808, 0.612552, -5.9e-5],
                           [1.7e-5,  -5.9e-5  , 16.534255]]) #KinematicTree requires COM inertia

linkCOM = pMid1 #if COM=0, gravity does not act on pendulum!


offsetsList = exu.Vector3DList([[0,0,0]])

rotList = exu.Matrix3DList([np.eye(3)])

linkCOMs=exu.Vector3DList([linkCOM])

linkInertiasCOM=exu.Matrix3DList([inertiaLinkCOM])



nGeneric = mbs.AddNode(NodeGenericODE2(referenceCoordinates=[0.],initialCoordinates=[0.],
                                       initialCoordinates_t=[0.],numberOfODE2Coordinates=1))


oKT = mbs.AddObject(ObjectKinematicTree(nodeNumber=nGeneric, jointTypes=[exu.JointType.RevoluteZ], linkParents=[-1],
                                  jointTransformations=rotList, jointOffsets=offsetsList, linkInertiasCOM=linkInertiasCOM,
                                  linkCOMs=linkCOMs, linkMasses=[93.26], gravity=g,
                                  visualization=VObjectKinematicTree(graphicsDataList = graphicsList)))



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

exu.SolveDynamic(mbs, simulationSettings=simulationSettings,
                 solverType=exu.DynamicSolverType.TrapezoidalIndex2)

from exudyn.interactive import SolutionViewer
SolutionViewer(mbs)
