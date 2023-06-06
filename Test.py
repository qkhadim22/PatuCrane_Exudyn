import exudyn as exu
from exudyn.itemInterface import *
SC = exu.SystemContainer()
mbs = SC.AddSystem()
#create simple system:
ground = mbs.AddObject(ObjectGround())
mbs.AddNode(NodePoint())
body = mbs.AddObject(MassPoint(physicsMass=1, nodeNumber=0))
m0 = mbs.AddMarker(MarkerBodyPosition(bodyNumber=ground))
m1 = mbs.AddMarker(MarkerBodyPosition(bodyNumber=body))
mbs.AddObject(CartesianSpringDamper(markerNumbers=[m0,m1], stiffness=[100,100,100]))
mbs.AddLoad(LoadForceVector(markerNumber=m1, loadVector=[10,10,10]))
#
mbs.Assemble()
simulationSettings = exu.SimulationSettings()
simulationSettings.timeIntegration.endTime = 10
success = exu.SolveDynamic(simulationSettings, storeSolver = True)
print("success =", success)
print("iterations = ", mbs.sys['dynamicSolver'].it)
print("pos=", mbs.GetObjectOutputBody(body,localPosition=[0,0,0],
      variableType=exu.OutputVariableType.Position))