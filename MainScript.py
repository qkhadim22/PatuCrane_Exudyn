#%%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#This file includes the simulation of hydraulically actuated rigid multibody crane (PATU655 crane) in Exudyn environment 
#using redundtant and minimal coordinates. The details about rigid bodies and hydraulics in the crane can be found in 
#PatuNotes.pdf.  

# Author            : Qasim Khadim
# Contributors      : Pauli Pöppönen (Graphics) MrValhe (Github)
# Contact           : qasim.khadim@outlook.com,qkhadim22 (Github)
# Dated             : 11-05-2023
# Organizations     : University of Oulu in the collaboration of LUT University and University of Innsbruck.
# Funding           : Business of Finland (SANTTU Project)
# Copyright         :

#%%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

FlexibleModel = True

if FlexibleModel:
    from Models.FlexibleSystem import FFRFHydraulics
    FFRFHydraulics(FineMesh=False, nModes=4, loadFromSavedNPY=True, ComputeModes =False, FreeModes = False, HCB= True, ModeAnimation = True,
                   StaticCase= True, Hydraulics=False, useFriction=True, Visualization=True, Plotting=False) 
    
else:
    from Models.RigidSystem import RigidMultibodyHydraulics
    RigidMultibodyHydraulics(RedundantCoordinates=True, Hydraulics=True, useFriction=True, 
                                                  Plotting=True)                          






# %%
