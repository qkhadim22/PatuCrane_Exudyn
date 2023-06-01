# PatuCrane_Exudyn
This file includes the simulation of hydraulically actuated rigid multibody crane (PATU655 crane) in Exudyn environment 
using redundtant and minimal coordinates. The details about rigid bodies and hydraulics in the crane can be found in PatuNotes.pdf.  

## Installation
Install Exudyn software for running PatuCrane_Exudyn.py using the instructions https://exudyn.readthedocs.io/en/latest/docs/RST/InstallationInstructions.html.

## Usage
Use 'RedundantCoordinates=False' to run Patu model with minimal coordinates. 'RedundantCoordinates=True' adds redundant coordinates. 'Hydraulics=True' and 'useFriction=True'
adds hydraulics and friction in the cylinders, respectively. Use 'Plotting=True' for making plots. Note the difference in the minimal coordinates w.r.t MATLAB is due to difference of definition of tilt boom angle in Exduyn and MATLAB environment. Otherwise, MATLAB and EXUDYN agrees very well.

## Contribution:
Contributors      : Pauli Pöppönen (Graphics) MrValhe (Github) 

## License:
This model is free to use.

## Contact:
Author            : Qasim Khadim
Contact           : qasim.khadim@outlook.com,qkhadim22 (Github)

## Acknowledgements:

The author would like to thanks Prof. Johannes Gerstmayr for hosting me at University of Innsbruck, Department of Mechatronics, Innsbruck, Austria. His supervision enabled me to understand and write this code in Exudyn environment. Special thanks to Stefan Holzinger for helping me in ANSYS mesh export and modeling FFRF in exudyn.  I would also like to thanks Prof. Aki Mikkola, Prof. Emil Kurvinen and Dr. Grzegorz Orzechowski. Special thanks to Prof. Emil Kurvinen for giving me this opportunity to work at University of Innsbruck.

## Funding:
SANTTU Project funded by Business Finland. 

## REFERENCES:
1. S. Holzinger, J. Schöberl, J. Gerstmayr. The equations of motion for a rigid body using non-redundant unified local velocity coordinates. Multibody System Dynamics, Vol. 48, pp. 283 – 309, 2020. 
2. M. Sereinig, P. Manzl, and J. Gerstmayr. Task Dependent Comfort Zone, a Base Placement Strategy for Autonomous Mobile Manipulators using Manipulability Measures, Robotics and Autonomous Systems, submitted.
3. Q. Khadim et al., "State Estimation in a Hydraulically Actuated Log Crane Using Unscented Kalman Filter," in IEEE Access, vol. 10, pp. 62863-62878, 2022, doi: 10.1109/ACCESS.2022.3179591 (Hydraulics).
<<<<<<< HEAD
4. Q. Khadim et al., "Experimental Investigation into the State Estimation of a Forestry Crane Using the Unscented Kalman Filter and a Multiphysics Model," in Mechanism and Machine Theory, Accepted, 2023, doi: https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4285871 (Hydraulics-Verification).
=======
4. Q. Khadim et al., "Experimental Investigation into the State Estimation of a Forestry Crane Using the Unscented Kalman Filter and a Multiphysics Model," in Mechanism and Machine Theory, Accepted, 2023, doi: https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4285871 (Hydraulics-Verification).
>>>>>>> 83a720af2265b144acd8ad8e536c8b10c5ed6f27
