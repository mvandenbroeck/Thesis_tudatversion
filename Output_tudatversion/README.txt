mainTSI output structure
------------------------

Output file columns (separated by ','):
 1st column:	Time [s],
 2nd column:	C [m/s] (first element of USM7 state vector)
 3rd column:	R_{f1} [m/s] (second element of USM7 state vector)
 4th column:	R_{f2} [m/s] (third element of USM7 state vector)
 5th column:	\epsilon_1 [-] (fourth element of USM7 state vector)
 6th column:	\epsilon_2 [-] (fifth element of USM7 state vector)
 7th column:	\epsilon_3 [-] (sixth element of USM7 state vector)
 8th column:	\eta [-] (seventh element of USM7 state vector)
 9th column:	Spacecraft mass [kg]
10th column:	v_{e1} [m/s]**
11th column:	v_{e2} [m/s]**
12th column:	\omega_3 [1/s]**
13th column:	Thrust acceleration in x direction of velocity frame* (m/s^2)
14th column:	Thrust acceleration in y direction of velocity frame* (m/s^2)
15th column:	Thrust acceleration in z direction of velocity frame* (m/s^2)

*As explained on page 16 of https://repository.tudelft.nl/islandora/object/uuid:2567c152-ab56-4323-bcfa-b076343664f9/datastream/OBJ/download,
the velocity frame consists of the folliwing axes:
x-axis		tangent to the trajectory, pointing in the direction of motion
y-axis		normal to the x-axis, in the orbital plane and in the direction of the positivie derivative of the tangential unit vector with respect to the arclength parameter of the trajectory
z-axis		perpendicular to the x-axis an y-axis

**See https://repository.tudelft.nl/islandora/object/uuid:2567c152-ab56-4323-bcfa-b076343664f9/datastream/OBJ/download for more information

mainRKusm output structure
--------------------------

Output file columns (separated by ','):
 1st column:	Time [s],
 2nd column:	C [m/s] (first element of USM7 state vector)
 3rd column:	R_{f1} [m/s] (second element of USM7 state vector)
 4th column:	R_{f2} [m/s] (third element of USM7 state vector)
 5th column:	\epsilon_1 [-] (fourth element of USM7 state vector)
 6th column:	\epsilon_2 [-] (fifth element of USM7 state vector)
 7th column:	\epsilon_3 [-] (sixth element of USM7 state vector)
 8th column:	\eta [-] (seventh element of USM7 state vector)
 9th column:	Spacecraft mass [kg]

The thrust accelerations are stored in a seperate file starting with "mainRKusmThrustAcc"
Output file columns (separated by ','):
 1st column:	Time [s],
 2nd column:	Thrust acceleration in x direction of velocity frame* (m/s^2)
 3nd column:	Thrust acceleration in y direction of velocity frame* (m/s^2)
 4nd column:	Thrust acceleration in z direction of velocity frame* (m/s^2)

*As explained on page 16 of https://repository.tudelft.nl/islandora/object/uuid:2567c152-ab56-4323-bcfa-b076343664f9/datastream/OBJ/download,
the velocity frame consists of the folliwing axes:
x-axis		tangent to the trajectory, pointing in the direction of motion
y-axis		normal to the x-axis, in the orbital plane and in the direction of the positivie derivative of the tangential unit vector with respect to the arclength parameter of the trajectory
z-axis		perpendicular to the x-axis an y-axis