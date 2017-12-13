------------------------
mainTSI_tudatversion.cpp
------------------------

THE PROGRAM

To propagate the state of a spacecraft expressed
 in the Unified State Model (USM) with the use of 
Taylor Series Integration (TSI), do the following:
1. Open mainTSI_tudatversion.cpp
2. Provide relevant inputs on the following lines:
(You can also find these inputs by using CTRL+F for "INPUT")
All inputs are in SI units, unless mentioned otherwise
	77: The output folder
	83: The gravitational parameter of your main body
	90: The specific impulse of the motor of your spacacraft
	94-101: The initial Keplerian state of your spacecraft
	103: The initial mass of the spacecraft
	111: The folder in which the thrust profile is specified
	127: The method used to interpolate the thrust profile (see list of options)
	135: The order of the TSI
	138: The integration start time
	139: The integration end time
	144-150: The absolute error tolerances on each of the seven USM7 elements
	153: Safety factor for the step-size, (used in combination with tolerances to enhance TSI accuracy
	156: Initial integration step-size
	159: Minimum integration step-size
	162: Maximum integration step-size
	165: Factor to interpolate between integration steps to refine the output
	320: Name of output file
3. Build application_mainTSI_tudatversion
4. Run application_mainTSI_tudatversion
5. The output is stored in a .txt file in the specified folder and file.
See the README.txt file in the default output folder to interprete the output.


THE ASSUMPTIONS

In the default setup, the following assumptions are made:
1. The spacecraft is subject to a central gravitational field.
2. Besides the gravitational force and the thrust force, no other forces act on the spacecraft.
3. The fourth and higher order derivatives of the thrust accelerations in the velocity frame are zero. (See pages 16, 51*)
4. The fourth and higher order derivatives of the thrust accelerations in the RTN frame are zero. (See pages 17, 56*)


THE VALIDATION

As explained in the Conclusion* the developed TSI-USM integrator was found to be 
sufficiently** accurate for cases in which the thrust force is constant.

*https://repository.tudelft.nl/islandora/object/uuid:2567c152-ab56-4323-bcfa-b076343664f9/datastream/OBJ/download

**A sufficient accuracy is defined as an absolute accuracy in position and velocity
of below 1 km and 1 m/s respectively for interplanetary low-thrust trajectories (semi-major axis > 1 AU) after an integration
period of 10 years.



--------------------------
mainRKusm_tudatversion.cpp
--------------------------

THE PROGRAM

To propagate the state of a spacecraft expressed
 in the Unified State Model (USM) with the use of 
the Runge-Kutta 8(7) integrator with 13 stages (RK8(7)13M),
 do the following:
1. Open mainRKusm_tudatversion.cpp
2. Provide relevant inputs on the following lines:
(You can also find these inputs by using CTRL+F for "INPUT")
All inputs are in SI units, unless mentioned otherwise
	116x: The gravitational parameter of your main body
	120: The specific impulse of the motor of your spacacraft
	125: The integration start time
	126: The integration end time
	127: Initial integration step-size
	128: Minimum integration step-size
	129: Maximum integration step-size
	131: The absolute and relative error tolerance limits
	140: The folder in which the thrust profile is specified
	145-153: The initial Keplerian state of your spacecraft
	155: The initial mass of the spacecraft
	289: Name of output file containing the spacecraft state history
	290: Output folder for the spacecraft state history
	298: Name of output file containing the thrust acceleration history
	299: Output folder for the thrust acceleration history
3. Build application_mainRKusm_tudatversion
4. Run application_mainRKusm_tudatversion
5. The output is stored in a .txt file in the specified folder and file.
See the README.txt file in the default output folder to interprete the output.


THE ASSUMPTIONS

In the default setup, the following assumptions are made:
1. The spacecraft is subject to a central gravitational field.
2. Besides the gravitational force and the thrust force, no other forces act on the spacecraft.


THE VALIDATION

The RK-USM integrator displays a sufficiently* accurate propagation for specific thrust cases (e.g. constant tangential thrust), but it requires
further testing of different cases to guarantee its validation.

*A sufficient accuracy is defined as an absolute accuracy in position and velocity
of below 1 km and 1 m/s respectively for interplanetary low-thrust trajectories (semi-major axis > 1 AU) after an integration
period of 10 years.