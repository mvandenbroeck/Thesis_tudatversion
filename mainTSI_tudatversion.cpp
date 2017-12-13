/*    Copyright (c) 2010-2015, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      161115    M. Van den Broeck Creation.
 *
 *    References
 *
 *    Notes
 *
 */

#include <Eigen/Core>
#include <math.h>
#include <iostream>

#include <boost/test/results_collector.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include <Tudat/InputOutput/matrixTextFileReader.h>
#include <Tudat/InputOutput/basicInputOutput.h>
#include <tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/applicationOutput.h>
#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
#include <tudat/Astrodynamics/BasicAstrodynamics/unifiedStateModelElementConversions.h>
#include <tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>

#include "Thesis/constants.h"
#include "Thesis/spacecraft.h"
#include "Thesis/stepSizeControlTSI.h"
#include "Thesis/integrationSettings.h"
#include "Thesis/frameTransformation.h"
#include "Thesis/taylorSeriesIntegratorGenMap.h"
#include "Thesis/discreteForceModelGen.h"
#include "Thesis/problemRecurrenceRelationsGen.h"

#include <cstdio>
#include <ctime>
#include <chrono>


int main()
{

    // Start measuring wall time and cpu time
    std::clock_t startCPUtime = std::clock( );
    auto startWALLtime = std::chrono::system_clock::now( );

    // Find filepath of this cpp file
    std::string cppFilePath( __FILE__ );

    // Find folder of this cpp file
    std::string cppFolder = cppFilePath.substr( 0 , cppFilePath.find_last_of("/\\")+1 );

    std::string outputFolder = cppFolder + "Output_tudatversion"; // (INPUT)

    /// Define all constants

    // Define constants (INPUT)
    double standardGravity = 9.80665;
    double centralBodyGravitationalParameter = 1.32712440018e20; // [m^3/s^2] Test case = orbit around Sun



    /// Define spacecraft parameters

    // Specific impulse (INPUT)
    double specificImpulse = 3000.0;

    // Keplerian state (INPUT)
    Eigen::Matrix< double, 6, 1 > initialStateKepler;
    initialStateKepler( 0 ) = 1.495960829265199e+11;
    initialStateKepler( 1 ) = 0.6;
    initialStateKepler( 2 ) = tudat::unit_conversions::convertDegreesToRadians( 169.0 );
    initialStateKepler( 3 )
            = tudat::unit_conversions::convertDegreesToRadians( 45.0 );
    initialStateKepler( 4 )
            = tudat::unit_conversions::convertDegreesToRadians( 80.0 );
    initialStateKepler( 5 ) = tudat::unit_conversions::convertDegreesToRadians( 15.0 );

    double initialMassOfSpacecraft = 2000.0;

    // Convert vehicle state from Keplerian elements to Cartesian elements.
    Eigen::VectorXd initialReducedStateUSM = tudat::orbital_element_conversions::convertKeplerianToUnifiedStateModelElements(
                initialStateKepler,
                centralBodyGravitationalParameter );

    // Load thrust file (INPUT)
    std::string fileName = cppFolder + "unitTests/thrustProfileExample.txt";

    std::cout << "Selected thrust profile: " << fileName << std::endl;

    /// Select interpolation type used to interpolate thrust vector:
    // L2 = 2nd order Lagrange Interpolation
    // L3 = 3rd order Lagrange Interpolation
    // CSI = Cubic Spline Interpolation
    // CSIcorrectGamma = CSI without gamma_dot = 0 assumption (MOST ACCURATE OPTION)
    // CSIimprovedMass = CSIcorrectGamma with improved mass propagation:
    //                   different time variable for the mass Taylor series.
    //                   h = t / dt so that (h - h_0)^k = 1 (for fixed step-sizes).
    // CSIseparateSum = CSIcorrectGamma with separated summation of the positive
    //                  and negative parts of the terms of the Taylor series in
    //                  order to prevent numerical errors to occur.

    std::string interpolationType = "CSIcorrectGamma"; // (INPUT)

    std::cout<< "Selected interpolation type: " << interpolationType << std::endl;


    /// Define integration settings

    // Order of Taylor Series Integration (INPUT)
    int order = 20;

    // Define initial and final time of integration (INPUT)
    double initialTime = 0.0; // [s]
    double finalTime = 10.0 * 365.25 * 86400.0; // [s]

    // Define tolerance four each state element (USM7)
    // Set the ABSOLUTE error tolerance (INPUT)
    Eigen::Matrix< double, 7, 1> errorTolerance;
    errorTolerance << 1.0e-15, //1.0e-16
            1.0e-15, //1.0e-16
            1.0e-15, //1.0e-16
            1.0e-18, //1.0e-19
            1.0e-18, //1.0e-19
            1.0e-18, //1.0e-19
            1.0e-18; //1.0e-19

    // Safety factor (INPUT)
    double safetyFactor = 0.001; // [-] Must be in interval (0.0, 1.0]

    // Initial step size (INPUT)
    double initialStepSize = 400.0;//0.01; // [s]

    // Minimum step size (set equal to maximumStepSize in case of fixed step size) (INPUT)
    double minimumStepSize = 400.0;//10; // [s]

    // Maximum step size (INPUT)
    double maximumStepSize = 400.0;//500; // [s]

    // Factor that determines how many interpolated points are added in between two steps of the integrator (INPUT)
    int refineFactor = 0; // [-]



    /// Make all the pointers

    // Constants
    ConstantsPointer constantsPointer =
            boost::make_shared< Constants >( centralBodyGravitationalParameter, standardGravity  );

    // Spacecraft
    SpacecraftPointer spacecraftPointer =
            boost::make_shared< Spacecraft >( initialMassOfSpacecraft,
                                              initialReducedStateUSM,
                                              specificImpulse,
                                              constantsPointer );
    spacecraftPointer->loadThrustForce( fileName );

    // Integration settings
    IntegrationSettingsPointer integrationSettingsPointer =
            boost::make_shared< IntegrationSettings >( initialTime, finalTime );
    integrationSettingsPointer->setErrorTolerance( errorTolerance ); // Set in pointer
    integrationSettingsPointer->setSafeftyFactor( safetyFactor ); // Set in pointer
    integrationSettingsPointer->setInitialStepSize( initialStepSize ); // Set in pointer
    integrationSettingsPointer->setMinimumStepSize( minimumStepSize ); // Set in pointer
    integrationSettingsPointer->setMaximumStepSize( maximumStepSize ); // Set in pointer
    integrationSettingsPointer->setOrderOfTaylorSeries( order ); // Set in pointer
    integrationSettingsPointer->setRefineFactor( refineFactor ); // Set in pointer

    // Step size control
    StepSizeControlTSIPointer stepSizeControlTSIPointer =
            boost::make_shared< StepSizeControlTSI >( integrationSettingsPointer );

    // Declare map's to save state and thrust history
    std::map< double, Eigen::MatrixXd > computedStateMap;
    std::map< double, Eigen::MatrixXd > computedThrustAccelerationMap;

    bool useCSI = false;
    bool useCorrectGamma = false;
    bool useImprovedMass = false;
    bool useSeparateSum = false;
    boost::shared_ptr< DiscreteForceModel > discreteForceModelPointer;
    if ( interpolationType == "L2" )
    {
        discreteForceModelPointer = boost::make_shared< DiscreteForceModelL2 >(
                    spacecraftPointer,
                    integrationSettingsPointer,
                    constantsPointer );
    }
    else if ( interpolationType == "L3" )
    {
        discreteForceModelPointer = boost::make_shared< DiscreteForceModelL3 >(
                    spacecraftPointer,
                    integrationSettingsPointer,
                    constantsPointer );
    }
    else if ( interpolationType == "CSI" )
    {
        discreteForceModelPointer = boost::make_shared< DiscreteForceModelCSI >(
                    spacecraftPointer,
                    integrationSettingsPointer,
                    constantsPointer );
        useCSI = true;
    }
    else if ( interpolationType == "CSIcorrectGamma" )
    {
        discreteForceModelPointer = boost::make_shared< DiscreteForceModelCSIcorrectGamma >(
                    spacecraftPointer,
                    integrationSettingsPointer,
                    constantsPointer );
        useCSI = true;
        useCorrectGamma = true;
    }
    else if ( interpolationType == "CSIimprovedMass" )
    {
        discreteForceModelPointer = boost::make_shared< DiscreteForceModelCSIimprovedMass >(
                    spacecraftPointer,
                    integrationSettingsPointer,
                    constantsPointer );
        useCSI = true;
        useCorrectGamma = true;
        useImprovedMass = true;
    }
    else if ( interpolationType == "CSIseparateSum")
    {
        discreteForceModelPointer = boost::make_shared< DiscreteForceModelCSIcorrectGamma >(
                    spacecraftPointer,
                    integrationSettingsPointer,
                    constantsPointer );
        useCSI = true;
        useCorrectGamma = true;
        useSeparateSum = true;
    }

    // Problem recurrence relations
    ProblemRecurrenceRelationsPointer problemRecurrenceRelationsPointer =
            boost::make_shared< ProblemRecurrenceRelations >(
                useCSI,
                useImprovedMass,
                useSeparateSum,
                order,
                spacecraftPointer,
                constantsPointer );

    // Make pointer of class TaylorSeriesIntegrator
    TaylorSeriesIntegratorPointer taylorSeriesIntegratorPointer =
            boost::make_shared< TaylorSeriesIntegrator >(
                useCSI,
                useCorrectGamma,
                problemRecurrenceRelationsPointer,
                discreteForceModelPointer,
                spacecraftPointer,
                integrationSettingsPointer,
                constantsPointer,
                stepSizeControlTSIPointer );

    /// Perform integration
    taylorSeriesIntegratorPointer->integrate();

    /// Get computed ouput
    computedStateMap = taylorSeriesIntegratorPointer->getStateMap();
    computedThrustAccelerationMap = taylorSeriesIntegratorPointer->getComputedThrustAccelerationMap();

    // Stop WALL and CPU time counter and print on screen
    double CPUduration = ( std::clock( ) - startCPUtime ) / (double)CLOCKS_PER_SEC;
    std::chrono::duration< double > WALLduration = ( std::chrono::system_clock::now() - startWALLtime );

    std::cout << std::endl << "-----" << std::endl;
    std::cout << "CPU time duration = " << CPUduration << " s." << std::endl;
    std::cout << "Wall time duration = " << WALLduration.count( ) << " s." << std::endl;

    // Store output in map
    std::map< double, Eigen::Matrix< double, 1, 14 > > mapOfOutput;

    // Put the computedStateMap and computedThrustAccelerationMap in one map with the same dependentVariable
    std::map< double, Eigen::MatrixXd >::iterator variableIterator1 = computedStateMap.begin( );
    std::map< double, Eigen::MatrixXd >::iterator variableIterator2 = computedThrustAccelerationMap.begin( );

    //??
    std::cout << "Size of computedStateMap is " << computedStateMap.size( ) << std::endl;

    for( unsigned int i = 0; i < computedStateMap.size( ) - 1; i++ )
    {
        mapOfOutput[ variableIterator1->first ] << variableIterator1->second.adjoint( ),
                variableIterator2->second;

        variableIterator1++;
        variableIterator2++;
    }
    // Add last element to the map
    mapOfOutput[ variableIterator1-> first ] << variableIterator1-> second.adjoint( ),
            NAN, NAN, NAN;

    // Write spacecraft propagation history to file. (INPUT)
    tudat::input_output::writeDataMapToTextFile( mapOfOutput,
                                                 "mainTSI_" + interpolationType + "_ExampleOutput.dat",
                                                 outputFolder,
                                                 "",
                                                 std::numeric_limits< double >::digits10,
                                                 std::numeric_limits< double >::digits10,
                                                 "," );
}
