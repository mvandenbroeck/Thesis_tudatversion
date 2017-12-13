/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <Tudat/SimulationSetup/tudatSimulationHeader.h> //inlcudes everything that has to be included for simulations
#include <tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/applicationOutput.h>

#include <tudat/Mathematics/NumericalIntegrators/rungeKuttaVariableStepSizeIntegrator.h>
#include <tudat/Mathematics/NumericalIntegrators/rungeKuttaCoefficients.h>
#include <tudat/Mathematics/NumericalIntegrators/numericalIntegrator.h>
#include <tudat/Mathematics/Interpolators/cubicSplineInterpolator.h>
#include <tudat/Mathematics/BasicMathematics/nearestNeighbourSearch.h>

#include <tudat/InputOutput/matrixTextFileReader.h>
#include <tudat/Astrodynamics/BasicAstrodynamics/unifiedStateModelElementConversions.h>

#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "constants.h"
#include "spacecraftGen.h"
#include "frameTransformation.h"
#include "USMstateDerivativeModel.h"

#include <cstdio>
#include <ctime>
#include <chrono>

/// This function computes the coefficients of all cubic polynomials of the cubic spline using
/// cubic spline interpolation for the thrust force in velocity frame in X, Y, Z direction and
/// resultant thrust force.
void computeCSIcoefficients( SpacecraftGenPointer spacecraftGenPointer )
{
    // Get thrustForceMatrix from spacecraftGenPointer
    Eigen::Matrix< double, Eigen::Dynamic, 4 > thrustForceMatrix =
            spacecraftGenPointer->getThrustForceMatrix();

    /// Interpolate thrust profile
    // Define time and thrust values in std::vector
    std::vector< double > timePointsToBeInterpolated( thrustForceMatrix.rows() );
    std::vector< double > thrustPointsToBeInterpolatedX( thrustForceMatrix.rows() );
    std::vector< double > thrustPointsToBeInterpolatedY( thrustForceMatrix.rows() );
    std::vector< double > thrustPointsToBeInterpolatedZ( thrustForceMatrix.rows() );
    std::vector< double > resultantThrustPointsToBeInterpolated( thrustForceMatrix.rows() );

    for ( int i = 0; i < thrustForceMatrix.rows(); i++ )
    {
        timePointsToBeInterpolated[ i ] = thrustForceMatrix( i, 0 );
        thrustPointsToBeInterpolatedX[ i ] = thrustForceMatrix( i, 1 );
        thrustPointsToBeInterpolatedY[ i ] = thrustForceMatrix( i, 2 );
        thrustPointsToBeInterpolatedZ[ i ] = thrustForceMatrix( i, 3 );
        resultantThrustPointsToBeInterpolated[ i ] = ( ( thrustForceMatrix.rightCols( 3 ) ).row( i ) ).norm();
    }

    // Calculate the four coefficients of all splines (number of splines = 1 - number of nodes)
    // of the cubic spline interpolation.
    tudat::interpolators::CubicSplineInterpolatorDoublePointer thrustInVelocicityFrameXcsiPointer =
            boost::make_shared< tudat::interpolators::CubicSplineInterpolatorDouble >(
                timePointsToBeInterpolated,
                thrustPointsToBeInterpolatedX );
    tudat::interpolators::CubicSplineInterpolatorDoublePointer thrustInVelocicityFrameYcsiPointer =
            boost::make_shared< tudat::interpolators::CubicSplineInterpolatorDouble >(
                timePointsToBeInterpolated,
                thrustPointsToBeInterpolatedY );
    tudat::interpolators::CubicSplineInterpolatorDoublePointer thrustInVelocicityFrameZcsiPointer =
            boost::make_shared< tudat::interpolators::CubicSplineInterpolatorDouble >(
                timePointsToBeInterpolated,
                thrustPointsToBeInterpolatedZ );
    tudat::interpolators::CubicSplineInterpolatorDoublePointer resultantThrustForceCSIpointer =
            boost::make_shared< tudat::interpolators::CubicSplineInterpolatorDouble >(
                timePointsToBeInterpolated,
                resultantThrustPointsToBeInterpolated );

    std::vector< std::vector< double > > CSIcoefficientsOfThrustForceInVelocityFrameX,
                                         CSIcoefficientsOfThrustForceInVelocityFrameY,
                                         CSIcoefficientsOfThrustForceInVelocityFrameZ,
                                         CSIcoefficientsOfResultantThrustForce;
    CSIcoefficientsOfThrustForceInVelocityFrameX = thrustInVelocicityFrameXcsiPointer->GetCoefficients();
    CSIcoefficientsOfThrustForceInVelocityFrameY = thrustInVelocicityFrameYcsiPointer->GetCoefficients();
    CSIcoefficientsOfThrustForceInVelocityFrameZ = thrustInVelocicityFrameZcsiPointer->GetCoefficients();
    CSIcoefficientsOfResultantThrustForce = resultantThrustForceCSIpointer->GetCoefficients();


    // Store coefficients in spacecraftGen pointer
    spacecraftGenPointer->setCSIcoefficientsOfThrustForceInVelocityFrameX( CSIcoefficientsOfThrustForceInVelocityFrameX );
    spacecraftGenPointer->setCSIcoefficientsOfThrustForceInVelocityFrameY( CSIcoefficientsOfThrustForceInVelocityFrameY );
    spacecraftGenPointer->setCSIcoefficientsOfThrustForceInVelocityFrameZ( CSIcoefficientsOfThrustForceInVelocityFrameZ );
    spacecraftGenPointer->setCSIcoefficientsOfResultantThrustForce( CSIcoefficientsOfResultantThrustForce );
}


//! Execute propagation of orbit of vehicle around the Earth. The vehicle is subject to a thrustforce, which is specified in
//! the thrustProfileExample.txt file. In that file, the first column is time in seconds, the last three columns give the x, y
//! and z components of the thrust force in the velocity frame (see page 16 of
//! https://repository.tudelft.nl/islandora/object/uuid:2567c152-ab56-4323-bcfa-b076343664f9/datastream/OBJ/download).
//! A thrust force in the x direction corresponds to thrust in the velocity direction.
int main()
{
    // Start measuring wall time and cpu time
    std::clock_t startCPUtime = std::clock( );
    auto startWALLtime = std::chrono::system_clock::now( );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        INTEGRATION SETTINGS           /////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define constants (INPUT)
    double centralBodyGravitationalParameter = 1.32712440018e20;
    double standardGravity = 9.80665;

    // Spacecraft settings (INPUT)
    double specificImpulse = 3000.0;

    // Step-size settings. (INPUT)
    // The minimum and maximum step-size are set such that the input data is fully accepted by the
    // integrator, to determine the steps to be taken.
    const double initialTime = 0.0;
    const double finalTime = 10.0 * 365.25 * 86400.0;
    const double initialStepSize = 400.0;
    const double minimumStepSize = 400.0;
    const double maximumStepSize = 400.0;
    const double relativeErrorTolerance = 1.0e12;//1.0e12;
    const double absoluteErrorTolerance = 1.0e12;//1.0e12;

    // Find filepath of this cpp file
    std::string cppFilePath( __FILE__ );

    // Find folder of this cpp file
    std::string cppFolder = cppFilePath.substr( 0 , cppFilePath.find_last_of("/\\")+1 );

    // Load thrust file (INPUT)
    std::string fileName = cppFolder + "unitTests/thrustProfileExample.txt";

    std::cout << "Selected thrust profile: " << fileName << std::endl;

    // Set Keplerian elements for vehicle. (INPUT)
    tudat::basic_mathematics::Vector6d initialStateKepler;
    initialStateKepler( 0 ) = 1.495960829265199e+11;
    initialStateKepler( 1 ) = 0.6;
    initialStateKepler( 2 ) = tudat::unit_conversions::convertDegreesToRadians( 169.0 );
    initialStateKepler( 3 )
            = tudat::unit_conversions::convertDegreesToRadians( 45.0 );
    initialStateKepler( 4 )
            = tudat::unit_conversions::convertDegreesToRadians( 80.0 );
    initialStateKepler( 5 ) = tudat::unit_conversions::convertDegreesToRadians( 15.0 );

    double initialMass = 2000.0; // (INPUT)

    // Convert vehicle state from Keplerian elements to Cartesian elements.
    Eigen::VectorXd initialReducedStateUSM = tudat::orbital_element_conversions::convertKeplerianToUnifiedStateModelElements(
                initialStateKepler,
                centralBodyGravitationalParameter );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        MAKE POINTERS          /////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize constants and spacecraftGen pointer
    // Constants
    ConstantsPointer constantsPointer =
            boost::make_shared< Constants >( centralBodyGravitationalParameter, standardGravity  );

    // Spacecraft
    SpacecraftGenPointer spacecraftGenPointer =
            boost::make_shared< SpacecraftGen >( initialMass,
                                              initialReducedStateUSM,
                                              specificImpulse,
                                              constantsPointer );
    spacecraftGenPointer->loadThrustForce( fileName );



    // Interpolate thrust profile and save coefficients of cubic spline of x, y, z
    // and resultant thrust in spacecraftGen pointer
    computeCSIcoefficients( spacecraftGenPointer );


    USMstateDerivativeModel usmStateDerivativeModel( spacecraftGenPointer,
                                                     constantsPointer );



    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PERFORM INTEGRATION           //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    Eigen::Matrix< double, 1, 8 > initialState = spacecraftGenPointer->getInitialState();

    boost::function< Eigen::VectorXd( const double, const Eigen::VectorXd& ) > stateDerivativeFunction =
        boost::bind( &USMstateDerivativeModel::computeUnifiedStateModelStateDerivative, &usmStateDerivativeModel, _1, _2 );

    tudat::numerical_integrators::RungeKuttaVariableStepSizeIntegratorXd integrator(
                         tudat::numerical_integrators::RungeKuttaCoefficients::get(
                     tudat::numerical_integrators::RungeKuttaCoefficients::rungeKutta87DormandPrince ),
                stateDerivativeFunction,
                initialTime,
                initialState,
                minimumStepSize,
                maximumStepSize,
                relativeErrorTolerance,
                absoluteErrorTolerance );

    // Set initial running time. This is updated after each step that the numerical integrator takes.
    double runningTime = 0.0;
    double stepSize = initialStepSize;

    // Initialize map to save state history
    std::map< double, Eigen::VectorXd > stateHistoryMap;
    std::map< double, Eigen::MatrixXd > thrustHistoryMap;

    // Compute current thrust accelerations at the initialTime and currentState (saved in spacecraftPointer already)
    usmStateDerivativeModel.computeCurrentAccelerations( initialTime );

    // Save initial state and initial thrust accelerations
    stateHistoryMap[ initialTime ] = initialState;
    thrustHistoryMap[ initialTime ] = spacecraftGenPointer->getCurrentAcceleration();

    // Make progress output
    double onePercentOfCPUtime = ( finalTime - initialTime ) / 100.0;
    double progress = 0.0;
    int percentage = 0;

    do
    {
        // Output progress
        if ( ( runningTime - progress ) >= onePercentOfCPUtime )
        {
            percentage++;
            progress = progress + onePercentOfCPUtime;
            std::cout << "[ " << percentage << "% ]" << std::endl;
        }


        // Make sure the integrator does not integrate beyond the end time.
        if ( std::fabs( finalTime - runningTime )
             <= std::fabs( stepSize ) * ( 1.0 + std::numeric_limits< double >::epsilon( ) ) )
        {
            stepSize = finalTime - integrator.getCurrentIndependentVariable( );
        }

        // Perform a single integration step.
        integrator.performIntegrationStep( stepSize );
        runningTime = integrator.getCurrentIndependentVariable( );
        stepSize = integrator.getNextStepSize( );

        // Store current state in spacecraftGenPointer
        spacecraftGenPointer->setCurrentState( integrator.getCurrentState( ) );

        // Compute accelerations at calculated new state and new time
        usmStateDerivativeModel.computeCurrentAccelerations( runningTime );

        // Store current acceleration and time in thrustHistoryMap
        thrustHistoryMap[ runningTime ] = spacecraftGenPointer->getCurrentAcceleration( );

        // Store current state and time in stateHistoryMap
        stateHistoryMap[ runningTime ] = integrator.getCurrentState( );

    } while( !( finalTime - runningTime <= std::numeric_limits< double >::epsilon( ) ) );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Stop WALL and CPU time counter and print on screen
    double CPUduration = ( std::clock( ) - startCPUtime ) / (double)CLOCKS_PER_SEC;
    std::chrono::duration< double > WALLduration = ( std::chrono::system_clock::now() - startWALLtime );

    std::cout << std::endl << "-----" << std::endl;
    std::cout << "CPU time duration = " << CPUduration << " s." << std::endl;
    std::cout << "Wall time duration = " << WALLduration.count( ) << " s." << std::endl;

    // Output size of stateHistoryMap and thrustHistoryMap
    std::cout << "Size of stateHistoryMap = " << stateHistoryMap.size() << std::endl;
    std::cout << "Size of thrustHistoryMap = " << thrustHistoryMap.size() << std::endl;

    // Write Apollo propagation history to file. (INPUT)
    tudat::input_output::writeDataMapToTextFile( stateHistoryMap,
                            "mainRKusm_ExampleOutput.dat",
                            cppFolder + "Output_tudatversion",
                            "",
                            std::numeric_limits< double >::digits10,
                            std::numeric_limits< double >::digits10,
                            "," );


    tudat::input_output::writeDataMapToTextFile( thrustHistoryMap,
                            "mainRKusmThrustAcc_ExampleOutput.dat",
                            cppFolder + "Output_tudatversion",
                            "",
                            std::numeric_limits< double >::digits10,
                            std::numeric_limits< double >::digits10,
                            "," );

}
