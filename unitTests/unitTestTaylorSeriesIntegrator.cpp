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
 *      161025    M. Van den Broeck Creation.
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <Eigen/Core>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/results_collector.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include <Tudat/InputOutput/matrixTextFileReader.h>
#include <Tudat/InputOutput/basicInputOutput.h>


#include "Thesis/constants.h"
#include "Thesis/spacecraft.h"
#include "Thesis/stepSizeControlTSI.h"
#include "Thesis/integrationSettings.h"
#include "Thesis/frameTransformation.h"
#include "Thesis/discreteForceModel.h"
#include "Thesis/taylorSeriesIntegrator.h"

namespace tudat
{
namespace unit_tests
{

inline bool current_test_passing()
{
  using namespace boost::unit_test;
  test_case::id_t id = framework::current_test_case().p_id;
  test_results rez = results_collector.results(id);
  return rez.passed();
}

BOOST_AUTO_TEST_SUITE( test_taylorSeriesIntegrator )

//! Test if future states are computed correctly.
BOOST_AUTO_TEST_CASE( test1_circular_zeroThrust )
{

    // Find filepath of this cpp file
    std::string cppFilePath( __FILE__ );

    // Find folder of this cpp file
    std::string cppFolder = cppFilePath.substr( 0 , cppFilePath.find_last_of("/\\")+1 );
    // std::cout << "cppFolder " << cppFolder << std::endl;


    /// Define all constants

    // Define constants
    double standardGravity = 9.81;
    double centralBodyGravitationalParameter = 3.986004418e14; // Test case = orbit around Earth



    /// Define spacecraft parameters

    // Specific impulse
    double specificImpulse = 3000.0;

    // USM state
    Vector7d initialReducedState;
    initialReducedState <<         // in USM elements
                 7433.800841666130509111099,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0;
    double initialMassOfSpacecraft = 2000.0;

    // Load thrust file
    std::string fileName = cppFolder + "zeroThrust.txt";



    /// Define integration settings

    // Order of Taylor Series Integration
    int order = 10;

    // Define initial and final time of integration
    double initialTime = 0.0; // [s]
    double finalTime = 1.828968111440366e5; // [s]

    // Define tolerance for each state element
    // Set the error tolerance
    Eigen::Matrix< double, 7, 1> errorTolerance;
    errorTolerance << 1.0e-13,
                      1.0e-13,
                      1.0e-13,
                      1.0e-16,
                      1.0e-16,
                      1.0e-16,
                      1.0e-16;

    // Safety factor
    double safetyFactor = 0.75; // [-] Must be in interval (0.0, 1.0]

    // Initial step size
    double initialStepSize = 0.01; // [s]

    // Minimum step size
    double minimumStepSize = 10; // [s]

    // Maximum step size
    double maximumStepSize = 500; // [s]

    // Factor that determines how many interpolated points are added in between two steps of the integrator
    int refineFactor = 0; // [-]



    /// Make all the pointers

    // Constants
    ConstantsPointer constantsPointer =
            boost::make_shared< Constants >( centralBodyGravitationalParameter, standardGravity  );

    // Spacecraft
    SpacecraftPointer spacecraftPointer =
            boost::make_shared< Spacecraft >( initialMassOfSpacecraft,
                                              initialReducedState,
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

    // Problem recurrence relations
    ProblemRecurrenceRelationsPointer problemRecurrenceRelationsPointer =
            boost::make_shared< ProblemRecurrenceRelations >(
                order,
                spacecraftPointer,
                constantsPointer );

    // Step size control
    StepSizeControlTSIPointer stepSizeControlTSIPointer =
            boost::make_shared< StepSizeControlTSI >( integrationSettingsPointer );

    // Make pointer of class TaylorSeriesIntegrator
    TaylorSeriesIntegratorPointer taylorSeriesIntegratorPointer =
            boost::make_shared< TaylorSeriesIntegrator >(
                problemRecurrenceRelationsPointer,
                spacecraftPointer,
                integrationSettingsPointer,
                constantsPointer,
                stepSizeControlTSIPointer );



    /// Perform integration
    taylorSeriesIntegratorPointer->integrate();


    /// Define expected ouput
    Eigen::MatrixXd expectedStateMatrix =
            input_output::readMatrixFromFile(
                cppFolder + "stateMatrix1.txt" , " \t", "#" );
    Eigen::MatrixXd expectedTimeVector =
            input_output::readMatrixFromFile(
                cppFolder + "timeVector1.txt", " \t", "#" );


    /// Get computed ouput

    Eigen::MatrixXd computedStateMatrix = taylorSeriesIntegratorPointer->getStateMatrix();
    Eigen::MatrixXd computedTimeVector = taylorSeriesIntegratorPointer->getTimeVector();

//    std::cout << "computedStateMatrix: Rows = " << computedStateMatrix.rows()
//              << ", Columns = " << computedStateMatrix.cols() << std::endl;

//    std::cout << "computedTimeVector: Rows = " << computedTimeVector.rows()
//              << ", Columns = " << computedTimeVector.cols() << std::endl;


    /// Verify output

    // Check if computed state matrix matches expected value.
    for ( int i = 0; i < expectedStateMatrix.rows(); i++ )
    {
        for ( int j = 0; j < expectedStateMatrix.cols(); j++ )
        {
            // std::cout << "i = " << i << ", j = " << j << std::endl;
            BOOST_CHECK_CLOSE_FRACTION( computedStateMatrix( i, j ), expectedStateMatrix( i, j ),
                                        1.0e-11 );
        }
    }

    // Check if computed time vector matches expected value.
    for ( int i = 0; i < expectedTimeVector.rows(); i++ )
    {
        // std::cout << "i = " << i << std::endl;
        BOOST_CHECK_CLOSE_FRACTION( computedTimeVector( i, 0 ), expectedTimeVector( i, 0 ),
                                    1.0e-15 );
    }


}

//! Test if future states are computed correctly.
BOOST_AUTO_TEST_CASE( test2_elliptic_zeroThrust )
{

    // Find filepath of this cpp file
    std::string cppFilePath( __FILE__ );

    // Find folder of this cpp file
    std::string cppFolder = cppFilePath.substr( 0 , cppFilePath.find_last_of("/\\")+1 );
    // std::cout << "cppFolder " << cppFolder << std::endl;

    /// Define all constants

    // Define constants
    double standardGravity = 9.81;
    double centralBodyGravitationalParameter = 3.986004418e14; // Test case = orbit around Earth



    /// Define spacecraft parameters

    // Specific impulse
    double specificImpulse = 3000.0;

    // USM state
    Vector7d initialReducedState;
    initialReducedState <<         // in USM elements
                    2.938467791467713e+03,
                   -1.444231139074881e+03,
                   -1.011261530497428e+03,
                    9.802738934708755e-01,
                    1.728487359030509e-01,
                    9.006554637692694e-02,
                    3.278117801412362e-02;
    double initialMassOfSpacecraft = 2000.0;

    // Load thrust file
    std::string fileName = cppFolder + "zeroThrust.txt";



    /// Define integration settings

    // Order of Taylor Series Integration
    int order = 10;

    // Define initial and final time of integration
    double initialTime = 0.0; // [s]
    double finalTime = 1927901.666656073415652; // [s]

    // Define tolerance for each state element
    // Set the error tolerance
    Eigen::Matrix< double, 7, 1> errorTolerance;
    errorTolerance << 1.0e-13,
                      1.0e-13,
                      1.0e-13,
                      1.0e-16,
                      1.0e-16,
                      1.0e-16,
                      1.0e-16;

    // Safety factor
    double safetyFactor = 0.75; // [-] Must be in interval (0.0, 1.0]

    // Initial step size
    double initialStepSize = 0.01; // [s]

    // Minimum step size
    double minimumStepSize = 10; // [s]

    // Maximum step size
    double maximumStepSize = 500; // [s]

    // Factor that determines how many interpolated points are added in between two steps of the integrator
    int refineFactor = 0; // [-]



    /// Make all the pointers

    // Constants
    ConstantsPointer constantsPointer =
            boost::make_shared< Constants >( centralBodyGravitationalParameter, standardGravity  );

    // Spacecraft
    SpacecraftPointer spacecraftPointer =
            boost::make_shared< Spacecraft >( initialMassOfSpacecraft,
                                              initialReducedState,
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

    // Problem recurrence relations
    ProblemRecurrenceRelationsPointer problemRecurrenceRelationsPointer =
            boost::make_shared< ProblemRecurrenceRelations >(
                order,
                spacecraftPointer,
                constantsPointer );

    // Step size control
    StepSizeControlTSIPointer stepSizeControlTSIPointer =
            boost::make_shared< StepSizeControlTSI >( integrationSettingsPointer );

    // Make pointer of class TaylorSeriesIntegrator
    TaylorSeriesIntegratorPointer taylorSeriesIntegratorPointer =
            boost::make_shared< TaylorSeriesIntegrator >(
                problemRecurrenceRelationsPointer,
                spacecraftPointer,
                integrationSettingsPointer,
                constantsPointer,
                stepSizeControlTSIPointer );



    /// Perform integration
    taylorSeriesIntegratorPointer->integrate();


    /// Define expected ouput
    Eigen::MatrixXd expectedStateMatrix =
            input_output::readMatrixFromFile(
                cppFolder + "stateMatrix2.txt", " \t", "#" );
    Eigen::MatrixXd expectedTimeVector =
            input_output::readMatrixFromFile(
                cppFolder + "timeVector2.txt", " \t", "#" );


    /// Get computed ouput

    Eigen::MatrixXd computedStateMatrix = taylorSeriesIntegratorPointer->getStateMatrix();
    Eigen::MatrixXd computedTimeVector = taylorSeriesIntegratorPointer->getTimeVector();

//    std::cout << "computedStateMatrix: Rows = " << computedStateMatrix.rows()
//              << ", Columns = " << computedStateMatrix.cols() << std::endl;

//    std::cout << "computedTimeVector: Rows = " << computedTimeVector.rows()
//              << ", Columns = " << computedTimeVector.cols() << std::endl;


    /// Verify output

    // Check if computed state matrix matches expected value.
    for ( int i = 0; i < expectedStateMatrix.rows(); i++ )
    {
        for ( int j = 0; j < expectedStateMatrix.cols(); j++ )
        {
            // std::cout << "i = " << i << ", j = " << j << std::endl;
            BOOST_CHECK_CLOSE_FRACTION( computedStateMatrix( i, j ), expectedStateMatrix( i, j ),
                                        1.3e-09 );
            BOOST_REQUIRE_MESSAGE( current_test_passing(), "location of failure: i = " << i << ", j = " << j );
        }
    }

    // Check if computed time vector matches expected value.
    for ( int i = 0; i < expectedTimeVector.rows(); i++ )
    {
        // std::cout << "i = " << i << std::endl;
        BOOST_CHECK_CLOSE_FRACTION( computedTimeVector( i, 0 ), expectedTimeVector( i, 0 ),
                                    1.0e-14 );
    }


}

//! Test if future states are computed correctly.
BOOST_AUTO_TEST_CASE( test3_circular_tangentialThrust )
{

    // Find filepath of this cpp file
    std::string cppFilePath( __FILE__ );

    // Find folder of this cpp file
    std::string cppFolder = cppFilePath.substr( 0 , cppFilePath.find_last_of("/\\")+1 );
    // std::cout << "cppFolder " << cppFolder << std::endl;

    /// Define all constants

    // Define constants
    double standardGravity = 9.81;
    double centralBodyGravitationalParameter = 3.986004418e14; // Test case = orbit around Earth



    /// Define spacecraft parameters

    // Specific impulse
    double specificImpulse = 3000.0;

    // USM state
    Vector7d initialReducedState;
    initialReducedState <<         // in USM elements
                   7433.800841666130509111099,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      1.0;
    double initialMassOfSpacecraft = 2000.0;

    // Load thrust file
    std::string fileName = cppFolder + "tangentialThrust.txt";



    /// Define integration settings

    // Order of Taylor Series Integration
    int order = 10;

    // Define initial and final time of integration
    double initialTime = 0.0; // [s]
    double finalTime = 1.828968111440366e5; // [s]

    // Define tolerance for each state element
    // Set the error tolerance
    Eigen::Matrix< double, 7, 1> errorTolerance;
    errorTolerance << 1.0e-13,
                      1.0e-13,
                      1.0e-13,
                      1.0e-16,
                      1.0e-16,
                      1.0e-16,
                      1.0e-16;

    // Safety factor
    double safetyFactor = 0.75; // [-] Must be in interval (0.0, 1.0]

    // Initial step size
    double initialStepSize = 0.01; // [s]

    // Minimum step size
    double minimumStepSize = 10; // [s]

    // Maximum step size
    double maximumStepSize = 500; // [s]

    // Factor that determines how many interpolated points are added in between two steps of the integrator
    int refineFactor = 0; // [-]



    /// Make all the pointers

    // Constants
    ConstantsPointer constantsPointer =
            boost::make_shared< Constants >( centralBodyGravitationalParameter, standardGravity  );

    // Spacecraft
    SpacecraftPointer spacecraftPointer =
            boost::make_shared< Spacecraft >( initialMassOfSpacecraft,
                                              initialReducedState,
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

    // Problem recurrence relations
    ProblemRecurrenceRelationsPointer problemRecurrenceRelationsPointer =
            boost::make_shared< ProblemRecurrenceRelations >(
                order,
                spacecraftPointer,
                constantsPointer );

    // Step size control
    StepSizeControlTSIPointer stepSizeControlTSIPointer =
            boost::make_shared< StepSizeControlTSI >( integrationSettingsPointer );

    // Make pointer of class TaylorSeriesIntegrator
    TaylorSeriesIntegratorPointer taylorSeriesIntegratorPointer =
            boost::make_shared< TaylorSeriesIntegrator >(
                problemRecurrenceRelationsPointer,
                spacecraftPointer,
                integrationSettingsPointer,
                constantsPointer,
                stepSizeControlTSIPointer );



    /// Perform integration
    taylorSeriesIntegratorPointer->integrate();


    /// Define expected ouput
    Eigen::MatrixXd expectedStateMatrix =
            input_output::readMatrixFromFile(
                cppFolder + "stateMatrix3.txt", " \t", "#" );
    Eigen::MatrixXd expectedTimeVector =
            input_output::readMatrixFromFile(
                cppFolder + "timeVector3.txt", " \t", "#" );


    /// Get computed ouput

    Eigen::MatrixXd computedStateMatrix = taylorSeriesIntegratorPointer->getStateMatrix();
    Eigen::MatrixXd computedTimeVector = taylorSeriesIntegratorPointer->getTimeVector();

//    std::cout << "expectedStateMatrix: Rows = " << expectedStateMatrix.rows()
//              << ", Columns = " << expectedStateMatrix.cols() << std::endl;

//    std::cout << "computedStateMatrix: Rows = " << computedStateMatrix.rows()
//              << ", Columns = " << computedStateMatrix.cols() << std::endl;

//    std::cout << "expectedTimeVector: Rows = " << expectedTimeVector.rows()
//              << ", Columns = " << expectedTimeVector.cols() << std::endl;

//    std::cout << "computedTimeVector: Rows = " << computedTimeVector.rows()
//              << ", Columns = " << computedTimeVector.cols() << std::endl;


    /// Verify output

    // Check if computed state matrix matches expected value.
    for ( int i = 0; i < expectedStateMatrix.rows(); i++ )
    {
        for ( int j = 0; j < expectedStateMatrix.cols(); j++ )
        {
            // std::cout << "i = " << i << ", j = " << j << std::endl;
            BOOST_CHECK_CLOSE_FRACTION( computedStateMatrix( i, j ), expectedStateMatrix( i, j ),
                                        5.0e-08 );
            BOOST_REQUIRE_MESSAGE( current_test_passing(), "location of failure: i = " << i << ", j = " << j );
        }
    }

    // Check if computed time vector matches expected value.
    for ( int i = 0; i < expectedTimeVector.rows(); i++ )
    {
        // std::cout << "i = " << i << std::endl;
        BOOST_CHECK_CLOSE_FRACTION( computedTimeVector( i, 0 ), expectedTimeVector( i, 0 ),
                                    1.0e-13 );
        // BOOST_REQUIRE_MESSAGE( current_test_passing(), "location of failure: i = " << i );
    }


}

//! Test if future states are computed correctly.
BOOST_AUTO_TEST_CASE( test4_circular_crosstrackThrust )
{

    // Find filepath of this cpp file
    std::string cppFilePath( __FILE__ );

    // Find folder of this cpp file
    std::string cppFolder = cppFilePath.substr( 0 , cppFilePath.find_last_of("/\\")+1 );
    // std::cout << "cppFolder " << cppFolder << std::endl;

    /// Define all constants

    // Define constants
    double standardGravity = 9.81;
    double centralBodyGravitationalParameter = 3.986004418e14; // Test case = orbit around Earth



    /// Define spacecraft parameters

    // Specific impulse
    double specificImpulse = 3000.0;

    // USM state
    Vector7d initialReducedState;
    initialReducedState <<         // in USM elements
                   7433.800841666130509111099,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      1.0;
    double initialMassOfSpacecraft = 2000.0;

    // Load thrust file
    std::string fileName = cppFolder + "crosstrackThrust.txt";



    /// Define integration settings

    // Order of Taylor Series Integration
    int order = 10;

    // Define initial and final time of integration
    double initialTime = 0.0; // [s]
    double finalTime = 1.828968111440366e5; // [s]

    // Define tolerance for each state element
    // Set the error tolerance
    Eigen::Matrix< double, 7, 1> errorTolerance;
    errorTolerance << 1.0e-13,
                      1.0e-13,
                      1.0e-13,
                      1.0e-16,
                      1.0e-16,
                      1.0e-16,
                      1.0e-16;

    // Safety factor
    double safetyFactor = 0.75; // [-] Must be in interval (0.0, 1.0]

    // Initial step size
    double initialStepSize = 0.01; // [s]

    // Minimum step size
    double minimumStepSize = 10; // [s]

    // Maximum step size
    double maximumStepSize = 500; // [s]

    // Factor that determines how many interpolated points are added in between two steps of the integrator
    int refineFactor = 0; // [-]



    /// Make all the pointers

    // Constants
    ConstantsPointer constantsPointer =
            boost::make_shared< Constants >( centralBodyGravitationalParameter, standardGravity  );

    // Spacecraft
    SpacecraftPointer spacecraftPointer =
            boost::make_shared< Spacecraft >( initialMassOfSpacecraft,
                                              initialReducedState,
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

    // Problem recurrence relations
    ProblemRecurrenceRelationsPointer problemRecurrenceRelationsPointer =
            boost::make_shared< ProblemRecurrenceRelations >(
                order,
                spacecraftPointer,
                constantsPointer );

    // Step size control
    StepSizeControlTSIPointer stepSizeControlTSIPointer =
            boost::make_shared< StepSizeControlTSI >( integrationSettingsPointer );

    // Make pointer of class TaylorSeriesIntegrator
    TaylorSeriesIntegratorPointer taylorSeriesIntegratorPointer =
            boost::make_shared< TaylorSeriesIntegrator >(
                problemRecurrenceRelationsPointer,
                spacecraftPointer,
                integrationSettingsPointer,
                constantsPointer,
                stepSizeControlTSIPointer );



    /// Perform integration
    taylorSeriesIntegratorPointer->integrate();


    /// Define expected ouput
    Eigen::MatrixXd expectedStateMatrix =
            input_output::readMatrixFromFile(
                cppFolder + "stateMatrix4.txt", " \t", "#" );
    Eigen::MatrixXd expectedTimeVector =
            input_output::readMatrixFromFile(
                cppFolder +  "timeVector4.txt", " \t", "#" );


    /// Get computed ouput

    Eigen::MatrixXd computedStateMatrix = taylorSeriesIntegratorPointer->getStateMatrix();
    Eigen::MatrixXd computedTimeVector = taylorSeriesIntegratorPointer->getTimeVector();

//    std::cout << "expectedStateMatrix: Rows = " << expectedStateMatrix.rows()
//              << ", Columns = " << expectedStateMatrix.cols() << std::endl;

//    std::cout << "computedStateMatrix: Rows = " << computedStateMatrix.rows()
//              << ", Columns = " << computedStateMatrix.cols() << std::endl;

//    std::cout << "expectedTimeVector: Rows = " << expectedTimeVector.rows()
//              << ", Columns = " << expectedTimeVector.cols() << std::endl;

//    std::cout << "computedTimeVector: Rows = " << computedTimeVector.rows()
//              << ", Columns = " << computedTimeVector.cols() << std::endl;


    /// Verify output

    // Check if computed state matrix matches expected value.
    for ( int i = 0; i < expectedStateMatrix.rows(); i++ )
    {
        for ( int j = 0; j < expectedStateMatrix.cols(); j++ )
        {
            // std::cout << "i = " << i << ", j = " << j << std::endl;
            BOOST_CHECK_CLOSE_FRACTION( computedStateMatrix( i, j ), expectedStateMatrix( i, j ),
                                        1.2e-07 );
            // BOOST_REQUIRE_MESSAGE( current_test_passing(), "location of failure: i = " << i << ", j = " << j );
        }
    }

    // Check if computed time vector matches expected value.
    for ( int i = 0; i < expectedTimeVector.rows(); i++ )
    {
        // std::cout << "i = " << i << std::endl;
        BOOST_CHECK_CLOSE_FRACTION( computedTimeVector( i, 0 ), expectedTimeVector( i, 0 ),
                                    1.0e-12 );
        // BOOST_REQUIRE_MESSAGE( current_test_passing(), "location of failure: i = " << i );
    }


}

//! Test if future states are computed correctly.
BOOST_AUTO_TEST_CASE( test5_circular_outofplaneThrust )
{

    // Find filepath of this cpp file
    std::string cppFilePath( __FILE__ );

    // Find folder of this cpp file
    std::string cppFolder = cppFilePath.substr( 0 , cppFilePath.find_last_of("/\\")+1 );
    // std::cout << "cppFolder " << cppFolder << std::endl;

    /// Define all constants

    // Define constants
    double standardGravity = 9.81;
    double centralBodyGravitationalParameter = 3.986004418e14; // Test case = orbit around Earth



    /// Define spacecraft parameters

    // Specific impulse
    double specificImpulse = 3000.0;

    // USM state
    Vector7d initialReducedState;
    initialReducedState <<         // in USM elements
                   7433.800841666130509111099,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      1.0;
    double initialMassOfSpacecraft = 2000.0;

    // Load thrust file
    std::string fileName = cppFolder + "outofplaneThrust.txt";



    /// Define integration settings

    // Order of Taylor Series Integration
    int order = 10;

    // Define initial and final time of integration
    double initialTime = 0.0; // [s]
    double finalTime = 1.828968111440366e5; // [s]

    // Define tolerance for each state element
    // Set the error tolerance
    Eigen::Matrix< double, 7, 1> errorTolerance;
    errorTolerance << 1.0e-13,
                      1.0e-13,
                      1.0e-13,
                      1.0e-16,
                      1.0e-16,
                      1.0e-16,
                      1.0e-16;

    // Safety factor
    double safetyFactor = 0.75; // [-] Must be in interval (0.0, 1.0]

    // Initial step size
    double initialStepSize = 0.01; // [s]

    // Minimum step size
    double minimumStepSize = 10; // [s]

    // Maximum step size
    double maximumStepSize = 500; // [s]

    // Factor that determines how many interpolated points are added in between two steps of the integrator
    int refineFactor = 0; // [-]



    /// Make all the pointers

    // Constants
    ConstantsPointer constantsPointer =
            boost::make_shared< Constants >( centralBodyGravitationalParameter, standardGravity  );

    // Spacecraft
    SpacecraftPointer spacecraftPointer =
            boost::make_shared< Spacecraft >( initialMassOfSpacecraft,
                                              initialReducedState,
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

    // Problem recurrence relations
    ProblemRecurrenceRelationsPointer problemRecurrenceRelationsPointer =
            boost::make_shared< ProblemRecurrenceRelations >(
                order,
                spacecraftPointer,
                constantsPointer );

    // Step size control
    StepSizeControlTSIPointer stepSizeControlTSIPointer =
            boost::make_shared< StepSizeControlTSI >( integrationSettingsPointer );

    // Make pointer of class TaylorSeriesIntegrator
    TaylorSeriesIntegratorPointer taylorSeriesIntegratorPointer =
            boost::make_shared< TaylorSeriesIntegrator >(
                problemRecurrenceRelationsPointer,
                spacecraftPointer,
                integrationSettingsPointer,
                constantsPointer,
                stepSizeControlTSIPointer );



    /// Perform integration
    taylorSeriesIntegratorPointer->integrate();


    /// Define expected ouput
    Eigen::MatrixXd expectedStateMatrix =
            input_output::readMatrixFromFile(
                cppFolder + "stateMatrix5.txt", " \t", "#" );
    Eigen::MatrixXd expectedTimeVector =
            input_output::readMatrixFromFile(
                cppFolder + "timeVector5.txt", " \t", "#" );


    /// Get computed ouput

    Eigen::MatrixXd computedStateMatrix = taylorSeriesIntegratorPointer->getStateMatrix();
    Eigen::MatrixXd computedTimeVector = taylorSeriesIntegratorPointer->getTimeVector();

//    std::cout << "expectedStateMatrix: Rows = " << expectedStateMatrix.rows()
//              << ", Columns = " << expectedStateMatrix.cols() << std::endl;

//    std::cout << "computedStateMatrix: Rows = " << computedStateMatrix.rows()
//              << ", Columns = " << computedStateMatrix.cols() << std::endl;

//    std::cout << "expectedTimeVector: Rows = " << expectedTimeVector.rows()
//              << ", Columns = " << expectedTimeVector.cols() << std::endl;

//    std::cout << "computedTimeVector: Rows = " << computedTimeVector.rows()
//              << ", Columns = " << computedTimeVector.cols() << std::endl;


    /// Verify output

    // Check if computed state matrix matches expected value.
    for ( int i = 0; i < expectedStateMatrix.rows(); i++ )
    {
        for ( int j = 0; j < expectedStateMatrix.cols(); j++ )
        {
            // std::cout << "i = " << i << ", j = " << j << std::endl;
            BOOST_CHECK_CLOSE_FRACTION( computedStateMatrix( i, j ), expectedStateMatrix( i, j ),
                                        5.0e-08 );
            // BOOST_REQUIRE_MESSAGE( current_test_passing(), "location of failure: i = " << i << ", j = " << j );
        }
    }

    // Check if computed time vector matches expected value.
    for ( int i = 0; i < expectedTimeVector.rows(); i++ )
    {
        // std::cout << "i = " << i << std::endl;
        BOOST_CHECK_CLOSE_FRACTION( computedTimeVector( i, 0 ), expectedTimeVector( i, 0 ),
                                    1.0e-12 );
        // BOOST_REQUIRE_MESSAGE( current_test_passing(), "location of failure: i = " << i );
    }


}

//! Test if future states are computed correctly.
BOOST_AUTO_TEST_CASE( test6_circular_nonconstantThrust )
{

    // Find filepath of this cpp file
    std::string cppFilePath( __FILE__ );

    // Find folder of this cpp file
    std::string cppFolder = cppFilePath.substr( 0 , cppFilePath.find_last_of("/\\")+1 );
    // std::cout << "cppFolder " << cppFolder << std::endl;

    /// Define all constants

    // Define constants
    double standardGravity = 9.81;
    double centralBodyGravitationalParameter = 3.986004418e14; // Test case = orbit around Earth



    /// Define spacecraft parameters

    // Specific impulse
    double specificImpulse = 3000.0;

    // USM state
    Vector7d initialReducedState;
    initialReducedState <<         // in USM elements
                   7433.800841666130509111099,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      1.0;
    double initialMassOfSpacecraft = 2000.0;

    // Load thrust file
    std::string fileName = cppFolder + "nonconstantThrust.txt";



    /// Define integration settings

    // Order of Taylor Series Integration
    int order = 10;

    // Define initial and final time of integration
    double initialTime = 0.0; // [s]
    double finalTime = 1.828968111440366e5; // [s]

    // Define tolerance for each state element
    // Set the error tolerance
    Eigen::Matrix< double, 7, 1> errorTolerance;
    errorTolerance << 1.0e-13,
                      1.0e-13,
                      1.0e-13,
                      1.0e-16,
                      1.0e-16,
                      1.0e-16,
                      1.0e-16;

    // Safety factor
    double safetyFactor = 0.75; // [-] Must be in interval (0.0, 1.0]

    // Initial step size
    double initialStepSize = 0.01; // [s]

    // Minimum step size
    double minimumStepSize = 10; // [s]

    // Maximum step size
    double maximumStepSize = 500; // [s]

    // Factor that determines how many interpolated points are added in between two steps of the integrator
    int refineFactor = 0; // [-]



    /// Make all the pointers

    // Constants
    ConstantsPointer constantsPointer =
            boost::make_shared< Constants >( centralBodyGravitationalParameter, standardGravity  );

    // Spacecraft
    SpacecraftPointer spacecraftPointer =
            boost::make_shared< Spacecraft >( initialMassOfSpacecraft,
                                              initialReducedState,
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

    // Problem recurrence relations
    ProblemRecurrenceRelationsPointer problemRecurrenceRelationsPointer =
            boost::make_shared< ProblemRecurrenceRelations >(
                order,
                spacecraftPointer,
                constantsPointer );

    // Step size control
    StepSizeControlTSIPointer stepSizeControlTSIPointer =
            boost::make_shared< StepSizeControlTSI >( integrationSettingsPointer );

    // Make pointer of class TaylorSeriesIntegrator
    TaylorSeriesIntegratorPointer taylorSeriesIntegratorPointer =
            boost::make_shared< TaylorSeriesIntegrator >(
                problemRecurrenceRelationsPointer,
                spacecraftPointer,
                integrationSettingsPointer,
                constantsPointer,
                stepSizeControlTSIPointer );



    /// Perform integration
    taylorSeriesIntegratorPointer->integrate();


    /// Define expected ouput
    Eigen::MatrixXd expectedStateMatrix =
            input_output::readMatrixFromFile(
                cppFolder + "stateMatrix6.txt", " \t", "#" );
    Eigen::MatrixXd expectedTimeVector =
            input_output::readMatrixFromFile(
                cppFolder + "timeVector6.txt", " \t", "#" );


    /// Get computed ouput

    Eigen::MatrixXd computedStateMatrix = taylorSeriesIntegratorPointer->getStateMatrix();
    Eigen::MatrixXd computedTimeVector = taylorSeriesIntegratorPointer->getTimeVector();

//    std::cout << "expectedStateMatrix: Rows = " << expectedStateMatrix.rows()
//              << ", Columns = " << expectedStateMatrix.cols() << std::endl;

//    std::cout << "computedStateMatrix: Rows = " << computedStateMatrix.rows()
//              << ", Columns = " << computedStateMatrix.cols() << std::endl;

//    std::cout << "expectedTimeVector: Rows = " << expectedTimeVector.rows()
//              << ", Columns = " << expectedTimeVector.cols() << std::endl;

//    std::cout << "computedTimeVector: Rows = " << computedTimeVector.rows()
//              << ", Columns = " << computedTimeVector.cols() << std::endl;


    /// Verify output

    // Check if computed state matrix matches expected value.
    for ( int i = 0; i < expectedStateMatrix.rows(); i++ )
    {
        for ( int j = 0; j < expectedStateMatrix.cols(); j++ )
        {
            // std::cout << "i = " << i << ", j = " << j << std::endl;
            BOOST_CHECK_CLOSE_FRACTION( computedStateMatrix( i, j ), expectedStateMatrix( i, j ),
                                        5.0e-07 );
            // BOOST_REQUIRE_MESSAGE( current_test_passing(), "location of failure: i = " << i << ", j = " << j );
        }
    }

    // Check if computed time vector matches expected value.
    for ( int i = 0; i < expectedTimeVector.rows(); i++ )
    {
        // std::cout << "i = " << i << std::endl;
        BOOST_CHECK_CLOSE_FRACTION( computedTimeVector( i, 0 ), expectedTimeVector( i, 0 ),
                                    1.0e-13 );
        // BOOST_REQUIRE_MESSAGE( current_test_passing(), "location of failure: i = " << i );
    }


}

//! Test if future states are computed correctly.
BOOST_AUTO_TEST_CASE( test7_elliptic_nonconstantThrust_refine )
{

    // Find filepath of this cpp file
    std::string cppFilePath( __FILE__ );

    // Find folder of this cpp file
    std::string cppFolder = cppFilePath.substr( 0 , cppFilePath.find_last_of("/\\")+1 );
    // std::cout << "cppFolder " << cppFolder << std::endl;

    /// Define all constants

    // Define constants
    double standardGravity = 9.81;
    double centralBodyGravitationalParameter = 3.986004418e14; // Test case = orbit around Earth



    /// Define spacecraft parameters

    // Specific impulse
    double specificImpulse = 3000.0;

    // USM state
    Vector7d initialReducedState;
    initialReducedState <<         // in USM elements
            2.938467791467713e+03,
           -1.444231139074881e+03,
           -1.011261530497428e+03,
            9.802738934708755e-01,
            1.728487359030509e-01,
            9.006554637692694e-02,
            3.278117801412362e-02;
    double initialMassOfSpacecraft = 2000.0;

    // Load thrust file
    std::string fileName = cppFolder + "nonconstantThrust.txt";



    /// Define integration settings

    // Order of Taylor Series Integration
    int order = 30;

    // Define initial and final time of integration
    double initialTime = 0.0; // [s]
    double finalTime = 963950.833328036707826; // [s]

    // Define tolerance four each state element
    // Set the error tolerance
    Eigen::Matrix< double, 7, 1> errorTolerance;
    errorTolerance << 1.0e-13,
                      1.0e-13,
                      1.0e-13,
                      1.0e-16,
                      1.0e-16,
                      1.0e-16,
                      1.0e-16;

    // Safety factor
    double safetyFactor = 0.75; // [-] Must be in interval (0.0, 1.0]

    // Initial step size
    double initialStepSize = 0.01; // [s]

    // Minimum step size
    double minimumStepSize = 10; // [s]

    // Maximum step size
    double maximumStepSize = 500; // [s]

    // Factor that determines how many interpolated points are added in between two steps of the integrator
    int refineFactor = 5; // [-]



    /// Make all the pointers

    // Constants
    ConstantsPointer constantsPointer =
            boost::make_shared< Constants >( centralBodyGravitationalParameter, standardGravity  );

    // Spacecraft
    SpacecraftPointer spacecraftPointer =
            boost::make_shared< Spacecraft >( initialMassOfSpacecraft,
                                              initialReducedState,
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

    // Problem recurrence relations
    ProblemRecurrenceRelationsPointer problemRecurrenceRelationsPointer =
            boost::make_shared< ProblemRecurrenceRelations >(
                order,
                spacecraftPointer,
                constantsPointer );

    // Step size control
    StepSizeControlTSIPointer stepSizeControlTSIPointer =
            boost::make_shared< StepSizeControlTSI >( integrationSettingsPointer );

    // Make pointer of class TaylorSeriesIntegrator
    TaylorSeriesIntegratorPointer taylorSeriesIntegratorPointer =
            boost::make_shared< TaylorSeriesIntegrator >(
                problemRecurrenceRelationsPointer,
                spacecraftPointer,
                integrationSettingsPointer,
                constantsPointer,
                stepSizeControlTSIPointer );



    /// Perform integration
    taylorSeriesIntegratorPointer->integrate();


    /// Define expected ouput
    Eigen::MatrixXd expectedStateMatrix =
            input_output::readMatrixFromFile(
                cppFolder + "stateMatrix7.txt", " \t", "#" );
    Eigen::MatrixXd expectedTimeVector =
            input_output::readMatrixFromFile(
                cppFolder + "timeVector7.txt", " \t", "#" );


    /// Get computed ouput

    Eigen::MatrixXd computedStateMatrix = taylorSeriesIntegratorPointer->getStateMatrix();
    Eigen::MatrixXd computedTimeVector = taylorSeriesIntegratorPointer->getTimeVector();

//    std::cout << "expectedStateMatrix: Rows = " << expectedStateMatrix.rows()
//              << ", Columns = " << expectedStateMatrix.cols() << std::endl;

//    std::cout << "computedStateMatrix: Rows = " << computedStateMatrix.rows()
//              << ", Columns = " << computedStateMatrix.cols() << std::endl;

//    std::cout << "expectedTimeVector: Rows = " << expectedTimeVector.rows()
//              << ", Columns = " << expectedTimeVector.cols() << std::endl;

//    std::cout << "computedTimeVector: Rows = " << computedTimeVector.rows()
//              << ", Columns = " << computedTimeVector.cols() << std::endl;


    /// Verify output

    // Check if computed state matrix matches expected value.
    for ( int i = 0; i < expectedStateMatrix.rows(); i++ )
    {
        for ( int j = 0; j < expectedStateMatrix.cols(); j++ )
        {
            // std::cout << "i = " << i << ", j = " << j << std::endl;
            BOOST_CHECK_CLOSE_FRACTION( computedStateMatrix( i, j ), expectedStateMatrix( i, j ),
                                        1.7e-09 );
            // BOOST_REQUIRE_MESSAGE( current_test_passing(), "location of failure: i = " << i << ", j = " << j );
        }
    }

    // Check if computed time vector matches expected value.
    for ( int i = 0; i < expectedTimeVector.rows(); i++ )
    {
        // std::cout << "i = " << i << std::endl;
        BOOST_CHECK_CLOSE_FRACTION( computedTimeVector( i, 0 ), expectedTimeVector( i, 0 ),
                                    2.3e-16 );
        // BOOST_REQUIRE_MESSAGE( current_test_passing(), "location of failure: i = " << i );
    }


}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
