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
 *      251017    M. Van den Broeck Creation.
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
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "Thesis/constants.h"
#include "Thesis/spacecraft.h"
#include "Thesis/stepSizeControlTSI.h"
#include "Thesis/integrationSettings.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_stepSizeControlTSI )

//! Test if next step size is computed correctly.
BOOST_AUTO_TEST_CASE( test1 )
{
    // Define initial and final time of integration
    double initialTime; // It does not matter what these numbers are for this test case.
    double finalTime;

    // Construct integrationSettings object
    IntegrationSettingsPointer integrationSettingsPointer =
            boost::make_shared< IntegrationSettings >( initialTime, finalTime );

    // Set the error tolerance
    Eigen::Matrix< double, 7, 1> errorTolerance;
    errorTolerance << 1.0e-16,
                      1.0e-16,
                      1.0e-16,
                      1.0e-19,
                      1.0e-19,
                      1.0e-19,
                      1.0e-19;
    integrationSettingsPointer->setErrorTolerance( errorTolerance ); // Set in pointer

    // Set safety factor
    double safetyFactor = 0.1;
    integrationSettingsPointer->setSafeftyFactor( safetyFactor );

    // Construct stepSizeControlTSI object
    StepSizeControlTSIPointer stepSizeControlTSIPointer =
           boost::make_shared< StepSizeControlTSI >( integrationSettingsPointer );

    // Define derivatives
    Eigen::MatrixXd derivatives( 5, 7);
    // Data from Matlab
    derivatives <<           7433.80084166613,                         0.0,                         0,        0.0,                     0.0,                        0.0,            1.0,
                            -0.005,                                    0.0,                      0.01,        0.0,                     0.0,        0.00051530575638889,       0.0,
                            -2.02848251913308e-07,        -5.1530575638889e-06,  4.09059521014254e-07,        0.0,                    0.0,        1.73298211562457e-10, -1.32770011283763e-07,
                            4.54929098599676e-11,         -2.8220962258952e-10,  -1.86097750199679e-09,       0.0,                    0.0,       -2.28009522242284e-11,  -8.93015659900336e-14,
                            -3.53057259164559e-16,        5.26135234416357e-13,  -1.09161632842322e-13,       0.0,                   0.0,        -8.51492361377252e-17,  2.93550896801263e-15;

    double previousStepSize = 1.0;

    // Define expected output
    double expectedNextStepSize = 0.00016368735504862;

    // Determine next step size
    double computedNextStepSize;
    computedNextStepSize = stepSizeControlTSIPointer->determineNextStepSize( derivatives, previousStepSize );

    // Check if computed next step size matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedNextStepSize, expectedNextStepSize,
                                3e-15 );

}

//! Test if next step size is computed correctly.
BOOST_AUTO_TEST_CASE( test2 )
{
    // Define initial and final time of integration
    double initialTime; // It does not matter what these numbers are for this test case.
    double finalTime;

    // Construct integrationSettings object
    IntegrationSettingsPointer integrationSettingsPointer =
            boost::make_shared< IntegrationSettings >( initialTime, finalTime );

    // Set the error tolerance
    Eigen::Matrix< double, 7, 1> errorTolerance;
    errorTolerance << 1.0e-13,
                      1.0e-13,
                      1.0e-13,
                      1.0e-16,
                      1.0e-16,
                      1.0e-16,
                      1.0e-16;
    integrationSettingsPointer->setErrorTolerance( errorTolerance ); // Set in pointer

    // Set safety factor
    double safetyFactor = 0.75;
    integrationSettingsPointer->setSafeftyFactor( safetyFactor );

    // Construct stepSizeControlTSI object
    StepSizeControlTSIPointer stepSizeControlTSIPointer =
           boost::make_shared< StepSizeControlTSI >( integrationSettingsPointer );

    // Define derivatives
    Eigen::MatrixXd derivatives( 4, 7);
    // Data from Matlab
    derivatives <<           2938.46779146771,           -1444.23113907488,       -1011.26153049743,  0.980273893470875,   0.172848735903051,      0.0900655463769269,  0.0327811780141236,
                            -0.0031502596668462,        -0.00559822074997089,   -0.00591062665231279, 1.37255898954843e-05,  -7.78416884377887e-05, 2.60309109790132e-06, -7.15193401226721e-06,
                            -1.51140185673749e-07,      2.24830773465602e-07,  -7.04564331544126e-07, -3.29758091411486e-09,  6.2871355543598e-10, -3.23209424510199e-10, 4.48158765743693e-12,
                            2.25970697002021e-11,      9.31862583082247e-11,   5.23652832140337e-11, 3.97906603876116e-14,  3.19329048620848e-13, -1.56611792466127e-15,  3.0946077076142e-14;

    double previousStepSize = 0.01;

    // Define expected output
    double expectedNextStepSize = 0.000130606113685119;

    // Determine next step size
    double computedNextStepSize;
    computedNextStepSize = stepSizeControlTSIPointer->determineNextStepSize( derivatives, previousStepSize );

    // Check if computed next step size matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedNextStepSize, expectedNextStepSize,
                                2.0e-15 );

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
